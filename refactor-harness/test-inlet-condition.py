#!/usr/bin/env python3
"""Check the inlet extraction and compare production bodies with controlled properties.

The directed fixture covers dispatch, call ordering, arithmetic, live state and
floating-point flags. It does not replace validation against real fluid models.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess


def function(source, name):
    match = re.search(r'(?m)^void ' + name + r'\(', source)
    if not match:
        raise AssertionError('Missing function: ' + name)
    opening = source.index('{', match.end())
    depth = 0
    for position in range(opening, len(source)):
        depth += (source[position] == '{') - (source[position] == '}')
        if depth == 0:
            return source[match.start():position+1]
    raise AssertionError('Unbalanced function: ' + name)


def tokens(source):
    pattern = r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/'
    source = re.sub(pattern, lambda match: '' if match[0].startswith(('//', '/*')) else match[0], source)
    return re.findall(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|\w+|[^\s]', source)


def verify(reference, current, old_header, new_header, old_internal, new_internal):
    names = ['refreshCompositionalInletQuality', 'refreshGasInletFractions']
    helpers = [function(current, name) for name in names]
    expanded = function(current, 'refreshInletCondition')
    mappings = [
        {'baseFluidGasMassFraction': 'titH', 'mixtureLiquidDensity': 'rlMix', 'gasMassPerLiquidVolume': 'val1'},
        {'standardGasDensity': 'rgST', 'oilDensity': 'liquidDensity', 'hydrocarbonGasMassFraction': 'quality',
         'standardComplementaryDensity': 'rcST', 'complementaryRatioTerm': 'val1',
         'liquidToGasVolumeRatio': 'val2', 'updatedInletQuality': 'titT'},
    ]
    for name, helper, mapping in zip(names, helpers, mappings):
        body = helper[helper.index('{')+1:helper.rindex('}')]
        alias = 'const Cel &inletCell = state.cells[0];'
        assert body.count(alias) == 1
        body = body.replace(alias, '').replace('inletCell.', 'state.cells[0].')
        body = re.sub(r'\bconst\s+double\b', 'double', body)
        for new, old in mapping.items(): body = re.sub(r'\b' + new + r'\b', old, body)
        call = name + '(state);'
        assert expanded.count(call) == 1
        expanded = expanded.replace(call, body)
    old_main = function(reference, 'refreshInletCondition')
    assert tokens(expanded) == tokens(old_main), 'Changed dispatcher, properties, writes or arithmetic'
    old_mini = function(reference, 'refreshFluidMiniTable')
    new_mini = function(current, 'refreshFluidMiniTable')
    assert old_mini.replace('        generateFluidMiniTables(state);', '    generateFluidMiniTables(state);') == new_mini
    unrelated = current
    for helper in helpers: unrelated = unrelated.replace(helper, '', 1)
    unrelated = unrelated.replace(function(current, 'refreshInletCondition'), old_main, 1)
    unrelated = unrelated.replace(new_mini, old_mini, 1)
    assert tokens(unrelated) == tokens(reference), 'Changed unrelated implementation'
    assert tokens(old_header) == tokens(new_header), 'Changed public declarations'
    for name in names:
        new_internal = new_internal.replace('void ' + name + '(const TransientStepState &state);', '', 1)
    assert tokens(old_internal) == tokens(new_internal), 'Changed unrelated private declarations'
    return {'inlet_body_equivalent_after_expansion': True, 'arithmetic_tokens_preserved': True,
            'mini_table_only_one_line_whitespace': True, 'unrelated_implementation_unchanged': True,
            'public_declarations_unchanged': True, 'new_private_helpers': names}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--reference-dir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    reference = args.reference_dir.resolve()
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    read = lambda base, path: (base / path).read_text()
    path = 'src/core/SisProdTransient.cpp'
    old_source, new_source = read(reference, path), read(root, path)
    proof = verify(old_source, new_source, read(reference, 'src/include/SisProdTransient.h'),
                   read(root, 'src/include/SisProdTransient.h'), read(reference, 'src/include/SisProdTransientInternal.h'),
                   read(root, 'src/include/SisProdTransientInternal.h'))
    (output / 'source-proof.json').write_text(json.dumps(proof, indent=2) + '\n')
    fixture = read(root, 'refactor-harness/inlet-condition-tests.cpp')
    compiler = os.environ.get('CXX', 'g++')
    env = dict(os.environ, ASAN_OPTIONS='detect_leaks=0:halt_on_error=1',
               UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    results = {'source_proof': proof, 'modes': {}, 'scope': 'production bodies with controlled property callbacks',
               'performance_measured': False}
    for mode, flags in [('sanitized', ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer']),
                        ('release', ['-O3', '-DNDEBUG', '-funroll-loops'])]:
        modes = {}
        for version, source in [('before', old_source), ('after', new_source)]:
            folder = output / mode / version; folder.mkdir(parents=True, exist_ok=True)
            names = [] if version == 'before' else ['refreshCompositionalInletQuality', 'refreshGasInletFractions']
            bodies = '\n\n'.join(function(source, name) for name in [*names, 'refreshInletCondition'])
            (folder / 'checks.cpp').write_text(fixture.replace('// @PRODUCTION_FUNCTIONS@', bodies))
            command = [compiler, '-std=c++20', '-ffp-contract=off', *flags, '-I' + str(root / 'src/include'),
                       str(folder / 'checks.cpp'), '-o', str(folder / 'checks')]
            with (folder / 'build.log').open('w') as log:
                log.write(shlex.join(command) + '\n'); log.flush()
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            with (folder / 'checks.log').open('w') as log:
                subprocess.run([str(folder / 'checks'), str(folder / 'state.bin')], stdout=log,
                               stderr=subprocess.STDOUT, env=env, check=True)
            data = (folder / 'state.bin').read_bytes()
            modes[version] = {'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
            print(mode, version, (folder / 'checks.log').read_text().strip(), flush=True)
        modes['bitwise_equal'] = modes['before'] == modes['after']
        assert modes['bitwise_equal'], mode + ' differs'
        results['modes'][mode] = modes
    (output / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
    print('PASS source proof and both binary comparisons')


if __name__ == '__main__': main()
