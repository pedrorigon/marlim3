"""Check synchronized outlet sweeps and their ordered numerical operations."""
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = 'af37be616068c9ccc760e15369023b051dace24d'
LEX = re.compile(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|\w+|\S')


def tokens(source):
    return [match[0] for match in LEX.finditer(source)
            if not match[0].startswith(('//', '/*'))]


def body_span(source):
    clean = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
                   lambda match: re.sub(r'[^\n]', ' ', match[0]), source)
    matches = list(re.finditer(r'(?m)^void evolveFractions\([^;{]*\)\s*\{', clean))
    assert len(matches) == 1, 'fraction evolution definition'
    start = clean.index('{', matches[0].start())
    end, depth = start + 1, 1
    while depth:
        depth += (clean[end] == '{') - (clean[end] == '}')
        end += 1
    return start + 1, end - 1


def bodies():
    source = subprocess.check_output(
        ['git', 'show', REFERENCE + ':src/core/SisProdTransient.cpp'], cwd=ROOT, text=True)
    start, end = body_span(source)
    original = source[start:end]
    updated = original
    for method, boundary, arguments in [
        ('avancalf', 'updateOutletVoidFraction', ', alfrev'),
        ('avancbet', 'updateOutletComplementaryFraction', ''),
    ]:
        pattern = (r'for \(int i = 0; i <= state\.lastCell; i\+\+\) \{\s*'
                   r'if \(i < state\.lastCell\) \{\s*state\.cells\[i\]\.' + method +
                   r'\(state\.restart, state\.lastCell\);\s*\}\s*'
                   r'if \(i == state\.lastCell\) \{\s*' + boundary +
                   r'\(state, i' + re.escape(arguments) + r'\);\s*\}\s*\}')
        replacement = (
            'for (int i = 0; i < state.lastCell; i++) {\n'
            '        state.cells[i].' + method + '(state.restart, state.lastCell);\n'
            '    }\n'
            '    if (state.lastCell >= 0) {\n'
            '        ' + boundary + '(state, state.lastCell' + arguments + ');\n'
            '    }')
        updated, count = re.subn(pattern, replacement, updated)
        assert count == 1, method
    return original, updated


def restore(source, path):
    if not path.endswith('SisProdTransient.cpp') or 'void captureValveAreaRatios(' not in source:
        return source
    start, end = body_span(source)
    original, updated = bodies()
    assert tokens(source[start:end]) == tokens(updated), (
        'fraction interior bounds, outlet synchronization, thread counts, pig order or restart flags')
    return source[:start] + original + source[end:]
