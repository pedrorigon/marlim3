"""Adapt verified context and valve names for historical calculation proofs."""
import json
from pathlib import Path
import re

FIELDS = json.loads(Path(__file__).with_name('simulation-context-fields.json').read_text())
INVERSE = {new: old for old, new in FIELDS.items() if old != new}
LEX = re.compile(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\S')
NAMES = {'SimulationContext':'varGlob1D', 'simulationContext':'vg1dSP',
         'networkContext':'vg1dRedeSimples','networkContexts':'vg1dRede',
         'convectionContext':'vg1dConvec2D','defaultContext':'vg1dTemp',
         'parametricContexts':'vg1dTramo','branchContext':'vg1dTramo'}

def tokens(source):
    return [m[0] for m in LEX.finditer(source) if not m[0].startswith(('//','/*'))]

def masked(source):
    return re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
                  lambda m:re.sub(r'[^\n]',' ',m[0]),source)

def restore(source, path):
    source = LEX.sub(lambda m:NAMES.get(m[0],m[0]),source)
    clean = masked(source)
    receiver = (r'(?:\(\s*\*\s*[\w.]*\b(?:globals|vg1dSP|DGP3D|vG1d)\s*\)\s*\.'
                r'|[\w.]*\b(?:globals|vg1dSP|DGP3D|vG1d)\s*->'
                r'|\bglobals\s*\.)\s*')
    edits=[]
    for hit in re.finditer(receiver+r'(?P<name>'+ '|'.join(INVERSE)+r')\b',clean):
        edits.append((hit.start('name'),hit.end('name'),INVERSE[hit['name']]))
    for start,end,text in reversed(edits):source=source[:start]+text+source[end:]
    if path.endswith('SisProd.cpp'):
        source=source.replace('transient.valveAreaRatiosBeforeUpdate','transient.masterValveRatiosInactive')
        source=source.replace('transient.valveAreaRatiosAfterUpdate','transient.masterValveRatiosActive')
    for new,old in [('valveAreaRatiosBeforeUpdate','masterRatio0'),('valveAreaRatiosAfterUpdate','masterRatio1')]:
        if path.endswith('SisProd.h'):
            old='masterValveRatiosInactive' if old=='masterRatio0' else 'masterValveRatiosActive'
        source=LEX.sub(lambda m:old if m[0]==new else m[0],source)
    if path.endswith('SisProd.h'):
        source=source.replace('void captureValveAreaRatios(sisprod::ValveSnapshotPhase phase);',
                              'void aberturaVal1();\n    void aberturaVal0();')
    if path.endswith('SisProdTransient.h'):
        source=source.replace('void captureValveAreaRatios(const TransientStepState &state, double *destination);',
                              'void valveOpeningLow(const TransientStepState &state);\nvoid valveOpeningHigh(const TransientStepState &state);')
    if path.endswith('SisProdTransient.cpp') and 'void captureValveAreaRatios(' in source:
        clean=masked(source)
        match=re.search(r'(?m)^void captureValveAreaRatios\(const TransientStepState &state, double \*destination\) \{',clean)
        assert match,'valve capture signature'
        start=clean.index('{',match.start());end=start+1;depth=1
        while depth:
            depth+=(clean[end]=='{')-(clean[end]=='}');end+=1
        body=source[start:end]
        expected='''{
            int celpos = state.input.master1.posic;
            destination[0] = state.cells[celpos].acsr.chk.AreaGarg / state.cells[celpos].duto.area;
            for (int i = 1; i <= state.input.nvalv; i++) {
                celpos = state.input.valv[i - 1].posicP;
                destination[i] = state.cells[celpos].acsr.chk.AreaGarg / state.cells[celpos].duto.area;
            }
        }'''
        assert tokens(body)==tokens(expected),'valve capture arithmetic or order'
        kernels=[]
        for name,field in [('valveOpeningLow','masterRatio0'),('valveOpeningHigh','masterRatio1')]:
            kernels.append('void '+name+'(const TransientStepState &state) '+body.replace('destination','state.'+field))
        source=source[:match.start()]+'\n\n'.join(kernels)+source[end:]
        source=source.replace('captureValveAreaRatios(state.step, state.step.masterRatio0)','valveOpeningLow(state.step)')
        source=source.replace('captureValveAreaRatios(state.step, state.step.masterRatio1)','valveOpeningHigh(state.step)')
    return source
