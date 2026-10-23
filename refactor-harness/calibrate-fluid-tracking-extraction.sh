#!/usr/bin/env bash
# Mutate temporary source copies to verify the fluid-tracking source proof.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 - "$script_dir" "${1:-}" <<'PY'
import hashlib
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile

harness=pathlib.Path(sys.argv[1])
root=harness.parent
cpp=root/'src/core/SisProdFluidTracking.cpp'
internal=root/'src/include/SisProdFluidTrackingInternal.h'
public=root/'src/include/SisProdFluidTracking.h'
proof=harness/'fluid-tracking-extraction-layout.py'
spec=importlib.util.spec_from_file_location('fluid_layout',proof)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
original=cpp.read_text();header=internal.read_text();contract=public.read_text()
original_hashes={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in [cpp,internal,public]}


def mutate_function(source,name,before,after):
    signature,body=module.parts(source,name)
    changed=module.replace_once(body,before,after)
    # Replace the body using the same brace parser, preserving the remaining definitions.
    import re
    match=re.search(r'(?m)^[\w:<>]+\s+'+re.escape(name)+r'\s*\(',module.masked(source))
    opening=source.index('{',match.end());end=module.closing(source,opening,'{','}')
    return source[:opening+1]+changed+source[end:]


mutations=[
 ('black-oil-water-cut','prepareBlackOilStockTankBalances',
  'balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw;',
  'balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw / cell.oilVolumeFactor;'),
 ('stock-tank-gas-threshold','prepareBlackOilAccessorySource','source.stockTankQuality > 1. - 1e-15','source.stockTankQuality > 1. - 1e-3'),
 ('balance-default','prepareBlackOilStockTankBalances','balance.MultPe = 0.;','balance.MultPe = 1.;'),
 ('source-default','prepareBlackOilCellSourceBasis','source.razdgl = 1.;','source.razdgl = 0.;'),
 ('liquid-flow-donor','upwindLiquidFractions',
  'prepareLiquidFaceTransport(right, state.cells[i].flui, state.cells[i].fluicol, state.cells[i].pres, state.cells[i].temp, state.cells[i].QL, state.cells[i].alf, fracMol1O, ncomp);',
  'prepareLiquidFaceTransport(right, state.cells[i].flui, state.cells[i].fluicol, state.cells[i].pres, state.cells[i].temp, state.cells[i + 1].QL, state.cells[i].alf, fracMol1O, ncomp);'),
 ('complementary-fluid-donor','upwindLiquidFractions',
  'prepareLiquidFaceTransport(right, state.cells[i + 1].flui, state.cells[i].fluicol, state.cells[i + 1].pres, state.cells[i + 1].temp, state.cells[i + 1].QL, state.cells[i + 1].alf, fracMol1O, ncomp);',
  'prepareLiquidFaceTransport(right, state.cells[i + 1].flui, state.cells[i + 1].fluicol, state.cells[i + 1].pres, state.cells[i + 1].temp, state.cells[i + 1].QL, state.cells[i + 1].alf, fracMol1O, ncomp);'),
 ('vapour-zero-flow','upwindVapourAndMixtureFractions','(state.cells[i].MC - state.cells[i].Mliqini) >= 0.','(state.cells[i].MC - state.cells[i].Mliqini) > 0.'),
 ('phase-vector-lifetime','prepareLiquidFaceTransport','vector<double> fracMolFase(ncomp);','vector<double> fracMolFase(ncomp + 1);'),
 ('mixture-weight-order','copyMixtureFaceComposition',
  'face.pesoMol += fluid.masMol[kfrac] * fluid.fracMol[kfrac];','face.pesoMol = fluid.masMol[kfrac] * fluid.fracMol[kfrac] + face.pesoMol;'),
 ('flash-unknown-seed','refreshPhaseSourceComposition','fluid.dCalculatedBeta < 0. || fluid.dCalculatedBeta > 1.','fluid.dCalculatedBeta <= 0. || fluid.dCalculatedBeta >= 1.'),
 ('leak-flash-order','readPhaseSourceFluid',
  'refreshPhaseSourceComposition(state, state.cells[i].acsr.fontechk.fluidoP, state.cells[i].pres, state.cells[i].temp); refreshPhaseSourceComposition(state, state.cells[i].acsr.fontechk.fluidoPamb, state.cells[i].pres, state.cells[i].temp);',
  'refreshPhaseSourceComposition(state, state.cells[i].acsr.fontechk.fluidoPamb, state.cells[i].pres, state.cells[i].temp); refreshPhaseSourceComposition(state, state.cells[i].acsr.fontechk.fluidoP, state.cells[i].pres, state.cells[i].temp);'),
 ('porous-flash-temperature','readPhaseSourceFluid',
  'refreshPhaseSourceComposition(state, state.cells[i].acsr.radialPoro.flup, state.cells[i].pres, tRes);',
  'refreshPhaseSourceComposition(state, state.cells[i].acsr.radialPoro.flup, state.cells[i].pres, state.cells[i].temp);'),
 ('mole-inventory-donor','updatePhaseMoleInventory',
  'fluC[i].MasEspLiq(state.cells[i].pres, state.cells[i].temp)',
  'state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp)'),
 ('phase-water-cut','transportPhaseLiquidProperties',
  'balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw / cell.oilVolumeFactor;',
  'balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw;'),
 ('denag-publication','transportBlackOilBalances','denag[i] = cell.denagini;','denag[i] = source.denag;'),
 ('property-evaluation','blackOilSourceCellFluid',
  '(void)(state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre));','(void)0;'),
 ('outlet-direction','updateBlackOilOutletProperties','state.cells[state.lastCell].Mliqini > -(*state.globals).localtiny','state.cells[state.lastCell].Mliqini >= -(*state.globals).localtiny'),
 ('omp-configuration','transportBlackOilProperties','#pragma omp parallel for num_threads((*state.globals).ntrd)','#pragma omp parallel for num_threads(1)'),
 ('buffer-construction','transportPhaseMolarFractions','Vcr<double> fracMolF(ncomp);','Vcr<double> fracMolF(ncomp + 1);'),
 ('gas-density-callback','publishBlackOilCellProperties','correctGasSpecificGravity(state.cells, state.input, state.steadyIteration, i);','correctGasSpecificGravity(state.cells, state.input, state.steadyIteration, i + 1);'),
]
records=[]
with tempfile.TemporaryDirectory(prefix='marlim-fluid-proof-') as temporary:
    directory=pathlib.Path(temporary)
    c,h,p=[directory/name for name in ['source.cpp','internal.h','public.h']]
    def check(source,private=header,api=contract):
        c.write_text(source);h.write_text(private);p.write_text(api)
        return subprocess.run([sys.executable,str(proof),str(c),str(h),str(p)],capture_output=True,text=True)
    for label,source in [('literal',original),('comments-and-spacing','// static control\n'+original.replace('    ','  '))]:
        result=check(source)
        assert result.returncode==0,(label,result.stderr)
        records.append({'name':label,'positive':True,'passed':True})
    for label,name,before,after in mutations:
        changed=mutate_function(original,name,before,after)
        assert changed!=original,label
        result=check(changed)
        assert result.returncode!=0,('undetected',label)
        assert 'AssertionError' in result.stderr,(label,result.stderr)
        print('caught',label)
        records.append({'name':label,'positive':False,'passed':True})
    for label,source,private,api in [
      ('helper-reference-contract',original.replace('const double &liquidVolumeFlow,','double liquidVolumeFlow,'),header,contract),
      ('private-layout',original,header.replace('    double stockTankQuality;','    float stockTankQuality;'),contract),
      ('public-copy-contract',original,header,contract.replace('ProFlu fluiRev','const ProFlu &fluiRev')),
    ]:
        assert (source,private,api)!=(original,header,contract),label
        result=check(source,private,api)
        assert result.returncode!=0 and 'AssertionError' in result.stderr,('undetected',label)
        print('caught',label)
        records.append({'name':label,'positive':False,'passed':True})
for path,digest in original_hashes.items():
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,path
if sys.argv[2]:
    pathlib.Path(sys.argv[2]).write_text(json.dumps(records,indent=2)+'\n')
print('FLUID TRACKING CALIBRATION PASSED: 23 mutations detected; 2 positive controls passed; production source unchanged')
PY
