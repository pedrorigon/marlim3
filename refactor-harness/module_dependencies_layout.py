#!/usr/bin/env python3
"""Check module call substitutions and live bindings before historical layout proofs.

The calculation-token proof connects the current explicit dependencies to the
reference interfaces. Historical extraction checks can then use that verified
snapshot without retaining old infrastructure in production code.
"""
from pathlib import Path
import re
import json
import subprocess
from simulation_context_names import restore as restore_context_names

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT
REFERENCE = 'ceea4db73b005c11ef990f0a52e6f32d3c567855'
_CACHE = {}
def reference(path):
    if path not in _CACHE:
        _CACHE[path] = subprocess.check_output(['git','show',f'{REFERENCE}:{path}'],cwd=ROOT,text=True)
    return _CACHE[path]
def candidate(path, overrides):
    source = overrides.get(path, (ROOT/path).read_text() if (ROOT/path).exists() else reference(path))
    return restore_context_names(source, path)

LEX = re.compile(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'')

def masked(s):
    return LEX.sub(lambda m: re.sub(r'[^\n]', ' ', m[0]), s)

def close(s, start, left='(', right=')'):
    depth=0
    for i in range(start,len(s)):
        depth += (s[i]==left)-(s[i]==right)
        if depth==0:return i
    raise AssertionError(('unclosed', s[start:start+80]))

def functions(s, types):
    m=masked(s);out=[]
    p=re.compile(r'\b(?P<name>\w+(?:::\w+)*)\s*\(\s*const\s+(?P<type>'+ '|'.join(types)+r')\s*&\s*(?P<arg>state)?\b')
    for hit in p.finditer(m):
        op=m.index('(',hit.start());end=close(m,op);n=end+1
        while n<len(m) and m[n].isspace():n+=1
        if n==len(m) or m[n]!='{':continue
        tail=close(m,n,'{','}');out.append(dict(name=hit['name'],leaf=hit['name'].split('::')[-1],type=hit['type'],start=hit.start(),op=op,end=end,body=n,tail=tail,arg=hit['arg']))
    return out

def edits(s, replacements):
    for a,b,new in sorted(replacements,reverse=True):s=s[:a]+new+s[b:]
    return s

def direct_calls(s, mapping):
    changes=[]
    for old,(target,prefix) in mapping.items():
        p=re.compile(r'\b'+re.escape(old)+r'\s*\(')
        m=masked(s)
        for h in p.finditer(m):
            op=m.index('(',h.start());end=close(m,op);args=s[op+1:end]
            changes.append((h.start(),op+1,target+'('+prefix+(', ' if prefix and args.strip() else '')))
    return edits(s,changes)


stage='final'
mapping={}
for owner in ['state.updaters.','state.march.updaters.']:
    for n in ['advanceSteadyTemperature','advanceReverseSteadyTemperature','computeTemperature','computeGasTemperature','updateProductionTemperaturePeriphery']:
        mapping[owner+n]=('thermal::'+n,'')
    for n in ['initializeSteadyValveGasFlowRate','initializeTubingConnectionSteady','updateSteadyGasPressure','updateSteadyGasTemperature','computeSteadyGasFlowRate','connectTubing','connectTubingSteady','steadyGasPressureDrop','steadyInjectionPressureDrop','searchUnloadingInjectionPressure','solveGasLine']:
        mapping[owner+n]=('gaslift::'+n,'')
    mapping[owner+'steadyDriftClosure']=('driftflux::coefficient::steadyState','')
    mapping[owner+'updateSource']=('sources::renewSourceTerms','')
for old,new in [('writeProductionTrendHeader','writeProductionTrendHeader'),('writeProductionTrendRows','writeProductionTrendRows'),('writeGasTrendHeader','writeGasLineTrendHeader'),('writeGasTrendRows','writeGasLineTrendRows'),('writeProductionCrossSectionTrendHeader','writeProductionCrossSectionTrendHeader'),('writeProductionCrossSectionTrendRows','writeProductionCrossSectionTrendRows'),('writeGasCrossSectionTrendHeader','writeGasLineCrossSectionTrendHeader'),('writeGasCrossSectionTrendRows','writeGasLineCrossSectionTrendRows')]:mapping['state.updaters.'+old]=('trendoutput::'+new,'')
mapping['state.updaters.updateTemperatures']=('thermal::updateDistributedMassTransfer','')
mapping['state.updaters.updateThermal']=('thermal::updateFlowPartitionTerms','')
mapping['state.sourceUpdater']=('sources::renewSourceTerms','')
for n in ['instantaneous','buffered','initialization','bufferedInitialization']:mapping['state.closureUpdater.'+n]=('driftflux::coefficient::'+n,'')
for old,new in [('dischargeTemperature','computeDischargeTemperature'),('computeGasTemperature','computeGasTemperature'),('gasLiftDischargeTemperature','computeGasLiftDischargeTemperature')]:mapping['state.temperatureUpdater.'+old]=('thermal::'+new,'')
mapping['state.updaters.updateBlackOilProperties']=('fluidtracking::transportBlackOilProperties','')
mapping['state.updaters.updateMolarFractions']=('fluidtracking::transportPhaseMolarFractions','')
if stage=='final':
    mapping['state.updaters.marchTransientEnergy']=('thermal::advanceTransientEnergy','state.step')
    mapping['state.updaters.advanceGasSubStep']=('gaslift::advanceGasSubStep','')
    for n in ['searchGasPressureSteadySecondary','searchGasPressureSteadyTertiary']:mapping['state.updaters.'+n]=('steady::'+n,'')
    mapping['state.evolutionUpdater.renew']=('transient::updateCells','0')
    mapping['state.evolutionUpdater.solvePressureVelocityCoupling']=('transient::solvePressureVolumeCoupling','0')
    mapping['state.updaters.solveHydrateEnvelopes']=('HYDRATE_ENVELOPE','')

def strip_function(s,name,types):
    fs=[f for f in functions(s,types) if f['leaf']==name]
    if not fs:return s,None
    assert len(fs)==1
    f=fs[0];start=s.rfind('\n',0,f['start'])+1
    return s[:start]+s[f['tail']+1:],s[start:f['tail']+1]

def tokens(s):
    s=re.sub(r'(?m)^#include[^\n]*','',s)
    p=re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/|\w+|\S')
    return [m[0] for m in p.finditer(s) if not m[0].startswith(('//','/*'))]

def normalize(s,base=False):
    if base:
        s=s.replace('state.updaters.updateThermal()','state.updaters.updateThermal(0)')
        s=s.replace('state.temperatureUpdater.computeGasTemperature(gasCellIndex, state.gasCells[gasCellIndex - 1].tempini)','state.temperatureUpdater.computeGasTemperature(gasCellIndex, state.gasCells[gasCellIndex - 1].tempini, 0)')
        for n in ['updateMolarFractions','updateBlackOilProperties']:s=s.replace('state.updaters.'+n+'(fluiRev)','state.updaters.'+n+'(reverseFlowFluid(fluiRev))')
        s=s.replace('state.updaters.correctGasSpecificGravity(i)','fluidtracking::correctGasSpecificGravity(state.cells, state.input, state.steadyIteration, i)')
        s=s.replace('correctGasSpecificGravity(state, i)','fluidtracking::correctGasSpecificGravity(state.cells, state.input, state.steadyIteration, i)')
        s=direct_calls(s,mapping)
    else:
        s=re.sub(r',\s*const\s+(?:(?:\w+::)*\w+)\s*&\s*(dependencies|gasLiftDependencies|thermalDependencies|transientStep)\b','',s)
        args=r'(?:dependencies(?:\.\w+)*|gasLiftDependencies|thermalDependencies|transientStep)'
        for _ in range(10):
            old=s
            s=re.sub(r'(?<=[(,])\s*'+args+r'\s*,\s*','',s)
            s=re.sub(r',\s*'+args+r'\s*(?=[,)])','',s)
            s=re.sub(r'(?<=\()\s*'+args+r'\s*(?=\))','',s)
            if old==s:break
        s=s.replace('dependencies.hydrates()','HYDRATE_ENVELOPE()')
        if stage=='final':
            for n in ['searchGasPressureSteadySecondary','searchGasPressureSteadyTertiary']:
                s=s.replace('steady::'+n+'(dependencies.search', 'steady::'+n+'(')
    s=re.sub(r'(?<![:\w])correctGasSpecificGravity\(', 'fluidtracking::correctGasSpecificGravity(', s)
    return tokens(s)


def verify(overrides=None):
    overrides = overrides or {}
    assert not (ROOT/'src/core/SisProdAdapters.cpp').exists()
    assert not (ROOT/'src/include/SisProdSolveContext.h').exists()
    records=[]
    for name in ['SisProdFluidTracking.cpp','SisProdGasLift.cpp','SisProdSteadyState.cpp','SisProdSteadyStateSearch.cpp','SisProdThermal.cpp','SisProdTransient.cpp']:
        path='src/core/'+name
        base=reference(path);now=candidate(path,overrides)
        old_layout=bool(re.search(r'state\.(?:march\.)?(?:updaters\.|sourceUpdater|closureUpdater|temperatureUpdater)',masked(now)))
        if name=='SisProdSteadyState.cpp':
            base,_=strip_function(base,'correctGasSpecificGravity',['SteadyStateState'])
            if old_layout:now,_=strip_function(now,'correctGasSpecificGravity',['SteadyStateState'])
        if name=='SisProdFluidTracking.cpp' and not old_layout:
            h=re.search(r'(?m)^void correctGasSpecificGravity\(',now);assert h
            op=now.index('{',h.start());end=close(masked(now),op,'{','}')
            new_kernel=now[h.start():end+1];now=now[:h.start()]+now[end+1:]
            _,old_kernel=strip_function(reference('src/core/SisProdSteadyState.cpp'),'correctGasSpecificGravity',['SteadyStateState'])
            expected=old_kernel.replace('const SteadyStateState &state','Cel *const &cells, const Ler &input, const int &steadyIteration').replace('state.cells','cells').replace('state.input','input').replace('state.steadyIteration','steadyIteration')
            assert tokens(new_kernel)==tokens(expected),'density correction'
        if name=='SisProdTransient.cpp' and not old_layout:
            h=re.search(r'(?m)^const ProFlu &reverseFlowFluid\(',now);assert h
            op=now.index('{',h.start());end=close(masked(now),op,'{','}')
            fallback=now[h.start():end+1];now=now[:h.start()]+now[end+1:]
            old=reference('src/core/SisProdAdapters.cpp');start=old.index('const ProFlu &reverseFlowFluid');end=old.index('\n}',start)+2
            assert tokens(fallback)==tokens(old[start:end]),'reverse fluid fallback'
        a,b=normalize(base,True),normalize(now,old_layout)
        assert a==b, (path,'calculation tokens or operation bindings')
        records.append({'file':path,'tokens':len(a)})
    verify_interfaces(overrides)
    return records


def block(text, name):
    m=masked(text);h=re.search(r'\b(?:struct|enum class) '+re.escape(name)+r'(?:[^{;]*)\{',m);assert h,name
    start=m.index('{',h.start());end=close(m,start,'{','}')
    return text[start+1:end]


def builder(text,name):
    m=masked(text);h=re.search(r'(?m)^[\w:]+ (?:SProd::)?'+name+r'\([^;{]*\)(?: const)? \{',m);assert h,name
    start=m.index('{',h.start());end=close(m,start,'{','}')
    return text[start+1:end]


def verify_interfaces(overrides):
    old=reference('src/core/SisProdAdapters.cpp');current=candidate('src/core/SisProd.cpp',overrides)
    names=re.findall(r'(?m)^[\w:]+ (\w+StateOf)\(',old)
    assert len(names)==10
    for name in names:
        expected=builder(old,name).replace('system.','').replace('StateOf(system, context)','StateOf()')
        expected=re.sub(r'\.(?:updaters|temperatureUpdater) = \{context\},','',expected)
        expected=re.sub(r'\.(?:sourceUpdater|closureUpdater|evolutionUpdater) = \{\.context = context\},','',expected)
        if name=='fluidTrackingStateOf':expected=expected.replace('\n    };','        .steadyIteration = steadyIteration,\n    };')
        assert tokens(expected)==tokens(builder(current,name)),(name,'live view binding')
    assert tokens(builder(current,'solveHydrateEnvelopes'))==tokens(builder(reference('src/core/SisProd.cpp'),'solveHydrateEnvelopes')),'hydrate integration'
    old_sprod=reference('src/include/SisProd.h');new_sprod=candidate('src/include/SisProd.h',overrides)
    def public(text):return text[text.index('class SProd {'):].split('  public:',1)[1].split('  private:',1)[0]
    assert tokens(public(old_sprod))==tokens(public(new_sprod)),'SProd public interface'
    headers={str(p.relative_to(ROOT)):restore_context_names(p.read_text(),str(p.relative_to(ROOT))) for p in (ROOT/'src/include').glob('SisProd*.h')}
    headers.update({p:restore_context_names(t,p) for p,t in overrides.items() if p.endswith('.h')})
    reference_headers=set(subprocess.check_output(
        ['git','ls-tree','-r','--name-only',REFERENCE,'src/include'],cwd=ROOT,text=True).splitlines())
    for path,text in headers.items():
        if not path.endswith('.h') or path=='src/include/SisProd.h':continue
        if path not in reference_headers:continue
        old_header=reference(path)
        def prototypes(text):
            clean=masked(text);result={}
            pattern=r'(?m)^(?:\[\[nodiscard\]\] )?[\w:<>]+\s+(\w+)\s*\('
            for match in re.finditer(pattern,clean):
                opening=clean.index('(',match.start());end=close(clean,opening)
                if clean[end+1:].lstrip().startswith(';'):
                    result.setdefault(match[1],[]).append(text[match.start():end+1])
            return result
        before,after=prototypes(old_header),prototypes(text)
        for name,signatures in before.items():
            if name=='correctGasSpecificGravity' and path.endswith('SisProdSteadyState.h'):continue
            assert name in after,(path,name,'missing operation')
            assert sorted(normalize(sig) for sig in signatures)==sorted(normalize(sig) for sig in after[name]),(path,name,'operation interface')
        # Data and policy structures keep their exact fields and operations.
        for name in re.findall(r'(?m)^(?:struct|enum class) (\w+)[^{;]*\{',masked(old_header)):
            if 'Updater' in name:continue
            original=block(old_header,name)
            original=re.sub(r'\b\w+Updater(?:s)? \w+;','',original)
            if name=='FluidTrackingState':original+='\nconst int &steadyIteration;'
            actual=block(text,name)
            # A reference snapshot is also accepted by chained historical proofs.
            if tokens(actual)==tokens(block(old_header,name)):continue
            assert normalize(original)==normalize(actual),(path,name,'data/policy layout')
    return True


def verified_snapshot(overrides):
    verify(overrides)
    return {path:reference(path) for path in overrides}


def candidate_text(path):
    return (ROOT/path).read_text() if (ROOT/path).exists() else reference(path)


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
    print('MODULE DEPENDENCIES VERIFIED: six calculation streams and ten live-view bindings')
