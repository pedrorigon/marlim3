// The runner inserts the production search; the march supplies controlled outcomes.
#include "SisProdConstants.h"
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sisprod;
void require(bool condition) { if (!condition) throw std::runtime_error("pressure search contract failed"); }
struct Failure { std::string message; };
[[noreturn]] void NumError(const std::string &message) { throw Failure{message}; }
struct Cell { struct { int tipo=kAccessoryLiquidInjection; struct {double QLiq=1.;} injl; struct {double QGas=2.;} injg; } acsr; };
struct Input { int AP=0; double buscaFC=.1; };
struct Globals { int networkMode=0, networkIteration=0; };
namespace sisprod::steady {
struct SteadyStateState {
    Cell *cells; Input &input; Globals *globals;
    mutable double convergenceMonitor=-1.; mutable int searchOrigin=0;
};
struct SteadyStateSearchState { SteadyStateState march; bool &reverseSteady; };
struct SteadyDependencies {};
namespace detail {
struct PressureMarchStop { int cell=-1; double pressure=0.; };
}
std::vector<double> events;
std::vector<double> residuals;
double fallbackResidual=0.;
int calls=0;
double controlledMarch(const SteadyStateState &state,double guess,detail::PressureMarchStop *stop) {
    events.insert(events.end(),{1.,guess,state.input.buscaFC,double(state.searchOrigin)});
    const double residual=calls<int(residuals.size())?residuals[calls]:fallbackResidual;
    ++calls;
    if(stop) *stop=residual==-kSteadyMarchLimit ? detail::PressureMarchStop{1+calls%3,.75} : detail::PressureMarchStop{};
    return residual;
}
double marchProductionPressureToPressure(const SteadyStateState &state,const SteadyDependencies &,double guess) {
    return controlledMarch(state,guess,nullptr);
}
namespace detail {
double marchProductionPressureToPressure(const SteadyStateState &state,const SteadyDependencies &,double guess,PressureMarchStop *stop) {
    return controlledMarch(state,guess,stop);
}
}
double estimateInletMixtureVelocity(const SteadyStateSearchState &,double guess) { events.insert(events.end(),{2.,guess}); return 1.; }
double solveSteadyRoot(const SteadyStateSearchState &,const SteadyDependencies &,double negative,double positive,int prod,int boundary) {
    events.insert(events.end(),{3.,negative,positive,double(prod),double(boundary)}); return 77.;
}
double searchReverseProductionPressureToPressure(const SteadyStateSearchState &,const SteadyDependencies &,double guess,double maximum,int count) {
    events.insert(events.end(),{4.,guess,maximum,double(count)}); return 99.;
}

// @PRESSURE_SEARCH@

}

int main(int argc,char **argv) {
    using namespace sisprod::steady;
    require(argc==2); auto *output=std::fopen(argv[1],"wb");require(output);
    const auto record=[&](double value){require(std::fwrite(&value,sizeof(value),1,output)==1);};
    int scenarios=0;
    for(int mode=0;mode<13;++mode) for(int network : {0,1,2}) {
        Cell cells[2]; Input input; Globals globals;
        bool reverse=true; SteadyStateSearchState state{{cells,input,&globals},reverse};
        globals.networkMode=network==0?0:1; globals.networkIteration=network==2?1:0;
        residuals.clear();events.clear();calls=0;
        double guess=10.,maximum=0.;int iterations=0,expectedCalls=0;
        switch(mode) {
        case 0: fallbackResidual=.0005;expectedCalls=1;break;
        case 1: residuals={5.,4.,-3.};expectedCalls=3;break;
        case 2: residuals={-5.,-4.,3.};expectedCalls=3;break;
        case 3: residuals={-kSteadyMarchLimit,-kSteadyMarchLimit,5.,5.,4.,-3.};expectedCalls=6;break;
        case 4: fallbackResidual=3.;expectedCalls=202;break;
        case 5: fallbackResidual=-kSteadyMarchLimit;guess=1e-4;expectedCalls=23;break;
        case 6: residuals={3.};fallbackResidual=-kSteadyMarchLimit;expectedCalls=202;break;
        case 7: fallbackResidual=-3.;expectedCalls=202;break;
        case 8: residuals={-3.};fallbackResidual=-kSteadyMarchLimit;expectedCalls=202;break;
        case 9: residuals={std::numeric_limits<double>::infinity()};fallbackResidual=-kSteadyMarchLimit;expectedCalls=202;break;
        case 10: fallbackResidual=std::numeric_limits<double>::quiet_NaN();expectedCalls=1;break;
        case 11: residuals={1.,2.,3.};expectedCalls=3;break;
        case 12: fallbackResidual=2.;maximum=1.;expectedCalls=network==0?202:2;break;
        }
        double result=0.;bool failed=false;
        try { result=searchProductionPressureToPressure(state,{},guess,maximum,iterations); }
        catch(const Failure &failure) {
            failed=true;
            require(network==0);
            require(failure.message.find("atingiu maximo de iteracoes")!=std::string::npos);
#if MARLIM_CORRECTED_BUILD
            require(failure.message.find("pressure march diagnostics:")!=std::string::npos);
            if(mode==5) {
                require(failure.message.find("preliminary=23, bracketing=0, completed=0, pressure-limited=23")!=std::string::npos);
                require(failure.message.find("last-completed-residual")==std::string::npos);
            }
            if(mode==6 || mode==8) {
                require(failure.message.find("preliminary=1, bracketing=201, completed=1, pressure-limited=201")!=std::string::npos);
                require(failure.message.find("stopping-cell=2, stopping-pressure[kgf/cm2]=0.75")!=std::string::npos);
            }
            if(mode==9) {
                require(failure.message.find("completed=0, pressure-limited=201, non-finite=1")!=std::string::npos);
                require(failure.message.find("last-completed-residual")==std::string::npos);
            }
            std::printf("diagnostic mode=%d %s\n",mode,failure.message.c_str());
#endif
        }
        require(calls==expectedCalls);
        if(mode==0) require(result==guess && !failed);
        if(mode==1 || mode==2 || mode==3 || mode==10) require(result==77. && !failed);
        if(mode==11) require(result==99. && !failed);
        if(mode==4 || mode==5 || mode==6 || mode==7 || mode==8 || mode==9 || (mode==12 && network==0)) {
            if(network==0) require(failed);
            else require(!failed && result==(network==1?kSteadyConvergenceFailure:-kSteadyConvergenceFailure));
        }
        if(mode==12 && network!=0) require(result==0. && cells[0].acsr.injl.QLiq==kNoFlowSolutionMarker);
        for(double v:{double(mode),double(network),double(calls),result,double(failed),double(reverse),
                       input.buscaFC,state.march.convergenceMonitor,double(state.march.searchOrigin),
                       cells[0].acsr.injl.QLiq,cells[0].acsr.injg.QGas,double(events.size())}) record(v);
        for(double v:events) record(v);
        ++scenarios;
    }
    std::fclose(output);std::printf("pressure search scenarios=%d\n",scenarios);
}
