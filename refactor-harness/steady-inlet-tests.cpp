// The runner inserts production bodies; fluid callbacks record their evaluation order.
#include "SisProdConstants.h"
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace sisprod;
std::vector<double> events;
void require(bool condition) { if (!condition) throw std::runtime_error("steady inlet contract failed"); }
bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
void record(std::FILE *file, double value) { require(std::fwrite(&value, sizeof(value), 1, file) == 1); }
void event(double id, double p = 0., double t = 0.) { events.insert(events.end(), {id,p,t}); }
struct Fluid {
    double id = 1., RGO = 700., BSW = .2, Deng = .7, yco2 = .03, API = 40., Denag = 1.;
    double dStockTankVaporMassFraction = .2, TempL = 20., TempH = 80., LVisL = 1., LVisH = 2.;
    double MasEspLiq(double p, double t) const { event(10.*id+1.,p,t); return 900.; }
    double MasEspGas(double p, double t) const { event(10.*id+2.,p,t); return 10.; }
    double BOFunc(double p, double t) const { event(10.*id+3.,p,t); return 1.2; }
    double BAFunc(double p, double t) const { event(10.*id+4.,p,t); return 1.01; }
    double VisOM(double t) const { event(10.*id+5.,t); return 5.; }
    void RenovaFluido() { event(10.*id+6.,RGO); }
};
struct ComplementaryFluid {
    double id = 1., TR = 5.;
    double MasEspFlu(double p, double t) const { event(100.+id,p,t); return 1050.+id; }
    double standardDensity() const { event(200.+id); return 1000.+id; }
};
struct LiquidSource { double QLiq = .01; Fluid FluidoPro{3.,1600.}; ComplementaryFluid fluidocol{3.}; };
struct Cell {
    double FW = .2, MC = 0., MComp = 0., Mliqini = 0., QL = .002;
    double VolLeveST = -1., alf = .3, bet = .1, fontemassCR = 0., fontemassLR = .001;
    double pres = 80., presaux = 79., temp = 80., dx = 5., dxL = 5.;
    Fluid flui;
    ComplementaryFluid fluicol;
    struct { LiquidSource injl; } acsr;
};
struct Input { int trackRGO = 0, flashCompleto = 0; };
struct Globals { double smallValueTolerance = 1e-5, maximumGasOilRatio = 9400000.; int networkMode = 0; };
struct SteadyStateState {
    Cell *cells; Input &input; Globals *globals;
    int productionFluidCount = 1, steadyIteration = 0;
};
namespace fluidtracking {
void correctGasSpecificGravity(Cell *, const Input &, int iteration, int i) { event(300.,iteration,i); }
}

// @LIQUID_SOURCE@

double totalStandardLiquidTrend(Cell *celula, int m) {
    const int i = m, trend = 0, linha = 0;
    int k = 0;
    struct { int qlstTot = 1; } trendp[1];
    double output = -1.; double *flut[] = {&output};
    // @TOTAL_LIQUID_TREND@
    require(k == 1);
    return output;
}

int main(int argc, char **argv) {
    require(argc == 3);
    auto *unchanged = std::fopen(argv[1],"wb"); auto *corrected = std::fopen(argv[2],"wb");
    require(unchanged && corrected);
    int sourceCases = 0, correctedSourceCases = 0, trendCases = 0, correctedTrendCases = 0;
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int i : {1,2}) for (int upstream : {0,1,2,3,4}) for (int track : {-1,0,1})
        for (double sourceFlow : {-1.,0.,1e-16,.01}) for (int phase : {0,1})
            for (double rate : {-1e-16,0.,1e-16,7.2585186596708266e-8,1e-5,std::nextafter(1e-5,infinity),.01,infinity,nan})
                for (double sourceGOR : {0.,1600.,-1.,infinity,nan}) {
        std::array<Cell,3> cells;
        for (int n=0;n<3;++n) { cells[n].flui.id=n+1.; cells[n].fluicol.id=n+1.; }
        Input input{track,0}; Globals globals;
        SteadyStateState state{cells.data(),input,&globals};
        auto &left = cells[i-1];
        left.MC = upstream == 1 ? 1. : 0.;
        left.Mliqini = upstream == 2 ? .001 : 0.;
        left.MComp = upstream == 3 ? .001 : 0.;
        left.MC = upstream == 4 ? nan : left.MC;
        left.acsr.injl.QLiq = sourceFlow;
        left.acsr.injl.FluidoPro.RGO = sourceGOR;
        left.acsr.injl.FluidoPro.dStockTankVaporMassFraction = phase == 0 ? .2 : 1.;
        const double sourceWaterFraction = .2*1.01/(1.2+1.01*.2-.2*1.2);
        left.fontemassLR = rate*1.2*900./(1.-sourceWaterFraction);
        const double actualSourceOil = phase == 0 ? left.fontemassLR*(1.-sourceWaterFraction)/(1.2*900.) : 0.;
        const double upstreamOil = left.Mliqini*(1.-left.FW)*(1.-left.bet)/(1.2*900.);
        const bool changed = i==1 && upstream==0 && track!=-1 && sourceFlow>1e-15 &&
            actualSourceOil>0. && std::fabs(actualSourceOil)<=globals.smallValueTolerance &&
            sourceGOR>=0. && std::isfinite(sourceGOR);
        double expected = 700.;
        if (std::fabs(sourceFlow)>1e-15) {
            if (track==-1) {
                const double light = (.3*10./(.7*kAirDensityAtStandardConditions)) + .7*.9*.8*100./1.2;
                expected = light*1.2/(.7*.9*.8);
                if (expected>globals.maximumGasOilRatio) expected=globals.maximumGasOilRatio;
            } else if (std::fabs(upstreamOil+actualSourceOil)>globals.smallValueTolerance) {
                expected=(upstreamOil*700.+actualSourceOil*sourceGOR)/(upstreamOil+actualSourceOil);
            } else if (changed && MARLIM_CORRECTED_BUILD) expected=sourceGOR;
            else expected=globals.maximumGasOilRatio;
        }
        events.clear(); double residence=-1.,bo=-2.,ba=-3.,fw=-4.;
        applySteadyMassLiquidInjection(state,i,20.,80.,1.2,1.01,100.,79.,residence,bo,ba,fw);
        require((std::isnan(expected) && std::isnan(cells[i].flui.RGO)) ||
                std::fabs(cells[i].flui.RGO-expected)<=std::fabs(expected)*1e-14 || same(expected,cells[i].flui.RGO));
        if (changed && MARLIM_CORRECTED_BUILD) require(same(cells[i].flui.RGO,sourceGOR));
        auto *out=changed ? corrected : unchanged;
        for (double v : {double(i),double(upstream),double(track),sourceFlow,double(phase),rate,sourceGOR,
                         cells[i].flui.RGO,cells[i].flui.Deng,cells[i].flui.yco2,cells[i].flui.API,
                         cells[i].flui.BSW,cells[i].flui.Denag,cells[i].bet,left.VolLeveST,residence,bo,ba,fw}) record(out,v);
        record(out,events.size()); for (double v:events) record(out,v);
        ++sourceCases; correctedSourceCases+=changed;
    }
    for (int m : {0,1,2}) for (double flow : {-.002,0.,.002})
        for (double beta : {0.,1e-16,.1,.7,1.}) for (int vapour : {0,1}) {
        std::array<Cell,3> cells;
        for (int n=0;n<3;++n) {
            cells[n].QL=flow; cells[n].bet=beta; cells[n].flui.BSW=.2;
            cells[n].flui.dStockTankVaporMassFraction=vapour==0?.2:1.;
            cells[n].fluicol.id=n+1.; cells[n].pres=80.-n; cells[n].temp=70.-n;
        }
        const bool changed=m==0 && beta>1e-15;
        if (changed && !MARLIM_CORRECTED_BUILD) continue;
        const double water=.2*1.01/(1.2+1.01*.2-.2*1.2);
        const double hydrocarbons=vapour==0?flow*(1.-water)*(1.-beta)*86400./1.2/.8:0.;
        const double id=(m==0?0:m-1)+1.;
        const double complement=beta>1e-15?86400.*flow*beta*(1050.+id)/(1000.+id):0.;
        events.clear(); const double actual=totalStandardLiquidTrend(cells.data(),m);
        const double expected=hydrocarbons+complement;
        require(same(actual,expected) || std::fabs(actual-expected)<=std::fabs(expected)*1e-14);
        auto *out=changed ? corrected : unchanged;
        record(out,actual);record(out,events.size());for(double v:events)record(out,v);
        ++trendCases;correctedTrendCases+=changed;
    }
    std::fclose(unchanged);std::fclose(corrected);
    std::printf("source cases=%d corrected=%d; trend cases=%d corrected=%d\n",sourceCases,correctedSourceCases,trendCases,correctedTrendCases);
}
