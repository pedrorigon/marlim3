#ifndef SISPRODFLUIDTRACKINGINTERNAL_H_
#define SISPRODFLUIDTRACKINGINTERNAL_H_

#include "SisProdFluidTracking.h"

#include <vector>

template <typename T>
class Vcr;

namespace sisprod::fluidtracking {
namespace {

/// What upwinding gives at one face of cell i: the liquid's black-oil properties
/// and the gas's density and composition, each taken from the side its phase
/// flows from.
struct BlackOilFace {
    double betI;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double solutionGasRatio;
    double dgO;
    double yco2O;
    double API;
    double BSW;
    double denag;
    double viscL;
    double viscH;
    double razdgd;
    double razdgl;
    double rhogST;
    double dgG;
    double yco2G;
};

/// Cell i itself: its liquid holdup, complementary-liquid fraction and in-situ black-oil
/// properties, and the fluid properties it holds before this step (the *ini
/// members).
struct BlackOilCell {
    double liquidHoldup;
    double complementaryFraction;
    double rholST;
    double rhog;
    double rhogST;
    double solutionGasRatio;
    double razdgd;
    double razdgl;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double dgini;
    double yco2ini;
    double APIini;
    double BSWini;
    double denagini;
    double viscLini;
    double viscHini;
};

/// The source of cell i, by accessory: the mass rates it brings and the
/// properties of the fluid it brings them in, which hold cell i's own values
/// when there is no source.
struct BlackOilSource {
    double dissolvedGas;
    double freeGas;
    double deadOil;
    double water;
    double API;
    double BSW;
    double denag;
    double dgO;
    double dgG;
    double yco2O;
    double yco2G;
    double rgo;
    double viscL;
    double viscH;
    double rholST;
    double rhogST;
    double razdgd;
    double razdgl;
    double stockTankQuality;
};

/// The balances of cell i: the face fluxes (Mult*), the volumes at the end of
/// the step (vol*Fim) and the residuals (residuo*).
struct BlackOilBalance {
    double MultOe;
    double MultOd;
    double MultGe;
    double MultGd;
    double volleveFim;
    double residuo;
    double volpesFim;
    double volaguaFim;
    double MultPe;
    double MultPd;
    double residuoP;
    double MultAe;
    double MultAd;
    double residuoA;
};

/// What upwinding gives at one face of cell i in the phase molar-fraction
/// transport: the liquid's properties, the phase mass flow rates and quality, and
/// the molar weights of the mixture and of each phase.
struct PhaseFace {
    double betI;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double BSW;
    double denag;
    double viscL;
    double viscH;
    double solutionGasRatio;
    double vazMasLiq;
    double vazMasGas;
    double titV;
    double pesoMol;
    double pesoMolO;
    double pesoMolG;
};

/// Cell i itself: its liquid holdup, complementary-liquid fraction and in-situ properties,
/// the properties it holds before this step (the *ini members), and its molar
/// weight.
struct PhaseCell {
    double liquidHoldup;
    double complementaryFraction;
    double solutionGasRatio;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double BSWini;
    double denagini;
    double viscLini;
    double viscHini;
    double pesoMol;
};

/// Thermodynamic working values for one mixture-balance or inventory evaluation.
struct PhaseTransported {
    double betIV;
    double solutionGasRatio;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bswV;
    double rhoOVol;
    double rhoWVol;
    double titVol;
    double tempMol;
};

/// The source of cell i, by accessory: the mass rates it brings and the
/// properties of the fluid it brings them in.
struct PhaseSource {
    double dissolvedGas;
    double freeGas;
    double deadOil;
    double water;
    double BSW;
    double denag;
    double viscL;
    double viscH;
    double rhoW;
    double rhoO;
    double waterCut;
    double waterVolumeFactor;
    double oilVolumeFactor;
    double tit;
    double pesoMol;
};

/// The balances of cell i: the face fluxes (Mult*), the volumes at the end of
/// the step (vol*Fim) and the residuals (residuo*).
struct PhaseBalance {
    double volpesFim;
    double volaguaFim;
    double MultPe;
    double MultPd;
    double residuoP;
    double MultAe;
    double MultAd;
    double residuoA;
};

/// Updates black-oil properties from the ordered stock-tank balances.
/// Gas, oil, and water transport retain their individual model and inflow gates.
/// The published water density retains the cell's previous value.
void transportBlackOilBalances(const FluidTrackingState &state, const BlackOilFace &left, const BlackOilFace &right, const BlackOilCell &cell, const BlackOilSource &source, BlackOilBalance &balance, double dx, double flowArea, int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double dt);

/// Prepares a source carried by a local or upstream choke or volumetric pump.
/// The inlet and upstream paths retain their individual property expressions.
void blackOilSourceCellFluid(const FluidTrackingState &state, BlackOilSource &source, int i, double temperatureHigh, double temperatureLow);

/// Prepares the 2D porous source using its transfer BSW and the selected
/// radialPoro.flup or cell fluid; selection follows the signed source mass rate.
void blackOilSourcePorous2D(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Prepares the radial porous source from its reservoir or cell fluid and BSW.
void blackOilSourceRadialPorous(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Prepares an annulus leak from the higher-pressure side using a local working
/// fluid. The ambient-gas rule controls the dissolved-gas source.
void blackOilSourceLeak(const FluidTrackingState &state, BlackOilSource &source, int i, double temperatureHigh, double temperatureLow);

/// Prepares an IPR source using the reservoir fluid only above cell pressure.
void blackOilSourceInflowPerformance(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Prepares a multiple source from the accessory fluid when its mass sum is positive.
void blackOilSourceMultipleSource(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Prepares a liquid injection from the accessory fluid only for positive injection.
void blackOilSourceLiquidInjection(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Prepares a wet gas injection with its selected fluid. The configured injection
/// fluid supplies the water-presence test and source water density.
void blackOilSourceWetGasInjection(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Selects left-face standard gas density and gas composition from the gas-flow donor.
/// Zero gas flow selects the cell fluid.
void upwindLeftFaceGasProperties(const FluidTrackingState &state, BlackOilFace &left, int i);

/// Replaces the right-face black-oil liquid properties for strictly negative liquid flow.
void upwindRightFaceLiquidProperties(const FluidTrackingState &state, BlackOilFace &right, int i, double temperatureHigh, double temperatureLow);

/// Prepares left-face liquid properties with independent inlet/pig fraction rules.
/// Solution gas ratio is converted after evaluating the formation volume factors.
void upwindLeftFaceBlackOilLiquid(const FluidTrackingState &state, BlackOilFace &left, int i, double temperatureHigh, double temperatureLow);

/// Prepares face and source properties, transports black-oil properties,
/// and publishes this cell's three stock-tank volumes.
void transportCellBlackOilProperties(const FluidTrackingState &state, int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double temperatureHigh, double temperatureLow, double dt);

/// Copies a phase composition and normalizes it only when its sequential sum is positive.
/// The minimum search does not shift the composition; the subtraction uses zero.
/// npseudo is positive and the destination already contains npseudo elements.
void normalizePhaseFractions(std::vector<double>& fracMolFase, double* fracMolOriginal, int npseudo);

/// Updates liquid stock-tank volumes, BSW, and dead-oil viscosities under the
/// compositional transport gates. Published water density retains its previous value.
void transportPhaseVolumesAndViscosities(const FluidTrackingState &state, const PhaseFace &left, const PhaseFace &right, const PhaseCell &cell, const PhaseSource &source, PhaseBalance &balance, double dx, double flowArea, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double dt);

/// Converts source oil/water mass rates to stock-tank volume rates and reads
/// accessory-specific liquid properties. Configured and selected fluid donors differ
/// for some accessories and retain their separate roles.
void phaseSourceStandardRatesByAccessory(const FluidTrackingState &state, PhaseSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow);

/// Advances mixture molar fractions through limCorrige ordered corrections,
/// shifts negative trial fractions, normalizes, and refreshes enabled flash properties.
/// The molecular-weight update uses the original cell composition.
void solveCellPhaseMolarFractions(const FluidTrackingState &state, const PhaseFace &left, const PhaseFace &right, PhaseCell &cell, PhaseTransported &transported, const PhaseSource &source, int limCorrige, ProFlu &fluF, int i, Vcr<double> &fracMol0, Vcr<double> &fracMol1, double dt, ProFlu *fluC, int ncomp);

/// Refreshes accessory flash state before selecting the working source fluid.
/// Both leak fluids are refreshed in pressure-side order; porous flashes use reservoir
/// temperature, while the selected liquid properties use cell temperature.
void readPhaseSourceFluid(const FluidTrackingState &state, PhaseSource &source, ProFlu &fluF, int i);

/// Prepares normalized vapour compositions and mixture compositions at both faces.
/// The mixture compositions and weights are consumed by the mixture balance;
/// phase-specific composition buffers remain separate working storage.
void upwindVapourAndMixtureFractions(const FluidTrackingState &state, PhaseFace &left, PhaseFace &right, int i, Vcr<double> &fracMol0, Vcr<double> &fracMol1, Vcr<double> &fracMol0G, Vcr<double> &fracMol1G, int ncomp);

/// Adjusts transported liquid mass for complementary fluid and water and prepares
/// normalized liquid compositions. Flow, complementary-fluid, and composition donors
/// retain their individual face and cell indices.
void upwindLiquidFractions(const FluidTrackingState &state, PhaseFace &left, PhaseFace &right, int i, Vcr<double> &fracMol0O, Vcr<double> &fracMol1O, int ncomp);

/// Prepares left-face liquid volume factors, water cut, and liquid properties
/// from the liquid-flow donor, preserving the independent inlet/pig fraction rule.
void upwindLeftFacePhaseLiquidProperties(const FluidTrackingState &state, PhaseFace &left, int i, double temperatureHigh, double temperatureLow);

/// Prepares phase and mixture donors, updates source flash state, advances trial
/// mixture composition and mole inventory, then transports the liquid volume properties.
void transportCellPhaseMolarFractions(const FluidTrackingState &state, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, Vcr<double> &fracMol0, Vcr<double> &fracMol1, Vcr<double> &fracMol0O, Vcr<double> &fracMol1O, Vcr<double> &fracMol0G, Vcr<double> &fracMol1G, double temperatureHigh, double temperatureLow, double dt, ProFlu *fluC, int ncomp);

}  // namespace
}  // namespace sisprod::fluidtracking

#endif  // SISPRODFLUIDTRACKINGINTERNAL_H_
