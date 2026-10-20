#include "SisProdGasLift.h"

#include "Leitura.h"
#include "SisProdConstants.h"
#include "Matriz.h"
#include "Vetor.h"
#include "celulaGas.h"
#include "celula3.h"
#include "chokegas.h"
#include "UnitConversions.h"

#include <cmath>
#include <numbers>

namespace sisprod::gaslift {

namespace {

/// Below this opening the injection choke is throttling and becomes the mass
/// source of the first cell instead of a free boundary.
inline constexpr double kThrottlingChokeOpening = 0.2;

/// Fixed gas viscosity used to prime the steady annulus connection, in cP.
inline constexpr double kInitialAnnulusGasViscosityCentipoise = 0.16;

/// Gas-lift valve geometry. The bore threshold and the two spring rates are
/// imperial, as the whole valve correlation is.
inline constexpr double kLargeBoreValveInches = 1.1;
inline constexpr double kLargeBoreSpringRate = 500.0;
inline constexpr double kSmallBoreSpringRate = 1950.0;

/// The unloading controller corrects the injection or the separator pressure when the completion-fluid flow
/// through the selected valve leaves the band from 1 - kUnloadingLowerBandMargin to 1 - kUnloadingUpperBandMargin
/// of the controlled flow rate.
inline constexpr double kUnloadingUpperBandMargin = 0.1;
inline constexpr double kUnloadingLowerBandMargin = 0.4;

/// The largest pressure correction of one step, as a fraction of the corrected pressure per second of the step.
inline constexpr double kMaximumUnloadingCorrectionPerSecond = 0.01;

/// The separator pressure counts as at the unloading's maximum or minimum within these factors.
inline constexpr double kReachedMaximumPressureFactor = 0.9999999;
inline constexpr double kReachedMinimumPressureFactor = 1.0000001;

/// presMaxDesc before the minimum over the IPR accessories: above any pressure the column reaches, in kgf/cm^2.
inline constexpr double kUnboundedUnloadingPressure = 100000.;

/// The pressure the unloading hydrostatics start the gas line from when nothing else gives it, in kgf/cm^2.
inline constexpr double kDefaultUnloadingStartPressure = 10.;

/// The weight of the step's own flow against the moving average, in the flow the unloading controller corrects.
inline constexpr double kCurrentFlowWeight = 0.5;

/// The fraction of the separator pressure the unloading takes off per second while the column's last cell but
/// one holds no completion fluid.
inline constexpr double kSeparatorPressureDecayPerSecond = 0.01;

/// The factors a pressure is relaxed by per step towards its target, after the unloading.
inline constexpr double kRelaxationDownFactor = 0.95;
inline constexpr double kRelaxationUpFactor = 1.05;

/// The ratio of its cell the interface must reach in the last cell but one before it is handed over to the end of
/// the gas line.
inline constexpr double kInterfaceHandoverRatio = 0.99;

/// The temperature of the gas line's last cell while the flow there runs backwards, in degrees Celsius.
inline constexpr double kReverseFlowOutletTemperatureCelsius = 20.;

/// The gas line's master2 valve throttles while its throat area lies between these fractions of the smaller pipe
/// area around it.
inline constexpr double kMaster2ThrottlingAreaFraction = 0.01;
inline constexpr double kMaster2ClosedAreaFraction = 1e-5;

/// The bellows of a calibrated valve: the temperature its dome charge is referred to, in degrees Fahrenheit, and
/// the fractional change of the dome pressure per degree Fahrenheit.
inline constexpr double kBellowsReferenceTemperatureFahrenheit = 80.;
inline constexpr double kBellowsPressureChangePerFahrenheit = 0.00215;

/// The steady gas temperature marches a cell only above this superficial gas velocity, in m/s, and between
/// diameters within this ratio of each other; elsewhere the cell takes the external temperature.
inline constexpr double kMinimumMarchGasVelocity = 1e-3;
inline constexpr double kMinimumMarchDiameterRatio = 0.5;

/// The steady gas temperature marches a cell in substeps when its length over its global thermal resistance
/// passes this value.
inline constexpr double kThermalSubstepLength = 1000.;

/// The diameter a Reynolds number is formed on: the bore when the segment is
/// bare, the hydraulic diameter when it is cased.
[[nodiscard]] double characteristicDiameter(int cased, double bore, double area,
                                            double perimeter) {
    return cased == 0 ? bore : 4 * area / perimeter;
}

/// Opening of a calibrated gas-lift valve, from the state already primed on the
/// choke. Imperial in, fraction out.
[[nodiscard]] double calibratedValveOpening(const ChokeGas &choke) {
    return calibratedValveArea(choke.pcalib * units::kPsiPerKgfPerCm2, choke.tcalib,
                               (choke.presEstag - units::kAtmosphereInKgfPerCm2) * units::kPsiPerKgfPerCm2,
                               (choke.presGarg - units::kAtmosphereInKgfPerCm2) * units::kPsiPerKgfPerCm2,
                               choke.dextern, choke.areagarg,
                               choke.areagarg / choke.areafole,
                               units::celsiusToFahrenheit(choke.tempEstag));
}

/// What flows through a gas-lift valve, and how its throat pressure is
/// recovered.
struct GasThroughValve {
    [[nodiscard]] static double recovery(const ChokeGas &choke) { return choke.frec; }
    [[nodiscard]] static double massFlowRate(ChokeGas &choke, const Ler &) { return choke.massica(); }
};
struct CompletionFluidThroughValve {
    [[nodiscard]] static double recovery(const ChokeGas &choke) { return choke.frecliq; }
    [[nodiscard]] static double massFlowRate(ChokeGas &choke, const Ler &input) { return choke.massica(1, input.salinDescarga); }
};
/// Gas-side recovery, completion-fluid mass flow: the combination the unloading
/// search uses.
struct UnloadingThroughValve {
    [[nodiscard]] static double recovery(const ChokeGas &choke) { return choke.frec; }
    [[nodiscard]] static double massFlowRate(ChokeGas &choke, const Ler &input) { return choke.massica(1, input.salinDescarga); }
};

/// Primes the choke from the two cells it spans and returns the mass source it
/// delivers, already scaled by the calibrated opening when the valve is one.
template <typename Phase>
[[nodiscard]] double primeValveMassSource(const GasLiftState &state, int valveIndex,
                                          const CelG &gasCell, const Cel &productionCell) {
    ChokeGas &choke = state.gasLiftChokes[valveIndex];
    const double recovery = Phase::recovery(choke);
    choke.presEstag = gasCell.pres;
    choke.presGarg = (productionCell.pres - choke.presEstag * recovery) / (1. - recovery);
    choke.tempEstag = gasCell.temp;
    double massSource = productionCell.pres < gasCell.pres ? Phase::massFlowRate(choke, state.input) : 0.;
    if (choke.tipo == kGasLiftCalibratedValve)
        massSource *= calibratedValveOpening(choke);
    return massSource;
}

/// The opening the injection choke presents to the gas line. Under a pressure
/// condition it is the choke's own, and a choke throttling the line also sets
/// the inlet mass source through it; under a flow-rate condition the source is
/// given and the choke counts as open.
[[nodiscard]] double injectionChokeOpening(const GasLiftState &state, InjectionPressureCondition) {
    double chokeOpeningFraction = state.injectionChoke.areagarg / state.gasCells[0].duto.area;
    if (chokeOpeningFraction < kThrottlingChokeOpening) {
        state.injectionChoke.presEstag = state.injectionPressure;
        state.injectionChoke.tempEstag = state.injectionTemperature;
        state.injectionChoke.presGarg = state.gasCells[0].pres;
        state.gasCells[0].massfonteCH = state.injectionChoke.massica();
    }
    return chokeOpeningFraction;
}

[[nodiscard]] double injectionChokeOpening(const GasLiftState &, InjectionFlowRateCondition) {
    return 1.;
}

/// Opens the injection choke if it is throttling, assembles the band system
/// from every gas cell and solves it.
void assembleAndSolveGasSystem(const GasLiftState &state) {
    double chokeOpeningFraction = withGasInletCondition(
        state.gasCells[0].tipoCC,
        [&](auto condition) { return injectionChokeOpening(state, condition); });
#pragma omp parallel for num_threads(state.globals->ntrd)
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        state.gasCells[gasCellIndex].GeraLocal(state.gasCellCount, state.injectionPressure, state.injectionTemperature, chokeOpeningFraction);
        for (int bandIndex = 0; bandIndex < 9; bandIndex++) {
            state.gasSystemMatrix[3 * gasCellIndex][bandIndex - 3] = state.gasCells[gasCellIndex].local[0][bandIndex];
            state.gasSystemMatrix[3 * gasCellIndex + 1][bandIndex - 4] = state.gasCells[gasCellIndex].local[1][bandIndex];
            state.gasSystemMatrix[3 * gasCellIndex + 2][bandIndex - 5] = state.gasCells[gasCellIndex].local[2][bandIndex];
        }
        state.gasRightHandSide[3 * gasCellIndex] = state.gasCells[gasCellIndex].TL[0];
        state.gasRightHandSide[3 * gasCellIndex + 1] = state.gasCells[gasCellIndex].TL[1];
        state.gasRightHandSide[3 * gasCellIndex + 2] = state.gasCells[gasCellIndex].TL[2];
    }
    state.gasSystemMatrix.GaussElimPP(state.gasRightHandSide);
}

/// Clamps a gas cell's temperature to the physical range and copies it to the
/// neighbours that face it.
void clampAndPropagateTemperature(const GasLiftState &state, int cellIndex) {
    if (state.gasCells[cellIndex].temp < kMinimumTemperatureCelsius)
        state.gasCells[cellIndex].temp = kMinimumTemperatureCelsius;
    if (state.gasCells[cellIndex].temp > kMaximumTemperatureCelsius)
        state.gasCells[cellIndex].temp = kMaximumTemperatureCelsius;

    if (cellIndex > 0)
        state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
    if (cellIndex < state.gasCellCount)
        state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
}

/// Walks a pressure towards a target by a fixed factor per step, without
/// overshooting it. Down at 5% a step, up at 5% a step -- the two rates are not
/// each other's inverse.
void relaxTowards(double &pressure, double target) {
    if (pressure > target) {
        pressure *= kRelaxationDownFactor;
        if (pressure < target)
            pressure = target;
    } else if (pressure < target) {
        pressure *= kRelaxationUpFactor;
        if (pressure > target)
            pressure = target;
    }
}

/// A pressure correction of the unloading control, limited to kMaximumUnloadingCorrectionPerSecond of the pressure
/// it corrects per second of the step, with its sign.
[[nodiscard]] double limitedCorrection(double correction, double pressure, double timeStep) {
    if (fabs(correction) > kMaximumUnloadingCorrectionPerSecond * pressure * timeStep)
        return (fabs(correction) / correction) * kMaximumUnloadingCorrectionPerSecond * pressure * timeStep;
    return correction;
}

/// The volume flow of a gas cell's source during the unloading: gas once the
/// interface is past the cell's middle, completion fluid before.
[[nodiscard]] double sourceVolumeFlowRate(const CelG &cell, double liquidDensity, double gasDensity) {
    return cell.razInter > 0.5 ? cell.massfonteCH / gasDensity : cell.massfonteCH / liquidDensity;
}

/// Publishes the transport properties of the gas line's fluid onto the production
/// cell that faces it, as that cell's external fluid.
struct ExternalFluidProperties {
    double conductivity;
    double specificHeat;
    double density;
    double viscosityCentipoise;
};

void publishExternalFluid(Cel &productionCell, const ExternalFluidProperties &fluid) {
    productionCell.calor.kextern1 = fluid.conductivity;
    productionCell.calor.cpextern1 = fluid.specificHeat;
    productionCell.calor.rhoextern1 = fluid.density;
    productionCell.calor.viscextern1 = fluid.viscosityCentipoise * units::kPascalSecondPerCentipoise;
}

/// The pressure the unloading hydrostatics start the gas line from when the
/// injection pressure does not give it: the unloading's own initial pressure when
/// the unloading is controlled, kDefaultUnloadingStartPressure otherwise.
[[nodiscard]] double defaultUnloadingStartPressure(const GasLiftState &state) {
    return state.input.controDesc == 1 ? state.input.presIniDescG : kDefaultUnloadingStartPressure;
}

/// The unloading starts from the configured injection pressure, unless the
/// unloading is controlled.
[[nodiscard]] double unloadingStartPressure(const GasLiftState &state, InjectionPressureCondition) {
    if (state.input.controDesc == 0)
        return state.input.gasinj.presinj[0];
    return defaultUnloadingStartPressure(state);
}

[[nodiscard]] double unloadingStartPressure(const GasLiftState &state, InjectionFlowRateCondition) {
    return defaultUnloadingStartPressure(state);
}

}  // namespace

void computeGasUnloadingHydrostatics(const GasLiftState &state) {
    // The configured condition, not the first gas cell's: the cells are being
    // given their initial state here.
    double meanPressure = withGasInletCondition(
        state.input.gasinj.tipoCC,
        [&](auto condition) { return unloadingStartPressure(state, condition); });
    if (state.input.controDesc == 1)
        state.injectionPressure = state.input.presIniDescG;

    state.gasCells[0].presL = meanPressure;
    state.gasCells[0].pres = meanPressure;
    state.gasCells[0].presini = meanPressure;
    state.gasCells[1].presL = meanPressure;
    double meanTemperature = state.gasCells[0].calor.Textern1;
    state.gasCells[0].tempL = meanTemperature;
    state.gasCells[0].temp = meanTemperature;
    state.gasCells[1].tempL = meanTemperature;
    double leftDensity = state.gasCells[0].flui.MasEspGas(meanPressure, meanTemperature);
    state.gasCells[0].u1L = state.gasCells[0].duto.area * leftDensity;
    state.gasCells[0].u1LL = state.gasCells[0].u1L;
    state.gasCells[1].u1LL = state.gasCells[0].u1L;
    state.gasCells[0].VGasL = 0;
    state.gasCells[0].VGasR = 0;
    state.gasCells[1].VGasL = 0;
    state.gasCells[0].massfonteCH = 0.;
    for (int gasCellIndex = 1; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double leftHalfLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        meanPressure -= leftDensity * kGravityUnloadingVariant * leftHalfLength * sin(state.gasCells[gasCellIndex - 1].duto.teta) / units::kPascalPerKgfPerCm2;
        meanTemperature = state.gasCells[gasCellIndex].calor.Textern1;
        const double localDensity = gasCellIndex < state.interfaceCell
                                ? state.gasCells[gasCellIndex].flui.MasEspGas(meanPressure, meanTemperature)
                                : state.gasCells[gasCellIndex].MasEspFlu(meanPressure, meanTemperature);
        meanPressure -= localDensity * kGravityUnloadingVariant * halfLocalLength * sin(state.gasCells[gasCellIndex].duto.teta) / units::kPascalPerKgfPerCm2;
        leftDensity = localDensity;

        state.gasCells[gasCellIndex].pres = meanPressure;
        state.gasCells[gasCellIndex].presini = meanPressure;
        state.gasCells[gasCellIndex - 1].presR = meanPressure;
        state.gasCells[gasCellIndex].temp = meanTemperature;
        state.gasCells[gasCellIndex - 1].tempR = meanTemperature;
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].duto.area * leftDensity;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        state.gasCells[gasCellIndex].VGasR = 0;
        state.gasCells[gasCellIndex - 1].VGasRR = 0;
        state.gasCells[gasCellIndex].massfonteCH = 0.;
        if (gasCellIndex < state.gasCellCount) {
            state.gasCells[gasCellIndex + 1].presL = meanPressure;
            state.gasCells[gasCellIndex + 1].tempL = meanTemperature;
            state.gasCells[gasCellIndex + 1].u1LL = state.gasCells[gasCellIndex].u1L;
            state.gasCells[gasCellIndex + 1].VGasL = 0;
        }
    }
}

void computeProductionUnloadingHydrostatics(const GasLiftState &state) {
    state.cells[0].massfonteCH = 0;
    state.cells[0].fontemassCL = 0;
    state.cells[0].fontemassLL = 0;
    state.cells[0].fontemassGL = 0;
    state.cells[0].fontemassCR = 0;
    state.cells[0].fontemassLR = 0;
    state.cells[0].fontemassGR = 0;
    const bool inletBelowInterface = 0 <= state.input.celdescargaP;
    state.cells[0].alf = inletBelowInterface ? 0. : 1.;
    state.cells[0].bet = inletBelowInterface ? 1. : 0.;
    state.cells[0].betI = state.cells[0].bet;
    state.cells[0].alfini = state.cells[0].alf;
    state.cells[0].betini = state.cells[0].bet;
    state.cells[0].alfPigD = state.cells[0].alf;
    state.cells[0].betPigD = state.cells[0].bet;
    state.cells[0].alfPigE = state.cells[0].alf;
    state.cells[0].betPigE = state.cells[0].bet;
    state.cells[1].alfL = state.cells[0].alf;
    state.cells[1].betL = state.cells[0].bet;
    state.cells[1].betLI = state.cells[0].bet;
    state.cells[1].alfLini = state.cells[0].alf;
    state.cells[1].betLini = state.cells[0].bet;
    double meanPressure = state.cells[0].acsr.ipr.Pres;
    state.cells[0].presL = meanPressure;
    state.cells[1].presL = meanPressure;
    state.cells[0].pres = meanPressure;
    state.cells[0].presini = meanPressure;
    state.cells[1].presLini = meanPressure;
    state.cells[0].presauxL = meanPressure;
    double meanTemperature = state.cells[0].calor.Textern1;
    state.cells[0].tempL = meanTemperature;
    state.cells[0].temp = meanTemperature;
    state.cells[0].tempini = meanTemperature;
    state.cells[1].tempL = meanTemperature;
    double leftDensity = inletBelowInterface ? state.cells[0].fluicol.MasEspFlu(meanPressure, meanTemperature)
                                      : state.cells[0].flui.MasEspGas(meanPressure, meanTemperature);
    double localDensity = state.cells[0].fluicol.MasEspFlu(meanPressure, meanTemperature);

    state.cells[0].FW = 0;
    state.cells[0].FWini = 0;
    state.cells[0].arranjo = 0;
    state.cells[0].QLL = 0;
    state.cells[1].QLL = state.cells[0].QLL;
    state.cells[0].QL = 0;
    state.cells[0].QG = 0;
    state.cells[0].rpL = localDensity;
    state.cells[0].rpC = localDensity;
    state.cells[1].rpL = state.cells[0].rpC;
    state.cells[0].rcL = localDensity;
    state.cells[0].rcC = localDensity;
    state.cells[1].rcL = state.cells[0].rcC;

    state.cells[0].MC = 0.;
    state.cells[0].ML = 0.;
    state.cells[1].ML = state.cells[0].MC;
    state.cells[0].Mliqini = state.cells[0].MC;
    state.cells[1].MliqiniL = state.cells[0].MC;
    state.cells[0].MliqiniL = state.cells[0].MC;

    state.cells[0].rpLi = localDensity;
    state.cells[0].rpCi = localDensity;
    state.cells[0].rcLi = localDensity;
    state.cells[0].rcCi = localDensity;

    for (int i = 1; i <= state.lastCell; i++) {
        const bool belowInterface = i <= state.input.celdescargaP;
        double leftHalfLength = 0.5 * state.cells[i].dxL;
        double localHalfLength = 0.5 * state.cells[i].dx;
        meanPressure -= leftDensity * kGravityUnloadingVariant * leftHalfLength * sin(state.cells[i - 1].duto.teta) / units::kPascalPerKgfPerCm2;
        meanTemperature = state.cells[i].calor.Textern1;
        double interpolatedTemperature = (leftHalfLength * state.cells[i - 1].temp + localHalfLength * meanTemperature) / (leftHalfLength + localHalfLength);
        state.cells[i].presaux = meanPressure;
        state.cells[i - 1].presauxR = state.cells[i].presaux;
        state.cells[i].presauxL = state.cells[i - 1].presaux;
        if (belowInterface)
            localDensity = state.cells[i].fluicol.MasEspFlu(meanPressure, meanTemperature);
        else
            localDensity = state.cells[i].flui.MasEspGas(meanPressure, meanTemperature);
        meanPressure -= localDensity * kGravityUnloadingVariant * localHalfLength * sin(state.cells[i].duto.teta) / units::kPascalPerKgfPerCm2;
        leftDensity = localDensity;

        state.cells[i].presL = state.cells[i - 1].pres;
        state.cells[i].pres = meanPressure;
        state.cells[i - 1].presR = meanPressure;
        state.cells[i].presini = meanPressure;
        state.cells[i].presLini = state.cells[i - 1].presini;
        state.cells[i].presauxL = state.cells[i - 1].presaux;
        state.cells[i].tempL = state.cells[i - 1].temp;
        state.cells[i].temp = meanTemperature;
        state.cells[i - 1].tempR = meanTemperature;
        state.cells[i].tempini = meanTemperature;

        state.cells[i].FW = 0;
        state.cells[i].FWini = 0;
        state.cells[i].arranjo = 0;
        state.cells[i].QL = 0;
        state.cells[i - 1].QLR = state.cells[i].QL;
        state.cells[i].QG = 0;
        state.cells[i].rpL = localDensity;
        state.cells[i].rpC = localDensity;
        state.cells[i - 1].rpR = state.cells[i].rpC;
        state.cells[i].rcL = localDensity;
        state.cells[i].rcC = localDensity;
        state.cells[i - 1].rcR = state.cells[i].rcC;

        const double interpolatedDensity = belowInterface ? state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, interpolatedTemperature)
                                             : state.cells[0].flui.MasEspGas(state.cells[i].presaux, interpolatedTemperature);
        state.cells[i].rpCi = interpolatedDensity;
        state.cells[i - 1].rpRi = interpolatedDensity;
        state.cells[i].rcCi = interpolatedDensity;
        state.cells[i - 1].rcRi = interpolatedDensity;

        state.cells[i].massfonteCH = 0;
        state.cells[i].fontemassCL = 0;
        state.cells[i].fontemassLL = 0;
        state.cells[i].fontemassGL = 0;
        state.cells[i].fontemassCR = 0;
        state.cells[i].fontemassLR = 0;
        state.cells[i].fontemassGR = 0;
        state.cells[i].alf = belowInterface ? 0. : 1.;
        state.cells[i].bet = belowInterface ? 1. : 0.;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].betPigE = state.cells[i].bet;

        if (i < state.lastCell) {
            state.cells[i + 1].tempL = meanTemperature;
            state.cells[i + 1].QLL = state.cells[i].QL;
            state.cells[i + 1].rpL = state.cells[i].rpC;
            state.cells[i + 1].rcL = state.cells[i].rcC;
            state.cells[i + 1].rpLi = interpolatedDensity;
            state.cells[i + 1].rcLi = interpolatedDensity;
            state.cells[i + 1].alfL = state.cells[i].alf;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
            state.cells[i + 1].alfLini = state.cells[i].alf;
            state.cells[i + 1].betLini = state.cells[i].bet;
        }
    }
}

double calibratedValveArea(double calibrationPressure, double calibrationTemperature, double valveOpeningPressure, double tubingPressure,
                           double externalDiameter, double throatArea, double valveRatio, double bottomHoleTemperatureFahrenheit) {
    // Imperial throughout: pressures in psi, temperature in Fahrenheit. The
    // valve opening pressure is the casing pressure; the ratio is of areas.

    double bellowsPressureAt80F = calibrationPressure * (1 - valveRatio);
    const double calibrationTemperatureRankine = units::celsiusToFahrenheit(calibrationTemperature) + units::kZeroFahrenheitInRankine;
    bellowsPressureAt80F = (bellowsPressureAt80F + units::kAtmosphereInPsi) * (kBellowsReferenceTemperatureFahrenheit + units::kZeroFahrenheitInRankine) / calibrationTemperatureRankine - units::kAtmosphereInPsi;
    double bellowsPressure = bellowsPressureAt80F * (1 + kBellowsPressureChangePerFahrenheit * (bottomHoleTemperatureFahrenheit - kBellowsReferenceTemperatureFahrenheit));
    // A closed valve is not forced to zero here: the opening fraction returned is
    // always the area ratio computed at the end.
    double bellowsArea = throatArea / valveRatio;
    // The spring rate depends only on whether the valve is a large-bore one.
    const double externalDiameterInches = externalDiameter * units::kCentimetrePerMetre / units::kCentimetrePerInch;
    const double bellowsSpringRate =
        (externalDiameterInches > kLargeBoreValveInches ? kLargeBoreSpringRate : kSmallBoreSpringRate) * bellowsArea;
    double stemTravel = ((valveOpeningPressure - bellowsPressure) * bellowsArea - (valveOpeningPressure - tubingPressure) * throatArea) / bellowsSpringRate;

    double throatDiameter = sqrt(throatArea * 4. / std::numbers::pi);
    double throatRadius = throatDiameter / 2.0;
    double bellowsRadius = sqrt(bellowsArea / std::numbers::pi);
    // The seat offset: the distance from the stem axis to where the bellows
    // meets the seat.
    const double seatOffset = sqrt(bellowsRadius * bellowsRadius - throatRadius * throatRadius);
    double openingArea = std::numbers::pi * throatRadius * stemTravel * (stemTravel + 2.0 * seatOffset);
    openingArea = openingArea / sqrt((stemTravel + seatOffset) * (stemTravel + seatOffset) + throatRadius * throatRadius);
    if (openingArea > throatArea)
        openingArea = throatArea;

    return openingArea / throatArea;
}

double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate, int valveIndex, double factor, int sign) {
    int gasValveCell = state.gasValveCellIndices[valveIndex];
    state.gasLiftChokes[valveIndex].presEstag = state.gasCells[gasValveCell].pres;
    state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[gasValveCell].temp;
    double completionFluidDensity = state.gasCells[gasValveCell].MasEspFlu(state.gasLiftChokes[valveIndex].presEstag, state.gasLiftChokes[valveIndex].tempEstag);
    double massFlowRate = (factor * state.input.vazDescControl - maximumFlowRate) * completionFluidDensity;
    const double pressureCorrection =
        pow(massFlowRate / state.gasLiftChokes[valveIndex].areagarg, 2.) / (2. * completionFluidDensity * units::kPascalPerKgfPerCm2);
    return sign * pressureCorrection;
}

void computeUnloadingValvePressure(const GasLiftState &state, double throatFlowRate, int valveIndex) {
    double meanPressure = state.separatorPressure;
    state.input.presMaxDesc = kUnboundedUnloadingPressure;
    for (int cellIndex = state.lastCell; cellIndex >= 0; cellIndex--) {
        double meanTemperature = state.cells[cellIndex].temp;
        double flowArea = state.cells[cellIndex].duto.area;
        double perimeter = state.cells[cellIndex].duto.peri;
        double cellLength = state.cells[cellIndex].dx;
        double gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        double productionLiquidDensity = state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature);
        double completionFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        double gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        double oilViscosity = state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature);
        double completionFluidViscosity = state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
        double voidFraction = state.cells[cellIndex].alf;
        double completionFluidFraction = state.cells[cellIndex].bet;
        double liquidDensity = completionFluidFraction * completionFluidDensity + (1. - completionFluidFraction) * productionLiquidDensity;
        double liquidViscosity = completionFluidFraction * completionFluidViscosity + (1. - completionFluidFraction) * oilViscosity;
        double mixtureDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;
        double mixtureViscosity = voidFraction * gasViscosity + (1. - voidFraction) * liquidViscosity;
        double mixtureVelocity = state.cells[cellIndex].QL / (flowArea * liquidDensity) + state.cells[cellIndex].QG / (flowArea * gasDensity);
        double frictionFactor = 0.;
        if (fabs(mixtureVelocity) > 1e-15) {
            const double reynolds = state.cells[cellIndex].Rey(
                characteristicDiameter(state.cells[cellIndex].duto.revest, state.cells[cellIndex].duto.a, flowArea, perimeter),
                mixtureVelocity, mixtureDensity, mixtureViscosity);
            frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / state.cells[cellIndex].duto.a);
        }
        double wallShearStress = frictionFactor * mixtureDensity * mixtureVelocity * fabs(mixtureVelocity) / 2.;
        meanPressure -= (-kGravity * mixtureDensity * sin(state.cells[cellIndex].duto.teta) - wallShearStress * perimeter / flowArea) * cellLength / units::kPascalPerKgfPerCm2;
        if (state.cells[cellIndex].acsr.tipo == kAccessoryInflowPerformance) {
            double candidateMaximumPressure = state.cells[cellIndex].acsr.ipr.Pres - (meanPressure - state.separatorPressure);
            if (candidateMaximumPressure < state.input.presMaxDesc)
                state.input.presMaxDesc = candidateMaximumPressure;
        }
    }
    if (state.separatorPressure >= state.input.presMaxDesc * kReachedMaximumPressureFactor && throatFlowRate > (1 - kUnloadingUpperBandMargin) * state.input.vazDescControl && state.gasCells[0].VGasR > 0.) {
        const double pressureCorrection = limitedCorrection(unloadingPressureCorrection(state, throatFlowRate, valveIndex, (1 - kUnloadingUpperBandMargin), -1), state.injectionPressure, state.cells[0].dt);
        state.injectionPressure += pressureCorrection;
        if (state.injectionPressure < state.input.presMinDescG) {
            state.injectionPressure = state.input.presMinDescG;
        }

    } else if (throatFlowRate <= (1 - kUnloadingLowerBandMargin) * state.input.vazDescControl && state.separatorPressure <= state.input.presMinDesc * kReachedMinimumPressureFactor) {
        const double pressureCorrection = limitedCorrection(unloadingPressureCorrection(state, throatFlowRate, valveIndex, 1 - kUnloadingLowerBandMargin, 1), state.injectionPressure, state.cells[0].dt);
        state.injectionPressure += pressureCorrection;
        if (state.injectionPressure > state.input.presMaxDescG) {
            state.injectionPressure = state.input.presMaxDescG;
        }
    }
}

double searchUnloadingInjectionPressure(const GasLiftState &state) {
    double maximumUnloadingFlowRate = 0;
    double instantaneousMaximumFlowRate = 0;

    if (state.input.descarga == 1) {
        double initialUnloadingPressure = state.cells[state.lastCell].pres;
        if (state.globals->lixo5 > state.input.tempoLatenciaDesc) {
            int valveCount = state.input.nvalvgas;
            int selectedValveIndex = 0;
            for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
                int gasValveCell = state.gasValveCellIndices[valveIndex];
                int productionValveCell = state.productionValveCellIndices[valveIndex];
                if (gasValveCell > state.interfaceCell) {
                    const CelG &valveCell = state.gasCells[gasValveCell];
                    double completionFluidDensity = valveCell.MasEspFlu(valveCell.pres, valveCell.temp);
                    const double massFlowRate = primeValveMassSource<UnloadingThroughValve>(
                        state, valveIndex, valveCell, state.cells[productionValveCell]);
                    double candidateFlowRate = massFlowRate / completionFluidDensity;
                    if (candidateFlowRate > maximumUnloadingFlowRate) {
                        maximumUnloadingFlowRate = candidateFlowRate;
                        selectedValveIndex = valveIndex;
                    }
                }
            }
            if (std::ssize(state.unloadingTimeStepSamples) > state.unloadingSampleCapacity || state.unloadingAveragingInterval > state.unloadingAveragingWindow) {
                state.unloadingAveragingInterval -= state.unloadingTimeStepSamples[0];
                state.unloadingTimeStepSamples.erase(state.unloadingTimeStepSamples.begin());
                state.meanUnloadingFlowRate -= state.maximumUnloadingFlowSamples[0];
                state.maximumUnloadingFlowSamples.erase(state.maximumUnloadingFlowSamples.begin());
            }

            state.maximumUnloadingFlowSamples.push_back(maximumUnloadingFlowRate * state.timeStep);
            state.meanUnloadingFlowRate += maximumUnloadingFlowRate * state.timeStep;
            state.unloadingTimeStepSamples.push_back(state.timeStep);
            state.unloadingAveragingInterval += state.timeStep;

            instantaneousMaximumFlowRate = maximumUnloadingFlowRate;
            maximumUnloadingFlowRate = kCurrentFlowWeight * maximumUnloadingFlowRate + (1. - kCurrentFlowWeight) * state.meanUnloadingFlowRate / state.unloadingAveragingInterval;
            computeUnloadingValvePressure(state, maximumUnloadingFlowRate, selectedValveIndex);
            if (state.gasValveCellIndices[selectedValveIndex] >= state.interfaceCell) {
                if (state.cells[state.lastCell - 1].MC > 0.0) {
                    if (maximumUnloadingFlowRate > (1. - kUnloadingUpperBandMargin) * state.input.vazDescControl) {
                        const double pressureCorrection = limitedCorrection(unloadingPressureCorrection(state, maximumUnloadingFlowRate, selectedValveIndex, 1. - kUnloadingUpperBandMargin, -1), state.separatorPressure, state.cells[0].dt);
                        state.separatorPressure -= pressureCorrection;
                        if (state.separatorPressure > state.input.presMaxDesc)
                            state.separatorPressure = state.input.presMaxDesc;
                    } else if (maximumUnloadingFlowRate <= (1 - kUnloadingLowerBandMargin) * state.input.vazDescControl) {
                        const double pressureCorrection = limitedCorrection(unloadingPressureCorrection(state, maximumUnloadingFlowRate, selectedValveIndex, 1. - kUnloadingLowerBandMargin, 1), state.separatorPressure, state.cells[0].dt);
                        state.separatorPressure -= pressureCorrection;
                        if (state.separatorPressure < state.input.presMinDesc)
                            state.separatorPressure = state.input.presMinDesc;
                    }
                } else {
                    state.separatorPressure *= (1 - kSeparatorPressureDecayPerSecond * state.timeStep);
                    if (state.separatorPressure < state.input.presMinDesc)
                        state.separatorPressure = state.input.presMinDesc;
                }
            }
        } else {
            if (state.globals->lixo5 < 0.5 * state.input.tempoLatenciaDesc)
                state.separatorPressure = initialUnloadingPressure - (initialUnloadingPressure - state.input.presMinDesc) * state.cells[0].dt / (0.5 * state.input.tempoLatenciaDesc - state.globals->lixo5);
            else
                state.separatorPressure = state.input.presMinDesc;
            state.injectionPressure = state.input.presIniDescG;
        }
    } else {
        state.separatorPressure *= kRelaxationDownFactor;
        if (state.separatorPressure < state.input.presMinDesc)
            state.separatorPressure = state.input.presMinDesc;
        relaxTowards(state.injectionPressure, state.input.presIniDesc);
    }
    return instantaneousMaximumFlowRate;
}

void solveUnloading(const GasLiftState &state) {
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double leftHalfLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        double localGasShare = 1.;
        if (state.gasCells[gasCellIndex].razInter <= 0.5)
            localGasShare = 2 * state.gasCells[gasCellIndex].razInter;
        double leftGasShare = 0.;
        if (state.gasCells[gasCellIndex - 1].razInter >= 0.5)
            leftGasShare = 2 * (state.gasCells[gasCellIndex - 1].razInter - 0.5);
        double leftGasLength = leftHalfLength * leftGasShare;
        double localGasLength = halfLocalLength * localGasShare;
        double leftLiquidLength = leftHalfLength - leftGasLength;
        double localLiquidLength = halfLocalLength - localGasLength;
        double leftTemperature = state.gasCells[gasCellIndex - 1].temp;
        double leftCellPressure = state.gasCells[gasCellIndex - 1].pres;
        double liquidDensity = state.gasCells[gasCellIndex].MasEspFlu(leftCellPressure, leftTemperature);
        double liquidViscosity = state.gasCells[gasCellIndex].VisFlu(leftCellPressure, leftTemperature);
        double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(leftCellPressure, leftTemperature);
        double gasViscosity = state.gasCells[gasCellIndex].flui.ViscGas(leftCellPressure, leftTemperature);
        const double carryingDensity =
            state.gasCells[gasCellIndex].razInter > state.globals->localtiny ? gasDensity : liquidDensity;
        double leftVelocity = state.gasCells[gasCellIndex].VGasL / (carryingDensity * state.gasCells[gasCellIndex - 1].duto.area);
        double localVelocity = state.gasCells[gasCellIndex].VGasL / (carryingDensity * state.gasCells[gasCellIndex].duto.area);
        const double leftDiameter = characteristicDiameter(state.gasCells[gasCellIndex - 1].duto.revest, state.gasCells[gasCellIndex - 1].duto.a, state.gasCells[gasCellIndex - 1].duto.area, state.gasCells[gasCellIndex - 1].duto.peri);
        const double localDiameter = characteristicDiameter(state.gasCells[gasCellIndex].duto.revest, state.gasCells[gasCellIndex].duto.a, state.gasCells[gasCellIndex].duto.area, state.gasCells[gasCellIndex].duto.peri);
        const double leftLiquidReynolds = state.gasCells[gasCellIndex - 1].Rey(leftDiameter, leftVelocity, liquidDensity, liquidViscosity);
        const double localLiquidReynolds = state.gasCells[gasCellIndex].Rey(localDiameter, localVelocity, liquidDensity, liquidViscosity);
        const double leftGasReynolds = state.gasCells[gasCellIndex - 1].Rey(leftDiameter, leftVelocity, gasDensity, gasViscosity);
        const double localGasReynolds = state.gasCells[gasCellIndex].Rey(localDiameter, localVelocity, gasDensity, gasViscosity);
        double leftLiquidFriction = state.gasCells[gasCellIndex - 1].fric(leftLiquidReynolds, state.gasCells[gasCellIndex - 1].duto.rug / state.gasCells[gasCellIndex - 1].duto.a) * leftLiquidLength;
        double localLiquidFriction = state.gasCells[gasCellIndex].fric(localLiquidReynolds, state.gasCells[gasCellIndex].duto.rug / state.gasCells[gasCellIndex].duto.a) * localLiquidLength;
        double leftLiquidHead = (kGravity * sin(state.gasCells[gasCellIndex - 1].duto.teta) * liquidDensity) * leftLiquidLength;
        double localLiquidHead = (kGravity * sin(state.gasCells[gasCellIndex].duto.teta) * liquidDensity) * localLiquidLength;
        double leftGasFriction = state.gasCells[gasCellIndex - 1].fric(leftGasReynolds, state.gasCells[gasCellIndex - 1].duto.rug / state.gasCells[gasCellIndex - 1].duto.a) * leftGasLength;
        double localGasFriction = state.gasCells[gasCellIndex].fric(localGasReynolds, state.gasCells[gasCellIndex].duto.rug / state.gasCells[gasCellIndex].duto.a) * localGasLength;
        double leftGasHead = (kGravity * sin(state.gasCells[gasCellIndex - 1].duto.teta) * gasDensity) * leftGasLength;
        double localGasHead = (kGravity * sin(state.gasCells[gasCellIndex].duto.teta) * gasDensity) * localGasLength;
        state.gasCells[gasCellIndex].pres = state.gasCells[gasCellIndex - 1].pres + (-0.5 * (leftLiquidFriction * liquidDensity + leftGasFriction * gasDensity) * leftVelocity * fabs(leftVelocity) * state.gasCells[gasCellIndex - 1].duto.peri / state.gasCells[gasCellIndex - 1].duto.area - 0.5 * (localLiquidFriction * liquidDensity + localGasFriction * gasDensity) * localVelocity * fabs(localVelocity) * state.gasCells[gasCellIndex].duto.peri / state.gasCells[gasCellIndex].duto.area - leftLiquidHead - localLiquidHead - leftGasHead - localGasHead) / units::kPascalPerKgfPerCm2;
        state.gasCells[gasCellIndex].presL = state.gasCells[gasCellIndex - 1].pres;
        state.gasCells[gasCellIndex - 1].presR = state.gasCells[gasCellIndex].pres;

        state.temperatureUpdater.dischargeTemperature(gasCellIndex);

        state.gasCells[gasCellIndex].u1L = ((1. - state.gasCells[gasCellIndex].razInter) * state.gasCells[gasCellIndex].MasEspFlu(state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp) + state.gasCells[gasCellIndex].razInter * state.gasCells[gasCellIndex].flui.MasEspGas(leftCellPressure, leftTemperature)) * state.gasCells[gasCellIndex].duto.area;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        state.gasCells[gasCellIndex].u1LL = state.gasCells[gasCellIndex - 1].u1L;
    }

    double totalSourceVolumeFlow = 0.;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        const double cellTemperature = state.gasCells[gasCellIndex].temp;
        const double cellPressure = state.gasCells[gasCellIndex].pres;
        const double liquidDensity = state.gasCells[gasCellIndex].MasEspFlu(cellPressure, cellTemperature);
        const double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(cellPressure, cellTemperature);
        totalSourceVolumeFlow += sourceVolumeFlowRate(state.gasCells[gasCellIndex], liquidDensity, gasDensity);
    }
    const double interfaceTemperature = state.gasCells[state.interfaceCell].temp;
    const double interfacePressure = state.gasCells[state.interfaceCell].pres;
    const double interfaceGasDensity = state.gasCells[state.interfaceCell].flui.MasEspGas(interfacePressure, interfaceTemperature);
    state.gasCells[state.interfaceCell].VGasL = totalSourceVolumeFlow * interfaceGasDensity;
    state.gasCells[state.interfaceCell - 1].VGasR = state.gasCells[state.interfaceCell].VGasL;
    state.gasCells[state.interfaceCell - 2].VGasRR = state.gasCells[state.interfaceCell].VGasL;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        const double cellTemperature = state.gasCells[gasCellIndex].temp;
        const double cellPressure = state.gasCells[gasCellIndex].pres;
        const double liquidDensity = state.gasCells[gasCellIndex].MasEspFlu(cellPressure, cellTemperature);
        const double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(cellPressure, cellTemperature);
        totalSourceVolumeFlow -= sourceVolumeFlowRate(state.gasCells[gasCellIndex], liquidDensity, gasDensity);
        if (gasCellIndex < state.gasCellCount) {
            state.gasCells[gasCellIndex + 1].VGasL = totalSourceVolumeFlow * liquidDensity;
            state.gasCells[gasCellIndex].VGasR = state.gasCells[gasCellIndex + 1].VGasL;
        } else
            state.gasCells[gasCellIndex].VGasR = 0.;
        state.gasCells[gasCellIndex - 1].VGasRR = state.gasCells[gasCellIndex].VGasL;
    }

    state.interfaceVelocity = state.gasCells[state.interfaceCell + 1].VGasL / (state.gasCells[state.interfaceCell + 1].MasEspFlu(state.gasCells[state.interfaceCell + 1].pres, state.gasCells[state.interfaceCell + 1].temp) * state.gasCells[state.interfaceCell].duto.area);
}

void advanceInterface(const GasLiftState &state) {
    state.gasCells[state.interfaceCell].razInter = (state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 + state.interfaceVelocity * state.timeStep) / state.gasCells[state.interfaceCell].dx0;
    if (state.interfaceCell == (state.gasCellCount - 1) && state.gasCells[state.interfaceCell].razInter >= kInterfaceHandoverRatio) {
        double meanPressure = state.gasCells[state.interfaceCell].pres;
        double meanTemperature = state.gasCells[state.interfaceCell].temp;
        double interfaceFlowArea = state.gasCells[state.interfaceCell].duto.area;
        double gasDensity = state.gasCells[state.interfaceCell].flui.MasEspGas(meanPressure, meanTemperature);
        double liquidDensity = state.gasCells[state.interfaceCell].MasEspFlu(meanPressure, meanTemperature);
        state.gasCells[state.interfaceCell].VGasR = state.gasCells[state.interfaceCell].VGasR * gasDensity / liquidDensity;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * gasDensity;
        state.gasCells[state.interfaceCell - 1].u1R = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell - 1].VGasRR = state.gasCells[state.interfaceCell].VGasR;
        state.gasCells[state.interfaceCell + 1].u1LL = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell + 1].VGasL = state.gasCells[state.interfaceCell].VGasR;

        state.gasCells[state.interfaceCell].razInter = 1.0;
        state.gasCells[state.interfaceCell + 1].razInter = 1.0;

        state.interfaceCell++;

        state.gasCells[state.interfaceCell].presini = state.gasCells[state.interfaceCell].pres;
        meanPressure = state.gasCells[state.interfaceCell].pres;
        meanTemperature = state.gasCells[state.interfaceCell].temp;
        interfaceFlowArea = state.gasCells[state.interfaceCell].duto.area;
        gasDensity = state.gasCells[state.interfaceCell].flui.MasEspGas(meanPressure, meanTemperature);
        liquidDensity = state.gasCells[state.interfaceCell].MasEspFlu(meanPressure, meanTemperature);
        state.gasCells[state.interfaceCell].VGasR = 0 * state.gasCells[state.interfaceCell].VGasR * gasDensity / liquidDensity;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * gasDensity;
        state.gasCells[state.interfaceCell].u1R = interfaceFlowArea * gasDensity;
        state.gasCells[state.interfaceCell - 1].u1R = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell - 1].VGasRR = state.gasCells[state.interfaceCell].VGasR;

        state.interfaceCell = kNoInterfaceCell;
        state.input.descarga = 0;
    }

    if (state.input.descarga == 1) {
        if (state.gasCells[state.interfaceCell].razInter <= state.globals->localtiny && state.gasCells[state.interfaceCell].razInter >= -state.globals->localtiny)
            state.gasCells[state.interfaceCell].razInter = 0;
        else if (state.gasCells[state.interfaceCell].razInter < -state.globals->localtiny) {
            const double interfaceCrossingTimeStep =
                -state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 / state.interfaceVelocity;
            if (interfaceCrossingTimeStep > state.globals->localtiny)
                state.interfaceTimeStep = interfaceCrossingTimeStep;
            state.gasCells[state.interfaceCell].razInter = 0.;
        } else if (state.gasCells[state.interfaceCell].razInter >= 1. - state.globals->localtiny && state.gasCells[state.interfaceCell].razInter <= 1. + state.globals->localtiny) {
            state.gasCells[state.interfaceCell].razInter = 1.;
        } else if (state.gasCells[state.interfaceCell].razInter > (1. + state.globals->localtiny)) {
            const double interfaceCrossingTimeStep =
                (1. - state.gasCells[state.interfaceCell].razInterIni) * state.gasCells[state.interfaceCell].dx0 / state.interfaceVelocity;
            if (interfaceCrossingTimeStep > state.globals->localtiny)
                state.interfaceTimeStep = interfaceCrossingTimeStep;
            state.gasCells[state.interfaceCell].razInter = 1.;
        }

        if (fabs(state.gasCells[state.interfaceCell].razInter) < state.globals->localtiny && state.interfaceVelocity < -fabs(state.globals->localtiny)) {
            state.gasCells[state.interfaceCell].razInter = 0.;
            if (state.interfaceCell > 0) {
                state.interfaceCell--;
                state.gasCells[state.interfaceCell].razInter = 1.;
            }
        } else if (fabs(state.gasCells[state.interfaceCell].razInter - 1.) < state.globals->localtiny && state.interfaceVelocity > state.globals->localtiny) {
            state.gasCells[state.interfaceCell].razInter = 1.;
            if (state.interfaceCell < state.gasCellCount) {
                state.interfaceCell++;
                state.gasCells[state.interfaceCell].razInter = 0.;
            }
        }
    }
}

void solveGasLine(const GasLiftState &state) {
    if (state.input.lingas <= 0)
        return;

    updateTransientGasValves(state);
    for (int valveIndex = 0; valveIndex < state.input.nvalvgas; valveIndex++) {
        int gasLiftProductionCell = state.productionValveCellIndices[valveIndex];
        int gasLiftGasCell = state.gasValveCellIndices[valveIndex];
        const bool valveDischargesGas =
            gasLiftGasCell < state.interfaceCell ||
            (gasLiftGasCell == state.interfaceCell && state.gasCells[gasLiftGasCell].razInter > 0.5);
        if (valveDischargesGas) {
            state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH * units::kSecondsPerDay / state.gasCells[gasLiftGasCell].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius);
            state.cells[gasLiftProductionCell].acsr.injg.tipoflu = kInjectedGas;
            const double gasTemperature = state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg
                                              ? state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex)
                                              : state.cells[gasLiftProductionCell].temp;
            state.cells[gasLiftProductionCell].acsr.injg.temp = gasTemperature;
            if (state.cells[gasLiftProductionCell].acsr.injg.temp < kMinimumTemperatureCelsius)
                state.cells[gasLiftProductionCell].acsr.injg.temp = kMinimumTemperatureCelsius;
        } else {
            state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH;
            state.cells[gasLiftProductionCell].acsr.injg.tipoflu = kInjectedCompletionFluid;
            state.cells[gasLiftProductionCell].acsr.injg.temp = state.gasLiftChokes[valveIndex].tempEstag;
        }
        state.gasCells[gasLiftGasCell].pEstag = state.gasLiftChokes[valveIndex].presEstag;
        state.gasCells[gasLiftGasCell].tEstag = state.gasLiftChokes[valveIndex].tempEstag;
        state.gasCells[gasLiftGasCell].pGarg = state.gasLiftChokes[valveIndex].presGarg;
        state.gasCells[gasLiftGasCell].tGarg = state.gasLiftChokes[valveIndex].tempGarg;
        state.gasCells[gasLiftGasCell].qGarg = state.gasLiftChokes[valveIndex].qGarg;
        state.gasCells[gasLiftGasCell].areaGarg = state.gasLiftChokes[valveIndex].areagarg;
    }
    advanceGasSubStep(state);
}

void updateTransientGasValves(const GasLiftState &state) {
    int valveCount = state.input.nvalvgas;
    int master2CellIndex = state.input.master2.posic;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].massfonteCH = 0.;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        const int valveCellIndex = state.gasValveCellIndices[valveIndex];
        CelG &valveCell = state.gasCells[valveCellIndex];
        const Cel &productionCell = state.cells[state.productionValveCellIndices[valveIndex]];
        // Behind the interface the valve discharges completion fluid, ahead of it gas.
        valveCell.massfonteCH =
            valveCellIndex >= state.interfaceCell && valveCell.razInter < 0.5
                ? primeValveMassSource<CompletionFluidThroughValve>(state, valveIndex, valveCell, productionCell)
                : primeValveMassSource<GasThroughValve>(state, valveIndex, valveCell, productionCell);
    }

    double master2FlowArea = state.gasCells[master2CellIndex].duto.area;
    if (master2FlowArea > state.gasCells[master2CellIndex].dutoR.area)
        master2FlowArea = state.gasCells[master2CellIndex].dutoR.area;
    if (state.gasCells[master2CellIndex].chkcell.areagarg <= kMaster2ThrottlingAreaFraction * master2FlowArea &&
        state.gasCells[master2CellIndex].chkcell.areagarg >= kMaster2ClosedAreaFraction * master2FlowArea) {
        state.gasCells[master2CellIndex].chkcell.tempEstag = state.gasCells[master2CellIndex].temp;
        double upstreamGasDensity = state.gasCells[master2CellIndex].flui.MasEspGas(state.gasCells[master2CellIndex].pres, state.gasCells[master2CellIndex].temp);
        double downstreamGasDensity = state.gasCells[master2CellIndex + 1].flui.MasEspGas(state.gasCells[master2CellIndex + 1].pres, state.gasCells[master2CellIndex + 1].temp);
        double upstreamHead = -0.5 * upstreamGasDensity * units::kStandardGravity * state.gasCells[master2CellIndex].dx0 * sin(state.gasCells[master2CellIndex].duto.teta) / units::kPascalPerKgfPerCm2;
        double downstreamHead = 0.5 * downstreamGasDensity * units::kStandardGravity * state.gasCells[master2CellIndex + 1].dx0 * sin(state.gasCells[master2CellIndex + 1].duto.teta) / units::kPascalPerKgfPerCm2;
        state.gasCells[master2CellIndex].chkcell.presEstag = state.gasCells[master2CellIndex].pres + upstreamHead;
        state.gasCells[master2CellIndex].chkcell.presGarg = state.gasCells[master2CellIndex + 1].pres + downstreamHead;
        if (state.gasCells[master2CellIndex + 1].pres < state.gasCells[master2CellIndex].pres)
            state.gasCells[master2CellIndex].fonteM2 = -state.gasCells[master2CellIndex].chkcell.massica();
        else {
            state.gasCells[master2CellIndex].chkcell.presEstag = state.gasCells[master2CellIndex + 1].pres;
            state.gasCells[master2CellIndex].chkcell.presGarg = state.gasCells[master2CellIndex].pres;
            state.gasCells[master2CellIndex].chkcell.tempEstag = state.gasCells[master2CellIndex + 1].temp;
            state.gasCells[master2CellIndex].fonteM2 = state.gasCells[master2CellIndex].chkcell.massica();
        }
        state.gasCells[master2CellIndex + 1].fonteM2 = -state.gasCells[master2CellIndex].fonteM2;
    } else {
        state.gasCells[master2CellIndex].fonteM2 = 0.;
        state.gasCells[master2CellIndex + 1].fonteM2 = 0.;
    }
}

void advanceGasSubStep(const GasLiftState &state) {
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].DeVoltaParaoFuturo();
    state.interfaceTimeStep = state.timeStep;
    if (state.input.descarga == 1) {
        state.previousInterfaceCell = state.interfaceCell;
        state.previousInterfaceTimeStep = state.interfaceTimeStep;
        state.previousInterfaceVelocity = state.interfaceVelocity;
        advanceInterface(state);
    }
    if (state.interfaceTimeStep < state.timeStep) {
        state.timeStep = state.interfaceTimeStep;
        for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
            state.cells[cellIndex].dt = state.timeStep;
            state.cells[cellIndex].dt2 = state.timeStep;
            state.cells[cellIndex].dtPig = state.timeStep;
        }
    }

    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        state.gasCells[gasCellIndex].dt = state.timeStep;
    }

    const int cellLimit = state.interfaceCell <= state.gasCellCount
                              ? state.interfaceCell - 1
                              : state.gasCellCount + 1;
    for (int gasCellIndex = 0; gasCellIndex < state.gasCellCount + 1; gasCellIndex++)
        state.gasCells[gasCellIndex].dTdt = 0;

    int maximumThermalCouplingCycles = state.input.cicloAcopTerm;
    for (int thermalCouplingCycle = 0; thermalCouplingCycle <= maximumThermalCouplingCycles; thermalCouplingCycle++) {
        assembleAndSolveGasSystem(state);
        updateGasLine(state);
        state.gasCells[0].temp = state.injectionTemperature;
        state.gasCells[0].dTdt = (state.gasCells[0].temp - state.gasCells[0].tempini) / state.timeStep;
        if (state.gasCells[state.gasCellCount].VGasR < 0.)
            state.gasCells[state.gasCellCount].temp = kReverseFlowOutletTemperatureCelsius;
#pragma omp parallel for num_threads(state.globals->ntrd)
        for (int gasCellIndex = 1; gasCellIndex < cellLimit; gasCellIndex++) {
            state.temperatureUpdater.computeGasTemperature(gasCellIndex, state.gasCells[gasCellIndex - 1].tempini);
            state.gasCells[gasCellIndex].dTdt = (state.gasCells[gasCellIndex].temp - state.gasCells[gasCellIndex].tempini) / state.timeStep;
        }
        if (thermalCouplingCycle < maximumThermalCouplingCycles)
            for (int renewedCellIndex = 0; renewedCellIndex <= state.gasCellCount; renewedCellIndex++)
                state.gasCells[renewedCellIndex].FeiticoDoTempo();
    }

#pragma omp parallel for num_threads(state.globals->ntrd)
    for (int gasCellIndex = 0; gasCellIndex < cellLimit; gasCellIndex++)
        state.gasCells[gasCellIndex].rg = state.gasCells[gasCellIndex].flui.MasEspGas(state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp);

    for (int gasCellIndex = 0; gasCellIndex < cellLimit; gasCellIndex++) {
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].rg * state.gasCells[gasCellIndex].duto.area;
        // The inlet feeds itself; every other cell is fed by the one upstream.
        state.gasCells[gasCellIndex].u1LL = gasCellIndex == 0
                                                ? state.gasCells[gasCellIndex].u1L
                                                : state.gasCells[gasCellIndex - 1].u1L;
        if (gasCellIndex > 0) {
            state.gasCells[gasCellIndex - 1].rgR = state.gasCells[gasCellIndex].rg;
            state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        }
        if (gasCellIndex == state.gasCellCount) {
            state.gasCells[gasCellIndex].u1R = state.gasCells[gasCellIndex].u1L;
            state.gasCells[state.gasCellCount].rgR = state.gasCells[state.gasCellCount].rg;
        }
    }

    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].presini = state.gasCells[gasCellIndex].pres;

    if (state.input.descarga == 1)
        solveUnloading(state);
}

void updateGasLine(const GasLiftState &state) {
    // The unknowns of a gas cell in the solved system: its pressure, then the gas mass flow through its right face.
    const auto pressureOf = [&](int cellIndex) { return state.gasRightHandSide[kGasLineUnknownsPerCell * cellIndex]; };
    const auto rightFaceFlowOf = [&](int cellIndex) {
        return state.gasRightHandSide[kGasLineUnknownsPerCell * cellIndex + 1];
    };
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        CelG &cell = state.gasCells[gasCellIndex];
        // The first cell counts as the first even when it is also the last.
        const bool first = gasCellIndex == 0;
        const bool last = !first && gasCellIndex == state.gasCellCount;
        cell.pres = pressureOf(gasCellIndex);
        cell.presL = first ? cell.pres : pressureOf(gasCellIndex - 1);
        cell.presR = last ? cell.pres : pressureOf(gasCellIndex + 1);
        cell.VGasR = rightFaceFlowOf(gasCellIndex);
        if (!first)
            cell.VGasL = rightFaceFlowOf(gasCellIndex - 1);
        if (!last)
            cell.VGasRR = rightFaceFlowOf(gasCellIndex + 1);
    }
}

void advanceBufferedGasSubStep(const GasLiftState &state) {
    for (int gasCellIndex = 0; gasCellIndex < state.gasCellCount + 1; gasCellIndex++)
        state.gasCells[gasCellIndex].dTdt = 0;
    assembleAndSolveGasSystem(state);
    updateBufferedGasLine(state);
}

void updateBufferedGasLine(const GasLiftState &state) {
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].VGasRBuf = state.gasRightHandSide[3 * gasCellIndex + 1];
}

void connectTubing(const GasLiftState &state) {
    for (int productionCellIndex = state.tubingCouplingStart; productionCellIndex >= state.tubingCouplingEnd; productionCellIndex--) {
        int gasCellIndex = state.tubingCouplingStart + state.annulusCouplingStart - productionCellIndex;
        state.cells[productionCellIndex].calor.Textern2 = state.gasCells[gasCellIndex].calor.Tcamada[0][0];
        state.cells[productionCellIndex].calor.betext = state.gasCells[gasCellIndex].calor.betint;
        int layerIndex = state.cells[productionCellIndex].calor.geom.ncamadas - 1;
        int discretizationIndex = state.cells[productionCellIndex].calor.ncamada[layerIndex] - 1;
        state.gasCells[gasCellIndex].calor.Tint2 = state.cells[productionCellIndex].calor.Tcamada[layerIndex][discretizationIndex];
        state.cells[productionCellIndex].calor.colunaDia = state.gasCells[gasCellIndex].duto.dia;
        state.cells[productionCellIndex].calor.geom.b = state.gasCells[gasCellIndex].calor.geom.a;

        state.cells[productionCellIndex].calor.Textern1 = state.gasCells[gasCellIndex].temp;
        state.cells[productionCellIndex].calor.Vextern1 = state.gasCells[gasCellIndex].VGasR / state.gasCells[gasCellIndex].u1L;
        // Ahead of the interface the annulus carries gas; behind it, completion
        // fluid. Same five fields either way, different fluid answering.
        const CelG &gasCell = state.gasCells[gasCellIndex];
        const double gasPressure = gasCell.pres;
        const double gasTemperature = gasCell.temp;
        ExternalFluidProperties fluid;
        if (gasCellIndex < state.interfaceCell) {
            const double specificHeat = gasCell.flui.CalorGas(gasPressure, gasTemperature);
            const double density = gasCell.flui.MasEspGas(gasPressure, gasTemperature);
            fluid = {gasCell.flui.CondGas(gasPressure, gasTemperature), specificHeat, density,
                     gasCell.flui.ViscGas(gasPressure, gasTemperature)};
        } else {
            const double specificHeat = gasCell.CalorLiq(gasPressure, gasTemperature);
            const double density = gasCell.MasEspFlu(gasPressure, gasTemperature);
            fluid = {gasCell.CondLiq(gasPressure, gasTemperature), specificHeat, density,
                     gasCell.VisFlu(gasPressure, gasTemperature)};
        }
        publishExternalFluid(state.cells[productionCellIndex], fluid);
    }
}

void connectTubingSteady(const GasLiftState &state) {
    for (int productionCellIndex = state.tubingCouplingStart; productionCellIndex >= state.tubingCouplingEnd; productionCellIndex--) {
        int gasCellIndex = state.tubingCouplingStart + state.annulusCouplingStart - productionCellIndex;
        state.cells[productionCellIndex].calor.Textern2 = state.gasCells[gasCellIndex].calor.Textern1;
        state.cells[productionCellIndex].calor.colunaDia = state.gasCells[gasCellIndex].duto.dia;

        const CelG &gasCell = state.gasCells[gasCellIndex];
        const double gasPressure = gasCell.pres;
        const double gasTemperature = gasCell.temp;
        const double gasSpecificHeat = gasCell.flui.CalorGas(gasPressure, gasTemperature);
        const double gasDensity = gasCell.flui.MasEspGas(gasPressure, gasTemperature);
        state.cells[productionCellIndex].calor.Textern1 = gasCell.calor.Textern1;
        state.cells[productionCellIndex].calor.Vextern1 = gasCell.VGasR / gasCell.u1L;
        publishExternalFluid(state.cells[productionCellIndex],
                             {gasCell.flui.CondGas(gasPressure, gasTemperature), gasSpecificHeat, gasDensity,
                              gasCell.flui.ViscGas(gasPressure, gasTemperature)});
    }
}

void initializeTubingConnectionSteady(const GasLiftState &state) {
    for (int productionCellIndex = state.tubingCouplingStart; productionCellIndex >= state.tubingCouplingEnd; productionCellIndex--) {
        int gasCellIndex = state.tubingCouplingStart + state.annulusCouplingStart - productionCellIndex;
        state.cells[productionCellIndex].calor.Textern2 = state.gasCells[gasCellIndex].calor.Textern2;
        state.cells[productionCellIndex].calor.colunaDia = state.gasCells[gasCellIndex].duto.dia;

        const CelG &gasCell = state.gasCells[gasCellIndex];
        // Not the gas cell's own state: the facing production cell's pressure and
        // the deck's external temperature.
        const double productionPressure = state.cells[productionCellIndex].pres;
        const double externalTemperature = state.input.celg[gasCellIndex].textern;
        const double gasSpecificHeat = gasCell.flui.CalorGas(productionPressure, externalTemperature);
        const double gasDensity = gasCell.flui.MasEspGas(productionPressure, externalTemperature);
        state.cells[productionCellIndex].calor.Textern1 = gasCell.calor.Textern1;
        state.cells[productionCellIndex].calor.Vextern1 = 1.;
        // The viscosity is a fixed figure here, not a correlation.
        publishExternalFluid(state.cells[productionCellIndex],
                             {gasCell.flui.CondGas(productionPressure, externalTemperature), gasSpecificHeat, gasDensity,
                              kInitialAnnulusGasViscosityCentipoise});
    }
}

double steadyGasPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.gasCells[cellIndex].dx0;
    double diameter = state.gasCells[cellIndex].duto.a;
    double area = 0.25 * std::numbers::pi * diameter * diameter;
    double perimeter = state.gasCells[cellIndex].duto.peri;
    double gasDensity = state.gasCells[cellIndex].flui.MasEspGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double leftFaceGasMassFlowRate = state.gasCells[cellIndex - 1].VGasR;

    double gasVelocity = leftFaceGasMassFlowRate / state.gasCells[cellIndex].u1L;

    double gasViscosity = state.gasCells[cellIndex].flui.ViscGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double reynolds = state.gasCells[cellIndex].Rey(
        characteristicDiameter(state.gasCells[cellIndex].duto.revest, state.gasCells[cellIndex].duto.a, area, perimeter),
        gasVelocity, gasDensity, gasViscosity);
    double frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    double frictionGradient = state.gasCells[cellIndex].dPdLFric * (0.5 * frictionFactor * gasDensity * (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area);
    double hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * (kGravity * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);

    double pressureDrop = (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    double meanPressure = state.gasCells[cellIndex].pres + pressureDrop;

    double meanTemperature = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].temp + state.gasCells[cellIndex].dxL * state.gasCells[cellIndex].tempL) / (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex].dxL);

    dx = 0.5 * state.gasCells[cellIndex].dxL;
    diameter = state.gasCells[cellIndex].dutoL.a;
    area = 0.25 * std::numbers::pi * diameter * diameter;
    perimeter = state.gasCells[cellIndex].dutoL.peri;
    gasDensity = state.gasCells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);

    gasVelocity = leftFaceGasMassFlowRate / (gasDensity * area);

    gasViscosity = state.gasCells[cellIndex - 1].flui.ViscGas(meanPressure, meanTemperature);

    reynolds = state.gasCells[cellIndex].Rey(
        characteristicDiameter(state.gasCells[cellIndex].dutoL.revest, state.gasCells[cellIndex].dutoL.a, area, perimeter),
        gasVelocity, gasDensity, gasViscosity);
    frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
    frictionGradient = state.gasCells[cellIndex - 1].dPdLFric * (0.5 * frictionFactor * gasDensity * (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area);
    hydrostaticGradient = state.gasCells[cellIndex - 1].dPdLHidro * (kGravity * sin(state.gasCells[cellIndex].dutoL.teta) * gasDensity * dx);
    pressureDrop += (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    return pressureDrop;
}

double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.cells[cellIndex].dx;
    double diameter = state.cells[cellIndex].duto.a;
    double area = 0.25 * std::numbers::pi * diameter * diameter;
    double perimeter = state.cells[cellIndex].duto.peri;
    double meanTemperature = (state.cells[cellIndex].dx * state.cells[cellIndex].temp + state.cells[cellIndex].dxL * state.cells[cellIndex - 1].temp) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
    double injectedFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presaux, meanTemperature);

    double fluidVelocity = state.cells[cellIndex - 1].QL / area;

    double fluidViscosity = state.cells[cellIndex].fluicol.VisFlu(state.cells[cellIndex].presaux, meanTemperature);

    double reynolds = state.cells[cellIndex].Rey(
        characteristicDiameter(state.cells[cellIndex].duto.revest, state.cells[cellIndex].duto.a, area, perimeter),
        fluidVelocity, injectedFluidDensity, fluidViscosity);
    double frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / diameter);
    double frictionGradient = state.cells[cellIndex].dPdLFric * 0.5 * frictionFactor * injectedFluidDensity * (fabs(fluidVelocity) * fluidVelocity) * perimeter * dx / area;
    double hydrostaticGradient = state.cells[cellIndex].dPdLHidro * kGravity * sin(state.cells[cellIndex].duto.teta) * injectedFluidDensity * dx;

    double pressureDrop = (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    double meanPressure = state.cells[cellIndex - 1].pres;

    dx = 0.5 * state.cells[cellIndex].dxL;
    diameter = state.cells[cellIndex].dutoL.a;
    area = 0.25 * std::numbers::pi * diameter * diameter;
    perimeter = state.cells[cellIndex].dutoL.peri;
    meanTemperature = state.cells[cellIndex - 1].temp;
    injectedFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);

    if (cellIndex > 1)
        fluidVelocity = state.cells[cellIndex - 2].QL / area;

    fluidViscosity = state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);

    reynolds = state.cells[cellIndex].Rey(
        characteristicDiameter(state.cells[cellIndex].dutoL.revest, state.cells[cellIndex].dutoL.a, area, perimeter),
        fluidVelocity, injectedFluidDensity, fluidViscosity);
    frictionFactor = state.cells[cellIndex - 1].fric(reynolds, state.cells[cellIndex].dutoL.rug / diameter);
    frictionGradient = state.cells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * injectedFluidDensity * (fabs(fluidVelocity) * fluidVelocity) * perimeter * dx / area;
    hydrostaticGradient = state.cells[cellIndex - 1].dPdLHidro * kGravity * sin(state.cells[cellIndex].dutoL.teta) * injectedFluidDensity * dx;
    pressureDrop += (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    return pressureDrop;
}

void updateSteadyGasPressure(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.gasCells[cellIndex].dxL;
    double diameter = state.gasCells[cellIndex].dutoL.a;
    double area = 0.25 * std::numbers::pi * diameter * diameter;
    double perimeter = state.gasCells[cellIndex].dutoL.peri;
    double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

    double leftFaceGasMassFlowRate = state.gasCells[cellIndex - 1].VGasR;
    double gasVelocity = leftFaceGasMassFlowRate / state.gasCells[cellIndex - 1].u1L;

    double gasViscosity = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

    double reynolds = state.gasCells[cellIndex - 1].Rey(
        characteristicDiameter(state.gasCells[cellIndex].dutoL.revest, state.gasCells[cellIndex].dutoL.a, area, perimeter),
        gasVelocity, gasDensity, gasViscosity);
    double frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
    double frictionGradient = state.gasCells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * gasDensity * (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area;
    double hydrostaticGradient = state.gasCells[cellIndex - 1].dPdLHidro * kGravity * sin(state.gasCells[cellIndex].dutoL.teta) * gasDensity * dx;
    state.gasCells[cellIndex - 1].termoFric = frictionGradient / dx;
    state.gasCells[cellIndex - 1].termoHidro = hydrostaticGradient / dx;

    double meanPressure = state.gasCells[cellIndex - 1].pres - (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;

    double meanTemperature = state.gasCells[cellIndex - 1].temp;

    dx = 0.5 * state.gasCells[cellIndex].dx0;
    diameter = state.gasCells[cellIndex].duto.a;
    area = 0.25 * std::numbers::pi * diameter * diameter;
    perimeter = state.gasCells[cellIndex].duto.peri;
    gasDensity = state.gasCells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
    double areaVariation = 1 / pow(state.gasCells[cellIndex].dutoL.area, 2.) - 1 / pow(state.gasCells[cellIndex].duto.area, 2.);
    double dynamicPressureTerm = 0.5 * leftFaceGasMassFlowRate * leftFaceGasMassFlowRate * areaVariation / gasDensity;

    gasVelocity = leftFaceGasMassFlowRate / (gasDensity * area);

    gasViscosity = state.gasCells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);

    reynolds = state.gasCells[cellIndex].Rey(
        characteristicDiameter(state.gasCells[cellIndex].duto.revest, state.gasCells[cellIndex].duto.a, area, perimeter),
        gasVelocity, gasDensity, gasViscosity);
    frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    frictionGradient = state.gasCells[cellIndex].dPdLFric * 0.5 * frictionFactor * gasDensity * (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area;
    hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * kGravity * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx;

    state.gasCells[cellIndex].pres = meanPressure - (frictionGradient + hydrostaticGradient - dynamicPressureTerm) / units::kPascalPerKgfPerCm2;
    state.gasCells[cellIndex - 1].presR = state.gasCells[cellIndex].pres;
    if (cellIndex < state.gasCellCount)
        state.gasCells[cellIndex + 1].presL = state.gasCells[cellIndex].pres;
    state.gasCells[cellIndex].presini = state.gasCells[cellIndex].pres;
}

void computeSteadyGasFlowRate(const GasLiftState &state, int cellIndex) {
    int valveCount = state.input.nvalvgas;
    bool matched = false;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        if (state.gasValveCellIndices[valveIndex] == cellIndex) {
            matched = true;
            CelG &valveCell = state.gasCells[state.gasValveCellIndices[valveIndex]];
            valveCell.massfonteCH = primeValveMassSource<GasThroughValve>(
                state, valveIndex, valveCell,
                state.cells[state.productionValveCellIndices[valveIndex]]);

            if (cellIndex > 0)
                state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR - state.gasCells[cellIndex].massfonteCH;
            else
                state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex].massfonteCH;

            int gasLiftProductionCell = state.productionValveCellIndices[valveIndex];
            int gasLiftGasCell = state.gasValveCellIndices[valveIndex];
            state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH * units::kSecondsPerDay /
                                            (kAirDensityAtStandardConditions * state.gasCells[gasLiftGasCell].flui.Deng);
            const double gasTemperature = state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg
                                              ? state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex)
                                              : state.cells[gasLiftProductionCell].temp;
            state.cells[gasLiftProductionCell].acsr.injg.temp = gasTemperature;
            if (state.cells[gasLiftProductionCell].acsr.injg.temp < kMinimumTemperatureCelsius)
                state.cells[gasLiftProductionCell].acsr.injg.temp = kMinimumTemperatureCelsius;
            state.cells[gasLiftProductionCell].fontemassGR = state.gasCells[gasLiftGasCell].massfonteCH;

            if (state.gasCells[cellIndex].VGasR <= 0. && cellIndex < state.gasCellCount) {
                state.gasCells[cellIndex].VGasR = 0.;
                for (int laterValveIndex = valveIndex + 1; laterValveIndex < valveCount; laterValveIndex++) {
                    state.gasCells[state.gasValveCellIndices[laterValveIndex]].massfonteCH = 0;
                    state.cells[state.productionValveCellIndices[laterValveIndex]].fontemassGR = 0.;
                }
            }

            if (cellIndex > 0)
                state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
            if (cellIndex < state.gasCellCount)
                state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
        }
    }
    if (!matched) {
        state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR - 0 * state.gasCells[cellIndex].massfonteCH;
        state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
        if (cellIndex < state.gasCellCount)
            state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
    }
}

namespace {

/// One valve's gas mass source at steady state. Under a flow-rate condition,
/// an equal share of the injection flow rate; under a pressure condition, what
/// the valve's choke passes from the casing pressure, if the casing is above
/// the production line.
void setSteadyValveGasMassSource(const GasLiftState &state, int valveIndex, int valveCount,
                                 double, InjectionFlowRateCondition) {
    double valveFlowRate = (state.input.gasinj.vazgas[0] * state.gasCells[0].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius) / units::kSecondsPerDay) / valveCount;
    state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = valveFlowRate;
}

void setSteadyValveGasMassSource(const GasLiftState &state, int valveIndex, int,
                                 double casingTemperature, InjectionPressureCondition) {
    int gasLiftGasCell = state.gasValveCellIndices[valveIndex];
    double casingPressure = state.gasCells[gasLiftGasCell].pres;
    state.gasLiftChokes[valveIndex].presEstag = casingPressure;
    state.gasLiftChokes[valveIndex].presGarg = state.cells[state.productionValveCellIndices[valveIndex]].pres;
    state.gasLiftChokes[valveIndex].tempEstag = casingTemperature;
    state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH =
        state.cells[state.productionValveCellIndices[valveIndex]].pres < casingPressure
            ? state.gasLiftChokes[valveIndex].massica()
            : 0.;
}

}  // namespace

void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex) {
    int valveCount = state.input.nvalvgas;
    double casingTemperature = state.gasCells[0].calor.Textern1;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        int gasLiftProductionCell = state.productionValveCellIndices[valveIndex];
        int gasLiftGasCell = state.gasValveCellIndices[valveIndex];
        if (cellIndex == gasLiftProductionCell) {
            // This site reads the flag the other way round from the rest: "1 is the flow
            // rate, anything else the pressure", where every other reads "0 is the
            // pressure, anything else the flow rate". The readings differ only for a
            // value outside {0, 1}, which the input does not reject.
            withGasInletCondition(
                state.gasCells[0].tipoCC == kGasInletInjectionFlowRate ? kGasInletInjectionFlowRate
                                                                       : kGasInletInjectionPressure,
                [&](auto condition) {
                    setSteadyValveGasMassSource(state, valveIndex, valveCount, casingTemperature, condition);
                });
            state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH * units::kSecondsPerDay /
                                            state.gasCells[gasLiftGasCell].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius);
            state.cells[gasLiftProductionCell].acsr.injg.temp = casingTemperature;
            state.cells[gasLiftProductionCell].fontemassGR = state.gasCells[gasLiftGasCell].massfonteCH;
        }
    }
}

void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex) {
    if (state.thermalSourceDisabled == 0) {
        if (cellIndex < state.annulusCouplingStart || cellIndex > state.annulusCouplingEnd) {
            double dx = state.gasCells[cellIndex].dx0;
            double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double meanTemperatureGradient = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].dTdLCor + state.gasCells[cellIndex - 1].dx0 * state.gasCells[cellIndex - 1].dTdLCor) / (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double area = state.gasCells[cellIndex].duto.area;
            double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double meanSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double gasJouleThomson = state.gasCells[cellIndex - 1].flui.JTG(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double hydrostaticTerm = (gasDensity * meanSuperficialGasVelocity) * area * kGravity * sin(state.gasCells[cellIndex].duto.teta);

            state.gasCells[cellIndex - 1].calor.Tint = state.gasCells[cellIndex - 1].temp;
            state.gasCells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity;
            state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex - 1].calor.rhoint = gasDensity;
            state.gasCells[cellIndex - 1].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp) * units::kPascalSecondPerCentipoise;
            state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
            double leftCellHeatFlux = state.gasCells[cellIndex - 1].fluxcal;
            double tubingHeatFlux = 0.;

            const double lengthRatio = dx / (dx + state.gasCells[cellIndex - 1].dx0);
            double temperatureMarchCoefficient = gasDensity * meanSuperficialGasVelocity * gasSpecificHeat * area;
            double pressureMarchCoefficient = gasDensity * meanSuperficialGasVelocity * gasJouleThomson * area;
            double pressureGradient = 2. * (state.gasCells[cellIndex].pres - ((1 - lengthRatio) * state.gasCells[cellIndex - 1].pres + lengthRatio * state.gasCells[cellIndex].pres)) * units::kPascalPerKgfPerCm2 / dx;

            double mixtureDensity = gasDensity;
            double leftCellLeftFaceSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasL / state.gasCells[cellIndex - 1].u1L;
            const double superficialVelocityGradient = (meanSuperficialGasVelocity - leftCellLeftFaceSuperficialGasVelocity) / dx;

            const double kineticTerm =
                state.input.nCompTotalUnidadesG / dx < 1e6
                    ? mixtureDensity * area * meanSuperficialGasVelocity * meanSuperficialGasVelocity * superficialVelocityGradient
                    : 0.;

            double gasMassSource = 0.;
            double liquidMassSource = 0.;

            if (meanSuperficialGasVelocity > kMinimumMarchGasVelocity && (state.gasCells[cellIndex].duto.a / state.gasCells[cellIndex - 1].duto.a > kMinimumMarchDiameterRatio &&
                                  state.gasCells[cellIndex - 1].duto.a / state.gasCells[cellIndex].duto.a > kMinimumMarchDiameterRatio)) {
                double workTerm = meanTemperatureGradient * (pressureMarchCoefficient * pressureGradient - kineticTerm - hydrostaticTerm + liquidMassSource + gasMassSource) / temperatureMarchCoefficient;
                double heatTerm = meanTemperatureGradient * (leftCellHeatFlux + tubingHeatFlux) / temperatureMarchCoefficient;
                int stepCount;
                double stepLength;
                if (meanCellLength / state.gasCells[cellIndex - 1].calor.resGlob < kThermalSubstepLength) {
                    stepCount = 0;
                    stepLength = meanCellLength;
                } else {
                    stepCount = static_cast<int>(2 * (meanCellLength / state.gasCells[cellIndex - 1].calor.resGlob) / kThermalSubstepLength + 1);
                    stepLength = meanCellLength / stepCount;
                }
                double stepTemperature = stepLength * (state.gasCells[cellIndex - 1].temp / stepLength + workTerm + heatTerm);
                for (int stepIndex = 1; stepIndex < stepCount; stepIndex++) {
                    state.gasCells[cellIndex - 1].calor.Tint = stepTemperature;
                    state.gasCells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity;
                    state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, stepTemperature);
                    state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
                    state.gasCells[cellIndex - 1].calor.rhoint = gasDensity;
                    state.gasCells[cellIndex - 1].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, stepTemperature) * units::kPascalSecondPerCentipoise;
                    state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
                    leftCellHeatFlux = state.gasCells[cellIndex - 1].fluxcal;
                    heatTerm = meanTemperatureGradient * (leftCellHeatFlux + tubingHeatFlux) / temperatureMarchCoefficient;
                    stepTemperature = stepLength * (stepTemperature / stepLength + workTerm + heatTerm);
                }

                state.gasCells[cellIndex].temp = stepTemperature;
            } else {
                state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
            }

            clampAndPropagateTemperature(state, cellIndex);
        } else {
            if (state.annulusThermalCoupling == 1) {
                int connectionCellIndex = cellIndex - state.annulusCouplingStart;
                int connectionIndex = state.tubingCouplingStart - connectionCellIndex;
                int secondConnectionIndex = connectionIndex;
                if (secondConnectionIndex == state.lastCell)
                    secondConnectionIndex -= 1;

                state.gasCells[cellIndex].temp = state.cells[secondConnectionIndex].calor.resGlob * state.cells[connectionIndex].calor.fluxFim + state.cells[connectionIndex].temp;

                if (state.input.correcaoContracorPerm == 1) {
                    double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
                    double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
                    double inverseHeatCapacity = 1. / (state.cells[connectionIndex].calor.cpint * state.cells[connectionIndex].MC) + 1. / (gasSpecificHeat * state.gasCells[cellIndex - 1].VGasR);
                    double logTemperature = -exp(-(1. / state.cells[secondConnectionIndex].calor.resGlob) * (1 + 0 * state.cells[connectionIndex].duto.peri) * meanCellLength * inverseHeatCapacity) *
                                         (state.cells[secondConnectionIndex + 1].temp - state.gasCells[cellIndex - 1].temp) +
                                     state.cells[secondConnectionIndex].temp;

                    if (logTemperature < state.gasCells[cellIndex].temp)
                        state.gasCells[cellIndex].temp = logTemperature;
                }
            }

            clampAndPropagateTemperature(state, cellIndex);

            double meanSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
            double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

            state.gasCells[cellIndex].calor.Tint = state.gasCells[cellIndex].temp;
            state.gasCells[cellIndex].calor.Vint = meanSuperficialGasVelocity;
            state.gasCells[cellIndex].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            state.gasCells[cellIndex].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex].calor.rhoint = gasDensity;
            state.gasCells[cellIndex].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp) * units::kPascalSecondPerCentipoise;
            state.gasCells[cellIndex].fluxcal = state.gasCells[cellIndex].calor.transperm();
        }
    } else {
        state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
    }
}

}  // namespace sisprod::gaslift
