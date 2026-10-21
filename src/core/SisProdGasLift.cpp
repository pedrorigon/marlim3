#include "SisProdGasLift.h"
#include "SisProdGasLiftInternal.h"

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

/// The diameter a Reynolds number is formed on: the bore when the segment is
/// bare, the hydraulic diameter when it is cased.
[[nodiscard]] double characteristicDiameter(int cased, double bore, double area,
                                            double perimeter) {
    return cased == 0 ? bore : 4 * area / perimeter;
}

/// The opening fraction of a calibrated gas-lift valve, from the state already primed on its choke, converted to
/// the units of the valve correlation.
[[nodiscard]] double calibratedValveOpening(const ChokeGas &choke) {
    return calibratedValveOpeningFraction(choke.pcalib * units::kPsiPerKgfPerCm2, choke.tcalib,
                                          (choke.presEstag - units::kAtmosphereInKgfPerCm2) * units::kPsiPerKgfPerCm2,
                                          (choke.presGarg - units::kAtmosphereInKgfPerCm2) * units::kPsiPerKgfPerCm2,
                                          choke.dextern, choke.areagarg, choke.areagarg / choke.areafole,
                                          units::celsiusToFahrenheit(choke.tempEstag));
}

/// Primes the choke from the two cells it spans and returns the mass source it
/// delivers, already scaled by the calibrated opening when the valve is one.
template <ValveFluids fluids>
[[nodiscard]] double primeValveMassSource(const GasLiftState &state, int valveIndex,
                                          const CelG &gasCell, const Cel &productionCell) {
    ChokeGas &choke = state.gasLiftChokes[valveIndex];
    const double recovery = fluids.recovery == kInjectedGas ? choke.frec : choke.frecliq;
    choke.presEstag = gasCell.pres;
    choke.presGarg = (productionCell.pres - choke.presEstag * recovery) / (1. - recovery);
    choke.tempEstag = gasCell.temp;
    double massSource = 0.;
    if (productionCell.pres < gasCell.pres)
        massSource = fluids.flow == kInjectedGas ? choke.massica()
                                                 : choke.massica(kInjectedCompletionFluid, state.input.salinDescarga);
    if (choke.tipo == kGasLiftCalibratedValve)
        massSource *= calibratedValveOpening(choke);
    return massSource;
}

/// Calls visit(productionCellIndex, gasCellIndex) for each pair of cells the column-annulus thermal coupling joins,
/// from tubingCouplingStart down to tubingCouplingEnd.
template <typename Visit>
void forEachCoupledCellPair(const GasLiftState &state, Visit &&visit) {
    for (int productionCellIndex = state.tubingCouplingStart; productionCellIndex >= state.tubingCouplingEnd;
         productionCellIndex--)
        visit(productionCellIndex, state.tubingCouplingStart + state.annulusCouplingStart - productionCellIndex);
}

/// Publishes the transport properties of the gas line's fluid onto the production
/// cell that faces it, as that cell's external fluid.
void publishExternalFluid(Cel &productionCell, const ExternalFluidProperties &fluid) {
    productionCell.calor.kextern1 = fluid.conductivity;
    productionCell.calor.cpextern1 = fluid.specificHeat;
    productionCell.calor.rhoextern1 = fluid.density;
    productionCell.calor.viscextern1 = fluid.viscosityCentipoise * units::kPascalSecondPerCentipoise;
}

}  // namespace

// Unloading the completion fluid through the gas-lift valves.

namespace {

/// Walks a pressure towards a target by kRelaxationDownFactor or kRelaxationUpFactor per step, without overshooting
/// it. The two factors are not each other's inverse.
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

/// Clears a production cell's mass sources.
void clearMassSources(Cel &cell) {
    cell.massfonteCH = 0;
    cell.fontemassCL = 0;
    cell.fontemassLL = 0;
    cell.fontemassGL = 0;
    cell.fontemassCR = 0;
    cell.fontemassLR = 0;
    cell.fontemassGR = 0;
}

/// Gives a production cell the phases of the unloading's initial state, completion fluid below the interface and
/// gas above it, with the copies the march and the pigs start from.
void setUnloadingPhase(Cel &cell, bool belowInterface) {
    cell.alf = belowInterface ? 0. : 1.;
    cell.bet = belowInterface ? 1. : 0.;
    cell.betI = cell.bet;
    cell.alfini = cell.alf;
    cell.betini = cell.bet;
    cell.alfPigD = cell.alf;
    cell.betPigD = cell.bet;
    cell.alfPigE = cell.alf;
    cell.betPigE = cell.bet;
}

/// Copies a production cell's phases to the face values its neighbours hold of it.
void copyPhaseToLeftNeighbour(const Cel &cell, Cel &left) {
    left.betR = cell.bet;
    left.betRini = cell.bet;
    left.alfR = cell.alf;
    left.alfRini = cell.alf;
}

void copyPhaseToRightNeighbour(const Cel &cell, Cel &right) {
    right.alfL = cell.alf;
    right.betL = cell.bet;
    right.betLI = cell.bet;
    right.alfLini = cell.alf;
    right.betLini = cell.bet;
}

/// The volume flow of a gas cell's source during the unloading: gas once the
/// interface is past the cell's middle, completion fluid before.
[[nodiscard]] double sourceVolumeFlowRate(const CelG &cell, double liquidDensity, double gasDensity) {
    return cell.razInter > kMidCellInterfaceRatio ? cell.massfonteCH / gasDensity : cell.massfonteCH / liquidDensity;
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
    for (int gasCellIndex = 1; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        double leftHalfLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        meanPressure -= leftDensity * kGravityUnloadingVariant * leftHalfLength *
                        sin(state.gasCells[gasCellIndex - 1].duto.teta) / units::kPascalPerKgfPerCm2;
        meanTemperature = state.gasCells[gasCellIndex].calor.Textern1;
        const double localDensity = gasCellIndex < state.interfaceCell
                                ? state.gasCells[gasCellIndex].flui.MasEspGas(meanPressure, meanTemperature)
                                : state.gasCells[gasCellIndex].MasEspFlu(meanPressure, meanTemperature);
        meanPressure -= localDensity * kGravityUnloadingVariant * halfLocalLength *
                        sin(state.gasCells[gasCellIndex].duto.teta) / units::kPascalPerKgfPerCm2;
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
        if (gasCellIndex < state.lastGasCell) {
            state.gasCells[gasCellIndex + 1].presL = meanPressure;
            state.gasCells[gasCellIndex + 1].tempL = meanTemperature;
            state.gasCells[gasCellIndex + 1].u1LL = state.gasCells[gasCellIndex].u1L;
            state.gasCells[gasCellIndex + 1].VGasL = 0;
        }
    }
}

void computeProductionUnloadingHydrostatics(const GasLiftState &state) {
    clearMassSources(state.cells[0]);
    const bool inletBelowInterface = 0 <= state.input.celdescargaP;
    setUnloadingPhase(state.cells[0], inletBelowInterface);
    copyPhaseToRightNeighbour(state.cells[0], state.cells[1]);
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

    for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) {
        const bool belowInterface = cellIndex <= state.input.celdescargaP;
        double leftHalfLength = 0.5 * state.cells[cellIndex].dxL;
        double localHalfLength = 0.5 * state.cells[cellIndex].dx;
        meanPressure -= leftDensity * kGravityUnloadingVariant * leftHalfLength *
                        sin(state.cells[cellIndex - 1].duto.teta) / units::kPascalPerKgfPerCm2;
        meanTemperature = state.cells[cellIndex].calor.Textern1;
        double interpolatedTemperature =
            (leftHalfLength * state.cells[cellIndex - 1].temp + localHalfLength * meanTemperature) /
            (leftHalfLength + localHalfLength);
        state.cells[cellIndex].presaux = meanPressure;
        state.cells[cellIndex - 1].presauxR = state.cells[cellIndex].presaux;
        state.cells[cellIndex].presauxL = state.cells[cellIndex - 1].presaux;
        if (belowInterface)
            localDensity = state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        else
            localDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        meanPressure -= localDensity * kGravityUnloadingVariant * localHalfLength *
                        sin(state.cells[cellIndex].duto.teta) / units::kPascalPerKgfPerCm2;
        leftDensity = localDensity;

        state.cells[cellIndex].presL = state.cells[cellIndex - 1].pres;
        state.cells[cellIndex].pres = meanPressure;
        state.cells[cellIndex - 1].presR = meanPressure;
        state.cells[cellIndex].presini = meanPressure;
        state.cells[cellIndex].presLini = state.cells[cellIndex - 1].presini;
        state.cells[cellIndex].presauxL = state.cells[cellIndex - 1].presaux;
        state.cells[cellIndex].tempL = state.cells[cellIndex - 1].temp;
        state.cells[cellIndex].temp = meanTemperature;
        state.cells[cellIndex - 1].tempR = meanTemperature;
        state.cells[cellIndex].tempini = meanTemperature;

        state.cells[cellIndex].FW = 0;
        state.cells[cellIndex].FWini = 0;
        state.cells[cellIndex].arranjo = 0;
        state.cells[cellIndex].QL = 0;
        state.cells[cellIndex - 1].QLR = state.cells[cellIndex].QL;
        state.cells[cellIndex].QG = 0;
        state.cells[cellIndex].rpL = localDensity;
        state.cells[cellIndex].rpC = localDensity;
        state.cells[cellIndex - 1].rpR = state.cells[cellIndex].rpC;
        state.cells[cellIndex].rcL = localDensity;
        state.cells[cellIndex].rcC = localDensity;
        state.cells[cellIndex - 1].rcR = state.cells[cellIndex].rcC;

        const double interpolatedDensity =
            belowInterface
                ? state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presaux, interpolatedTemperature)
                : state.cells[0].flui.MasEspGas(state.cells[cellIndex].presaux, interpolatedTemperature);
        state.cells[cellIndex].rpCi = interpolatedDensity;
        state.cells[cellIndex - 1].rpRi = interpolatedDensity;
        state.cells[cellIndex].rcCi = interpolatedDensity;
        state.cells[cellIndex - 1].rcRi = interpolatedDensity;

        clearMassSources(state.cells[cellIndex]);
        setUnloadingPhase(state.cells[cellIndex], belowInterface);
        copyPhaseToLeftNeighbour(state.cells[cellIndex], state.cells[cellIndex - 1]);

        if (cellIndex < state.lastCell) {
            state.cells[cellIndex + 1].tempL = meanTemperature;
            state.cells[cellIndex + 1].QLL = state.cells[cellIndex].QL;
            state.cells[cellIndex + 1].rpL = state.cells[cellIndex].rpC;
            state.cells[cellIndex + 1].rcL = state.cells[cellIndex].rcC;
            state.cells[cellIndex + 1].rpLi = interpolatedDensity;
            state.cells[cellIndex + 1].rcLi = interpolatedDensity;
            copyPhaseToRightNeighbour(state.cells[cellIndex], state.cells[cellIndex + 1]);
        }
    }
}

double calibratedValveOpeningFraction(double calibrationPressure, double calibrationTemperature,
                                      double valveOpeningPressure, double tubingPressure, double externalDiameter,
                                      double throatArea, double valveRatio, double bottomHoleTemperatureFahrenheit) {
    double domePressure = calibrationPressure * (1 - valveRatio);
    const double calibrationTemperatureRankine =
        units::celsiusToFahrenheit(calibrationTemperature) + units::kZeroFahrenheitInRankine;
    domePressure = (domePressure + units::kAtmosphereInPsi) *
                       (kBellowsReferenceTemperatureFahrenheit + units::kZeroFahrenheitInRankine) /
                       calibrationTemperatureRankine -
                   units::kAtmosphereInPsi;
    const double bellowsPressure =
        domePressure * (1 + kBellowsPressureChangePerFahrenheit *
                                (bottomHoleTemperatureFahrenheit - kBellowsReferenceTemperatureFahrenheit));
    double bellowsArea = throatArea / valveRatio;
    // The spring rate depends only on whether the valve is a large-bore one.
    const double externalDiameterInches = externalDiameter * units::kCentimetrePerMetre / units::kCentimetrePerInch;
    const double bellowsSpringRate =
        (externalDiameterInches > kLargeBoreValveInches ? kLargeBoreSpringRate : kSmallBoreSpringRate) * bellowsArea;
    double stemTravel = ((valveOpeningPressure - bellowsPressure) * bellowsArea -
                         (valveOpeningPressure - tubingPressure) * throatArea) /
                        bellowsSpringRate;

    double throatDiameter = sqrt(throatArea * 4. / std::numbers::pi);
    double throatRadius = throatDiameter / 2.0;
    double bellowsRadius = sqrt(bellowsArea / std::numbers::pi);
    // The seat offset: the distance from the stem axis to where the bellows
    // meets the seat.
    const double seatOffset = sqrt(bellowsRadius * bellowsRadius - throatRadius * throatRadius);
    double openingArea = std::numbers::pi * throatRadius * stemTravel * (stemTravel + 2.0 * seatOffset);
    openingArea =
        openingArea / sqrt((stemTravel + seatOffset) * (stemTravel + seatOffset) + throatRadius * throatRadius);
    if (openingArea > throatArea)
        openingArea = throatArea;

    return openingArea / throatArea;
}

double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate, int valveIndex, double factor,
                                   int sign) {
    int gasValveCell = state.gasValveCellIndices[valveIndex];
    state.gasLiftChokes[valveIndex].presEstag = state.gasCells[gasValveCell].pres;
    state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[gasValveCell].temp;
    double completionFluidDensity = state.gasCells[gasValveCell].MasEspFlu(state.gasLiftChokes[valveIndex].presEstag,
                                                                           state.gasLiftChokes[valveIndex].tempEstag);
    double massFlowRate = (factor * state.input.vazDescControl - maximumFlowRate) * completionFluidDensity;
    const double throatMassFlux = massFlowRate / state.gasLiftChokes[valveIndex].areagarg;
    const double pressureCorrection =
        throatMassFlux * throatMassFlux / (2. * completionFluidDensity * units::kPascalPerKgfPerCm2);
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
        double liquidDensity =
            completionFluidFraction * completionFluidDensity + (1. - completionFluidFraction) * productionLiquidDensity;
        double liquidViscosity =
            completionFluidFraction * completionFluidViscosity + (1. - completionFluidFraction) * oilViscosity;
        double mixtureDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;
        double mixtureViscosity = voidFraction * gasViscosity + (1. - voidFraction) * liquidViscosity;
        double mixtureVelocity = state.cells[cellIndex].QL / (flowArea * liquidDensity) +
                                 state.cells[cellIndex].QG / (flowArea * gasDensity);
        double frictionFactor = 0.;
        if (fabs(mixtureVelocity) > kStillMixtureVelocity) {
            const double reynolds =
                state.cells[cellIndex].Rey(characteristicDiameter(state.cells[cellIndex].duto.revest,
                                                                  state.cells[cellIndex].duto.a, flowArea, perimeter),
                                           mixtureVelocity, mixtureDensity, mixtureViscosity);
            frictionFactor =
                state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / state.cells[cellIndex].duto.a);
        }
        double wallShearStress = frictionFactor * mixtureDensity * mixtureVelocity * fabs(mixtureVelocity) / 2.;
        meanPressure -= (-kGravity * mixtureDensity * sin(state.cells[cellIndex].duto.teta) -
                         wallShearStress * perimeter / flowArea) *
                        cellLength / units::kPascalPerKgfPerCm2;
        if (state.cells[cellIndex].acsr.tipo == kAccessoryInflowPerformance) {
            double candidateMaximumPressure =
                state.cells[cellIndex].acsr.ipr.Pres - (meanPressure - state.separatorPressure);
            if (candidateMaximumPressure < state.input.presMaxDesc)
                state.input.presMaxDesc = candidateMaximumPressure;
        }
    }
    if (state.separatorPressure >= state.input.presMaxDesc * kReachedMaximumPressureFactor &&
        throatFlowRate > kUnloadingUpperBandEdge * state.input.vazDescControl && state.gasCells[0].VGasR > 0.) {
        const double pressureCorrection = limitedCorrection(
            unloadingPressureCorrection(state, throatFlowRate, valveIndex, kUnloadingUpperBandEdge, -1),
            state.injectionPressure, state.cells[0].dt);
        state.injectionPressure += pressureCorrection;
        if (state.injectionPressure < state.input.presMinDescG) {
            state.injectionPressure = state.input.presMinDescG;
        }

    } else if (throatFlowRate <= kUnloadingLowerBandEdge * state.input.vazDescControl &&
               state.separatorPressure <= state.input.presMinDesc * kReachedMinimumPressureFactor) {
        const double pressureCorrection = limitedCorrection(
            unloadingPressureCorrection(state, throatFlowRate, valveIndex, kUnloadingLowerBandEdge, 1),
            state.injectionPressure, state.cells[0].dt);
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
                    const double massFlowRate = primeValveMassSource<kUnloadingThroughValve>(
                        state, valveIndex, valveCell, state.cells[productionValveCell]);
                    double candidateFlowRate = massFlowRate / completionFluidDensity;
                    if (candidateFlowRate > maximumUnloadingFlowRate) {
                        maximumUnloadingFlowRate = candidateFlowRate;
                        selectedValveIndex = valveIndex;
                    }
                }
            }
            if (std::ssize(state.unloadingTimeStepSamples) > state.unloadingSampleCapacity ||
                state.unloadingAveragingInterval > state.unloadingAveragingWindow) {
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
            maximumUnloadingFlowRate =
                kCurrentFlowWeight * maximumUnloadingFlowRate +
                (1. - kCurrentFlowWeight) * state.meanUnloadingFlowRate / state.unloadingAveragingInterval;
            computeUnloadingValvePressure(state, maximumUnloadingFlowRate, selectedValveIndex);
            if (state.gasValveCellIndices[selectedValveIndex] >= state.interfaceCell) {
                if (state.cells[state.lastCell - 1].MC > 0.0) {
                    if (maximumUnloadingFlowRate > kUnloadingUpperBandEdge * state.input.vazDescControl) {
                        const double pressureCorrection = limitedCorrection(
                            unloadingPressureCorrection(state, maximumUnloadingFlowRate, selectedValveIndex,
                                                        kUnloadingUpperBandEdge, -1),
                            state.separatorPressure, state.cells[0].dt);
                        state.separatorPressure -= pressureCorrection;
                        if (state.separatorPressure > state.input.presMaxDesc)
                            state.separatorPressure = state.input.presMaxDesc;
                    } else if (maximumUnloadingFlowRate <= kUnloadingLowerBandEdge * state.input.vazDescControl) {
                        const double pressureCorrection = limitedCorrection(
                            unloadingPressureCorrection(state, maximumUnloadingFlowRate, selectedValveIndex,
                                                        kUnloadingLowerBandEdge, 1),
                            state.separatorPressure, state.cells[0].dt);
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
            if (state.globals->lixo5 < kSeparatorRampLatencyFraction * state.input.tempoLatenciaDesc)
                state.separatorPressure =
                    initialUnloadingPressure -
                    (initialUnloadingPressure - state.input.presMinDesc) * state.cells[0].dt /
                        (kSeparatorRampLatencyFraction * state.input.tempoLatenciaDesc - state.globals->lixo5);
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
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        double leftHalfLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        double localGasShare = 1.;
        if (state.gasCells[gasCellIndex].razInter <= kMidCellInterfaceRatio)
            localGasShare = 2 * state.gasCells[gasCellIndex].razInter;
        double leftGasShare = 0.;
        if (state.gasCells[gasCellIndex - 1].razInter >= kMidCellInterfaceRatio)
            leftGasShare = 2 * (state.gasCells[gasCellIndex - 1].razInter - kMidCellInterfaceRatio);
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
        double leftVelocity =
            state.gasCells[gasCellIndex].VGasL / (carryingDensity * state.gasCells[gasCellIndex - 1].duto.area);
        double localVelocity =
            state.gasCells[gasCellIndex].VGasL / (carryingDensity * state.gasCells[gasCellIndex].duto.area);
        const double leftDiameter = characteristicDiameter(
            state.gasCells[gasCellIndex - 1].duto.revest, state.gasCells[gasCellIndex - 1].duto.a,
            state.gasCells[gasCellIndex - 1].duto.area, state.gasCells[gasCellIndex - 1].duto.peri);
        const double localDiameter =
            characteristicDiameter(state.gasCells[gasCellIndex].duto.revest, state.gasCells[gasCellIndex].duto.a,
                                   state.gasCells[gasCellIndex].duto.area, state.gasCells[gasCellIndex].duto.peri);
        const double leftLiquidReynolds =
            state.gasCells[gasCellIndex - 1].Rey(leftDiameter, leftVelocity, liquidDensity, liquidViscosity);
        const double localLiquidReynolds =
            state.gasCells[gasCellIndex].Rey(localDiameter, localVelocity, liquidDensity, liquidViscosity);
        const double leftGasReynolds =
            state.gasCells[gasCellIndex - 1].Rey(leftDiameter, leftVelocity, gasDensity, gasViscosity);
        const double localGasReynolds =
            state.gasCells[gasCellIndex].Rey(localDiameter, localVelocity, gasDensity, gasViscosity);
        double leftLiquidFriction =
            state.gasCells[gasCellIndex - 1].fric(leftLiquidReynolds, state.gasCells[gasCellIndex - 1].duto.rug /
                                                                          state.gasCells[gasCellIndex - 1].duto.a) *
            leftLiquidLength;
        double localLiquidFriction =
            state.gasCells[gasCellIndex].fric(localLiquidReynolds, state.gasCells[gasCellIndex].duto.rug /
                                                                       state.gasCells[gasCellIndex].duto.a) *
            localLiquidLength;
        double leftLiquidHead =
            (kGravity * sin(state.gasCells[gasCellIndex - 1].duto.teta) * liquidDensity) * leftLiquidLength;
        double localLiquidHead =
            (kGravity * sin(state.gasCells[gasCellIndex].duto.teta) * liquidDensity) * localLiquidLength;
        double leftGasFriction =
            state.gasCells[gasCellIndex - 1].fric(leftGasReynolds, state.gasCells[gasCellIndex - 1].duto.rug /
                                                                       state.gasCells[gasCellIndex - 1].duto.a) *
            leftGasLength;
        double localGasFriction =
            state.gasCells[gasCellIndex].fric(localGasReynolds, state.gasCells[gasCellIndex].duto.rug /
                                                                    state.gasCells[gasCellIndex].duto.a) *
            localGasLength;
        double leftGasHead = (kGravity * sin(state.gasCells[gasCellIndex - 1].duto.teta) * gasDensity) * leftGasLength;
        double localGasHead = (kGravity * sin(state.gasCells[gasCellIndex].duto.teta) * gasDensity) * localGasLength;
        state.gasCells[gasCellIndex].pres =
            state.gasCells[gasCellIndex - 1].pres +
            (-0.5 * (leftLiquidFriction * liquidDensity + leftGasFriction * gasDensity) * leftVelocity *
                 fabs(leftVelocity) * state.gasCells[gasCellIndex - 1].duto.peri /
                 state.gasCells[gasCellIndex - 1].duto.area -
             0.5 * (localLiquidFriction * liquidDensity + localGasFriction * gasDensity) * localVelocity *
                 fabs(localVelocity) * state.gasCells[gasCellIndex].duto.peri / state.gasCells[gasCellIndex].duto.area -
             leftLiquidHead - localLiquidHead - leftGasHead - localGasHead) /
                units::kPascalPerKgfPerCm2;
        state.gasCells[gasCellIndex].presL = state.gasCells[gasCellIndex - 1].pres;
        state.gasCells[gasCellIndex - 1].presR = state.gasCells[gasCellIndex].pres;

        state.temperatureUpdater.dischargeTemperature(gasCellIndex);

        state.gasCells[gasCellIndex].u1L =
            ((1. - state.gasCells[gasCellIndex].razInter) *
                 state.gasCells[gasCellIndex].MasEspFlu(state.gasCells[gasCellIndex].pres,
                                                        state.gasCells[gasCellIndex].temp) +
             state.gasCells[gasCellIndex].razInter *
                 state.gasCells[gasCellIndex].flui.MasEspGas(leftCellPressure, leftTemperature)) *
            state.gasCells[gasCellIndex].duto.area;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        state.gasCells[gasCellIndex].u1LL = state.gasCells[gasCellIndex - 1].u1L;
    }

    double totalSourceVolumeFlow = 0.;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        const double cellTemperature = state.gasCells[gasCellIndex].temp;
        const double cellPressure = state.gasCells[gasCellIndex].pres;
        const double liquidDensity = state.gasCells[gasCellIndex].MasEspFlu(cellPressure, cellTemperature);
        const double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(cellPressure, cellTemperature);
        totalSourceVolumeFlow += sourceVolumeFlowRate(state.gasCells[gasCellIndex], liquidDensity, gasDensity);
    }
    const double interfaceTemperature = state.gasCells[state.interfaceCell].temp;
    const double interfacePressure = state.gasCells[state.interfaceCell].pres;
    const double interfaceGasDensity =
        state.gasCells[state.interfaceCell].flui.MasEspGas(interfacePressure, interfaceTemperature);
    state.gasCells[state.interfaceCell].VGasL = totalSourceVolumeFlow * interfaceGasDensity;
    state.gasCells[state.interfaceCell - 1].VGasR = state.gasCells[state.interfaceCell].VGasL;
    state.gasCells[state.interfaceCell - 2].VGasRR = state.gasCells[state.interfaceCell].VGasL;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        const double cellTemperature = state.gasCells[gasCellIndex].temp;
        const double cellPressure = state.gasCells[gasCellIndex].pres;
        const double liquidDensity = state.gasCells[gasCellIndex].MasEspFlu(cellPressure, cellTemperature);
        const double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(cellPressure, cellTemperature);
        totalSourceVolumeFlow -= sourceVolumeFlowRate(state.gasCells[gasCellIndex], liquidDensity, gasDensity);
        if (gasCellIndex < state.lastGasCell) {
            state.gasCells[gasCellIndex + 1].VGasL = totalSourceVolumeFlow * liquidDensity;
            state.gasCells[gasCellIndex].VGasR = state.gasCells[gasCellIndex + 1].VGasL;
        } else
            state.gasCells[gasCellIndex].VGasR = 0.;
        state.gasCells[gasCellIndex - 1].VGasRR = state.gasCells[gasCellIndex].VGasL;
    }

    state.interfaceVelocity =
        state.gasCells[state.interfaceCell + 1].VGasL /
        (state.gasCells[state.interfaceCell + 1].MasEspFlu(state.gasCells[state.interfaceCell + 1].pres,
                                                           state.gasCells[state.interfaceCell + 1].temp) *
         state.gasCells[state.interfaceCell].duto.area);
}

void advanceInterface(const GasLiftState &state) {
    state.gasCells[state.interfaceCell].razInter =
        (state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 +
         state.interfaceVelocity * state.timeStep) /
        state.gasCells[state.interfaceCell].dx0;
    if (state.interfaceCell == (state.lastGasCell - 1) &&
        state.gasCells[state.interfaceCell].razInter >= kInterfaceHandoverRatio) {
        double meanPressure = state.gasCells[state.interfaceCell].pres;
        double meanTemperature = state.gasCells[state.interfaceCell].temp;
        double interfaceFlowArea = state.gasCells[state.interfaceCell].duto.area;
        double gasDensity = state.gasCells[state.interfaceCell].flui.MasEspGas(meanPressure, meanTemperature);
        double liquidDensity = state.gasCells[state.interfaceCell].MasEspFlu(meanPressure, meanTemperature);
        state.gasCells[state.interfaceCell].VGasR =
            state.gasCells[state.interfaceCell].VGasR * gasDensity / liquidDensity;
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
        state.gasCells[state.interfaceCell].VGasR =
            0 * state.gasCells[state.interfaceCell].VGasR * gasDensity / liquidDensity;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * gasDensity;
        state.gasCells[state.interfaceCell].u1R = interfaceFlowArea * gasDensity;
        state.gasCells[state.interfaceCell - 1].u1R = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell - 1].VGasRR = state.gasCells[state.interfaceCell].VGasR;

        state.interfaceCell = kNoInterfaceCell;
        state.input.descarga = 0;
    }

    if (state.input.descarga == 1) {
        const double tiny = state.globals->localtiny;
        CelG &cell = state.gasCells[state.interfaceCell];
        if (cell.razInter <= tiny && cell.razInter >= -tiny)
            cell.razInter = 0;
        else if (cell.razInter < -tiny) {
            const double interfaceCrossingTimeStep = -cell.razInterIni * cell.dx0 / state.interfaceVelocity;
            if (interfaceCrossingTimeStep > tiny)
                state.interfaceTimeStep = interfaceCrossingTimeStep;
            cell.razInter = 0.;
        } else if (cell.razInter >= 1. - tiny && cell.razInter <= 1. + tiny) {
            cell.razInter = 1.;
        } else if (cell.razInter > 1. + tiny) {
            const double interfaceCrossingTimeStep = (1. - cell.razInterIni) * cell.dx0 / state.interfaceVelocity;
            if (interfaceCrossingTimeStep > tiny)
                state.interfaceTimeStep = interfaceCrossingTimeStep;
            cell.razInter = 1.;
        }

        if (fabs(cell.razInter) < tiny && state.interfaceVelocity < -fabs(tiny)) {
            cell.razInter = 0.;
            if (state.interfaceCell > 0) {
                state.interfaceCell--;
                state.gasCells[state.interfaceCell].razInter = 1.;
            }
        } else if (fabs(cell.razInter - 1.) < tiny && state.interfaceVelocity > tiny) {
            cell.razInter = 1.;
            if (state.interfaceCell < state.lastGasCell) {
                state.interfaceCell++;
                state.gasCells[state.interfaceCell].razInter = 0.;
            }
        }
    }
}

// The gas line in the transient.

namespace {

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

}  // namespace

void solveGasLine(const GasLiftState &state) {
    if (state.input.lingas <= 0)
        return;

    updateTransientGasValves(state);
    for (int valveIndex = 0; valveIndex < state.input.nvalvgas; valveIndex++) {
        int productionValveCell = state.productionValveCellIndices[valveIndex];
        int gasValveCell = state.gasValveCellIndices[valveIndex];
        const bool valveDischargesGas =
            gasValveCell < state.interfaceCell ||
            (gasValveCell == state.interfaceCell && state.gasCells[gasValveCell].razInter > kMidCellInterfaceRatio);
        if (valveDischargesGas) {
            state.cells[productionValveCell].acsr.injg.QGas =
                state.gasCells[gasValveCell].massfonteCH * units::kSecondsPerDay /
                state.gasCells[gasValveCell].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius);
            state.cells[productionValveCell].acsr.injg.tipoflu = kInjectedGas;
            const double gasTemperature =
                state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg
                    ? state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex)
                    : state.cells[productionValveCell].temp;
            state.cells[productionValveCell].acsr.injg.temp = gasTemperature;
            if (state.cells[productionValveCell].acsr.injg.temp < kMinimumTemperatureCelsius)
                state.cells[productionValveCell].acsr.injg.temp = kMinimumTemperatureCelsius;
        } else {
            state.cells[productionValveCell].acsr.injg.QGas = state.gasCells[gasValveCell].massfonteCH;
            state.cells[productionValveCell].acsr.injg.tipoflu = kInjectedCompletionFluid;
            state.cells[productionValveCell].acsr.injg.temp = state.gasLiftChokes[valveIndex].tempEstag;
        }
        state.gasCells[gasValveCell].pEstag = state.gasLiftChokes[valveIndex].presEstag;
        state.gasCells[gasValveCell].tEstag = state.gasLiftChokes[valveIndex].tempEstag;
        state.gasCells[gasValveCell].pGarg = state.gasLiftChokes[valveIndex].presGarg;
        state.gasCells[gasValveCell].tGarg = state.gasLiftChokes[valveIndex].tempGarg;
        state.gasCells[gasValveCell].qGarg = state.gasLiftChokes[valveIndex].qGarg;
        state.gasCells[gasValveCell].areaGarg = state.gasLiftChokes[valveIndex].areagarg;
    }
    advanceGasSubStep(state);
}

void updateTransientGasValves(const GasLiftState &state) {
    int valveCount = state.input.nvalvgas;
    int master2CellIndex = state.input.master2.posic;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.lastGasCell; gasCellIndex++)
        state.gasCells[gasCellIndex].massfonteCH = 0.;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        const int gasValveCell = state.gasValveCellIndices[valveIndex];
        CelG &valveCell = state.gasCells[gasValveCell];
        const Cel &productionCell = state.cells[state.productionValveCellIndices[valveIndex]];
        // Behind the interface the valve discharges completion fluid, ahead of it gas.
        valveCell.massfonteCH =
            gasValveCell >= state.interfaceCell && valveCell.razInter < kMidCellInterfaceRatio
                ? primeValveMassSource<kCompletionFluidThroughValve>(state, valveIndex, valveCell, productionCell)
                : primeValveMassSource<kGasThroughValve>(state, valveIndex, valveCell, productionCell);
    }

    double master2FlowArea = state.gasCells[master2CellIndex].duto.area;
    if (master2FlowArea > state.gasCells[master2CellIndex].dutoR.area)
        master2FlowArea = state.gasCells[master2CellIndex].dutoR.area;
    if (state.gasCells[master2CellIndex].chkcell.areagarg <= kMaster2ThrottlingAreaFraction * master2FlowArea &&
        state.gasCells[master2CellIndex].chkcell.areagarg >= kMaster2ClosedAreaFraction * master2FlowArea) {
        state.gasCells[master2CellIndex].chkcell.tempEstag = state.gasCells[master2CellIndex].temp;
        double upstreamGasDensity = state.gasCells[master2CellIndex].flui.MasEspGas(
            state.gasCells[master2CellIndex].pres, state.gasCells[master2CellIndex].temp);
        double downstreamGasDensity = state.gasCells[master2CellIndex + 1].flui.MasEspGas(
            state.gasCells[master2CellIndex + 1].pres, state.gasCells[master2CellIndex + 1].temp);
        double upstreamHead = -0.5 * upstreamGasDensity * units::kStandardGravity *
                              state.gasCells[master2CellIndex].dx0 * sin(state.gasCells[master2CellIndex].duto.teta) /
                              units::kPascalPerKgfPerCm2;
        double downstreamHead = 0.5 * downstreamGasDensity * units::kStandardGravity *
                                state.gasCells[master2CellIndex + 1].dx0 *
                                sin(state.gasCells[master2CellIndex + 1].duto.teta) / units::kPascalPerKgfPerCm2;
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
    for (int gasCellIndex = 0; gasCellIndex <= state.lastGasCell; gasCellIndex++)
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

    for (int gasCellIndex = 0; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        state.gasCells[gasCellIndex].dt = state.timeStep;
    }

    // The temperature, density and flow loops below stop before this cell: past the last one without an interface,
    // at the interface cell's left neighbour otherwise.
    const int gasLoopEnd = state.interfaceCell <= state.lastGasCell ? state.interfaceCell - 1 : state.lastGasCell + 1;
    for (int gasCellIndex = 0; gasCellIndex <= state.lastGasCell; gasCellIndex++)
        state.gasCells[gasCellIndex].dTdt = 0;

    int maximumThermalCouplingCycles = state.input.cicloAcopTerm;
    for (int thermalCouplingCycle = 0; thermalCouplingCycle <= maximumThermalCouplingCycles; thermalCouplingCycle++) {
        assembleAndSolveGasSystem(state);
        updateGasLine(state);
        state.gasCells[0].temp = state.injectionTemperature;
        state.gasCells[0].dTdt = (state.gasCells[0].temp - state.gasCells[0].tempini) / state.timeStep;
        if (state.gasCells[state.lastGasCell].VGasR < 0.)
            state.gasCells[state.lastGasCell].temp = kReverseFlowOutletTemperatureCelsius;
#pragma omp parallel for num_threads(state.globals->ntrd)
        for (int gasCellIndex = 1; gasCellIndex < gasLoopEnd; gasCellIndex++) {
            state.temperatureUpdater.computeGasTemperature(gasCellIndex, state.gasCells[gasCellIndex - 1].tempini);
            state.gasCells[gasCellIndex].dTdt =
                (state.gasCells[gasCellIndex].temp - state.gasCells[gasCellIndex].tempini) / state.timeStep;
        }
        if (thermalCouplingCycle < maximumThermalCouplingCycles)
            for (int renewedCellIndex = 0; renewedCellIndex <= state.lastGasCell; renewedCellIndex++)
                state.gasCells[renewedCellIndex].FeiticoDoTempo();
    }

#pragma omp parallel for num_threads(state.globals->ntrd)
    for (int gasCellIndex = 0; gasCellIndex < gasLoopEnd; gasCellIndex++)
        state.gasCells[gasCellIndex].rg = state.gasCells[gasCellIndex].flui.MasEspGas(
            state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp);

    for (int gasCellIndex = 0; gasCellIndex < gasLoopEnd; gasCellIndex++) {
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].rg * state.gasCells[gasCellIndex].duto.area;
        // The inlet feeds itself; every other cell is fed by the one upstream.
        state.gasCells[gasCellIndex].u1LL = gasCellIndex == 0
                                                ? state.gasCells[gasCellIndex].u1L
                                                : state.gasCells[gasCellIndex - 1].u1L;
        if (gasCellIndex > 0) {
            state.gasCells[gasCellIndex - 1].rgR = state.gasCells[gasCellIndex].rg;
            state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        }
        if (gasCellIndex == state.lastGasCell) {
            state.gasCells[gasCellIndex].u1R = state.gasCells[gasCellIndex].u1L;
            state.gasCells[state.lastGasCell].rgR = state.gasCells[state.lastGasCell].rg;
        }
    }

    for (int gasCellIndex = 0; gasCellIndex <= state.lastGasCell; gasCellIndex++)
        state.gasCells[gasCellIndex].presini = state.gasCells[gasCellIndex].pres;

    if (state.input.descarga == 1)
        solveUnloading(state);
}

void assembleAndSolveGasSystem(const GasLiftState &state) {
    double chokeOpeningFraction = withGasInletCondition(
        state.gasCells[0].tipoCC,
        [&](auto condition) { return injectionChokeOpening(state, condition); });
#pragma omp parallel for num_threads(state.globals->ntrd)
    for (int gasCellIndex = 0; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        state.gasCells[gasCellIndex].GeraLocal(state.lastGasCell, state.injectionPressure, state.injectionTemperature,
                                               chokeOpeningFraction);
        const int firstRow = kGasLineUnknownsPerCell * gasCellIndex;
        for (int bandIndex = 0; bandIndex < kGasLineStencilCoefficients; bandIndex++) {
            state.gasSystemMatrix[firstRow][bandIndex - kGasLineUnknownsPerCell] =
                state.gasCells[gasCellIndex].local[0][bandIndex];
            state.gasSystemMatrix[firstRow + 1][bandIndex - kGasLineUnknownsPerCell - 1] =
                state.gasCells[gasCellIndex].local[1][bandIndex];
            state.gasSystemMatrix[firstRow + 2][bandIndex - kGasLineUnknownsPerCell - 2] =
                state.gasCells[gasCellIndex].local[2][bandIndex];
        }
        state.gasRightHandSide[firstRow] = state.gasCells[gasCellIndex].TL[0];
        state.gasRightHandSide[firstRow + 1] = state.gasCells[gasCellIndex].TL[1];
        state.gasRightHandSide[firstRow + 2] = state.gasCells[gasCellIndex].TL[2];
    }
    state.gasSystemMatrix.GaussElimPP(state.gasRightHandSide);
}

void updateGasLine(const GasLiftState &state) {
    // The unknowns of a gas cell in the solved system: its pressure, then the gas mass flow through its right face.
    const auto pressureOf = [&](int cellIndex) { return state.gasRightHandSide[kGasLineUnknownsPerCell * cellIndex]; };
    const auto rightFaceFlowOf = [&](int cellIndex) {
        return state.gasRightHandSide[kGasLineUnknownsPerCell * cellIndex + 1];
    };
    for (int gasCellIndex = 0; gasCellIndex <= state.lastGasCell; gasCellIndex++) {
        CelG &cell = state.gasCells[gasCellIndex];
        // The first cell counts as the first even when it is also the last.
        const bool first = gasCellIndex == 0;
        const bool last = !first && gasCellIndex == state.lastGasCell;
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

void connectTubing(const GasLiftState &state) {
    forEachCoupledCellPair(state, [&](int productionCellIndex, int gasCellIndex) {
        state.cells[productionCellIndex].calor.Textern2 = state.gasCells[gasCellIndex].calor.Tcamada[0][0];
        state.cells[productionCellIndex].calor.betext = state.gasCells[gasCellIndex].calor.betint;
        int layerIndex = state.cells[productionCellIndex].calor.geom.ncamadas - 1;
        int discretizationIndex = state.cells[productionCellIndex].calor.ncamada[layerIndex] - 1;
        state.gasCells[gasCellIndex].calor.Tint2 =
            state.cells[productionCellIndex].calor.Tcamada[layerIndex][discretizationIndex];
        state.cells[productionCellIndex].calor.colunaDia = state.gasCells[gasCellIndex].duto.dia;
        state.cells[productionCellIndex].calor.geom.b = state.gasCells[gasCellIndex].calor.geom.a;

        state.cells[productionCellIndex].calor.Textern1 = state.gasCells[gasCellIndex].temp;
        state.cells[productionCellIndex].calor.Vextern1 =
            state.gasCells[gasCellIndex].VGasR / state.gasCells[gasCellIndex].u1L;
        // Ahead of the interface the annulus carries gas; behind it, completion fluid.
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
    });
}

// The gas line in the steady state.

namespace {

/// The flow area of a bore of this diameter.
[[nodiscard]] double boreArea(double diameter) {
    return 0.25 * std::numbers::pi * diameter * diameter;
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
    if (cellIndex < state.lastGasCell)
        state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
}

/// One valve's gas mass source at steady state. Under a flow-rate condition,
/// an equal share of the injection flow rate; under a pressure condition, what
/// the valve's choke passes from the casing pressure, if the casing is above
/// the production line.
void setSteadyValveGasMassSource(const GasLiftState &state, int valveIndex, int valveCount,
                                 double, InjectionFlowRateCondition) {
    double valveFlowRate = (state.input.gasinj.vazgas[0] *
                            state.gasCells[0].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius) /
                            units::kSecondsPerDay) /
                           valveCount;
    state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = valveFlowRate;
}

void setSteadyValveGasMassSource(const GasLiftState &state, int valveIndex, int,
                                 double casingTemperature, InjectionPressureCondition) {
    int gasValveCell = state.gasValveCellIndices[valveIndex];
    double casingPressure = state.gasCells[gasValveCell].pres;
    state.gasLiftChokes[valveIndex].presEstag = casingPressure;
    state.gasLiftChokes[valveIndex].presGarg = state.cells[state.productionValveCellIndices[valveIndex]].pres;
    state.gasLiftChokes[valveIndex].tempEstag = casingTemperature;
    state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH =
        state.cells[state.productionValveCellIndices[valveIndex]].pres < casingPressure
            ? state.gasLiftChokes[valveIndex].massica()
            : 0.;
}

/// The steady gas temperature of a cell outside the column-annulus coupling, marched from its left neighbour.
void marchSteadyGasTemperature(const GasLiftState &state, int cellIndex) {
    double dx = state.gasCells[cellIndex].dx0;
    double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
    double meanTemperatureGradient = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].dTdLCor +
                                      state.gasCells[cellIndex - 1].dx0 * state.gasCells[cellIndex - 1].dTdLCor) /
                                     (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
    double area = state.gasCells[cellIndex].duto.area;
    double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres,
                                                                     state.gasCells[cellIndex - 1].temp);
    double meanSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
    double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres,
                                                                         state.gasCells[cellIndex - 1].temp);
    double gasJouleThomson =
        state.gasCells[cellIndex - 1].flui.JTG(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
    double hydrostaticTerm =
        (gasDensity * meanSuperficialGasVelocity) * area * kGravity * sin(state.gasCells[cellIndex].duto.teta);

    state.gasCells[cellIndex - 1].calor.Tint = state.gasCells[cellIndex - 1].temp;
    state.gasCells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity;
    state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(
        state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
    state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
    state.gasCells[cellIndex - 1].calor.rhoint = gasDensity;
    state.gasCells[cellIndex - 1].calor.viscint =
        state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres,
                                                   state.gasCells[cellIndex - 1].temp) *
        units::kPascalSecondPerCentipoise;
    state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
    double leftCellHeatFlux = state.gasCells[cellIndex - 1].fluxcal;
    double tubingHeatFlux = 0.;

    const double lengthRatio = dx / (dx + state.gasCells[cellIndex - 1].dx0);
    double temperatureMarchCoefficient = gasDensity * meanSuperficialGasVelocity * gasSpecificHeat * area;
    double pressureMarchCoefficient = gasDensity * meanSuperficialGasVelocity * gasJouleThomson * area;
    double pressureGradient =
        2. *
        (state.gasCells[cellIndex].pres -
         ((1 - lengthRatio) * state.gasCells[cellIndex - 1].pres + lengthRatio * state.gasCells[cellIndex].pres)) *
        units::kPascalPerKgfPerCm2 / dx;

    double mixtureDensity = gasDensity;
    double leftCellLeftFaceSuperficialGasVelocity =
        state.gasCells[cellIndex - 1].VGasL / state.gasCells[cellIndex - 1].u1L;
    const double superficialVelocityGradient =
        (meanSuperficialGasVelocity - leftCellLeftFaceSuperficialGasVelocity) / dx;

    const double kineticTerm = state.input.nCompTotalUnidadesG / dx < kKineticTermMaximumLengthRatio
                                   ? mixtureDensity * area * meanSuperficialGasVelocity * meanSuperficialGasVelocity *
                                         superficialVelocityGradient
                                   : 0.;

    double gasMassSource = 0.;
    double liquidMassSource = 0.;

    if (meanSuperficialGasVelocity > kMinimumMarchGasVelocity &&
        (state.gasCells[cellIndex].duto.a / state.gasCells[cellIndex - 1].duto.a > kMinimumMarchDiameterRatio &&
         state.gasCells[cellIndex - 1].duto.a / state.gasCells[cellIndex].duto.a > kMinimumMarchDiameterRatio)) {
        double workTerm = meanTemperatureGradient *
                          (pressureMarchCoefficient * pressureGradient - kineticTerm - hydrostaticTerm +
                           liquidMassSource + gasMassSource) /
                          temperatureMarchCoefficient;
        double heatTerm = meanTemperatureGradient * (leftCellHeatFlux + tubingHeatFlux) / temperatureMarchCoefficient;
        int stepCount;
        double stepLength;
        if (meanCellLength / state.gasCells[cellIndex - 1].calor.resGlob < kThermalSubstepLength) {
            stepCount = 0;
            stepLength = meanCellLength;
        } else {
            stepCount = static_cast<int>(
                2 * (meanCellLength / state.gasCells[cellIndex - 1].calor.resGlob) / kThermalSubstepLength + 1);
            stepLength = meanCellLength / stepCount;
        }
        double stepTemperature = stepLength * (state.gasCells[cellIndex - 1].temp / stepLength + workTerm + heatTerm);
        for (int stepIndex = 1; stepIndex < stepCount; stepIndex++) {
            state.gasCells[cellIndex - 1].calor.Tint = stepTemperature;
            state.gasCells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity;
            state.gasCells[cellIndex - 1].calor.kint =
                state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, stepTemperature);
            state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex - 1].calor.rhoint = gasDensity;
            state.gasCells[cellIndex - 1].calor.viscint =
                state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, stepTemperature) *
                units::kPascalSecondPerCentipoise;
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
}

/// The steady gas temperature of a cell inside the column-annulus coupling, and the properties its heat
/// transfer is computed with.
void coupledSteadyGasTemperature(const GasLiftState &state, int cellIndex) {
    if (state.annulusThermalCoupling) {
        int connectionCellIndex = cellIndex - state.annulusCouplingStart;
        int connectionIndex = state.tubingCouplingStart - connectionCellIndex;
        int secondConnectionIndex = connectionIndex;
        if (secondConnectionIndex == state.lastCell)
            secondConnectionIndex -= 1;

        state.gasCells[cellIndex].temp =
            state.cells[secondConnectionIndex].calor.resGlob * state.cells[connectionIndex].calor.fluxFim +
            state.cells[connectionIndex].temp;

        if (state.input.correcaoContracorPerm == 1) {
            double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres,
                                                                                 state.gasCells[cellIndex - 1].temp);
            double inverseHeatCapacity =
                1. / (state.cells[connectionIndex].calor.cpint * state.cells[connectionIndex].MC) +
                1. / (gasSpecificHeat * state.gasCells[cellIndex - 1].VGasR);
            double logTemperature =
                -exp(-(1. / state.cells[secondConnectionIndex].calor.resGlob) *
                     (1 + 0 * state.cells[connectionIndex].duto.peri) * meanCellLength * inverseHeatCapacity) *
                    (state.cells[secondConnectionIndex + 1].temp - state.gasCells[cellIndex - 1].temp) +
                state.cells[secondConnectionIndex].temp;

            if (logTemperature < state.gasCells[cellIndex].temp)
                state.gasCells[cellIndex].temp = logTemperature;
        }
    }

    clampAndPropagateTemperature(state, cellIndex);

    double meanSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
    double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres,
                                                                     state.gasCells[cellIndex - 1].temp);
    double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres,
                                                                         state.gasCells[cellIndex - 1].temp);

    state.gasCells[cellIndex].calor.Tint = state.gasCells[cellIndex].temp;
    state.gasCells[cellIndex].calor.Vint = meanSuperficialGasVelocity;
    state.gasCells[cellIndex].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(
        state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
    state.gasCells[cellIndex].calor.cpint = gasSpecificHeat;
    state.gasCells[cellIndex].calor.rhoint = gasDensity;
    state.gasCells[cellIndex].calor.viscint =
        state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres,
                                                   state.gasCells[cellIndex - 1].temp) *
        units::kPascalSecondPerCentipoise;
    state.gasCells[cellIndex].fluxcal = state.gasCells[cellIndex].calor.transperm();
}

}  // namespace

void connectTubingSteady(const GasLiftState &state) {
    forEachCoupledCellPair(state, [&](int productionCellIndex, int gasCellIndex) {
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
    });
}

void initializeTubingConnectionSteady(const GasLiftState &state) {
    forEachCoupledCellPair(state, [&](int productionCellIndex, int gasCellIndex) {
        state.cells[productionCellIndex].calor.Textern2 = state.gasCells[gasCellIndex].calor.Textern2;
        state.cells[productionCellIndex].calor.colunaDia = state.gasCells[gasCellIndex].duto.dia;

        const CelG &gasCell = state.gasCells[gasCellIndex];
        const double productionPressure = state.cells[productionCellIndex].pres;
        const double externalTemperature = state.input.celg[gasCellIndex].textern;
        const double gasSpecificHeat = gasCell.flui.CalorGas(productionPressure, externalTemperature);
        const double gasDensity = gasCell.flui.MasEspGas(productionPressure, externalTemperature);
        state.cells[productionCellIndex].calor.Textern1 = gasCell.calor.Textern1;
        state.cells[productionCellIndex].calor.Vextern1 = 1.;
        publishExternalFluid(state.cells[productionCellIndex],
                             {gasCell.flui.CondGas(productionPressure, externalTemperature), gasSpecificHeat,
                              gasDensity, kInitialAnnulusGasViscosityCentipoise});
    });
}

double steadyGasPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.gasCells[cellIndex].dx0;
    double diameter = state.gasCells[cellIndex].duto.a;
    double area = boreArea(diameter);
    double perimeter = state.gasCells[cellIndex].duto.peri;
    double gasDensity =
        state.gasCells[cellIndex].flui.MasEspGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double leftFaceGasMassFlowRate = state.gasCells[cellIndex - 1].VGasR;

    double gasVelocity = leftFaceGasMassFlowRate / state.gasCells[cellIndex].u1L;

    double gasViscosity =
        state.gasCells[cellIndex].flui.ViscGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double reynolds =
        state.gasCells[cellIndex].Rey(characteristicDiameter(state.gasCells[cellIndex].duto.revest,
                                                             state.gasCells[cellIndex].duto.a, area, perimeter),
                                      gasVelocity, gasDensity, gasViscosity);
    double frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    double frictionGradient =
        state.gasCells[cellIndex].dPdLFric *
        (0.5 * frictionFactor * gasDensity * (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area);
    double hydrostaticGradient =
        state.gasCells[cellIndex].dPdLHidro * (kGravity * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);

    double pressureDrop = (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    double meanPressure = state.gasCells[cellIndex].pres + pressureDrop;

    double meanTemperature = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].temp +
                              state.gasCells[cellIndex].dxL * state.gasCells[cellIndex].tempL) /
                             (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex].dxL);

    dx = 0.5 * state.gasCells[cellIndex].dxL;
    diameter = state.gasCells[cellIndex].dutoL.a;
    area = boreArea(diameter);
    perimeter = state.gasCells[cellIndex].dutoL.peri;
    gasDensity = state.gasCells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);

    gasVelocity = leftFaceGasMassFlowRate / (gasDensity * area);

    gasViscosity = state.gasCells[cellIndex - 1].flui.ViscGas(meanPressure, meanTemperature);

    reynolds = state.gasCells[cellIndex].Rey(characteristicDiameter(state.gasCells[cellIndex].dutoL.revest,
                                                                    state.gasCells[cellIndex].dutoL.a, area, perimeter),
                                             gasVelocity, gasDensity, gasViscosity);
    frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
    frictionGradient = state.gasCells[cellIndex - 1].dPdLFric *
                       (0.5 * frictionFactor * gasDensity * (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area);
    hydrostaticGradient = state.gasCells[cellIndex - 1].dPdLHidro *
                          (kGravity * sin(state.gasCells[cellIndex].dutoL.teta) * gasDensity * dx);
    pressureDrop += (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    return pressureDrop;
}

double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.cells[cellIndex].dx;
    double diameter = state.cells[cellIndex].duto.a;
    double area = boreArea(diameter);
    double perimeter = state.cells[cellIndex].duto.peri;
    double meanTemperature = (state.cells[cellIndex].dx * state.cells[cellIndex].temp +
                              state.cells[cellIndex].dxL * state.cells[cellIndex - 1].temp) /
                             (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
    double injectedFluidDensity =
        state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presaux, meanTemperature);

    double fluidVelocity = state.cells[cellIndex - 1].QL / area;

    double fluidViscosity = state.cells[cellIndex].fluicol.VisFlu(state.cells[cellIndex].presaux, meanTemperature);

    double reynolds = state.cells[cellIndex].Rey(
        characteristicDiameter(state.cells[cellIndex].duto.revest, state.cells[cellIndex].duto.a, area, perimeter),
        fluidVelocity, injectedFluidDensity, fluidViscosity);
    double frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / diameter);
    double frictionGradient = state.cells[cellIndex].dPdLFric * 0.5 * frictionFactor * injectedFluidDensity *
                              (fabs(fluidVelocity) * fluidVelocity) * perimeter * dx / area;
    double hydrostaticGradient =
        state.cells[cellIndex].dPdLHidro * kGravity * sin(state.cells[cellIndex].duto.teta) * injectedFluidDensity * dx;

    double pressureDrop = (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    double meanPressure = state.cells[cellIndex - 1].pres;

    dx = 0.5 * state.cells[cellIndex].dxL;
    diameter = state.cells[cellIndex].dutoL.a;
    area = boreArea(diameter);
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
    frictionGradient = state.cells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * injectedFluidDensity *
                       (fabs(fluidVelocity) * fluidVelocity) * perimeter * dx / area;
    hydrostaticGradient = state.cells[cellIndex - 1].dPdLHidro * kGravity * sin(state.cells[cellIndex].dutoL.teta) *
                          injectedFluidDensity * dx;
    pressureDrop += (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;
    return pressureDrop;
}

void updateSteadyGasPressure(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.gasCells[cellIndex].dxL;
    double diameter = state.gasCells[cellIndex].dutoL.a;
    double area = boreArea(diameter);
    double perimeter = state.gasCells[cellIndex].dutoL.peri;
    double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres,
                                                                     state.gasCells[cellIndex - 1].temp);

    double leftFaceGasMassFlowRate = state.gasCells[cellIndex - 1].VGasR;
    double gasVelocity = leftFaceGasMassFlowRate / state.gasCells[cellIndex - 1].u1L;

    double gasViscosity = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres,
                                                                     state.gasCells[cellIndex - 1].temp);

    double reynolds =
        state.gasCells[cellIndex - 1].Rey(characteristicDiameter(state.gasCells[cellIndex].dutoL.revest,
                                                                 state.gasCells[cellIndex].dutoL.a, area, perimeter),
                                          gasVelocity, gasDensity, gasViscosity);
    double frictionFactor =
        state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
    double frictionGradient = state.gasCells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * gasDensity *
                              (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area;
    double hydrostaticGradient = state.gasCells[cellIndex - 1].dPdLHidro * kGravity *
                                 sin(state.gasCells[cellIndex].dutoL.teta) * gasDensity * dx;
    state.gasCells[cellIndex - 1].termoFric = frictionGradient / dx;
    state.gasCells[cellIndex - 1].termoHidro = hydrostaticGradient / dx;

    double meanPressure =
        state.gasCells[cellIndex - 1].pres - (frictionGradient + hydrostaticGradient) / units::kPascalPerKgfPerCm2;

    double meanTemperature = state.gasCells[cellIndex - 1].temp;

    dx = 0.5 * state.gasCells[cellIndex].dx0;
    diameter = state.gasCells[cellIndex].duto.a;
    area = boreArea(diameter);
    perimeter = state.gasCells[cellIndex].duto.peri;
    gasDensity = state.gasCells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
    const double leftFlowArea = state.gasCells[cellIndex].dutoL.area;
    const double flowArea = state.gasCells[cellIndex].duto.area;
    double areaVariation = 1 / (leftFlowArea * leftFlowArea) - 1 / (flowArea * flowArea);
    double dynamicPressureTerm = 0.5 * leftFaceGasMassFlowRate * leftFaceGasMassFlowRate * areaVariation / gasDensity;

    gasVelocity = leftFaceGasMassFlowRate / (gasDensity * area);

    gasViscosity = state.gasCells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);

    reynolds = state.gasCells[cellIndex].Rey(characteristicDiameter(state.gasCells[cellIndex].duto.revest,
                                                                    state.gasCells[cellIndex].duto.a, area, perimeter),
                                             gasVelocity, gasDensity, gasViscosity);
    frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    frictionGradient = state.gasCells[cellIndex].dPdLFric * 0.5 * frictionFactor * gasDensity *
                       (fabs(gasVelocity) * gasVelocity) * perimeter * dx / area;
    hydrostaticGradient =
        state.gasCells[cellIndex].dPdLHidro * kGravity * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx;

    state.gasCells[cellIndex].pres =
        meanPressure - (frictionGradient + hydrostaticGradient - dynamicPressureTerm) / units::kPascalPerKgfPerCm2;
    state.gasCells[cellIndex - 1].presR = state.gasCells[cellIndex].pres;
    if (cellIndex < state.lastGasCell)
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
            valveCell.massfonteCH = primeValveMassSource<kGasThroughValve>(
                state, valveIndex, valveCell,
                state.cells[state.productionValveCellIndices[valveIndex]]);

            if (cellIndex > 0)
                state.gasCells[cellIndex].VGasR =
                    state.gasCells[cellIndex - 1].VGasR - state.gasCells[cellIndex].massfonteCH;
            else
                state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex].massfonteCH;

            int productionValveCell = state.productionValveCellIndices[valveIndex];
            int gasValveCell = state.gasValveCellIndices[valveIndex];
            state.cells[productionValveCell].acsr.injg.QGas =
                state.gasCells[gasValveCell].massfonteCH * units::kSecondsPerDay /
                (kAirDensityAtStandardConditions * state.gasCells[gasValveCell].flui.Deng);
            const double gasTemperature =
                state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg
                    ? state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex)
                    : state.cells[productionValveCell].temp;
            state.cells[productionValveCell].acsr.injg.temp = gasTemperature;
            if (state.cells[productionValveCell].acsr.injg.temp < kMinimumTemperatureCelsius)
                state.cells[productionValveCell].acsr.injg.temp = kMinimumTemperatureCelsius;
            state.cells[productionValveCell].fontemassGR = state.gasCells[gasValveCell].massfonteCH;

            if (state.gasCells[cellIndex].VGasR <= 0. && cellIndex < state.lastGasCell) {
                state.gasCells[cellIndex].VGasR = 0.;
                for (int laterValveIndex = valveIndex + 1; laterValveIndex < valveCount; laterValveIndex++) {
                    state.gasCells[state.gasValveCellIndices[laterValveIndex]].massfonteCH = 0;
                    state.cells[state.productionValveCellIndices[laterValveIndex]].fontemassGR = 0.;
                }
            }

            if (cellIndex > 0)
                state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
            if (cellIndex < state.lastGasCell)
                state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
        }
    }
    if (!matched) {
        state.gasCells[cellIndex].VGasR =
            state.gasCells[cellIndex - 1].VGasR - 0 * state.gasCells[cellIndex].massfonteCH;
        state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
        if (cellIndex < state.lastGasCell)
            state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
    }
}

void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex) {
    int valveCount = state.input.nvalvgas;
    double casingTemperature = state.gasCells[0].calor.Textern1;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        int productionValveCell = state.productionValveCellIndices[valveIndex];
        int gasValveCell = state.gasValveCellIndices[valveIndex];
        if (cellIndex == productionValveCell) {
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
            state.cells[productionValveCell].acsr.injg.QGas =
                state.gasCells[gasValveCell].massfonteCH * units::kSecondsPerDay /
                state.gasCells[gasValveCell].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius);
            state.cells[productionValveCell].acsr.injg.temp = casingTemperature;
            state.cells[productionValveCell].fontemassGR = state.gasCells[gasValveCell].massfonteCH;
        }
    }
}

void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex) {
    if (!state.thermalSourceDisabled) {
        if (cellIndex < state.annulusCouplingStart || cellIndex > state.annulusCouplingEnd)
            marchSteadyGasTemperature(state, cellIndex);
        else
            coupledSteadyGasTemperature(state, cellIndex);
    } else {
        state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
    }
}

}  // namespace sisprod::gaslift
