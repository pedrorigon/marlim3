#include "SisProdGasLift.h"

#include "Leitura.h"
#include "SisProdConstants.h"
#include "Matriz.h"
#include "Vetor.h"
#include "celulaGas.h"
#include "celula3.h"
#include "chokegas.h"

#include <math.h>

namespace sisprod::gaslift {

namespace {

/// The diameter a Reynolds number is formed on: the bore when the segment is
/// bare, the hydraulic diameter when it is cased.
///
/// This replaces eleven if/else pairs. Several of them also chose which CELL to
/// call Rey on -- and that choice never mattered: CelG::Rey and Cel::Rey are
///
///     return dia * fabs(vel) * rho / (vis * 1e-3);
///
/// with no member access at all, so the receiver is irrelevant and two branches
/// that differed only in it computed the same number. One of those pairs even
/// selected OPPOSITE cells under the same predicate in its two halves, which
/// looks like a bug and is simply indifferent.
/// Below this opening the injection choke is throttling and becomes the mass
/// source of the first cell instead of a free boundary.
inline constexpr double kThrottlingChokeOpening = 0.2;

/// Fixed gas viscosity used to prime the steady annulus connection, in cP.
inline constexpr double kInitialAnnulusGasViscosityCentipoise = 0.16;

/// Gas-lift valve geometry. The bore threshold and the two spring rates are
/// imperial, as the whole valve correlation is.
inline constexpr double kCentimetrePerMetre = 100.;
inline constexpr double kCentimetrePerInch = 2.54;
inline constexpr double kLargeBoreValveInches = 1.1;
inline constexpr double kLargeBoreSpringRate = 500.0;
inline constexpr double kSmallBoreSpringRate = 1950.0;

[[nodiscard]] double characteristicDiameter(int cased, double bore, double area,
                                            double perimeter) {
    return cased == 0 ? bore : 4 * area / perimeter;
}

/// Opening of a calibrated gas-lift valve, from the state already primed on the
/// choke. Imperial in, fraction out.
[[nodiscard]] double calibratedValveOpening(const ChokeGas &choke) {
    return calibratedValveArea(choke.pcalib * kPsiPerKgfPerCm2, choke.tcalib,
                               (choke.presEstag - kAtmosphereInKgfPerCm2) * kPsiPerKgfPerCm2,
                               (choke.presGarg - kAtmosphereInKgfPerCm2) * kPsiPerKgfPerCm2,
                               choke.dextern, choke.areagarg,
                               choke.areagarg / choke.areafole,
                               celsiusToFahrenheit(choke.tempEstag));
}

/// What flows through a gas-lift valve, and how its throat pressure is
/// recovered. Four call sites shared one twelve-line body and differed in
/// exactly these two things, so they are what the policy carries.
struct GasThroughValve {
    static double recovery(const ChokeGas &choke) { return choke.frec; }
    static double massFlow(ChokeGas &choke, const Ler &) { return choke.massica(); }
};
struct CompletionFluidThroughValve {
    static double recovery(const ChokeGas &choke) { return choke.frecliq; }
    static double massFlow(ChokeGas &choke, const Ler &input) { return choke.massica(1, input.salinDescarga); }
};
/// Gas-side recovery, completion-fluid mass flow. Not a combination anyone
/// would invent; it is what the unloading search does and it is preserved.
struct UnloadingThroughValve {
    static double recovery(const ChokeGas &choke) { return choke.frec; }
    static double massFlow(ChokeGas &choke, const Ler &input) { return choke.massica(1, input.salinDescarga); }
};

/// Primes the choke from the two cells it spans and returns the mass source it
/// delivers, already scaled by the calibrated opening when the valve is one.
///
/// The four sites this replaces each wrote the result into a different place --
/// three into the gas cell's source term, one into a local -- so the value is
/// returned rather than stored.
template <typename Phase>
[[nodiscard]] double primeValveMassSource(const GasLiftState &state, int valveIndex,
                                          const CelG &gasCell, const Cel &productionCell) {
    ChokeGas &choke = state.gasLiftChokes[valveIndex];
    const double recovery = Phase::recovery(choke);
    choke.presEstag = gasCell.pres;
    choke.presGarg = (productionCell.pres - choke.presEstag * recovery) / (1. - recovery);
    choke.tempEstag = gasCell.temp;
    double massSource = productionCell.pres < gasCell.pres ? Phase::massFlow(choke, state.input) : 0.;
    if (choke.tipo == 1)
        massSource *= calibratedValveOpening(choke);
    return massSource;
}

/// Opens the injection choke if it is throttling, assembles the band system
/// from every gas cell and solves it.
///
/// advanceGasSubStep and advanceBufferedGasSubStep carried these twenty-one
/// lines byte for byte. The two differ only in what they do with the solution
/// afterwards, and in that one runs the assembly inside a thermal-coupling
/// cycle while the other runs it once.
void assembleAndSolveGasSystem(const GasLiftState &state) {
    double chokeOpeningFraction = 1.;
    if (state.gasCells[0].tipoCC == 0) {
        chokeOpeningFraction = state.injectionChoke.areagarg / state.gasCells[0].duto.area;
        if (chokeOpeningFraction < kThrottlingChokeOpening) {
            state.injectionChoke.presEstag = state.initialGasPressure;
            state.injectionChoke.tempEstag = state.initialGasTemperature;
            state.injectionChoke.presGarg = state.gasCells[0].pres;
            state.gasCells[0].massfonteCH = state.injectionChoke.massica();
        }
    }
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {

        state.gasCells[gasCellIndex].GeraLocal(state.gasCellCount, state.initialGasPressure, state.initialGasTemperature, chokeOpeningFraction);
        for (int bandIndex = 0; bandIndex < 9; bandIndex++) {
            state.gasSystemMatrix[3 * gasCellIndex][bandIndex - 3] = state.gasCells[gasCellIndex].local[0][bandIndex];
            state.gasSystemMatrix[3 * gasCellIndex + 1][bandIndex - 4] = state.gasCells[gasCellIndex].local[1][bandIndex];
            state.gasSystemMatrix[3 * gasCellIndex + 2][bandIndex - 5] = state.gasCells[gasCellIndex].local[2][bandIndex];
        }
        state.gasFreeTerms[3 * gasCellIndex] = state.gasCells[gasCellIndex].TL[0];
        state.gasFreeTerms[3 * gasCellIndex + 1] = state.gasCells[gasCellIndex].TL[1];
        state.gasFreeTerms[3 * gasCellIndex + 2] = state.gasCells[gasCellIndex].TL[2];
    }
    state.gasSystemMatrix.GaussElimPP(state.gasFreeTerms);
}

/// Clamps a gas cell's temperature to the physical range and copies it to the
/// neighbours that face it. Both arms of updateSteadyGasTemperature carried
/// these four statements verbatim.
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
/// each other's inverse and are preserved as written.
void relaxTowards(double &pressure, double target) {
    if (pressure > target) {
        pressure *= 0.95;
        if (pressure < target)
            pressure = target;
    } else if (pressure < target) {
        pressure *= 1.05;
        if (pressure > target)
            pressure = target;
    }
}

/// Publishes one fluid's transport properties onto the annulus cell that faces
/// a tubing cell. connectTubing, connectTubingSteady and
/// initializeTubingConnectionSteady each wrote these five fields in the same
/// order; only which fluid answers, and at which pressure and temperature,
/// differed.
struct ExternalFluidProperties {
    double conductivity;
    double specificHeat;
    double density;
    double viscosityCentipoise;
};

void publishExternalFluid(Cel &annulusCell, const ExternalFluidProperties &fluid) {
    annulusCell.calor.kextern1 = fluid.conductivity;
    annulusCell.calor.cpextern1 = fluid.specificHeat;
    annulusCell.calor.rhoextern1 = fluid.density;
    annulusCell.calor.viscextern1 = fluid.viscosityCentipoise * kPascalSecondPerCentipoise;
}

}  // namespace

void computeGasUnloadingHydrostatics(const GasLiftState &state) {
    double meanPressure;
    double meanTemperature;

    if (state.input.gasinj.tipoCC == 0 && state.input.controDesc == 0)
        meanPressure = state.input.gasinj.presinj[0];
    else {
        meanPressure = 10.;
        if (state.input.controDesc == 1) {
            meanPressure = state.input.presIniDescG;
            state.initialGasPressure = state.input.presIniDescG;
        }
    }

    state.gasCells[0].presL = meanPressure;
    state.gasCells[0].pres = meanPressure;
    state.gasCells[0].presini = meanPressure;
    state.gasCells[1].presL = meanPressure;
    meanTemperature = state.gasCells[0].calor.Textern1;
    state.gasCells[0].tempL = meanTemperature;
    state.gasCells[0].temp = meanTemperature;
    state.gasCells[1].tempL = meanTemperature;
    double rho0 = state.gasCells[0].flui.MasEspGas(meanPressure, meanTemperature);
    double rho1;
    state.gasCells[0].u1L = state.gasCells[0].duto.area * rho0;
    state.gasCells[0].u1LL = state.gasCells[0].u1L;
    state.gasCells[1].u1LL = state.gasCells[0].u1L;
    state.gasCells[0].VGasL = 0;
    state.gasCells[0].VGasR = 0;
    state.gasCells[1].VGasL = 0;
    state.gasCells[0].massfonteCH = 0.;
    for (int gasCellIndex = 1; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double halfUpstreamLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        meanPressure -= rho0 * kGravityUnloadingVariant * halfUpstreamLength * sin(state.gasCells[gasCellIndex - 1].duto.teta) / kPascalPerKgfPerCm2Variant;
        meanTemperature = state.gasCells[gasCellIndex].calor.Textern1;
        if (gasCellIndex < state.interfaceCell)
            rho1 = state.gasCells[gasCellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        else
            rho1 = state.gasCells[gasCellIndex].MasEspFlu(meanPressure, meanTemperature);
        meanPressure -= rho1 * kGravityUnloadingVariant * halfLocalLength * sin(state.gasCells[gasCellIndex].duto.teta) / kPascalPerKgfPerCm2Variant;
        rho0 = rho1;

        state.gasCells[gasCellIndex].pres = meanPressure;
        state.gasCells[gasCellIndex].presini = meanPressure;
        state.gasCells[gasCellIndex - 1].presR = meanPressure;
        state.gasCells[gasCellIndex].temp = meanTemperature;
        state.gasCells[gasCellIndex - 1].tempR = meanTemperature;
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].duto.area * rho0;
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

void updateGasLine(const GasLiftState &state) {
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        if (gasCellIndex != 0 && gasCellIndex != state.gasCellCount) {
            state.gasCells[gasCellIndex].pres = state.gasFreeTerms[3 * gasCellIndex];
            state.gasCells[gasCellIndex].presL = state.gasFreeTerms[3 * gasCellIndex - 3];
            state.gasCells[gasCellIndex].presR = state.gasFreeTerms[3 * gasCellIndex + 3];
            state.gasCells[gasCellIndex].VGasR = state.gasFreeTerms[3 * gasCellIndex + 1];
            state.gasCells[gasCellIndex].VGasL = state.gasFreeTerms[3 * gasCellIndex - 2];
            state.gasCells[gasCellIndex].VGasRR = state.gasFreeTerms[3 * gasCellIndex + 4];
        } else if (gasCellIndex == 0) {
            state.gasCells[gasCellIndex].pres = state.gasFreeTerms[3 * gasCellIndex];
            state.gasCells[gasCellIndex].presL = state.gasCells[gasCellIndex].pres;
            state.gasCells[gasCellIndex].presR = state.gasFreeTerms[3 * gasCellIndex + 3];
            state.gasCells[gasCellIndex].VGasR = state.gasFreeTerms[3 * gasCellIndex + 1];
            state.gasCells[gasCellIndex].VGasRR = state.gasFreeTerms[3 * gasCellIndex + 4];
        } else {
            state.gasCells[gasCellIndex].pres = state.gasFreeTerms[3 * gasCellIndex];
            state.gasCells[gasCellIndex].presL = state.gasFreeTerms[3 * gasCellIndex - 3];
            state.gasCells[gasCellIndex].VGasR = state.gasFreeTerms[3 * gasCellIndex + 1];
            state.gasCells[gasCellIndex].VGasL = state.gasFreeTerms[3 * gasCellIndex - 2];
            state.gasCells[gasCellIndex].presR = state.gasCells[gasCellIndex].pres;
        }
    }
}

void updateBufferedGasLine(const GasLiftState &state) {
    // This had three branches -- interior cell, inlet, outlet -- carrying the
    // same assignment byte for byte. The conditions are comparisons on an int
    // and a const int&, so they had no effect beyond choosing which copy of one
    // statement to run.
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].VGasRBuf = state.gasFreeTerms[3 * gasCellIndex + 1];
}

double calibratedValveArea(double calibrationPressure, double calibrationTemperature, double valveOpeningPressure, double tubingPressure,
                           double externalDiameter, double throatArea, double valveRatio, double bottomHoleTemperatureFahrenheit) {
    // Imperial throughout: pressures in psi, temperature in Fahrenheit. The
    // valve opening pressure is the casing pressure; the ratio is of areas.

    double bellowsPressureAt80F = calibrationPressure * (1 - valveRatio);
    bellowsPressureAt80F = (bellowsPressureAt80F + 14.6959488) * (80 + 460.67) / (calibrationTemperature * 1.8 + 491.67) - 14.6959488;
    double bellowsPressure = bellowsPressureAt80F * (1 + 0.00215 * (bottomHoleTemperatureFahrenheit - 80));
    // The two assignments to the opening fraction that used to stand here --
    // one on `openingCriterion > bellowsPressure`, one on the closed-valve rule
    // below -- were both overwritten unconditionally by the area ratio at the
    // end of the function, so neither ever reached the return. The criterion
    // itself fed nothing else and is gone with them.
    //
    // This is reported, not silently repaired: it means the documented rule
    // "IF THE VALVE IS CLOSED, QG = 0" has no effect today. Restoring it would
    // change results and is an engineering decision, not a refactoring one.
    double bellowsArea = throatArea / valveRatio;
    // The spring rate depends only on whether the valve is a large-bore one.
    const double externalDiameterInches = externalDiameter * kCentimetrePerMetre / kCentimetrePerInch;
    const double bellowsSpringRate =
        (externalDiameterInches > kLargeBoreValveInches ? kLargeBoreSpringRate : kSmallBoreSpringRate) * bellowsArea;
    double stemTravel = ((valveOpeningPressure - bellowsPressure) * bellowsArea - (valveOpeningPressure - tubingPressure) * throatArea) / bellowsSpringRate;

    double throatDiameter = sqrt(throatArea * 4. / M_PI);
    double throatRadius = throatDiameter / 2.0;
    double bellowsRadius = sqrt(bellowsArea / M_PI);
    // Named once instead of evaluated three times: the seat offset, the
    // distance from the stem axis to where the bellows meets the seat.
    const double seatOffset = sqrt(bellowsRadius * bellowsRadius - throatRadius * throatRadius);
    double openingArea = M_PI * throatRadius * stemTravel * (stemTravel + 2.0 * seatOffset);
    openingArea = openingArea / sqrt((stemTravel + seatOffset) * (stemTravel + seatOffset) + throatRadius * throatRadius);
    if (openingArea > throatArea)
        openingArea = throatArea;

    return openingArea / throatArea;
}

double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate, int valveIndex, double factor, int sign) {
    int gasValveCell = state.gasValveCellIndices[valveIndex];
    state.gasLiftChokes[valveIndex].presEstag = state.gasCells[gasValveCell].pres;
    state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[gasValveCell].temp;
    double rho0 = state.gasCells[gasValveCell].MasEspFlu(state.gasLiftChokes[valveIndex].presEstag, state.gasLiftChokes[valveIndex].tempEstag);
    double massFlowRate = (factor * state.input.vazDescControl - maximumFlowRate) * rho0;
    const double pressureCorrection =
        pow(massFlowRate / state.gasLiftChokes[valveIndex].areagarg, 2.) / (2. * rho0 * kPascalPerKgfPerCm2Variant);
    return sign * pressureCorrection;
}

double computeUnloadingValvePressure(const GasLiftState &state, double throatFlowRate, int valveIndex) {

    double maximumVelocity = 0;
    double upperBandMargin = 0.1;
    double lowerBandMargin = 0.4;

    double meanPressure;
    meanPressure = state.gasSurfacePressure;
    state.input.presMaxDesc = 100000.;
    for (int cellIndex = state.lastCell; cellIndex >= 0; cellIndex--) {

        double meanTemperature = state.cells[cellIndex].temp;
        double flowArea = state.cells[cellIndex].duto.area;
        double perimeter = state.cells[cellIndex].duto.peri;
        double cellLength = state.cells[cellIndex].dx;
        double gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        double rhoP = state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature);
        double rhoC = state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        double viscG = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        double viscP = state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature);
        double viscC = state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
        double voidFraction = state.cells[cellIndex].alf;
        double composition = state.cells[cellIndex].bet;
        double liquidDensity = composition * rhoC + (1. - composition) * rhoP;
        double liquidViscosity = composition * viscC + (1. - composition) * viscP;
        double mixtureDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;
        double mixtureViscosity = voidFraction * viscG + (1. - voidFraction) * liquidViscosity;
        double vel1 = state.cells[cellIndex].QL / (flowArea * liquidDensity) + state.cells[cellIndex].QG / (flowArea * gasDensity);
        // One guard where there were two testing the same expression. The
        // first left `reynolds` uninitialised on its false path, which was safe
        // only because the second never read it there.
        double frictionFactor = 0.;
        if (fabs(vel1) > 1e-15) {
            const double reynolds = state.cells[cellIndex].Rey(
                characteristicDiameter(state.cells[cellIndex].duto.revest, state.cells[cellIndex].duto.a, flowArea, perimeter),
                vel1, mixtureDensity, mixtureViscosity);
            frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / state.cells[cellIndex].duto.a);
        }
        double tens1 = frictionFactor * mixtureDensity * vel1 * fabs(vel1) / 2.;
        meanPressure -= (-kGravity * mixtureDensity * sin(state.cells[cellIndex].duto.teta) - tens1 * perimeter / flowArea) * cellLength / kPascalPerKgfPerCm2;
        if (state.cells[cellIndex].acsr.tipo == 3) {
            double candidateMaximumPressure = state.cells[cellIndex].acsr.ipr.Pres - (meanPressure - state.gasSurfacePressure);
            if (candidateMaximumPressure < state.input.presMaxDesc)
                state.input.presMaxDesc = candidateMaximumPressure;
        }
    }
    if (state.gasSurfacePressure >= state.input.presMaxDesc * 0.9999999 && throatFlowRate > (1 - upperBandMargin) * state.input.vazDescControl && state.gasCells[0].VGasR > 0.) {
        double pressureCorrection = unloadingPressureCorrection(state, throatFlowRate, valveIndex, (1 - upperBandMargin), -1);
        if (fabs(pressureCorrection) > 0.01 * state.initialGasPressure * state.cells[0].dt)
            pressureCorrection = (fabs(pressureCorrection) / pressureCorrection) * 0.01 * state.initialGasPressure * state.cells[0].dt;
        state.initialGasPressure += pressureCorrection;
        if (state.initialGasPressure < state.input.presMinDescG) {
            state.initialGasPressure = state.input.presMinDescG;
        }

    } else if (throatFlowRate <= (1 - lowerBandMargin) * state.input.vazDescControl && state.gasSurfacePressure <= state.input.presMinDesc * 1.0000001) {
        double pressureCorrection = unloadingPressureCorrection(state, throatFlowRate, valveIndex, 1 - lowerBandMargin, 1);
        if (fabs(pressureCorrection) > 0.01 * state.initialGasPressure * state.cells[0].dt)
            pressureCorrection = (fabs(pressureCorrection) / pressureCorrection) * 0.01 * state.initialGasPressure * state.cells[0].dt;
        state.initialGasPressure += pressureCorrection;
        if (state.initialGasPressure > state.input.presMaxDescG) {
            state.initialGasPressure = state.input.presMaxDescG;
        }
    }

    return maximumVelocity;
}

void solveUnloading(const GasLiftState &state) {
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double halfUpstreamLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        double RgasR = 1.;
        if (state.gasCells[gasCellIndex].razInter <= 0.5)
            RgasR = 2 * state.gasCells[gasCellIndex].razInter;
        double RgasL = 0.;
        if (state.gasCells[gasCellIndex - 1].razInter >= 0.5)
            RgasL = 2 * (state.gasCells[gasCellIndex - 1].razInter - 0.5);
        double upstreamGasLength = halfUpstreamLength * RgasL;
        double localGasLength = halfLocalLength * RgasR;
        double upstreamLiquidLength = halfUpstreamLength - upstreamGasLength;
        double localLiquidLength = halfLocalLength - localGasLength;
        double temp = state.gasCells[gasCellIndex - 1].temp;
        double upstreamPressure = state.gasCells[gasCellIndex - 1].pres;
        double rhoL = state.gasCells[gasCellIndex].MasEspFlu(upstreamPressure, temp);
        double viscL = state.gasCells[gasCellIndex].VisFlu(upstreamPressure, temp);
        double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(upstreamPressure, temp);
        double viscG = state.gasCells[gasCellIndex].flui.ViscGas(upstreamPressure, temp);
        // Which density carries the velocity, decided once. Written as an
        // assignment followed by a conditional overwrite, this did two
        // divisions per velocity and discarded the first.
        const double carryingDensity =
            state.gasCells[gasCellIndex].razInter > (*state.globals).localtiny ? gasDensity : rhoL;
        double vel1 = state.gasCells[gasCellIndex].VGasL / (carryingDensity * state.gasCells[gasCellIndex - 1].duto.area);
        double vel2 = state.gasCells[gasCellIndex].VGasL / (carryingDensity * state.gasCells[gasCellIndex].duto.area);
        double re1G;
        double re1L;
        double re2G;
        double re2L;
        re1L = state.gasCells[gasCellIndex - 1].Rey(characteristicDiameter(state.gasCells[gasCellIndex - 1].duto.revest, state.gasCells[gasCellIndex - 1].duto.a, state.gasCells[gasCellIndex - 1].duto.area, state.gasCells[gasCellIndex - 1].duto.peri),
                                 vel1, rhoL, viscL);
        re2L = state.gasCells[gasCellIndex].Rey(characteristicDiameter(state.gasCells[gasCellIndex].duto.revest, state.gasCells[gasCellIndex].duto.a, state.gasCells[gasCellIndex].duto.area, state.gasCells[gasCellIndex].duto.peri),
                                 vel2, rhoL, viscL);
        re1G = state.gasCells[gasCellIndex - 1].Rey(characteristicDiameter(state.gasCells[gasCellIndex - 1].duto.revest, state.gasCells[gasCellIndex - 1].duto.a, state.gasCells[gasCellIndex - 1].duto.area, state.gasCells[gasCellIndex - 1].duto.peri),
                                 vel1, gasDensity, viscG);
        re2G = state.gasCells[gasCellIndex].Rey(characteristicDiameter(state.gasCells[gasCellIndex].duto.revest, state.gasCells[gasCellIndex].duto.a, state.gasCells[gasCellIndex].duto.area, state.gasCells[gasCellIndex].duto.peri),
                                 vel2, gasDensity, viscG);
        double upstreamLiquidFriction = state.gasCells[gasCellIndex - 1].fric(re1L, state.gasCells[gasCellIndex - 1].duto.rug / state.gasCells[gasCellIndex - 1].duto.a) * upstreamLiquidLength;
        double localLiquidFriction = state.gasCells[gasCellIndex].fric(re2L, state.gasCells[gasCellIndex].duto.rug / state.gasCells[gasCellIndex].duto.a) * localLiquidLength;
        double hidro1L = (kGravity * sin(state.gasCells[gasCellIndex - 1].duto.teta) * rhoL) * upstreamLiquidLength;
        double hidro2L = (kGravity * sin(state.gasCells[gasCellIndex].duto.teta) * rhoL) * localLiquidLength;
        double upstreamGasFriction = state.gasCells[gasCellIndex - 1].fric(re1G, state.gasCells[gasCellIndex - 1].duto.rug / state.gasCells[gasCellIndex - 1].duto.a) * upstreamGasLength;
        double localGasFriction = state.gasCells[gasCellIndex].fric(re2G, state.gasCells[gasCellIndex].duto.rug / state.gasCells[gasCellIndex].duto.a) * localGasLength;
        double hidro1G = (kGravity * sin(state.gasCells[gasCellIndex - 1].duto.teta) * gasDensity) * upstreamGasLength;
        double hidro2G = (kGravity * sin(state.gasCells[gasCellIndex].duto.teta) * gasDensity) * localGasLength;
        state.gasCells[gasCellIndex].pres = state.gasCells[gasCellIndex - 1].pres + (-0.5 * (upstreamLiquidFriction * rhoL + upstreamGasFriction * gasDensity) * vel1 * fabs(vel1) * state.gasCells[gasCellIndex - 1].duto.peri / state.gasCells[gasCellIndex - 1].duto.area - 0.5 * (localLiquidFriction * rhoL + localGasFriction * gasDensity) * vel2 * fabs(vel2) * state.gasCells[gasCellIndex].duto.peri / state.gasCells[gasCellIndex].duto.area - hidro1L - hidro2L - hidro1G - hidro2G) / kPascalPerKgfPerCm2Variant;
        state.gasCells[gasCellIndex].presL = state.gasCells[gasCellIndex - 1].pres;
        state.gasCells[gasCellIndex - 1].presR = state.gasCells[gasCellIndex].pres;

        state.temperatureUpdater.dischargeTemperature(gasCellIndex);

        state.gasCells[gasCellIndex].u1L = ((1. - state.gasCells[gasCellIndex].razInter) * state.gasCells[gasCellIndex].MasEspFlu(state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp) + state.gasCells[gasCellIndex].razInter * state.gasCells[gasCellIndex].flui.MasEspGas(upstreamPressure, temp)) * state.gasCells[gasCellIndex].duto.area;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        state.gasCells[gasCellIndex].u1LL = state.gasCells[gasCellIndex - 1].u1L;
    }

    double Qtotal = 0.;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double temp = state.gasCells[gasCellIndex].temp;
        double upstreamPressure = state.gasCells[gasCellIndex].pres;
        double rhoL = state.gasCells[gasCellIndex].MasEspFlu(upstreamPressure, temp);
        double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(upstreamPressure, temp);
        double sourceVolumeFlow = state.gasCells[gasCellIndex].massfonteCH / rhoL;
        if (state.gasCells[gasCellIndex].razInter > 0.5)
            sourceVolumeFlow = state.gasCells[gasCellIndex].massfonteCH / gasDensity;
        Qtotal += sourceVolumeFlow;
    }
    double temp = state.gasCells[state.interfaceCell].temp;
    double upstreamPressure = state.gasCells[state.interfaceCell].pres;
    double gasDensity = state.gasCells[state.interfaceCell].flui.MasEspGas(upstreamPressure, temp);
    state.gasCells[state.interfaceCell].VGasL = Qtotal * gasDensity;
    state.gasCells[state.interfaceCell - 1].VGasR = state.gasCells[state.interfaceCell].VGasL;
    state.gasCells[state.interfaceCell - 2].VGasRR = state.gasCells[state.interfaceCell].VGasL;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double temp = state.gasCells[gasCellIndex].temp;
        double upstreamPressure = state.gasCells[gasCellIndex].pres;
        double rhoL = state.gasCells[gasCellIndex].MasEspFlu(upstreamPressure, temp);
        double gasDensity = state.gasCells[gasCellIndex].flui.MasEspGas(upstreamPressure, temp);
        double sourceVolumeFlow = state.gasCells[gasCellIndex].massfonteCH / rhoL;
        if (state.gasCells[gasCellIndex].razInter > 0.5)
            sourceVolumeFlow = state.gasCells[gasCellIndex].massfonteCH / gasDensity;
        Qtotal -= sourceVolumeFlow;
        if (gasCellIndex < state.gasCellCount) {
            state.gasCells[gasCellIndex + 1].VGasL = (Qtotal)*rhoL;
            state.gasCells[gasCellIndex].VGasR = state.gasCells[gasCellIndex + 1].VGasL;
        } else
            state.gasCells[gasCellIndex].VGasR = 0.;
        state.gasCells[gasCellIndex - 1].VGasRR = state.gasCells[gasCellIndex].VGasL;
    }

    state.interfaceVelocity = state.gasCells[state.interfaceCell + 1].VGasL / (state.gasCells[state.interfaceCell + 1].MasEspFlu(state.gasCells[state.interfaceCell + 1].pres, state.gasCells[state.interfaceCell + 1].temp) * state.gasCells[state.interfaceCell].duto.area);
}

void advanceInterface(const GasLiftState &state) {

    state.gasCells[state.interfaceCell].razInter = (state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 + state.interfaceVelocity * state.timeStep) / state.gasCells[state.interfaceCell].dx0;
    if (state.interfaceCell == (state.gasCellCount - 1) && state.gasCells[state.interfaceCell].razInter >= 0.99) {

        double meanPressure = state.gasCells[state.interfaceCell].pres;
        double meanTemperature = state.gasCells[state.interfaceCell].temp;
        double interfaceFlowArea = state.gasCells[state.interfaceCell].duto.area;
        double rho1 = state.gasCells[state.interfaceCell].flui.MasEspGas(meanPressure, meanTemperature);
        double rhoL = state.gasCells[state.interfaceCell].MasEspFlu(meanPressure, meanTemperature);
        state.gasCells[state.interfaceCell].VGasR = state.gasCells[state.interfaceCell].VGasR * rho1 / rhoL;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * rho1;
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
        rho1 = state.gasCells[state.interfaceCell].flui.MasEspGas(meanPressure, meanTemperature);
        rhoL = state.gasCells[state.interfaceCell].MasEspFlu(meanPressure, meanTemperature);
        state.gasCells[state.interfaceCell].VGasR = 0 * state.gasCells[state.interfaceCell].VGasR * rho1 / rhoL;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * rho1;
        state.gasCells[state.interfaceCell].u1R = interfaceFlowArea * rho1;
        state.gasCells[state.interfaceCell - 1].u1R = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell - 1].VGasRR = state.gasCells[state.interfaceCell].VGasR;

        state.interfaceCell = 1e7;
        state.input.descarga = 0;
    }

    if (state.input.descarga == 1) {
        if (((state.gasCells[state.interfaceCell].razInter <= (*state.globals).localtiny) && (state.gasCells[state.interfaceCell].razInter >= -(*state.globals).localtiny)))
            state.gasCells[state.interfaceCell].razInter = 0;
        else if (state.gasCells[state.interfaceCell].razInter < -(*state.globals).localtiny) {
            // Both arms of the branch that used to stand here set the ratio to
            // the same value; only the time-step side effect was conditional.
            const double interfaceCrossingTimeStep =
                -state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 / state.interfaceVelocity;
            if (interfaceCrossingTimeStep > (*state.globals).localtiny)
                state.interfaceTimeStep = interfaceCrossingTimeStep;
            state.gasCells[state.interfaceCell].razInter = 0.;
        } else if ((state.gasCells[state.interfaceCell].razInter >= (1. - (*state.globals).localtiny) && state.gasCells[state.interfaceCell].razInter <= (1. + (*state.globals).localtiny))) {
            state.gasCells[state.interfaceCell].razInter = 1.;
        } else if (state.gasCells[state.interfaceCell].razInter > (1. + (*state.globals).localtiny)) {
            const double interfaceCrossingTimeStep =
                (1. - state.gasCells[state.interfaceCell].razInterIni) * state.gasCells[state.interfaceCell].dx0 / state.interfaceVelocity;
            if (interfaceCrossingTimeStep > (*state.globals).localtiny)
                state.interfaceTimeStep = interfaceCrossingTimeStep;
            state.gasCells[state.interfaceCell].razInter = 1.;
        }

        if (fabs(state.gasCells[state.interfaceCell].razInter) < (*state.globals).localtiny && state.interfaceVelocity < -fabs((*state.globals).localtiny)) {
            state.gasCells[state.interfaceCell].razInter = 0.;
            if (state.interfaceCell > 0) {
                state.interfaceCell--;
                state.gasCells[state.interfaceCell].razInter = 1.;
            }
        } else if (fabs(state.gasCells[state.interfaceCell].razInter - 1.) < (*state.globals).localtiny && state.interfaceVelocity > (*state.globals).localtiny) {
            state.gasCells[state.interfaceCell].razInter = 1.;
            if (state.interfaceCell < state.gasCellCount) {
                state.interfaceCell++;
                state.gasCells[state.interfaceCell].razInter = 0.;
            }
        }
    }
}

void updateTransientGasValves(const GasLiftState &state) {
    int valveCount = state.input.nvalvgas;
    int master2CellIndex = state.input.master2.posic;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        CelG &gasCell = state.gasCells[state.gasValveCellIndices[valveIndex]];
        const Cel &productionCell = state.cells[state.productionValveCellIndices[valveIndex]];
        if (state.gasValveCellIndices[valveIndex] < state.interfaceCell)
            gasCell.massfonteCH =
                primeValveMassSource<GasThroughValve>(state, valveIndex, gasCell, productionCell);
    }
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {

        state.gasCells[gasCellIndex].massfonteCH = 0.;
        // Behind the interface the valve discharges completion fluid, ahead of
        // it gas. That was the ONLY difference between the two loops that used
        // to stand here; everything else in them was identical.
        const bool valveBehindInterface = state.gasCells[gasCellIndex].razInter < 0.5;
        for (int candidateValveIndex = 0; candidateValveIndex < valveCount; candidateValveIndex++) {
            if (gasCellIndex == state.gasValveCellIndices[candidateValveIndex]) {
                CelG &valveCell = state.gasCells[state.gasValveCellIndices[candidateValveIndex]];
                const Cel &productionCell = state.cells[state.productionValveCellIndices[candidateValveIndex]];
                valveCell.massfonteCH =
                    valveBehindInterface
                        ? primeValveMassSource<CompletionFluidThroughValve>(state, candidateValveIndex, valveCell, productionCell)
                        : primeValveMassSource<GasThroughValve>(state, candidateValveIndex, valveCell, productionCell);
            }
        }
    }

    double master2FlowArea = state.gasCells[master2CellIndex].duto.area;
    if (master2FlowArea > state.gasCells[master2CellIndex].dutoR.area)
        master2FlowArea = state.gasCells[master2CellIndex].dutoR.area;
    if (state.gasCells[master2CellIndex].chkcell.areagarg <= 0.01 * master2FlowArea &&
        state.gasCells[master2CellIndex].chkcell.areagarg >= 1e-5 * master2FlowArea) {
        state.gasCells[master2CellIndex].chkcell.tempEstag = state.gasCells[master2CellIndex].temp;
        double rhoM = state.gasCells[master2CellIndex].flui.MasEspGas(state.gasCells[master2CellIndex].pres, state.gasCells[master2CellIndex].temp);
        double rhoJ = state.gasCells[master2CellIndex + 1].flui.MasEspGas(state.gasCells[master2CellIndex + 1].pres, state.gasCells[master2CellIndex + 1].temp);
        double hidroM = -0.5 * rhoM * state.gasCells[master2CellIndex].dx0 * sin(state.gasCells[master2CellIndex].duto.teta) / kPascalPerKgfPerCm2Variant;
        double hidroJ = 0.5 * rhoJ * state.gasCells[master2CellIndex + 1].dx0 * sin(state.gasCells[master2CellIndex + 1].duto.teta) / kPascalPerKgfPerCm2Variant;
        state.gasCells[master2CellIndex].chkcell.presEstag = state.gasCells[master2CellIndex].pres + hidroM;
        state.gasCells[master2CellIndex].chkcell.presGarg = state.gasCells[master2CellIndex + 1].pres + hidroJ;
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

double searchUnloadingInjectionPressure(const GasLiftState &state) {
    double maximumUnloadingFlowRate = 0;
    double instantaneousMaximumFlowRate = 0;

    double upperBandMargin = 0.1;
    double lowerBandMargin = 0.4;
    if (state.input.descarga == 1) {
        double initialUnloadingPressure = state.cells[state.lastCell].pres;
        if ((*state.globals).lixo5 > state.input.tempoLatenciaDesc) {
            int valveCount = state.input.nvalvgas;
            int selectedValveIndex = 0;
            for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
                int gasValveCell = state.gasValveCellIndices[valveIndex];
                int productionValveCell = state.productionValveCellIndices[valveIndex];
                if (gasValveCell > state.interfaceCell) {
                    const CelG &valveCell = state.gasCells[gasValveCell];
                    double rho1 = valveCell.MasEspFlu(valveCell.pres, valveCell.temp);
                    const double massFlowRate = primeValveMassSource<UnloadingThroughValve>(
                        state, valveIndex, valveCell, state.cells[productionValveCell]);
                    double candidateFlowRate = massFlowRate / (rho1);
                    if (candidateFlowRate > maximumUnloadingFlowRate) {
                        maximumUnloadingFlowRate = candidateFlowRate;
                        selectedValveIndex = valveIndex;
                    }
                }
            }
            if (state.unloadingTimeSteps.size() > state.maximumContinuousUnloadingCount || state.meanUnloadingTemperature > state.continuousMeanUnloadingTemperature) {
                state.meanUnloadingTemperature -= state.unloadingTimeSteps[0];
                state.unloadingTimeSteps.erase(state.unloadingTimeSteps.begin());
                state.meanUnloadingFlowRate -= state.maximumMeanUnloadingFlowRates[0];
                state.maximumMeanUnloadingFlowRates.erase(state.maximumMeanUnloadingFlowRates.begin());
            }

            state.maximumMeanUnloadingFlowRates.push_back(maximumUnloadingFlowRate * state.timeStep);
            state.meanUnloadingFlowRate += maximumUnloadingFlowRate * state.timeStep;
            state.unloadingTimeSteps.push_back(state.timeStep);
            state.meanUnloadingTemperature += state.timeStep;

            double weight = 0.5;
            instantaneousMaximumFlowRate = maximumUnloadingFlowRate;
            maximumUnloadingFlowRate = weight * maximumUnloadingFlowRate + (1. - weight) * state.meanUnloadingFlowRate / state.meanUnloadingTemperature;
            computeUnloadingValvePressure(state, maximumUnloadingFlowRate, selectedValveIndex);
            if (state.gasValveCellIndices[selectedValveIndex] >= state.interfaceCell) {
                if (state.cells[state.lastCell - 1].MC > 0.0) {
                    if (maximumUnloadingFlowRate > (1. - upperBandMargin) * state.input.vazDescControl) {

                        double pressureCorrection = 0.;
                        pressureCorrection = unloadingPressureCorrection(state, maximumUnloadingFlowRate, selectedValveIndex, 1. - upperBandMargin, -1);
                        if (fabs(pressureCorrection) > 0.01 * state.gasSurfacePressure * state.cells[0].dt)
                            pressureCorrection = (fabs(pressureCorrection) / pressureCorrection) * 0.01 * state.gasSurfacePressure * state.cells[0].dt;
                        state.gasSurfacePressure -= pressureCorrection;
                        if (state.gasSurfacePressure > state.input.presMaxDesc)
                            state.gasSurfacePressure = state.input.presMaxDesc;
                    } else if (maximumUnloadingFlowRate <= (1 - lowerBandMargin) * state.input.vazDescControl) {

                        double pressureCorrection = 0.;
                        pressureCorrection = unloadingPressureCorrection(state, maximumUnloadingFlowRate, selectedValveIndex, 1. - lowerBandMargin, 1);
                        if (fabs(pressureCorrection) > 0.01 * state.gasSurfacePressure * state.cells[0].dt)
                            pressureCorrection = (fabs(pressureCorrection) / pressureCorrection) * 0.01 * state.gasSurfacePressure * state.cells[0].dt;
                        state.gasSurfacePressure -= pressureCorrection;
                        if (state.gasSurfacePressure < state.input.presMinDesc)
                            state.gasSurfacePressure = state.input.presMinDesc;
                    }
                } else {
                    state.gasSurfacePressure *= (1 - 0.01 * state.timeStep);
                    if (state.gasSurfacePressure < state.input.presMinDesc)
                        state.gasSurfacePressure = state.input.presMinDesc;
                }
            }
        } else {
            if ((*state.globals).lixo5 < 0.5 * state.input.tempoLatenciaDesc)
                state.gasSurfacePressure = initialUnloadingPressure - (initialUnloadingPressure - state.input.presMinDesc) * state.cells[0].dt / (0.5 * state.input.tempoLatenciaDesc - (*state.globals).lixo5);
            else
                state.gasSurfacePressure = state.input.presMinDesc;
            state.initialGasPressure = state.input.presIniDescG;
        }
    } else {
        state.gasSurfacePressure *= 0.95;
        if (state.gasSurfacePressure < state.input.presMinDesc)
            state.gasSurfacePressure = state.input.presMinDesc;
        relaxTowards(state.initialGasPressure, state.input.presIniDesc);
    }
    return instantaneousMaximumFlowRate;
}

void advanceGasSubStep(const GasLiftState &state) {

    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].DeVoltaParaoFuturo();
    state.interfaceTimeStep = state.timeStep;
    if (state.input.descarga == 1) {
        state.initialInterfaceCell = state.interfaceCell;
        state.initialInterfaceTimeStep = state.interfaceTimeStep;
        state.initialInterfaceVelocity = state.interfaceVelocity;
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
        state.gasCells[0].temp = state.initialGasTemperature;
        state.gasCells[0].dTdt = (state.gasCells[0].temp - state.gasCells[0].tempini) / state.timeStep;
        if (state.gasCells[state.gasCellCount].VGasR < 0.)
            state.gasCells[state.gasCellCount].temp = 20.;
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int gasCellIndex = 1; gasCellIndex < cellLimit; gasCellIndex++) {
            state.temperatureUpdater.gasTemperature(gasCellIndex, state.gasCells[gasCellIndex - 1].tempini);
            state.gasCells[gasCellIndex].dTdt = (state.gasCells[gasCellIndex].temp - state.gasCells[gasCellIndex].tempini) / state.timeStep;
        }
        if (thermalCouplingCycle < maximumThermalCouplingCycles)
            for (int renewedCellIndex = 0; renewedCellIndex <= state.gasCellCount; renewedCellIndex++)
                state.gasCells[renewedCellIndex].FeiticoDoTempo();
    }

#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int gasCellIndex = 0; gasCellIndex < cellLimit; gasCellIndex++)
        state.gasCells[gasCellIndex].rg = state.gasCells[gasCellIndex].flui.MasEspGas(state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp);

    for (int gasCellIndex = 0; gasCellIndex < cellLimit; gasCellIndex++) {
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].rg * state.gasCells[gasCellIndex].duto.area;
        // The inlet feeds itself; every other cell is fed by the one upstream.
        state.gasCells[gasCellIndex].u1LL = gasCellIndex == 0
                                                ? state.gasCells[gasCellIndex].u1L
                                                : state.gasCells[gasCellIndex - 1].u1L;
        // Two passes over the same predicate became one. Moving the density
        // hand-back after u1L is safe: nothing between them reads rgR, and
        // u1LL above reads the UPSTREAM cell's u1L, set an iteration ago.
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

void advanceBufferedGasSubStep(const GasLiftState &state) {

    for (int gasCellIndex = 0; gasCellIndex < state.gasCellCount + 1; gasCellIndex++)
        state.gasCells[gasCellIndex].dTdt = 0;
    assembleAndSolveGasSystem(state);
    updateBufferedGasLine(state);
}

void connectTubing(const GasLiftState &state) {
    for (int annulusCellIndex = state.annulusTubingStart; annulusCellIndex >= state.annulusTubingEnd; annulusCellIndex--) {
        int tubingCellIndex = state.annulusTubingStart + state.tubingAnnulusStart - annulusCellIndex;
        state.cells[annulusCellIndex].calor.Textern2 = state.gasCells[tubingCellIndex].calor.Tcamada[0][0];
        state.cells[annulusCellIndex].calor.betext = state.gasCells[tubingCellIndex].calor.betint;
        int layerIndex = state.cells[annulusCellIndex].calor.geom.ncamadas - 1;
        int discretizationIndex = state.cells[annulusCellIndex].calor.ncamada[layerIndex] - 1;
        state.gasCells[tubingCellIndex].calor.Tint2 = state.cells[annulusCellIndex].calor.Tcamada[layerIndex][discretizationIndex];
        state.cells[annulusCellIndex].calor.colunaDia = state.gasCells[tubingCellIndex].duto.dia;
        state.cells[annulusCellIndex].calor.geom.b = state.gasCells[tubingCellIndex].calor.geom.a;

        state.cells[annulusCellIndex].calor.Textern1 = state.gasCells[tubingCellIndex].temp;
        state.cells[annulusCellIndex].calor.Vextern1 = state.gasCells[tubingCellIndex].VGasR / state.gasCells[tubingCellIndex].u1L;
        // Ahead of the interface the annulus carries gas; behind it, completion
        // fluid. Same five fields either way, different fluid answering.
        const CelG &tubingCell = state.gasCells[tubingCellIndex];
        const double pres = tubingCell.pres;
        const double temp = tubingCell.temp;
        // The four property calls stay in the order the original made them --
        // specific heat, density, conductivity, viscosity -- rather than the
        // order the struct lists its fields. A braced initialiser would have
        // reordered them, and nothing here proves these correlations are free
        // of side effects.
        ExternalFluidProperties fluid;
        if (tubingCellIndex < state.interfaceCell) {
            const double specificHeat = tubingCell.flui.CalorGas(pres, temp);
            const double density = tubingCell.flui.MasEspGas(pres, temp);
            fluid = {tubingCell.flui.CondGas(pres, temp), specificHeat, density,
                     tubingCell.flui.ViscGas(pres, temp)};
        } else {
            const double specificHeat = tubingCell.CalorLiq(pres, temp);
            const double density = tubingCell.MasEspFlu(pres, temp);
            fluid = {tubingCell.CondLiq(pres, temp), specificHeat, density,
                     tubingCell.VisFlu(pres, temp)};
        }
        publishExternalFluid(state.cells[annulusCellIndex], fluid);
    }
}

void solveGasLine(const GasLiftState &state) {
    if (state.input.lingas <= 0)
        return;

    {
        double gasTemperature;
        // dt,celula,ColunaAnulaIni,ColunaAnulaFim,
        updateTransientGasValves(state);
        for (int valveIndex = 0; valveIndex < state.input.nvalvgas; valveIndex++) {
            int gasLiftProductionCell = state.productionValveCellIndices[valveIndex];
            int gasLiftGasCell = state.gasValveCellIndices[valveIndex];
            // The same three-term test decided both blocks below. Named once.
            const bool valveDischargesGas =
                gasLiftGasCell < state.interfaceCell ||
                (gasLiftGasCell == state.interfaceCell && state.gasCells[gasLiftGasCell].razInter > 0.5);
            if (valveDischargesGas) {
                state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH * kSecondsPerDay / state.gasCells[gasLiftGasCell].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius);
                state.cells[gasLiftProductionCell].acsr.injg.tipoflu = 0;
            } else {
                state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH;
                state.cells[gasLiftProductionCell].acsr.injg.tipoflu = 1;
            }
            if (valveDischargesGas) {
                if (state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg) {
                    gasTemperature = state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex);
                } else
                    gasTemperature = state.cells[gasLiftProductionCell].temp;
                state.cells[gasLiftProductionCell].acsr.injg.temp = gasTemperature;
                if (state.cells[gasLiftProductionCell].acsr.injg.temp < kMinimumTemperatureCelsius)
                    state.cells[gasLiftProductionCell].acsr.injg.temp = kMinimumTemperatureCelsius;
            } else
                state.cells[gasLiftProductionCell].acsr.injg.temp = state.gasLiftChokes[valveIndex].tempEstag;
            state.gasCells[gasLiftGasCell].pEstag = state.gasLiftChokes[valveIndex].presEstag;
            state.gasCells[gasLiftGasCell].tEstag = state.gasLiftChokes[valveIndex].tempEstag;
            state.gasCells[gasLiftGasCell].pGarg = state.gasLiftChokes[valveIndex].presGarg;
            state.gasCells[gasLiftGasCell].tGarg = state.gasLiftChokes[valveIndex].tempGarg;
            state.gasCells[gasLiftGasCell].qGarg = state.gasLiftChokes[valveIndex].qGarg;
            state.gasCells[gasLiftGasCell].areaGarg = state.gasLiftChokes[valveIndex].areagarg;
        }
        advanceGasSubStep(state);
    }
}

void connectTubingSteady(const GasLiftState &state) {
    for (int annulusCellIndex = state.annulusTubingStart; annulusCellIndex >= state.annulusTubingEnd; annulusCellIndex--) {
        int tubingCellIndex = state.annulusTubingStart + state.tubingAnnulusStart - annulusCellIndex;

        state.cells[annulusCellIndex].calor.Textern2 = state.gasCells[tubingCellIndex].calor.Textern1;
        state.cells[annulusCellIndex].calor.colunaDia = state.gasCells[tubingCellIndex].duto.dia;

        const CelG &tubingCell = state.gasCells[tubingCellIndex];
        const double pres = tubingCell.pres;
        const double temp = tubingCell.temp;
        const double gasSpecificHeat = tubingCell.flui.CalorGas(pres, temp);
        const double gasDensity = tubingCell.flui.MasEspGas(pres, temp);
        state.cells[annulusCellIndex].calor.Textern1 = tubingCell.calor.Textern1;
        state.cells[annulusCellIndex].calor.Vextern1 = tubingCell.VGasR / tubingCell.u1L;
        publishExternalFluid(state.cells[annulusCellIndex],
                             {tubingCell.flui.CondGas(pres, temp), gasSpecificHeat, gasDensity,
                              tubingCell.flui.ViscGas(pres, temp)});
    }
}

void initializeTubingConnectionSteady(const GasLiftState &state) {
    for (int annulusCellIndex = state.annulusTubingStart; annulusCellIndex >= state.annulusTubingEnd; annulusCellIndex--) {
        int tubingCellIndex = state.annulusTubingStart + state.tubingAnnulusStart - annulusCellIndex;
        state.cells[annulusCellIndex].calor.Textern2 = state.gasCells[tubingCellIndex].calor.Textern2;
        state.cells[annulusCellIndex].calor.colunaDia = state.gasCells[tubingCellIndex].duto.dia;

        const CelG &tubingCell = state.gasCells[tubingCellIndex];
        // Not the tubing cell's own state: the annulus pressure and the deck's
        // external temperature. That difference is why this cannot share the
        // steady connection above.
        const double pres = state.cells[annulusCellIndex].pres;
        const double temp = state.input.celg[tubingCellIndex].textern;
        const double gasSpecificHeat = tubingCell.flui.CalorGas(pres, temp);
        const double gasDensity = tubingCell.flui.MasEspGas(pres, temp);
        state.cells[annulusCellIndex].calor.Textern1 = tubingCell.calor.Textern1;
        state.cells[annulusCellIndex].calor.Vextern1 = 1.;
        // The viscosity is a fixed figure here, not a correlation.
        publishExternalFluid(state.cells[annulusCellIndex],
                             {tubingCell.flui.CondGas(pres, temp), gasSpecificHeat, gasDensity,
                              kInitialAnnulusGasViscosityCentipoise});
    }
}

double steadyGasPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.gasCells[cellIndex].dx0;
    double diameter = state.gasCells[cellIndex].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.gasCells[cellIndex].duto.peri;
    double gasDensity = state.gasCells[cellIndex].flui.MasEspGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double meanLocalGasFlow;
    meanLocalGasFlow = state.gasCells[cellIndex - 1].VGasR;

    double vel1 = meanLocalGasFlow / state.gasCells[cellIndex].u1L;

    double visc = state.gasCells[cellIndex].flui.ViscGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double reynolds;
    reynolds = state.gasCells[cellIndex].Rey(
        characteristicDiameter(state.gasCells[cellIndex].duto.revest, state.gasCells[cellIndex].duto.a, area, perimeter),
        vel1, gasDensity, visc);
    double frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    double frictionGradient = state.gasCells[cellIndex].dPdLFric * (0.5 * frictionFactor * gasDensity * (fabs(vel1) * vel1) * perimeter * dx / area);
    double hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * (kGravity * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);

    double pressureDrop = (frictionGradient + hydrostaticGradient) / kPascalPerKgfPerCm2;
    double meanPressure = state.gasCells[cellIndex].pres + pressureDrop;

    double meanTemperature;
    meanTemperature = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].temp + state.gasCells[cellIndex].dxL * state.gasCells[cellIndex].tempL) / (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex].dxL);

    dx = 0.5 * state.gasCells[cellIndex].dxL;
    diameter = state.gasCells[cellIndex].dutoL.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.gasCells[cellIndex].dutoL.peri;
    gasDensity = state.gasCells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);

    meanLocalGasFlow = state.gasCells[cellIndex - 1].VGasR;
    vel1 = meanLocalGasFlow / (gasDensity * area);

    if (cellIndex > 0)
        visc = state.gasCells[cellIndex - 1].flui.ViscGas(meanPressure, meanTemperature);
    else
        visc = state.gasCells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);

    reynolds = state.gasCells[cellIndex].Rey(
        characteristicDiameter(state.gasCells[cellIndex].dutoL.revest, state.gasCells[cellIndex].dutoL.a, area, perimeter),
        vel1, gasDensity, visc);
    double dpFLoc;
    double dpHLoc;
    if (cellIndex > 0) {
        frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
        dpFLoc = state.gasCells[cellIndex - 1].dPdLFric;
        dpHLoc = state.gasCells[cellIndex - 1].dPdLHidro;
    } else {
        frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
        dpFLoc = state.gasCells[cellIndex].dPdLFric;
        dpHLoc = state.gasCells[cellIndex].dPdLHidro;
    }
    frictionGradient = dpFLoc * (0.5 * frictionFactor * gasDensity * (fabs(vel1) * vel1) * perimeter * dx / area);
    hydrostaticGradient = dpHLoc * (kGravity * sin(state.gasCells[cellIndex].dutoL.teta) * gasDensity * dx);
    pressureDrop += (frictionGradient + hydrostaticGradient) / kPascalPerKgfPerCm2;
    return pressureDrop;
}

double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.cells[cellIndex].dx;
    double diameter = state.cells[cellIndex].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.cells[cellIndex].duto.peri;
    double meanTemperature;
    meanTemperature = (state.cells[cellIndex].dx * state.cells[cellIndex].temp + state.cells[cellIndex].dxL * state.cells[cellIndex - 1].temp) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
    double completionFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presaux, meanTemperature);

    double vel1 = state.cells[cellIndex - 1].QL / (area);

    double visc = state.cells[cellIndex].fluicol.VisFlu(state.cells[cellIndex].presaux, meanTemperature);

    double reynolds;
    reynolds = state.cells[cellIndex].Rey(
        characteristicDiameter(state.cells[cellIndex].duto.revest, state.cells[cellIndex].duto.a, area, perimeter),
        vel1, completionFluidDensity, visc);
    double frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / diameter);
    double frictionGradient = state.cells[cellIndex].dPdLFric * 0.5 * frictionFactor * completionFluidDensity * (fabs(vel1) * vel1) * perimeter * dx / area;
    double hydrostaticGradient = state.cells[cellIndex].dPdLHidro * kGravity * sin(state.cells[cellIndex].duto.teta) * completionFluidDensity * dx;

    double pressureDrop = (frictionGradient + hydrostaticGradient) / kPascalPerKgfPerCm2;
    double meanPressure = state.cells[cellIndex - 1].pres;

    dx = 0.5 * state.cells[cellIndex].dxL;
    diameter = state.cells[cellIndex].dutoL.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.cells[cellIndex].dutoL.peri;
    meanTemperature = state.cells[cellIndex - 1].temp;
    completionFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);

    if (cellIndex > 1)
        vel1 = state.cells[cellIndex - 2].QL / (area);

    visc = state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);

    reynolds = state.cells[cellIndex].Rey(
        characteristicDiameter(state.cells[cellIndex].dutoL.revest, state.cells[cellIndex].dutoL.a, area, perimeter),
        vel1, completionFluidDensity, visc);
    frictionFactor = state.cells[cellIndex - 1].fric(reynolds, state.cells[cellIndex].dutoL.rug / diameter);
    frictionGradient = state.cells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * completionFluidDensity * (fabs(vel1) * vel1) * perimeter * dx / area;
    hydrostaticGradient = state.cells[cellIndex - 1].dPdLHidro * kGravity * sin(state.cells[cellIndex].dutoL.teta) * completionFluidDensity * dx;
    pressureDrop += (frictionGradient + hydrostaticGradient) / kPascalPerKgfPerCm2;
    return pressureDrop;
}

void updateSteadyGasPressure(const GasLiftState &state, int cellIndex) {

    double dx = 0.5 * state.gasCells[cellIndex].dxL;
    double diameter = state.gasCells[cellIndex].dutoL.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.gasCells[cellIndex].dutoL.peri;
    double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

    double meanUpstreamGasFlow;
    meanUpstreamGasFlow = state.gasCells[cellIndex - 1].VGasR;
    double vel1 = meanUpstreamGasFlow / state.gasCells[cellIndex - 1].u1L;

    double visc = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

    double reynolds;
    reynolds = state.gasCells[cellIndex - 1].Rey(
        characteristicDiameter(state.gasCells[cellIndex].dutoL.revest, state.gasCells[cellIndex].dutoL.a, area, perimeter),
        vel1, gasDensity, visc);
    double frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
    double frictionGradient = state.gasCells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * gasDensity * (fabs(vel1) * vel1) * perimeter * dx / area;
    double hydrostaticGradient = state.gasCells[cellIndex - 1].dPdLHidro * kGravity * sin(state.gasCells[cellIndex].dutoL.teta) * gasDensity * dx;
    state.gasCells[cellIndex - 1].termoFric = frictionGradient / dx;
    state.gasCells[cellIndex - 1].termoHidro = hydrostaticGradient / dx;

    double meanPressure = state.gasCells[cellIndex - 1].pres - (frictionGradient + hydrostaticGradient) / kPascalPerKgfPerCm2;

    // The branch on steadyIteration that used to stand here assigned a weighted
    // mean in one arm and the upstream temperature in the other, and the next
    // line then overwrote it with the upstream temperature unconditionally.
    // Neither arm reached anything.
    double meanTemperature = state.gasCells[cellIndex - 1].temp;

    dx = 0.5 * state.gasCells[cellIndex].dx0;
    diameter = state.gasCells[cellIndex].duto.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.gasCells[cellIndex].duto.peri;
    gasDensity = state.gasCells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
    double areaVariation = 1 / pow(state.gasCells[cellIndex].dutoL.area, 2.) - 1 / pow(state.gasCells[cellIndex].duto.area, 2.);
    double dynamicPressureTerm = 0.5 * meanUpstreamGasFlow * meanUpstreamGasFlow * areaVariation / gasDensity;

    meanUpstreamGasFlow = state.gasCells[cellIndex - 1].VGasR;
    vel1 = meanUpstreamGasFlow / (gasDensity * area);

    visc = state.gasCells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);

    reynolds = state.gasCells[cellIndex].Rey(
        characteristicDiameter(state.gasCells[cellIndex].duto.revest, state.gasCells[cellIndex].duto.a, area, perimeter),
        vel1, gasDensity, visc);
    frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    frictionGradient = state.gasCells[cellIndex].dPdLFric * 0.5 * frictionFactor * gasDensity * (fabs(vel1) * vel1) * perimeter * dx / area;
    hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * kGravity * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx;

    state.gasCells[cellIndex].pres = meanPressure - (frictionGradient + hydrostaticGradient - dynamicPressureTerm) / kPascalPerKgfPerCm2;
    state.gasCells[cellIndex - 1].presR = state.gasCells[cellIndex].pres;
    if (cellIndex < state.gasCellCount)
        state.gasCells[cellIndex + 1].presL = state.gasCells[cellIndex].pres;
    state.gasCells[cellIndex].presini = state.gasCells[cellIndex].pres;
}

void computeSteadyGasFlowRate(const GasLiftState &state, int cellIndex) {

    int valveCount = state.input.nvalvgas;
    int match = 0;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        if (state.gasValveCellIndices[valveIndex] == cellIndex) {
            match = 1;
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
            double gasTemperature;
            state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH * kSecondsPerDay /
                                            (kAirDensityAtStandardConditions * state.gasCells[gasLiftGasCell].flui.Deng); // celulaG[gasLiftGasCell].flui.MasEspGas(1.,15.);
            if (state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg) {
                gasTemperature = state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex);
            } else
                gasTemperature = state.cells[gasLiftProductionCell].temp;
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

            state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
            if (cellIndex < state.gasCellCount)
                state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
        }
    }
    if (match == 0) {
        if (cellIndex > 0)
            state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR - 0 * state.gasCells[cellIndex].massfonteCH;
        else
            state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR + 0 * state.gasCells[cellIndex].massfonteCH;
        state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
        if (cellIndex < state.gasCellCount)
            state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
    }
}

void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex) {
    int valveCount = state.input.nvalvgas;
    double valveFlowRate;
    double casingPressure = state.initialGasPressure;
    double casingTemperature = state.gasCells[0].calor.Textern1;
    for (int valveIndex = 0; valveIndex < valveCount; valveIndex++) {
        int gasLiftProductionCell = state.productionValveCellIndices[valveIndex];
        int gasLiftGasCell = state.gasValveCellIndices[valveIndex];
        if (cellIndex == gasLiftProductionCell) {
            // Only how the mass source is obtained differs. The three
            // statements that publish it were identical in both arms.
            if (state.gasCells[0].tipoCC == 1) {
                valveFlowRate = (state.input.gasinj.vazgas[0] * state.gasCells[0].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius) / kSecondsPerDay) / valveCount;
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = valveFlowRate;
            } else {
                casingPressure = state.gasCells[gasLiftGasCell].pres;
                state.gasLiftChokes[valveIndex].presEstag = casingPressure;
                state.gasLiftChokes[valveIndex].presGarg = state.cells[state.productionValveCellIndices[valveIndex]].pres;
                state.gasLiftChokes[valveIndex].tempEstag = casingTemperature;
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH =
                    state.cells[state.productionValveCellIndices[valveIndex]].pres < casingPressure
                        ? state.gasLiftChokes[valveIndex].massica()
                        : 0.;
            }
            state.cells[gasLiftProductionCell].acsr.injg.QGas = state.gasCells[gasLiftGasCell].massfonteCH * kSecondsPerDay /
                                            state.gasCells[gasLiftGasCell].flui.MasEspGas(kStandardPressureKgfPerCm2, kStandardTemperatureCelsius);
            state.cells[gasLiftProductionCell].acsr.injg.temp = casingTemperature;
            state.cells[gasLiftProductionCell].fontemassGR = state.gasCells[gasLiftGasCell].massfonteCH;
        }
    }
}

void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex) {
    if (state.thermalSourceDisabled == 0) {
        if (((cellIndex) < state.tubingAnnulusStart || (cellIndex) > state.tubingAnnulusEnd)) {
            double dx = state.gasCells[cellIndex].dx0;
            double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double meanTemperatureGradient = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].dTdLCor + state.gasCells[cellIndex - 1].dx0 * state.gasCells[cellIndex - 1].dTdLCor) / (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double area = state.gasCells[cellIndex].duto.area;
            double meanSuperficialGasVelocity;
            double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            meanSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double gasJouleThomson = state.gasCells[cellIndex - 1].flui.JTG(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double hidro = (gasDensity * meanSuperficialGasVelocity) * area * kGravity * sin(state.gasCells[cellIndex].duto.teta);

            state.gasCells[cellIndex - 1].calor.Tint = state.gasCells[cellIndex - 1].temp;
            state.gasCells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity;
            state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex - 1].calor.rhoint = gasDensity;
            state.gasCells[cellIndex - 1].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp) * kPascalSecondPerCentipoise;
            state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
            double upstreamHeatFlux = state.gasCells[cellIndex - 1].fluxcal;
            double tubingHeatFlux = 0.;

            double lengthRatio;
            if (cellIndex > 0)
                lengthRatio = dx / (dx + state.gasCells[cellIndex - 1].dx0);
            else
                lengthRatio = 0.5;
            double temperatureMarchCoefficient = gasDensity * meanSuperficialGasVelocity * gasSpecificHeat * area;
            double pressureMarchCoefficient = gasDensity * meanSuperficialGasVelocity * gasJouleThomson * area;
            double dpdx;
            dpdx = 2. * (state.gasCells[cellIndex].pres - ((1 - lengthRatio) * state.gasCells[cellIndex - 1].pres + lengthRatio * state.gasCells[cellIndex].pres)) * 98600. / dx;

            double kineticTerm;
            double superficialVelocityGradient = 0.;
            double mixtureDensity = gasDensity;
            double upstreamSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasL / state.gasCells[cellIndex - 1].u1L;
            superficialVelocityGradient = (meanSuperficialGasVelocity - upstreamSuperficialGasVelocity) / dx;

            if (state.input.nCompTotalUnidadesG / dx < 1e6)
                kineticTerm = mixtureDensity * area * meanSuperficialGasVelocity * meanSuperficialGasVelocity * superficialVelocityGradient;
            else
                kineticTerm = 0.;

            double gasMassSource = 0.;
            double liquidMassSource = 0.;


            if (meanSuperficialGasVelocity > 1e-3 && (state.gasCells[cellIndex].duto.a / state.gasCells[cellIndex - 1].duto.a > 0.5 &&
                                  state.gasCells[cellIndex - 1].duto.a / state.gasCells[cellIndex].duto.a > 0.5)) {
                double energyTerm1 = meanTemperatureGradient * (pressureMarchCoefficient * dpdx - kineticTerm - hidro + liquidMassSource + gasMassSource) / temperatureMarchCoefficient;
                double energyTerm2 = meanTemperatureGradient * (upstreamHeatFlux + tubingHeatFlux) / temperatureMarchCoefficient;
                int stepCount;
                double stepLength;
                if (meanCellLength / state.gasCells[cellIndex - 1].calor.resGlob < 1000.) {
                    stepCount = 0.;
                    stepLength = meanCellLength;
                } else {
                    stepCount = 2 * (meanCellLength / state.gasCells[cellIndex - 1].calor.resGlob) / 1000 + 1;
                    stepLength = meanCellLength / stepCount;
                }
                double stepTemperature = state.gasCells[cellIndex - 1].temp;

                stepTemperature = stepLength * (state.gasCells[cellIndex - 1].temp / stepLength + energyTerm1 + energyTerm2);
                for (int stepIndex = 1; stepIndex < stepCount; stepIndex++) {
                    state.gasCells[cellIndex - 1].calor.Tint = stepTemperature;
                    state.gasCells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity;
                    state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, stepTemperature);
                    state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
                    state.gasCells[cellIndex - 1].calor.rhoint = gasDensity;
                    state.gasCells[cellIndex - 1].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, stepTemperature) * kPascalSecondPerCentipoise;
                    state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
                    upstreamHeatFlux = state.gasCells[cellIndex - 1].fluxcal;
                    energyTerm2 = meanTemperatureGradient * (upstreamHeatFlux + tubingHeatFlux) / temperatureMarchCoefficient;
                    stepTemperature = stepLength * (stepTemperature / stepLength + energyTerm1 + energyTerm2);
                }

                state.gasCells[cellIndex].temp = stepTemperature;
            } else {
                state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
            }

            clampAndPropagateTemperature(state, cellIndex);
        } else {
            if (state.networkCoupled == 1) {
                int connectionCellIndex = (cellIndex)-state.tubingAnnulusStart;
                int connectionIndex = state.annulusTubingStart - connectionCellIndex;
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

            double meanSuperficialGasVelocity;
            meanSuperficialGasVelocity = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
            double gasDensity = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

            state.gasCells[cellIndex].calor.Tint = state.gasCells[cellIndex].temp;
            state.gasCells[cellIndex].calor.Vint = meanSuperficialGasVelocity;
            state.gasCells[cellIndex].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            state.gasCells[cellIndex].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex].calor.rhoint = gasDensity;
            state.gasCells[cellIndex].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp) * kPascalSecondPerCentipoise;
            state.gasCells[cellIndex].fluxcal = state.gasCells[cellIndex].calor.transperm();
        }
    } else {
        state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
    }
}

}  // namespace sisprod::gaslift
