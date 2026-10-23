/// The connection between SProd and the modules: the state view each module
/// receives (sisprod::adapters), the solve context that builds the views of a call,
/// and the callbacks through which a module reaches another module's code with them.
#include "SisProd.h"
#include "SisProdSolveContext.h"
#include "SisProdFluidTracking.h"
#include "SisProdGasLift.h"
#include "SisProdSources.h"
#include "SisProdSteadyState.h"
#include "SisProdSteadyStateSearch.h"
#include "SisProdThermal.h"
#include "SisProdTransient.h"
#include "SisProdTrendOutput.h"







// ------------------------------------------- steady-state callbacks ----
//
// The steady march reaches the drift closures, the source terms, and the gas-line
// and thermal steps through these, with the views of the solve context: its own
// view does not carry their state.
















// The two searches a march calls; their results are discarded at the call
// site.
void sisprod::steady::SteadyStateUpdaters::searchGasPressureSteadySecondary() const {
    const sisprod::thermal::ThermalDependencies thermalDependencies{context.sources(), context.closure()};
    sisprod::steady::searchGasPressureSteadySecondary(context.search(), {context.sources(), context.closure(), context.thermal(), thermalDependencies, context.gasLift()});
}
void sisprod::steady::SteadyStateUpdaters::searchGasPressureSteadyTertiary() const {
    const sisprod::thermal::ThermalDependencies thermalDependencies{context.sources(), context.closure()};
    sisprod::steady::searchGasPressureSteadyTertiary(context.search(), {context.sources(), context.closure(), context.thermal(), thermalDependencies, context.gasLift()});
}











void sisprod::thermal::ThermalEvolutionUpdater::solvePressureVelocityCoupling(
    int cycle) const {
    sisprod::transient::solvePressureVolumeCoupling(context.transientStep(), 0, cycle);
}

void sisprod::thermal::ThermalEvolutionUpdater::renew() const {
    sisprod::transient::updateCells(context.transientStep(), 0);
}

void sisprod::transient::TransientStepUpdaters::advanceGasSubStep() const {
    sisprod::gaslift::advanceGasSubStep(context.gasLift(), context.thermal());
}





void sisprod::transient::TransientSolveUpdaters::solveHydrateEnvelopes() const {
    context.system().solveHydrateEnvelopes();
}










void sisprod::transient::TransientSolveUpdaters::marchTransientEnergy(int ciclo, int ciclomax) const {
    sisprod::thermal::advanceTransientEnergy(context.thermal(), ciclo, ciclomax);
}
