#ifndef SISPRODSOLVECONTEXT_H_
#define SISPRODSOLVECONTEXT_H_

#include "DriftFluxClosure.h"
#include "SisProdComposition.h"
#include "SisProdGasLift.h"
#include "SisProdSources.h"
#include "SisProdSteadyState.h"
#include "SisProdSteadyStateSearch.h"
#include "SisProdThermal.h"
#include "SisProdTransient.h"
#include "SisProdTrendOutput.h"

class SProd;

namespace sisprod {

/// What a call into the modules reaches one system through: the system, and its views.
///
/// The views are built once, with the context, before the call reaches any parallel loop, and
/// are only read afterwards, so the threads of a loop share them. A view holds references to the
/// system's members, never copies, so it is of the system as it is when read. The updaters of
/// the views hold the context, so a module that calls back reaches the system and its views
/// through it.
class SolveContext {
  public:
    explicit SolveContext(SProd &system);
    SolveContext(const SolveContext &) = delete;
    SolveContext &operator=(const SolveContext &) = delete;

    /// The system the views are of.
    SProd &system() const { return system_; }

    const sisprod::gaslift::GasLiftState &gasLift() const { return gasLift_; }
    const sisprod::steady::SteadyStateState &steady() const { return steady_; }
    const sisprod::steady::SteadyStateSearchState &search() const { return search_; }
    const sisprod::transient::TransientStepState &transientStep() const { return transientStep_; }
    const sisprod::composition::CompositionState &composition() const { return composition_; }
    const sisprod::transient::TransientSolveState &transientSolve() const { return transientSolve_; }
    const driftflux::coefficient::ClosureState &closure() const { return closure_; }
    const sisprod::sources::SourceState &sources() const { return sources_; }
    const sisprod::thermal::ThermalState &thermal() const { return thermal_; }
    const trendoutput::TrendState &trends() const { return trends_; }

  private:
    SProd &system_;
    sisprod::gaslift::GasLiftState gasLift_;
    sisprod::steady::SteadyStateState steady_;
    sisprod::steady::SteadyStateSearchState search_;
    sisprod::transient::TransientStepState transientStep_;
    sisprod::composition::CompositionState composition_;
    sisprod::transient::TransientSolveState transientSolve_;
    driftflux::coefficient::ClosureState closure_;
    sisprod::sources::SourceState sources_;
    sisprod::thermal::ThermalState thermal_;
    trendoutput::TrendState trends_;
};

}  // namespace sisprod

#endif  // SISPRODSOLVECONTEXT_H_
