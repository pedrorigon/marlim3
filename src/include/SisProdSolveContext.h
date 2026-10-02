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

/// What a call into the modules reaches one system through: the system, and its views. Each
/// accessor builds the view it returns, so the view is of the system as it is when asked for, and
/// threads that ask at once get one each. The updaters of the views it builds hold the context, so
/// a module that calls back reaches the system and its views through it.
class SolveContext {
  public:
    explicit SolveContext(SProd &system) : system_(system) {}
    SolveContext(const SolveContext &) = delete;
    SolveContext &operator=(const SolveContext &) = delete;

    /// The system the views are of.
    SProd &system() const { return system_; }

    sisprod::gaslift::GasLiftState gasLift();
    sisprod::steady::SteadyStateState steady();
    sisprod::steady::SteadyStateSearchState search();
    sisprod::transient::TransientStepState transientStep();
    sisprod::composition::CompositionState composition();
    sisprod::transient::TransientSolveState transientSolve();
    driftflux::coefficient::ClosureState closure();
    sisprod::sources::SourceState sources();
    sisprod::thermal::ThermalState thermal();
    trendoutput::TrendState trends();

  private:
    SProd &system_;
};

}  // namespace sisprod

#endif  // SISPRODSOLVECONTEXT_H_
