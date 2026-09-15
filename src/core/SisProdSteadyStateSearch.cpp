#include "SisProdSteadyStateSearch.h"

#include "Leitura.h"
#include "celula3.h"
#include "celulaGas.h"
#include "chokegas.h"
#include "variaveisGlobais1D.h"

#include <math.h>

namespace sisprod::steady {

namespace {

/// The nine rows of the steady-state march dispatch table.
///
/// Named after the methods they select, deliberately: a dispatch table is read
/// by checking that each row goes where it says, and `case
/// SteadyMarch::marchaProdPresPres1Rev: return marchaProdPresPres1Rev(chute);`
/// makes a mis-wired row visible without cross-referencing anything.
enum class SteadyMarch {
    marchaInjPerm1,
    marchaGasPerm2,
    marchaGasPerm3,
    marchaProdPerm1,
    marchaProdPerm1Rev,
    marchaProdPerm2,
    marchaProdPresPres1,
    marchaProdPresPres1Rev,
    marchaProdPresPres2,
};

/// Resolves the four selectors to one row of the table.
///
/// This was a conditional chain nested four deep inside multMarcha. Pulled out,
/// it is a pure function of five values -- which is what makes the table
/// verifiable: refactor-harness/verify-dispatch.py sweeps every combination of
/// the selectors against the original chain carved out of the baseline commit.
/// Nothing else can check this. Measured over the 360 zriddr calls the corpus
/// makes, only three (prod, tipoCC) pairs ever occur -- (0,0), (1,0) and (1,1) --
/// so the corpus reaches at most FOUR of the nine rows. The injection row, the
/// second gas-line row and all three pressure-pressure rows never run, and a
/// mis-wired row among those five would leave every gate green.
///
/// Documented row by row in evidencia/tabela-despacho.md.
/// productionChokeOpening is the array, not the value, so that the subscript
/// happens only in the branch that needs it -- as it did when this was a nested
/// chain. Ler::copia_chokeSup leaves chokep.abertura null when parserie is not
/// positive, and reading it on every dispatch would turn a conditional
/// dereference into an unconditional one. No corpus model takes that path, so
/// no gate would have said anything.
SteadyMarch selectSteadyMarch(int injectorWell, int prod, int tipoCC, int reverseMarch,
                              const double *productionChokeOpening) {
    if (injectorWell != 0)
        return SteadyMarch::marchaInjPerm1;
    if (prod == 0)
        return tipoCC == 0 ? SteadyMarch::marchaGasPerm2 : SteadyMarch::marchaGasPerm3;
    if (prod == 1) {
        if (tipoCC != 0)
            return SteadyMarch::marchaProdPerm2;
        return reverseMarch == 0 ? SteadyMarch::marchaProdPerm1
                                 : SteadyMarch::marchaProdPerm1Rev;
    }
    if (tipoCC != 0) {
        // A2-01. Both arms select the same march, and that is correct, not a
        // bug -- which took measuring to establish. The obvious reading is that
        // the else should reach marchaProdPresPres3, since that function exists,
        // has no caller, and Num4Main splits on this very condition to choose
        // buscaProdPresPresPerm3. It should not. marchaProdPresPres2 already
        // reduces to marchaProdPresPres3 when the choke is shut: the throat area
        // is `abertura[0] * area`, so vazmaxSachd and vazmassSachd both go to
        // zero with it, and the residual becomes the same `0. - MR` that
        // marchaProdPresPres3 returns literally. Routing here to
        // marchaProdPresPres3 would swap a guarded, general march for a narrower
        // ancestor -- a regression wearing the shape of a fix.
        //
        // The branch is kept rather than collapsed because it is the only
        // surviving record that the two cases were once distinct. See A2-01 in
        // evidencia/anomalias.md for the measurement.
        return productionChokeOpening[0] > 1e-15 ? SteadyMarch::marchaProdPresPres2
                                                 : SteadyMarch::marchaProdPresPres2;
    }
    return reverseMarch == 0 ? SteadyMarch::marchaProdPresPres1
                             : SteadyMarch::marchaProdPresPres1Rev;
}

}  // namespace

double dispatchMarch(const SteadyStateSearchState &state, double chute, int prod, int tipoCC) {
    switch (selectSteadyMarch(state.march.input.pocinjec, prod, tipoCC, state.reverseSteady,
                              state.march.input.chokep.abertura)) {
    case SteadyMarch::marchaGasPerm2:         return marchGasSteadySecondary(state.march, chute);
    case SteadyMarch::marchaGasPerm3:         return marchGasSteadyTertiary(state.march, chute);
    case SteadyMarch::marchaProdPerm1:        return marchProductionSteady(state.march, chute);
    case SteadyMarch::marchaProdPerm1Rev:     return marchReverseProductionSteady(state.march, chute);
    case SteadyMarch::marchaProdPerm2:        return marchProductionSteadySecondary(state.march, chute);
    case SteadyMarch::marchaProdPresPres1:    return marchProductionPressureToPressure(state.march, chute);
    case SteadyMarch::marchaProdPresPres1Rev: return marchReverseProductionPressureToPressure(state.march, chute);
    case SteadyMarch::marchaProdPresPres2:    return marchProductionPressureToPressureSecondary(state.march, chute);
    // Falls out of the switch rather than returning inside it. An exhaustive
    // switch over a scoped enum still trips -Wreturn-type on GCC 11, and gate 1
    // admits no new warnings; a default: label would add an unreachable path
    // the original did not have.
    case SteadyMarch::marchaInjPerm1:         break;
    }
    return marchInjectionSteady(state.march, chute);
}

double solveSteadyRoot(const SteadyStateSearchState &state, double x1, double x2, int prod, int tipoCC) {
    // Hoisted out of the solver: arq is input-deck configuration, and a
    // generic root finder has no business reading it. minit gates three early
    // returns inside the solver; see A2-05 in evidencia/anomalias.md for what
    // that gating reaches.
    int minit=0;
    if(state.march.input.acopColAnulPermForte == 1)minit=10;
    return rootfinding::zriddr(
        x1, x2,
        [&](double guess) { return dispatchMarch(state, guess, prod, tipoCC); },
        // Domain feedback, kept on the domain side: the solver composes it as
        // monitor(objective(x)), which is the order the original evaluated in.
        [&](double residual) {
            double normalized = residual / state.march.baseConvergenceMonitor;
            if (prod != 0) {
                state.march.convergenceMonitor = fabs(normalized);
            }
            return normalized;
        },
        state.reverseSteady, minit);
}

}  // namespace sisprod::steady
