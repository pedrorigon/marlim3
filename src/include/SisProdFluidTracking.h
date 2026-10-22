#ifndef SISPRODFLUIDTRACKING_H_
#define SISPRODFLUIDTRACKING_H_

// The state holds references; its builders include the referenced definitions.
class Cel;
class Ler;
class ProFlu;
struct varGlob1D;
class SProd;
namespace sisprod { class SolveContext; }

namespace sisprod::fluidtracking {

/// Gas specific gravity correction through the steady-state view of the same system.
struct FluidTrackingUpdaters {
    SolveContext &context;

    void correctGasSpecificGravity(int i) const;
};

/// State shared by black-oil property tracking and phase molar-fraction transport.
///
/// Each transport reads its time step from the first cell.
///
/// Scalars are held by reference so each transport reads them only when its
/// selected fluid model and tracking options require them.
struct FluidTrackingState {
    Cel *&cells;                          // celula
    int &lastCell;                        // ncel
    Ler &input;                           // arq
    varGlob1D *&globals;                  // vg1dSP
    int &endNode;                         // noextremo
    double &inletPressure;                // presE
    double &inletTemperature;             // tempE
    double &inletComplementaryFraction;   // betaE
    int &trackGasOilRatio;                // trackRGO
    int &trackGasGravity;                 // trackDeng
    int &compositionalRefreshCounter;     // kontaRenovaComp

    FluidTrackingUpdaters updaters;
};

/// Transports gas-oil ratio, gas specific gravity, CO2 fraction, API gravity,
/// water content and density, and dead-oil viscosities as enabled by the fluid
/// model and tracking options. fluiRev is the fluid entering through the outlet
/// under reverse flow, taken by value.
void transportBlackOilProperties(const FluidTrackingState &state, ProFlu fluiRev);

/// Transports liquid and vapour molar fractions separately, together with the
/// water content and density and dead-oil viscosities. Refreshes the cell fluids
/// by flash. fluiRev is the outlet's reverse-flow fluid, taken by value.
void transportPhaseMolarFractions(const FluidTrackingState &state, ProFlu fluiRev);

}  // namespace sisprod::fluidtracking

#endif  // SISPRODFLUIDTRACKING_H_
