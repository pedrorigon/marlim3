# General

Global simulation settings define the physics contract for the entire model: which physical effects are active, which correlations describe multiphase interactions, and how the numerical solver behaves.

All settings in this section live inside a global configuration object.

> **JSON key:** `initialConfig` (EN) · `configuracaoInicial` (PT)

## System Type

Every simulation starts by declaring whether the system is a **producer** (multiphase flow from reservoir toward a separator) or an **injector** (flow from surface into the formation).

> **JSON key:** `system` (EN) · `sistema` (PT)
> Values: `"PROD"` / `"INJ"` (EN) · `"MULTIFASICO"` / `"INJETOR"` (PT)

---

## Interphase Slip

In multiphase flow, gas and liquid travel at different velocities. This **slip** determines holdup, pressure gradient, and flow-pattern transitions. Marlim3 lets you enable or disable slip treatment independently for steady-state and transient solvers.

- When slip is active, the simulator uses a drift-flux formulation to compute the velocity difference between phases.
- When slip is off, gas and liquid share the same velocity (homogeneous model), which is simpler but less physical.

> **JSON keys:**
>
> - Steady-state slip: `steadyStateSlip` (EN) · `escorregamentoPermanente` (PT) — default `true`
> - Transient slip: `transientSlip` (EN) · `escorregamentoTransiente` (PT) — default `true`

## Drift-Flux vs. Black-Box Correlations

For steady-state calculations, two modeling philosophies are available:

1. **Drift-flux model** — Mechanistic approach that resolves slip using physical parameters (distribution coefficient, drift velocity). More general and physically grounded.
2. **Black-box correlations** — Empirical correlations that directly predict pressure gradient and holdup without resolving the underlying velocity profile. Less general, but calibrated for specific conditions.

In transient mode, the drift-flux model is always used regardless of this setting.

> **JSON key:** `driftModel` (EN) · `modeloDrift` (PT) — default `true`
> `true` = drift-flux; `false` = black-box correlations (steady-state only)

## Flow-Pattern Map

The flow pattern (stratified, slug, bubble, annular, churn) governs which closure laws apply at each point in the pipe. Marlim3 uses the Barnea flow-pattern map in two variants:

- **Simplified Barnea (0):** Computationally lighter, adequate for most engineering applications.
- **Full Barnea (1):** More detailed transition criteria, useful when flow-regime sensitivity is important.

> **JSON key:** `flowPatternMap` (EN) · `mapaArranjo` (PT) — default `0`

## Slip Correlations by Flow Pattern

Once the flow pattern is identified, specific correlations compute the slip (drift velocity and distribution parameter). Different correlations can be assigned to each regime family:

| Flow-pattern family | Available correlations |
|--------------------|-----------------------|
| Stratified | `0`: Choi et al · `1`: Bhagwat & Ghajar (default) · `2`: Franca & Lahey · `4`: modified Bhagwat & Ghajar |
| Slug / Bubble | `0`: Choi et al · `1`: Bhagwat & Ghajar (default) · `4`: modified Bhagwat & Ghajar |
| Annular / Churn | `0`: Choi et al · `1`: Bhagwat & Ghajar (default) · `3`: Hibiki & Ishii · `4`: modified Bhagwat & Ghajar |

> **JSON key:** `correlationsByPattern` (EN) · `correlacoesPorArranjo` (PT)
> Sub-keys: `stratified` / `slugBubble` / `annularChurn`

## Interphase Mass Transfer

Gas and liquid exchange mass via evaporation, condensation, and flashing. The mass-transfer model determines how quickly and accurately this exchange is computed:

| Option | Behavior |
|--------|----------|
| `0` | Full model, implicit method — most stable (recommended) |
| `1` | Full model, explicit method — same physics, less stable |
| `2` | Simplified isothermal — neglects transient mass-transfer terms |
| `3` | No mass transfer — phases are frozen |

> **JSON key:** `massTransfer` (EN) · `transferenciaMassa` (PT) — default `0`

---

## Geometry Reference Direction

The **geometry direction** defines whether measured-length zero corresponds to the pipe inlet (where flow enters) or the pipe outlet. This affects how pipe segments, accessories, and sources are positioned.

- When geometry follows flow, zero is the inlet and length increases downstream.
- When reversed, zero is the outlet. Angles are still defined with respect to the flow direction for backward compatibility.

> **JSON key:** `geometryFollowsFlow` (EN) · `sentidoGeometriaSegueEscoamento` (PT) — default `true`

## Coordinate-Based Inclination (XY Mode)

Instead of specifying inclination angles directly for each pipe segment, you can provide endpoint XY coordinates and let Marlim3 compute the inclinations. This is natural when working from well survey data or pipeline route coordinates.

When XY mode is active, you must also provide the starting coordinates of the production line (and service line, if present).

> **JSON keys:**
>
> - Enable: `xyMode` (EN) · `modoXY` (PT) — default `false`
> - Production start: `xProdStart` / `yProdStart` (EN) · `xProdInicio` / `yProdInicio` (PT)
> - Service start: `xServiceStart` / `yServiceStart` (EN) · `xServInicio` / `yServInicio` (PT)

## Service (Gas) Line

Some production systems include a parallel service line for gas injection (gas-lift, gas circulation). Enabling this feature tells the simulator to expect service-pipe geometry and gas-injection objects.

> **JSON key:** `gasLine` (EN) · `linhaGas` (PT) — default `false`

---

## Thermal Initialization

When starting from a user-defined initial condition, the simulator needs a strategy for distributing temperature across the pipe wall layers:

- **Thermal equilibrium (true):** Linearly interpolates between ambient and fluid temperature across the wall cross-section.
- **Cold start (false):** Uses ambient temperature throughout the entire wall.

> **JSON key:** `thermalEquilibrium` (EN) · `equilibrioTermico` (PT) — default `true`

## Condensation Latent Heat

In the energy equation, condensation releases latent heat. In rare cases of extreme retrograde condensation (exotic compositional fluids), this term can cause numerical instability and can be disabled.

> **JSON key:** `latentHeatCond` (EN) · `condlatente` (PT) — default `true`

## Reverse-Flow Temperature

During transient simulations, flow may temporarily reverse at the outlet boundary (e.g., pressure surges). When gas enters from the outlet, its temperature must be defined. If not specified, ambient temperature is assumed.

> **JSON key:** `reverseTemp` (EN) · `tempReves` (PT) — unit: °C

---

## Steady-State Numerical Order

The steady-state solver integrates ordinary differential equations along the pipe. Two methods are available:

- **First order (1):** Simpler, faster, adequate for most cases.
- **Second order (2):** Runge-Kutta method, better accuracy for steep gradients.

In transient mode, the method is always first order.

> **JSON key:** `steadyStateOrder` (EN) · `ordemperm` (PT) — default `1`

## Outlet Check Valve

A check valve at the production-system outlet prevents reverse gas inflow when the downstream pressure exceeds the computed outlet pressure. Without it, gas can enter the system from the outlet boundary.

> **JSON key:** `checkValve` (EN) · `CheckValve` (PT) — default `0`
> `0` = no check valve (reverse flow allowed); `1` = check valve active

## Steady-State Initial Guess

The steady-state solver uses a shooting method that needs an initial estimate (either flow rate or pressure, depending on boundary conditions). If not provided or set to a negative value, the simulator auto-computes a guess from hydrostatic balance.

> **JSON key:** `steadyGuess` (EN) · `chutePerm` (PT) — default `-1` (auto)

## Sensitivity Analysis Parallelization

When running parametric sweeps (sensitivity analysis), individual cases can be executed in parallel using multiple threads.

> **JSON key:** `parallelizeSA` (EN) · `paralelizaAS` (PT) — default `false`

---

## Performance Acceleration

### Precomputed Property Tables

Computing fluid compressibilities and derivatives at every iteration can dominate transient runtime. Pre-building interpolation tables before the simulation starts can significantly accelerate execution:

- **Production fluid table:** Pre-computes compressibility for all declared black-oil fluids using reduced pressure and temperature.
- **Gas table:** Same concept applied to the service-line gas.

> **JSON keys:**
>
> - Production: `pressureTable` (EN) · `tabP` (PT) — default `false`
> - Service gas: `gasTable` (EN) · `tabG` (PT) — default `false`

### Dynamic Compositional Table

For compositional steady-state simulations (especially networks), a posterior table can be built from an initial black-oil pass. This avoids repeated flash calculations at each iteration by reusing the P-T mapping from the first solution.

> **JSON key:** `dynamicTableModel` (EN) · `modeloTabelaDinamica` (PT) — default `false`

### Solubility-Ratio Table

Pre-building the solution-gas-ratio (RS) lookup table can improve performance when using computationally expensive RS models (e.g., Livia Fulchignoni's model).

> **JSON key:** `srBpTable` (EN) · `tabelaRSPB` (PT) — default `false`

### Stream Mixing Tracking

When multiple fluid streams merge inside the system (at sources), mixture properties must be updated using mixing rules. Two tracking options control which variables are updated:

- **GOR tracking:** Updates gas-oil ratio, API, and BSW at merging points.
- **Gas-density tracking:** Updates gas density and CO₂ fraction at merging points.

> **JSON keys:**
>
> - GOR: `trackGOR` (EN) · `trackRgo` (PT) — default `true`
> - Gas density: `trackGasDensity` (EN) · `trackDensidadeGas` (PT) — default `true`

---

## Advanced Expert Controls

These parameters fine-tune solver robustness and rarely need modification. They live inside a nested advanced configuration object.

> **JSON key:** `initialConfig.advanced` (EN) · `configuracaoInicial.avancado` (PT)

### Single-Phase Threshold

Very small void fractions (near-zero gas content) can cause numerical problems. Below a threshold, the solver treats the flow as pure single-phase liquid, avoiding instability.

> **JSON key:** `monophasicCriterion` (EN) · default `1e-5`

### Condensation Threshold

Similarly, mass transfer is disabled below a minimum void fraction to prevent unphysical condensation behavior that would violate holdup limits.

> **JSON key:** `condensationCriterion` (EN) · default `0.001`

### Cell-Boundary Pressure

Two methods exist for computing pressure at cell boundaries:

- **Simple (true):** Average of adjacent cell-center pressures. Fast and stable.
- **Full marching (false):** Hydrostatic + friction marching from cell center to boundary. More accurate but can introduce oscillations.

> **JSON key:** `simplePressureFrontier` (EN) · default `true`

### Mass-Transfer Rate Limiter

When interphase mass transfer (condensation or evaporation) exceeds a given rate, latent-heat terms in the energy equation are disabled to prevent numerical blow-up in extreme retrograde scenarios.

> **JSON key:** `massTransferLimit` (EN) · default `10 kg/(s·m)`

### Steady-State Convergence Accelerator

A simplification in the steady-state solver that accelerates convergence and improves stability. This is the simulator default and rarely needs to be disabled.

> **JSON key:** `accelerateSteadyConvergence` (EN) · default `true`

### Boundary-Cell Slip Control

The last control volume in the production system can exhibit oscillations due to rapid flow-pattern changes at low pressure. Disabling slip in that cell can stabilize the solution.

> **JSON key:** `slipBoundaryCell` (EN) · default `true`

### Time-Step Penalization Controls

Several mechanisms penalize (restrict) time-step growth to maintain stability:

- **Choke oscillation penalization:** Restricts growth when liquid flow oscillates through the surface choke.
- **Valve event control:** Limits step size during valve opening/closing ramps.
- **Holdup oscillation penalization:** Restricts step size during abrupt holdup changes (e.g., shutdown segregation). Can be disabled for cases where the penalization is overly conservative.

> **JSON keys:**
>
> - Choke: `relaxChokeTimestep` (EN) — default `false`
> - Valve: `valveTimestepControl` (EN) — default `false`
> - Disable holdup penalization: `disablePenalizeTimestep` (EN) — default `false`

### Sonic Time-Stepping

Pressure waves propagate at sonic velocity. To capture these events, the simulator can temporarily force time steps small enough to satisfy the sonic CFL condition. This is controlled by time windows:

> **JSON keys:** `sonicTime` (EN) + `sonicFlag` (EN)

### Full Compressibility Model Windows

The full model (including time derivatives of liquid density) can be selectively activated in time windows where strong depressurization occurs. Outside these windows, a simplified one-step evolution is used.

> **JSON keys:** `compModelCorrectionTime` (EN) + `compModelCorrectionFlag` (EN)
> Threshold for auto-selection: `despressRate` (EN) — default `0.01 kgf/(cm²·s)`

### Parallelism

Simulation threads and matrix-solve threads can be configured independently.

> **JSON keys:** `threads` (EN) — default `1`; `matrixThreads` (EN) — default `1`

### Flash Mini-Table Refresh

For transient compositional simulations, local mini-tables of P-T flash data can be periodically regenerated to accelerate property lookups:

> **JSON keys:** `dynTableMinDelay` (interval in time steps), `dynTableMinDp` (pressure increment, kgf/cm²), `dynTableMinDt` (temperature increment, °C)

---

## Example JSON

```json
{
  "system": "PROD",
  "initialConfig": {
    "steadyStateSlip": true,
    "transientSlip": true,
    "driftModel": true,
    "flowPatternMap": 0,
    "massTransfer": 0,
    "steadyStateOrder": 1,
    "checkValve": 0,
    "geometryFollowsFlow": true,
    "gasLine": false,
    "pressureTable": true,
    "correlationsByPattern": {
      "stratified": 1,
      "slugBubble": 1,
      "annularChurn": 1
    },
    "advanced": {
      "threads": 1,
      "massTransferLimit": 10.0,
      "simplePressureFrontier": true,
      "accelerateSteadyConvergence": true
    }
  }
}
```

!!! tip
    Keep advanced options at their defaults until a baseline case is validated. Then tune one control at a time, comparing results against the validated baseline.
