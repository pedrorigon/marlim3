# General

This section covers global simulation settings defined in the `initialConfig` object. These control which physical models are active, how equations are solved, and how the simulation behaves.

## Concept

`initialConfig` is the simulator's "physics and numerics contract". It defines which physical effects are active, which correlations are used, and how the governing equations are solved. All other objects (fluids, pipes, accessories) are interpreted in the context of these global settings.

## System Definition

| Field | Values | Meaning |
|-------|--------|---------|
| `system` | `PROD` / `INJ` | Producer or injector system configuration |

A producer system (`PROD`) solves multiphase flow from reservoir or source towards a separator. An injector system (`INJ`) solves single-phase or multiphase flow from surface towards the formation.

## Flow-Model Choices

| Field | Type | Default | Purpose |
|-------|------|---------|---------|
| `steadyStateSlip` | bool | `true` | Enables interphase-slip (drift) treatment in steady-state solver |
| `transientSlip` | bool | `true` | Enables interphase-slip treatment in transient solver |
| `driftModel` | bool | `true` | Uses drift-flux model; if `false`, uses black-box correlations (steady-state only — transient always uses drift-flux) |
| `flowPatternMap` | int | `0` | Flow-regime map: `0` = simplified Barnea; `1` = full Barnea |
| `massTransfer` | int | `0` | Interphase mass-transfer model: `0` = full implicit (most stable); `1` = full explicit; `2` = simplified isothermal; `3` = none |
| `correlationsByPattern` | object | — | Slip correlation family by flow regime (stratified, slug/bubble, annular/churn) |

### Slip Correlations by Flow Pattern

When drift-flux is active, different correlations can be selected per regime:

| Flow pattern | Options |
|-------------|---------|
| Stratified | `0`: Choi et al; `1`: Bhagwat & Ghajar (default); `2`: Franca & Lahey; `4`: modified Bhagwat & Ghajar |
| Slug/Bubble | `0`: Choi et al; `1`: Bhagwat & Ghajar (default); `4`: modified Bhagwat & Ghajar |
| Annular/Churn | `0`: Choi et al; `1`: Bhagwat & Ghajar (default); `3`: Hibiki & Ishii; `4`: modified Bhagwat & Ghajar |

## Thermal and Geometry Controls

| Field | Type | Default | Purpose |
|-------|------|---------|---------|
| `thermalEquilibrium` | bool | `true` | Wall temperature initialization: `true` = linear interpolation between ambient and fluid; `false` = ambient throughout wall |
| `latentHeatCond` | bool | `true` | When `false`, disables condensation latent heat in energy equation (rare retrograde-condensation cases) |
| `reverseTemp` | number | ambient | Return-gas temperature at outlet boundary during reverse flow (°C) |
| `geometryFollowsFlow` | bool | `true` | If `true`, measured length zero is the pipe inlet; if `false`, zero is the pipe outlet |
| `xyMode` | bool | `false` | Enables XY-coordinate-based inclination reconstruction |
| `xProdStart`, `yProdStart` | number | `0` | Start coordinates for production line (when `xyMode = true`) |
| `xServiceStart`, `yServiceStart` | number | `0` | Start coordinates for service line (when `xyMode = true` and `gasLine = true`) |
| `gasLine` | bool | `false` | Indicates presence of a service/injection gas line |

## Numerical and Performance Controls

| Field | Type | Default | Purpose |
|-------|------|---------|---------|
| `steadyStateOrder` | int | `1` | Steady-state ODE method: `1` = first order; `2` = second-order Runge-Kutta |
| `checkValve` | int | `0` | `1` = adds check valve at production outlet (prevents reverse flow); `0` = reverse gas inflow allowed |
| `steadyGuess` | number | `-1` | Initial guess for steady-state closure (flow rate or pressure). Negative = auto-compute from hydrostatics |
| `parallelizeSA` | bool | `false` | Enables parallel execution of sensitivity-analysis cases |

## Performance Optimization Options

| Field | Type | Default | Purpose |
|-------|------|---------|---------|
| `pressureTable` | bool | `false` | Pre-builds compressibility table for production fluids (reduces transient runtime) |
| `gasTable` | bool | `false` | Pre-builds compressibility table for service-line gas |
| `dynamicTableModel` | bool | `false` | Post-processes compositional table from initial BO result (steady-state networks) |
| `srBpTable` | bool | `false` | Pre-builds solubility-ratio table (useful for expensive RS models) |
| `trackGOR` | bool | `true` | Updates GOR, API, BSW along pipeline at stream merging points |
| `trackGasDensity` | bool | `true` | Updates gas density and CO₂ fraction at merging points |

## Advanced Controls (`initialConfig.advanced`)

Expert tuning parameters for robustness and performance. These rarely need modification:

| Field | Default | Purpose |
|-------|---------|---------|
| `threads` | `1` | Number of execution threads |
| `monophasicCriterion` | `1e-5` | Minimum void fraction below which flow is treated as single-phase |
| `condensationCriterion` | `0.001` | Minimum void fraction for interphase mass-transfer activation |
| `simplePressureFrontier` | `true` | Cell-boundary pressure: `true` = average of neighbors; `false` = full hydrostatic/friction marching |
| `massTransferLimit` | `10 kg/(s·m)` | Disables latent heat above this mass-transfer rate |
| `accelerateSteadyConvergence` | `true` | Simplification that accelerates and stabilizes steady-state solver |
| `slipBoundaryCell` | `true` | Controls slip in last control volume (oscillation mitigation) |
| `relaxChokeTimestep` | `false` | Penalizes time-step growth during choke oscillations |
| `valveTimestepControl` | `false` | Restricts time-step increments during valve opening/closing |
| `disablePenalizeTimestep` | `false` | Disables holdup-oscillation time-step penalization (shutdown scenarios) |
| `sonicTime` / `sonicFlag` | — | Time windows for sonic CFL-based time-stepping (pressure-wave capture) |
| `compModelCorrectionTime` / `compModelCorrectionFlag` | — | Time windows for full compressibility model (liquid density time derivatives) |
| `despressRate` | `0.01 kgf/(cm²·s)` | Threshold for simplified vs. full two-step evolution |
| `dynTableMinDelay`, `dynTableMinDp`, `dynTableMinDt` | `0` | Flash mini-table refresh for compositional transient runs |

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
