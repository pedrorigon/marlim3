# General

Global settings: system type, slip models, correlations, and numerical parameters.

## Overview

The General tab contains simulation-wide configuration that affects all aspects of the solver: the system type, flow correlations, heat transfer options, geometry conventions, and numerical solver tuning.

## System Type

| Value | Description |
|-------|-------------|
| **MULTIFASICO** | Production system (multiphase flow) |
| **INJETOR** | Injection well system |

For branch/network JSON files, the equivalent values in `schema_branch.json` are:

| Branch schema value | Meaning |
|---------------------|---------|
| **PROD** | Producer system |
| **INJ** | Injector system |

## Basic Settings

| Parameter | Description |
|-----------|-------------|
| **Sensitivity analysis** | If true, performs parameter sensitivity sweeps |
| **Classic output** | Controls the format of simulation output |
| **Detailed screen output** | Controls verbosity during simulation |

## Slip & Correlations

### Slip Models

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Slip in steady-state** | true | Consider phase velocity differences in steady-state |
| **Slip in transient** | true | Consider slip in transient solver |
| **Drift-type model** | true | Use drift-flux model (vs. black-box correlations) |

### Flow Pattern Map

| Code | Description |
|------|-------------|
| 0 | Simplified Barnea |
| 1 | Complete Barnea |

### Mass Transfer Model

| Code | Description |
|------|-------------|
| 0 | Complete (implicit) |
| 1 | Complete (explicit) |
| 2 | Simplified isothermal |
| 3 | No mass transfer |

### Slip Correlations by Flow Pattern

| Flow Pattern | Options |
|--------------|---------|
| **Stratified** | `0` Choi et al, `1` Bhagwat & Ghajar, `2` França & Lahey, `4` B&G modified |
| **Bubble/Slug** | `0` Choi et al, `1` Bhagwat & Ghajar, `4` B&G modified |
| **Annular/Churn** | `0` Choi et al, `1` Bhagwat & Ghajar, `3` Hibiki & Ishii, `4` B&G modified |

## Heat Transfer

| Parameter | Description |
|-----------|-------------|
| **Thermal equilibrium** | Strategy for pipe wall temperature distribution |
| **Latent heat in condensation** | If false, disables latent heat in condensation |
| **Mass transfer limit** | Maximum allowed mass transfer rate [kg/(s·m)] |
| **Reverse temperature** | Gas return temperature for reverse flow [°C] |
| **Col-annular coupling steps** | Strong coupling iterations between column and annular |

## Geometry & Coordinates

| Parameter | Description |
|-----------|-------------|
| **Geometry follows flow direction** | If true, filling direction matches flow direction |
| **X/Y prod start** | Initial coordinates of production line (XY mode) |
| **X/Y service start** | Initial coordinates of service line (XY mode) |

## Numerical Aspects

| Parameter | Description |
|-----------|-------------|
| **Track GOR** | Update GOR along pipe using mixing rules |
| **Track gas density** | Update gas density along pipe |
| **Gas density BO correction** | Consider in-situ deviation for free gas density |
| **RS/PB table** | Pre-build RS table using black oil before simulation |
| **Parallelize SA** | Parallelize sensitivity analysis runs |

### Additional Branch-Schema Controls

The branch schema exposes additional controls that may not all appear in the current GUI:

| Parameter | Description |
|-----------|-------------|
| **checkValve** | Enables check valve at production outlet to block reverse flow |
| **steadyStateOrder** | Steady-state numerical order: 1st order or 2nd order RK |
| **steadyGuess** | User initial guess for steady-state solver objective |
| **fluidProp** | Legacy hybrid black-oil/PVTSim property mode |
| **diffusion3dMode** | Enables 3D thermal diffusion model |
| **diffusion3dThreads** | Thread count for 3D diffusion |
| **diffusion3dJson** | Input file for 3D thermal diffusion model |
| **waxMode** | Enables wax deposition modeling |

### Advanced Block Highlights

In `initialConfig.advanced`, branch schema also defines expert numerical controls such as:

- `monophasicCriterion` and `condensationCriterion`
- `regulaFalsiSearchCriterion`
- `simplePressureFrontier`
- `relaxChokeTimestep`, `disablePenalizeTimestep`, `valveTimestepControl`
- `threads` and `matrixThreads`
- Dynamic flash mini-table controls: `dynTableMinDelay`, `dynTableMinDp`, `dynTableMinDt`
- Sonic-event controls: `sonicTime`, `sonicFlag`

Use these options only when reproducing legacy cases or tuning challenging simulations.

## JSON Structure

```json
{
  "system": "PROD",
  "versaoJSON": "3.0.0",
  "initialConfig": {
    "sensitivityAnalysis": false,
    "steadyStateSlip": true,
    "transientSlip": true,
    "driftModel": true,
    "flowPatternMap": 0,
    "massTransfer": 0,
    "checkValve": 0,
    "steadyStateOrder": 1,
    "steadyGuess": -1,
    "correlationsByPattern": {
      "stratified": 1,
      "slugBubble": 1,
      "annularChurn": 1
    },
    "thermalEquilibrium": false,
    "geometryFollowsFlow": true,
    "advanced": {
      "threads": 1,
      "massTransferLimit": 10.0,
      "simplePressureFrontier": true
    }
  }
}
```
