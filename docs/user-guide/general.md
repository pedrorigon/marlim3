# General

General configuration defines model scope, multiphase physics options, geometry conventions, and solver controls that apply to the entire simulation.

> **JSON key:** `initialConfig` (EN) · `configuracaoInicial` (PT)

---

## System Declaration

At root level, choose producer or injector simulation:

> **JSON key:** `system` (EN) · `sistema` (PT)
> Values: `PROD` / `INJ` (EN) · `MULTIFASICO` / `INJETOR` (PT)

---

## Core Multiphase Model Controls

### Slip Activation

- `steadyStateSlip` · `escorregamentoPermanente`
- `transientSlip` · `escorregamentoTransiente`

Both default to `true`.

### Drift vs. Black-Box in Steady State

`driftModel` controls the steady-state closure family:

- `true` = drift-based
- `false` = black-box correlations

> **JSON key:** `driftModel` (EN) · `tipoModeloDrift` (PT)

In transient mode, drift model is enforced by the engine.

### Flow-Pattern Map and Correlations

- `flowPatternMap` · `mapaArranjo` (0 simplified Barnea, 1 full Barnea)
- `correlationsByPattern` · `correlacoesPorArranjo`
  - `stratified` · `estratificado`
  - `slugBubble` · `bolhaGolfada`
  - `annularChurn` · `anularChurn`

### Mass Transfer

`massTransfer` · `transferenciaMassa`:

- `0` full implicit
- `1` full explicit
- `2` simplified isothermal
- `3` disabled

---

## Geometry and Coordinate Conventions

- `geometryFollowsFlow` · `sentidoGeometriaSegueEscoamento`
- `xyMode` · `modoXY`
- `xProdStart`, `yProdStart` · `xProdInicio`, `yProdInicio`
- `xServiceStart`, `yServiceStart` · `xServInicio`, `yServInicio`
- `gasLine` · `linhaGas`

When `xyMode` is active, segment inclinations are inferred from coordinates instead of direct angles.

---

## Thermal/Fluid Auxiliary Controls

- `thermalEquilibrium` · `equilibrioTermico`
- `latentHeatCond` · `condlatente`
- `reverseTemp` · `tempReves`
- `steadyStateOrder` · `ordemperm`
- `steadyGuess` · `chutePerm`

`steadyStateOrder = 2` applies only to steady-state; transient uses first-order marching.

---

## Table and Tracking Performance Options

- `pressureTable` · `tabP`
- `gasTable` · `tabG`
- `dynamicTableModel` · `modeloTabelaDinamica`
- `srBpTable` · `tabelaRSPB`
- `trackGOR` · `trackRgo`
- `trackGasDensity` · `trackDensidadeGas`

These options trade pre-processing and memory for runtime speed or richer stream-property tracking.

---

## Check Valve and SA Parallelization

- `checkValve` · `CheckValve`
- `parallelizeSA` · `paralelizaAS`

Schema default for `checkValve` is `0`. Parser has legacy behavior that may assume check valve active if key is omitted, so set it explicitly for reproducibility.

---

## Advanced Block

Advanced numerical controls are nested in:

> **JSON key:** `initialConfig.advanced` (EN) · `configuracaoInicial.Avancado` (PT)

Commonly tuned fields:

- `monophasicCriterion` · `CriterioMonofasico`
- `condensationCriterion` · `CriterioCondensacao`
- `simplePressureFrontier` · `MedSimpPresFront`
- `massTransferLimit` · `limTransMass`
- `accelerateSteadyConvergence` · `AceleraConvergPerm`
- `slipBoundaryCell` · `escorregamentoCelulaContorno`
- `relaxChokeTimestep` · `RelaxaDTChoke`
- `valveTimestepControl` · `controleDTvalv`
- `disablePenalizeTimestep` · `desligaPenalizaDT`
- `compModelCorrectionTime` · `TcorrecaoModComp`
- `compModelCorrectionFlag` · `correcaoModComp`
- `despressRate` · `taxaDespre`
- `sonicTime` · `Tsonico`
- `sonicFlag` · `sonico`
- `threads` · `nthrd`
- `matrixThreads` · `nthrdMatriz`

---

## Example

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
  Keep advanced options at defaults until a baseline case is validated. Then tune one control at a time and compare against that baseline.
