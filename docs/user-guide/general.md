# General Configuration

The general configuration block governs the overarching setup of your simulation — from physical modeling choices and flow conventions to core numerical settings and performance tweaks.
These parameters define simulation modes, solver behaviors, and the relationship between geometry, physics, and available fluid and thermal models.

> **JSON key:** `initialConfig` (EN) · `configuracaoInicial` (PT) — top-level object

---

## 1. System Declaration

At the root level, declare the system role:

- `system` (EN) · `sistema` (PT):  
  - `"PROD"` / `"INJ"` (EN)  
  - `"MULTIFASICO"` / `"INJETOR"` (PT)

Defines producer vs injector top-level logic.

---

## 2. Simulation Geometry & Flow Direction

- `geometryFollowsFlow` · `sentidoGeometriaSegueEscoamento`:  
  If `true`, indices increase in the flow direction (`0` = pipe inlet, default).  
  If `false`, indices increase against the flow; angles are still measured with respect to the flow for compatibility.

- `gasLine` · `linhaGas`:  
  If `true`, includes a gas lift/service line.

- `xyMode` · `modoXY`:  
  If `true`, pipe inclinations are inferred from XY coordinates;  
  Provide:  
    - `xProdStart` · `xProdInicio`, `yProdStart` · `yProdInicio` (for production)  
    - `xServiceStart` · `xServInicio`, `yServiceStart` · `yServInicio` (for service/gas lines with gasLine = true)  
  If `false`, inclinations are specified explicitly in pipe segment definitions.

---

## 3. Multiphase Model & Solver Controls

### Core Physics Flags

- `steadyStateSlip` · `escorregamentoPermanente`
- `transientSlip` · `escorregamentoTransiente`

Both control whether to enable slip (phase velocity difference) modeling in steady and transient solver, respectively.  
**Default:** `true` for both.

- `driftModel` · `tipoModeloDrift`  
  `true` = Drift-flux closure (default and always enforced in transient),  
  `false` = Black-box correlations (steady state only).

### Mass Transfer Model

- `massTransfer` · `transferenciaMassa`:

  - `0` : Full implicit, complete mass transfer
  - `1` : Full explicit, complete mass transfer
  - `2` : Simplified, isothermal, neglect transients
  - `3` : Disabled (no mass transfer)

---

### Flow-Pattern Map and Correlation Selection

- `flowPatternMap` · `mapaArranjo`  
  `0` = Barnea simplified  
  `1` = Barnea complete

- `correlationsByPattern` · `correlacoesPorArranjo`  
    - `stratified` · `estratificado`
    - `slugBubble` · `bolhaGolfada`
    - `annularChurn` · `anularChurn`

Select integer code per flow pattern; see schema for details.

---

## 4. Thermal and Fluid Properties Controls

- `thermalEquilibrium` · `equilibrioTermico`  
  Linear interpolation between ambient and fluid (default: `true`).  
  If `false`, ambient temperature is used throughout the wall cross-section.

- `latentHeatCond` · `condlatente`  
  Controls inclusion of latent heat during phase changes (default: `true`).

- `steadyStateOrder` · `ordemperm`  
  1 = First order (default, robust)  
  2 = Runge-Kutta second order (for increased accuracy, steady state only).

- `steadyGuess` · `chutePerm`  
  Initial guess for steady solver (pressure or flow, as per boundary).

- `reverseTemp` · `tempReves`  
  Temperature to use for reverse gas flow at boundary (default: ambient).

---

## 5. Table Lookup, Tracking & Performance

These keywords enable data-driven performance improvements:

- `pressureTable` · `tabP`:  
  Precomputes compressibility tables for `black oil` estimation, indexed by reduced P/T (improves speed in long transient runs).

- `gasTable` · `tabG`:  
  Same as above but for gas in the service line.

- `dynamicTableModel` · `modeloTabelaDinamica`:  
  Generates flash property tables post-hoc from black-oil+compositional runs, for high-fidelity, multi-segment simulations.

- `srBpTable` · `tabelaRSPB`:  
  Solubility-ratio tables for advanced black-oil models.

- `trackGOR` · `trackRgo`:  
  Enables tracking GOR (and associated stream properties) with user mixing rules (default: `true`).

- `trackGasDensity` · `trackDensidadeGas`:  
  Tracks gas density along the system with mixing (default: `true`).

---

## 6. Valve/Boundary Controls & Misc

- `checkValve` · `CheckValve`:  
  Include (`1`) or omit (`0`, default) a check valve at outlet, preventing backflow.

- `parallelizeSA` · `paralelizaAS`:  
  If `true`, parallelizes sensitivity analysis.

---

## 7. Fluid Modeling Flags

- `blackOilModel` (traditional)
- `hybridFluidProps` / `propFluido` — combine black oil with PVTSim for some properties
- `flashTableModel` / `modeloFluidoTabelaFlash` — PVT flash table
- `compositionalModel` / `modeloFluidoComposicional` — true for full compositional modeling

**Best practice:**  
Start with pure black oil (`propFluido: 0`); only enable hybrid, flash or compositional models once baseline validation is reached.

---

## 8. Advanced Block

Fine tuning for convergence, solver stability, rare or edge-case configurations:

> **JSON key:** `initialConfig.advanced` (EN) · `configuracaoInicial.Avancado` (PT)

Notable controls:

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

Only tune after baseline results are established.

---

## 9. Formations & Lithology

For simulations involving heat exchange with surrounding rock, use:

> **JSON key:** `initialConfig.Formacao` (EN/PT)

- List properties (thermal conductivity, heat capacity, density — all per rock ID)
- `productionTime` · `TempoProducao`: controls heated radius calculation

---

## 10. Boundary Objects

**Pressure and flow boundary conditions** are formally specified as:

- `pressureCondition` · `condicaoPressao`
- `flowPressureCondition` · `condicaoVazPres`

Populate each with relevant vectors:  
`tempo` (time array), `pressao`, `temperatura`, etc.

!!! note
    If neither boundary object is declared, the system is closed, and fluids must enter through internal mass sources (IPR, fonteMassa, etc).

---

## 11. Output & Display

- `classicalOutput` · `saidaClassica`:  
  Cosmetic flag for style of final text output (no impact on results).

- `screenOutput` · `saidaTela`:  
  If `true`, shows detailed simulation progress; if `false` (default), only the percentage is shown.

---

## 12. Special/Developer Features

- `HISEP`, `modoDifus3D`, `modoDifus3DJson`, `modoParafina`:  
  Advanced, currently experimental, or developer-focused controls (see schema for further detail).

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

!!! tip
Keep advanced and optional features at defaults until your first baseline results are validated. When tuning for performance, stability or physics fidelity, only change one setting at a time — always compare against your previous working result.

---

## Practical Guidance

- Performance: Use precomputed tables (pressureTable, gasTable) for most transient cases. Only use full compositional/PVT/fancy models after baseline results with black oil models are checked.
- Slip & Drift: Leave slip enabled and use drift model unless you have robust justification otherwise.
- Geometry: Only use XY inference when direct segment angles are not available or reliable.
- Boundary conditions: Always declare at least one inlet/outlet boundary object or set up internal source(s).
- Sanity checks: Don’t deactivate penalizations or stability features unless you truly understand the trade-offs.