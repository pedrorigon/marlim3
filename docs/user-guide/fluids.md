# Fluids

Fluid modeling is the thermodynamic foundation of the simulation. Nearly every governing equation — momentum, energy, mass transfer — depends on fluid properties: density, viscosity, enthalpy, compressibility, and phase split. If fluid data is inconsistent, pressure and temperature predictions degrade quickly.

## Thermodynamic Model Family

Marlim3 supports three approaches to compute fluid properties, each with different fidelity and cost:

### Black-Oil Model

The simplest and fastest approach. Fluid behavior is characterized by correlations that depend on API gravity, GOR, gas density, and BSW. Phase split is determined by a solution-gas-ratio (RS) correlation. Adequate for most engineering studies.

### Compositional Model

A rigorous thermodynamic model where the fluid is described by pseudocomponent molar fractions and equilibrium is computed via flash calculations. Required for rich-gas/condensate systems or when composition changes significantly along the flow path. More expensive computationally.

### Flash-Table Model

Precomputed PVT properties from external software (PVTSim) stored in `.tab` or `.ctm` files. The simulator interpolates properties at each P-T condition from this table. Accurate within the table range but depends entirely on the quality and coverage of the external data.

> **JSON keys (inside global config):**
>
> - Flash-table mode: `flashTableFluidModel` (EN) · `modeloFluidoTabelaFlash` (PT) — default `false`
> - Compositional mode: `compositionalFluidModel` (EN) · `modeloFluidoComposicional` (PT) — default `false`
> - PVT file: `pvtFile` (EN) · `pvtsimArq` (PT)

When both flags are `false`, the black-oil model is used.

!!! note
    In transient simulations, black-oil is the lightest option. Consider enabling precomputed property tables to reduce runtime further.

## Fluid-Type Classification

The fluid type determines how interphase mass transfer is modeled:

- **Liquid-dominated (0):** Mass transfer follows black-oil logic (gas dissolves in/comes out of oil based on solubility ratio).
- **Gas-dominated (1):** Mass transfer uses quality/void-fraction relations (condensate drops out of gas).

> **JSON key:** `fluidType` (EN) · `tipoFluido` (PT) — default `0`

---

## Production Fluid

Each simulation can define one or more production fluids. Sources, initial conditions, and mixing rules reference fluids by their integer identifier.

> **JSON key:** `productionFluid` (EN) · `fluidosProducao` (PT) — array of objects

### Essential Properties

Every production fluid requires at minimum:

- **Oil API gravity** — characterizes oil density and correlations. *(key: `api`)*
- **Gas-oil ratio** — volume of dissolved gas per volume of oil at standard conditions. *(key: `gor` / `rgo`, unit: Sm³/Sm³)*
- **Gas relative density** — gas density divided by air density, both at standard conditions. *(key: `gasDensity` / `densidadeGas`)*
- **Water cut (BSW)** — volumetric fraction of water in the liquid stream. *(key: `bsw`, unit: m³/m³, default: 0)*

### Dead-Oil Viscosity

Dead-oil viscosity is the foundation of the viscosity model — it represents oil viscosity at atmospheric pressure with no dissolved gas. Several empirical correlations are available:

| Option | Correlation | Notes |
|--------|-------------|-------|
| 0 | ASTM | Requires two temperature-viscosity points |
| 1 | Beggs & Robinson | — |
| 2 | Modified Beggs & Robinson | — |
| 3 | Glaso | Default |
| 4 | Kartoatmodjo-Schmidt | — |
| 5 | Petrosky-Farshad | — |
| 6 | Beal | — |
| 7 | User table | Requires temperature and viscosity arrays |

For the ASTM model (0), provide two calibration points: temperature and viscosity at each point.

> **JSON keys:** `deadOilModel` (EN) · `modeloOleoMorto` (PT)
> ASTM points: `temp1`, `visc1`, `temp2`, `visc2`
> Table: `deadOilTemp` / `deadOilVisc`

### Live-Oil and Undersaturated-Oil Viscosity

Live-oil viscosity modifies dead-oil viscosity to account for dissolved gas. Undersaturated-oil viscosity further corrects for pressures above bubble point.

**Live-oil models:**

| Option | Correlation |
|--------|-------------|
| 0 | Beggs & Robinson (default) |
| 1 | Kartoatmodjo-Schmidt |
| 2 | Petrosky-Farshad |

**Undersaturated-oil models:**

| Option | Correlation |
|--------|-------------|
| 0 | Vasquez & Beggs (default) |
| 1 | Kartoatmodjo-Schmidt |
| 2 | Petrosky-Farshad |
| 3 | Beal |
| 4 | Khan |

> **JSON keys:**
>
> - Live oil: `liveOilModel` (EN) · `modeloOleoVivo` (PT)
> - Undersaturated: `undersaturatedOilModel` (EN) · `modeloOleoSubSaturado` (PT)

### Solution-Gas Ratio (Bubble-Point Behavior)

The solution-gas ratio (RS) determines how much gas is dissolved in oil at a given pressure and temperature. At the bubble point, RS equals the specified GOR.

| Option | Correlation |
|--------|-------------|
| 0 | Vazquez & Beggs (default) |
| 1 | Lasater |
| 2 | Standing |
| 3 | Glaso |
| 4 | Livia Fulchignoni |

> **JSON key:** `srBpModel` (EN) · `modeloRsPb` (PT)

### Gas Critical Properties

Gas compressibility calculations require pseudo-critical temperature and pressure. Different correlations suit different gas compositions, particularly when CO₂ content is significant.

| Option | Correlation | Notes |
|--------|-------------|-------|
| 0 | Standard | General purpose |
| 1 | Brown et al | Better for CO₂-rich gases |
| 2 | Piper et al | Better for CO₂-rich gases |

> **JSON key:** `criticalCorrelation` (EN) · `correlacaoCritica` (PT) — default `1`

### Emulsion Effects

When water and oil flow together, emulsions can dramatically increase apparent viscosity. The emulsion model modifies oil viscosity based on water content:

| Option | Model | Additional input required |
|--------|-------|--------------------------|
| 0 | Linear water-quality weighting | None (default) |
| 1 | Weak Woelflin | — |
| 2 | Medium Woelflin | — |
| 3 | Strong Woelflin | — |
| 4 | Exponential | Coefficients A and B |
| 5 | Pal-Rhodes | Multiplier φ₁₀₀ |
| 6 | User table | BSW and multiplier vectors |
| 7 | Below-saturation BSW model | — |

> **JSON key:** `emulsionType` (EN) · `tipoEmulsao` (PT)
> Additional: `emulsionCoefA` / `emulsionCoefB`, `phi100`, `bswVec` / `emulVec`, `bswInversion`

### Compositional Fluid Options

When compositional mode is active, the fluid can be defined by molar fractions of pseudocomponents (same order as in the `.ctm` file). Optionally, compositions can be corrected to match a specified GOR.

> **JSON keys:**
>
> - User provides composition: `userMolarFraction` (EN) · `fracaoMolarUsuario` (PT) — default `false`
> - Molar fractions: `molarFraction` (EN) · `fracaoMolar` (PT) — array
> - Correct composition to match GOR: `userGORComp` (EN) · `rgoCompUsuario` (PT) — default `false`

---

## Complementary Fluid

Some systems carry a **third liquid phase** that is inert, does not dissolve in gas, and does not slip relative to the hydrocarbon mixture. Examples include:

- Glycol or ethanol cushions (hydrate prevention before master valve opening)
- Chemical inhibitors
- Completion fluids (brine)

> **JSON key:** `complementaryFluid` (EN) · `fluidoComplementar` (PT)

### Fluid-Type Options

The type determines how properties are computed:

- **Generic (0):** User provides all physical properties (density, compressibility, viscosity, etc.)
- **Water-based (1):** Only salinity is required — correlations handle the rest.
- **Friction-reducer (2):** Like generic, but with internal friction-reduction treatment.

> **JSON key:** `complementaryFluidType` (EN) · `tipoFluidoComplementar` (PT)

For the generic type, the following properties define the fluid:

- **Density** at standard conditions [kg/m³] *(key: `density` / `densidade`)*
- **Compressibility** [1/Pa] *(key: `compressibility` / `compressibilidade`)*
- **Thermal expansivity** [1/K] *(key: `thermalExpansivity` / `expansividadeTermica`)*
- **Surface tension** [N/m] *(key: `surfaceTension` / `tensaoSuperficial`)*
- **Specific heat** [J/(kg·K)] *(key: `specificHeat` / `calorEspecifico`)*
- **Conductivity** [W/(m·K)] *(key: `conductivity` / `condutividade`)*
- **Viscosity** via two ASTM points *(keys: `temp1`, `visc1`, `temp2`, `visc2`)*

For water-based type, provide only:

- **Salinity** [g/(kg water)] *(key: `salinity` / `salinidade`)*

---

## Gas Fluid (Service Line)

When a service/injection gas line is present, the gas properties must be defined separately. This represents dry gas (not associated with production fluids).

> **JSON key:** `gasFluid` (EN) · `fluidoGas` (PT)

Key concepts:

- **Gas relative density** — same convention as production fluid gas. *(key: `gasDensity` / `densidadeGas`)*
- **CO₂ fraction** — molar fraction of CO₂ in the gas. *(key: `CO2Fraction` / `fracaoCO2`, default: 0)*
- **Critical-property correlation** — `1`: Brown et al; `2`: Piper et al. *(key: `criticalCorrelation` / `correlacaoCritica`)*
- **Flash-table usage** — whether to use the flash table for service-line gas when flash-table mode is active globally. *(key: `useFlashTable` / `usarTabelaFlash`, default: `false`)*
- **Composition** — explicit molar fractions for compositional runs. *(keys: `userMolarFraction`, `molarFraction`)*

---

## Specific-Heat Model

The specific heat (Cp) controls how temperature responds to pressure changes and heat exchange. Two options:

- **Black-oil correlation (0):** Internal model based on fluid parameters.
- **PVT table interpolation (1):** Reads Cp from the provided PVT file.

> **JSON key:** `cpModel` (EN) · `modeloCp` (PT) — default `0`

## Liquid Joule-Thomson Derivative

The Joule-Thomson effect describes temperature changes during isenthalpic expansion. By default the black-oil model computes the liquid density derivative with temperature (dρ_l/dT) internally. Alternatively, this value can be read from a PVTSim table for higher accuracy.

> **JSON key:** `jtlModel` (EN) · `modeloJTL` (PT) — default `0` (black-oil); `1` = use table

## Enthalpy from PVT Table

For more accurate thermal predictions, liquid and gas enthalpies can be read from the PVT file instead of using black-oil estimates.

> **JSON key:** `latentHeat` (EN) · `latente` (PT) — default `false`

## Free-Gas Density Correction

By default, in-situ gas density is computed from standard-condition values. When this correction is enabled, higher in-situ pressure is accounted for, recognizing that lighter hydrocarbons remain in solution under those conditions.

> **JSON key:** `freeGasDensityCorrectionBO` (EN) · `correcaoDenGasLivreBlackOil` (PT) — default `false`

---

## Precomputed Property Tables

Computing compressibility and its derivatives at every iteration can be expensive in transient simulations. Marlim3 can pre-build interpolation tables spanning a user-defined pressure-temperature grid before the simulation starts.

The grid is defined separately:

> **JSON key:** `compTable` (EN) · `tabelaComp` (PT)

Required grid parameters:

- **Number of points** in each dimension *(key: `numPoints` / `numPontos`)*
- **Pressure range** [kgf/cm²] *(keys: `minPressure` / `maxPressure` · `pressaoMin` / `pressaoMax`)*
- **Temperature range** [°C] *(keys: `minTemperature` / `maxTemperature` · `temperaturaMin` / `temperaturaMax`)*

!!! warning
    The grid must fully cover the expected operating envelope. Extrapolation outside the table bounds is unreliable and can cause numerical errors.

---

## Example JSON

```json
{
  "initialConfig": {
    "flashTableFluidModel": false,
    "compositionalFluidModel": false,
    "fluidType": 0,
    "cpModel": 0,
    "jtlModel": 0,
    "pressureTable": true,
    "pvtFile": "PVT_CASE.ctm"
  },
  "productionFluid": [
    {
      "id": 0,
      "active": true,
      "api": 28.0,
      "gasDensity": 0.65,
      "gor": 150.0,
      "bsw": 0.10,
      "deadOilModel": 3,
      "liveOilModel": 0,
      "criticalCorrelation": 1,
      "emulsionType": 0
    }
  ],
  "compTable": {
    "active": true,
    "numPoints": 50,
    "minPressure": 1.0,
    "maxPressure": 500.0,
    "minTemperature": 4.0,
    "maxTemperature": 120.0
  }
}
```

!!! tip
    Start with one fluid and validate a simple steady-state case before adding advanced options (compositional corrections, emulsions, multiple fluids) incrementally.
