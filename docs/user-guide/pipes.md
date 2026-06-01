# Pipes

Pipes define the axial geometry of the system — the path along which the 1D conservation equations are solved. Each pipe is composed of segments, and each segment is discretized into computational cells where mass, momentum, and energy balances are evaluated.

---

## Production and Service Lines

Two parallel line families can exist in a model:

- **Production line** — the main multiphase flow path (wellbore, riser, flowline, pipeline). Always present.
- **Service line** — a secondary line for gas injection or gas circulation (gas-lift systems). Present only when the service-line feature is enabled in global settings.

> **JSON keys:**
>
> - Production: `productionPipe` (EN) · `dutosProducao` (PT) — array of segments
> - Service: `servicePipe` (EN) · `dutosServico` (PT) — array of segments

---

## Pipe Inclination and Gravity

Inclination determines the hydrostatic pressure component — the dominant pressure contribution in vertical wells and risers. It also affects flow-pattern transitions (stratified flow only occurs at near-horizontal inclinations).

Two approaches exist:

### Direct Angle Specification

Each segment declares its inclination directly. Positive angles mean upward flow; negative means downward flow. Zero is horizontal.

> **JSON key:** `inclination` (EN) · `inclinacao` (PT) — unit: degrees (from horizontal)
> Python API alternative: `angle` — unit: radians

### Coordinate-Based Reconstruction (XY Mode)

Instead of specifying angles, you provide endpoint coordinates for each segment. The simulator computes inclinations from the geometry. This is natural when working from well-survey or pipeline-route data.

> See [General](general.md) for XY-mode configuration.

---

## Radial Structure (Cross-Section Reference)

Each pipe segment references a cross section (by integer ID) that defines its internal diameter, roughness, and concentric wall layers. Multiple segments can share the same cross section.

> **JSON key:** `crossSectionId` (EN) · `idSecaoTransversal` (PT)

---

## Thermal Environment

Each pipe segment has an external thermal boundary condition defined by:

### External Temperature

The temperature of the medium surrounding the pipe at this segment's position. For subsea flowlines this is seawater temperature; for wellbores it may follow a geothermal gradient.

> **JSON key:** `externalTemperature` (EN) · `temperaturaExterna` (PT) — unit: °C

### Environment Type

The nature of the external medium affects the heat-transfer coefficient at the outer boundary:

- **Buried / formation (0)** — heat conducts into surrounding rock (uses formation model).
- **Water (1)** — external convection by seawater or river water.
- **Atmosphere (2)** — external convection by air.

> **JSON key:** `environment` (EN) · `ambiente` (PT)
> Values: `0` = buried, `1` = water, `2` = atmosphere

### Formation Reference

For buried segments, a formation ID links to the formation thermal properties (conductivity, density, specific heat) that control long-term heat exchange.

> **JSON key:** `formationId` (EN) · `idFormacao` (PT)

### Variable External Conditions

When ambient temperature or external convection velocity varies along a segment, initial-condition profiles can be specified at normalized positions (0 to 1):

- **Ambient temperature profile** — external temperature distribution along the segment.
- **External velocity** — convection velocity of the surrounding medium (seawater current, wind speed).

> **JSON keys (inside `initialConditions`):**
>
> - Position: `measuredPosition` (EN) · `posicaoMedida` (PT)
> - Temperature: `ambientTemp` (EN) · `temperaturaAmbiente` (PT) — unit: °C
> - Velocity: `ambientVel` (EN) · `velocidadeAmbiente` (PT) — unit: m/s

---

## Discretization (Mesh)

Each pipe segment is subdivided into computational cells. The cell length determines numerical resolution:

- **Fine mesh** (small cells): Better resolution of thermal/pressure gradients and wave propagation. Higher computational cost.
- **Coarse mesh** (large cells): Faster execution, adequate for screening and steady-state validation. Less local detail in transient mode.

### Grouped Cells (Recommended)

Define blocks of uniform cells with a given count and individual cell length. Multiple blocks can compose a segment with variable resolution.

> **JSON keys:** `grouping: true`, then `discretization` array with objects containing:
>
> - `numCells` (EN) · `numCelulas` (PT) — number of cells in this block
> - `length` (EN) · `comprimento` (PT) — individual cell length [m]

### Explicit Cells

For non-uniform meshes, provide explicit cell-length arrays directly.

> **JSON key:** `grouping: false`, with explicit `dx` vectors

**Rules of thumb:**

- Flowlines: 10–20 m cells for general studies.
- Risers and wellbores: 5–10 m cells where gradients are steeper.
- Near accessories (valves, pumps, ESPs): consider local refinement.

---

## Example: Simple Horizontal Pipeline

```json
{
  "initialConfig": {
    "gasLine": false,
    "xyMode": false
  },
  "productionPipe": [
    {
      "id": 0,
      "active": true,
      "crossSectionId": 0,
      "formationId": 0,
      "inclination": 0.0,
      "externalTemperature": 4.0,
      "environment": 1,
      "grouping": true,
      "discretization": [
        {"numCells": 50, "length": 20.0}
      ]
    }
  ]
}
```

## Example: Multi-Segment Well (Vertical + Horizontal)

```json
{
  "productionPipe": [
    {
      "id": 0,
      "active": true,
      "crossSectionId": 0,
      "formationId": 0,
      "inclination": -90.0,
      "externalTemperature": 80.0,
      "environment": 0,
      "grouping": true,
      "discretization": [
        {"numCells": 30, "length": 10.0}
      ]
    },
    {
      "id": 1,
      "active": true,
      "crossSectionId": 1,
      "formationId": 0,
      "inclination": 0.0,
      "externalTemperature": 4.0,
      "environment": 1,
      "grouping": true,
      "discretization": [
        {"numCells": 100, "length": 25.0}
      ]
    }
  ]
}
```

!!! tip
    Start coarse (20–50 m cells) to validate boundary setup and overall behavior, then refine mesh in critical regions: near valves, pumps, strong thermal gradients, and elevation changes.
