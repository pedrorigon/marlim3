# Pipes

This section defines the axial geometry of the production and service lines solved by Marlim3.

## Concept

Each pipe is discretized into a sequence of axial segments. For each segment, the solver needs:

- **Geometry:** Radial structure (via `crossSectionId`), inclination, and length.
- **Thermal environment:** External temperature and formation properties.
- **Discretization:** Cell count and cell lengths that control numerical resolution.

Two line families exist in the model:

| Line | Object | Purpose |
|------|--------|---------|
| Production line | `productionPipe` | Main multiphase flow path (wellbore, riser, flowline) |
| Service line | `servicePipe` | Gas injection/service line (when `initialConfig.gasLine = true`) |

## Segment-Level Fields

| Field | Unit | Purpose |
|-------|------|---------|
| `id` | — | Segment identifier |
| `active` | — | Whether this segment is part of the simulation |
| `crossSectionId` | — | Links to radial geometry and roughness |
| `formationId` | — | Links to surrounding formation thermal properties |
| `inclination` | deg | Pipe inclination from horizontal (positive = upward flow); used when `xyMode = false` |
| `angle` | rad | Alternative inclination in radians (Python API convention) |
| `externalTemperature` | °C | Ambient/sea/formation temperature at this segment |
| `environment` | int | External environment type (`0`: buried/formation; `1`: water; `2`: atmosphere) |

## Coordinate Strategies

Marlim3 supports two equivalent geometry definitions:

### 1. Direct Inclination (`xyMode = false`, default)

Each segment specifies its inclination angle directly. This is the simplest approach for conceptual models.

### 2. XY-Based Reconstruction (`xyMode = true`)

Inclinations are computed from endpoint coordinates. Requires:

- `initialConfig.xProdStart`, `initialConfig.yProdStart` (production line start)
- `initialConfig.xServiceStart`, `initialConfig.yServiceStart` (service line, if present)
- Each segment provides endpoint XY coordinates

Use XY mode when trajectory data comes from well surveys or pipeline route coordinates.

## Discretization

Discretization controls numerical resolution, gradient capture, and runtime cost:

| Approach | Fields | When to use |
|----------|--------|-------------|
| Grouped cells | `grouping: true`, array of `{numCells, length}` | Uniform segments (most common) |
| Explicit cells | `grouping: false`, explicit `dx` arrays | Non-uniform mesh refinement |

**Trade-offs:**

- **Fine mesh** (small cells): Better thermal/pressure gradient resolution, more accurate wave propagation, higher runtime.
- **Coarse mesh** (large cells): Faster execution, less local detail, adequate for screening studies.

**Rule of thumb:** Start with 10–20 m cell length for flowlines, 5–10 m for risers and wellbores where gradients are steeper.

## Initial Conditions Along the Pipe

Each pipe segment can carry initial-condition profiles:

| Field | Description |
|-------|-------------|
| `initialConditions.measuredPosition` | Normalized positions (0–1) along segment |
| `initialConditions.ambientTemp` | External temperature profile [°C] |
| `initialConditions.ambientVel` | External convection velocity [m/s] |

## Example JSON

### Simple Horizontal Pipeline

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

### Multi-Segment Well (Vertical + Horizontal)

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
