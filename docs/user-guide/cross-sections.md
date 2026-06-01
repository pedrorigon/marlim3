# Cross Sections

Cross sections define the radial geometry of the pipe: the flow area, friction surface, and concentric layers through which heat is conducted.

## Concept

The axial (1D) conservation equations depend on radial definitions to compute:

- **Hydraulics:** Flow area, hydraulic diameter, and wall roughness determine velocity, friction factor, and pressure losses.
- **Heat transfer:** Layer thicknesses, discretization, and material conductivities determine radial temperature profiles and heat flux.

Each axial pipe segment (`productionPipe` or `servicePipe`) references a `crossSection` via `crossSectionId`.

## Core Parameters

| Parameter | Unit | Role |
|-----------|------|------|
| `id` | — | Integer identifier referenced by pipe segments |
| `active` | — | Whether this cross section is available for use |
| `innerDiameter` | m | Sets the flow area ($A = \pi d^2 / 4$) and velocity scale |
| `roughness` | m | Absolute wall roughness affecting friction factor (Colebrook/Moody) |
| `annular` | bool | Indicates annular-flow geometry (production in annulus around tubing) |

## Layer Model

Each cross section contains an array of concentric `layers` defined from the inner wall outward to the external boundary. Layers control radial heat-transfer resolution.

| Field | Type | Meaning |
|-------|------|---------|
| `layerMeasurementType` | string | `"THICKNESS"` or `"OUTER_DIAMETER"` |
| `thickness` | m | Layer thickness (when type is `THICKNESS`) |
| `outerDiameter` | m | Layer outer diameter (when type is `OUTER_DIAMETER`) |
| `discretization` | int | Number of radial nodes within this layer (default `1`) |
| `materialId` | int | Links to the corresponding `material` entry |

### How Layers Map to Materials

```
┌─────────────── Outer boundary (formation or ambient)
│  Layer N  ← materialId → cement, rock contact
│  Layer 2  ← materialId → insulation
│  Layer 1  ← materialId → steel pipe wall
│  ═════════  Inner surface (roughness applies here)
│  Fluid flow area (innerDiameter)
└─────────────── Pipe center
```

## Practical Guidance

- **Radial discretization:** Use at least 2–3 nodes in thick insulation layers to resolve thermal gradients during transient cooldown. Thin steel walls can use 1 node.
- **Diameter consistency:** Layers must nest concentrically — each layer's outer boundary should equal or exceed the previous layer's outer boundary.
- **Roughness selection:** Typical values range from 1.5×10⁻⁶ m (polished tubing) to 4.6×10⁻⁴ m (corroded old pipe). The README example uses 1.83×10⁻⁴ m.
- **Annular geometry:** When `annular = true`, the flow area is the annulus between tubing OD and casing ID. The `innerDiameter` then refers to the hydraulic equivalent.

## Example JSON

### Simple Pipeline (Steel + Insulation)

```json
{
  "crossSection": [
    {
      "id": 0,
      "active": true,
      "annular": false,
      "innerDiameter": 0.15,
      "roughness": 5e-5,
      "layers": [
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.012,
          "discretization": 1,
          "materialId": 0
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.050,
          "discretization": 3,
          "materialId": 1
        }
      ]
    }
  ]
}
```

### Well Completion (Steel + Completion Fluid + Casing + Cement)

```json
{
  "crossSection": [
    {
      "id": 1,
      "active": true,
      "annular": false,
      "innerDiameter": 0.10,
      "roughness": 1.83e-4,
      "layers": [
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.008,
          "discretization": 1,
          "materialId": 0
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.020,
          "discretization": 2,
          "materialId": 1
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.010,
          "discretization": 1,
          "materialId": 0
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.040,
          "discretization": 2,
          "materialId": 2
        }
      ]
    }
  ]
}
```

Where material IDs might map to: `0` = steel, `1` = completion fluid, `2` = cement.

!!! tip
    When in doubt about discretization, start with 1 node per layer for steady-state validation, then increase to 2–3 in insulation/cement layers for transient cooldown studies.
