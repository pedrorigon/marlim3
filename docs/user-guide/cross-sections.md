# Cross Sections

Define pipe cross-sections with material layers (inner to outer).

## Overview

Cross-sections describe the radial structure of the pipe wall — from the internal flow surface outward through multiple material layers to the external boundary. Each pipe segment references a cross-section by ID.

## Parameters

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Inner Diameter** | m | Internal diameter of the flow conduit |
| **Roughness** | m | Pipe wall roughness (for friction calculation) |
| **Annular** | — | If true, the cross-section has an annular geometry |

## Layers

Each cross-section contains one or more concentric layers defined from inside to outside:

| Layer Parameter | Description |
|-----------------|-------------|
| **Layer measurement type** | `THICKNESS` or `OUTER_DIAMETER` |
| **Thickness** | Layer thickness [m] (when type = THICKNESS) |
| **Outer diameter** | Layer outer diameter [m] (when type = OUTER_DIAMETER) |
| **Discretization** | Number of radial cells in this layer (for heat transfer) |
| **Material ID** | Reference to a material definition |

## Example: Typical Subsea Flowline

```
┌─────────────────────────────────┐
│  Layer 3: Concrete (60 mm)      │  materialId: 2
├─────────────────────────────────┤
│  Layer 2: PU Insulation (50 mm) │  materialId: 1
├─────────────────────────────────┤
│  Layer 1: Steel (12 mm)         │  materialId: 0
├─────────────────────────────────┤
│         Flow (D = 0.15 m)       │
└─────────────────────────────────┘
```

## JSON Structure

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
          "thickness": 0.05,
          "discretization": 3,
          "materialId": 1
        }
      ]
    }
  ]
}
```

!!! tip
    Increase layer discretization for thick insulation layers where the radial temperature gradient is significant.
