# Pipes

Define production and service line geometry, segments, and initial conditions.

## Overview

The Pipes tab configures the physical layout of the pipeline system — including the production pipe and an optional gas service line. Each pipe is discretized into cells for numerical computation.

## Pipe Settings

| Parameter | Description |
|-----------|-------------|
| **Gas service line** | Enable the gas-lift / gas service line |
| **XY coordinate mode** | Infer pipe inclinations from XY coordinates instead of explicit angles |

## Production Pipe

The production pipe is defined as a sequence of segments (ducts), each with:

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Cross Section ID** | — | Reference to a cross-section definition |
| **Formation ID** | — | Reference to a rock formation |
| **Inclination** | ° | Pipe inclination from horizontal (positive = upward) |
| **External temperature** | °C | Ambient temperature outside the pipe |
| **Elevation** | m | Vertical position (when using XY mode) |

## Discretization

Each pipe segment is discretized into computational cells. Two modes are available:

### Grouped Discretization

Define groups of uniform cells:

| Parameter | Description |
|-----------|-------------|
| **Number of cells** | How many cells in this group |
| **Cell length** | Length of each cell [m] |

### Cell-by-Cell (ungrouped)

Specify individual cell lengths (array of dx values).

## Service Pipe

When gas service line is enabled, a separate pipe geometry is defined for the gas injection path. It uses the same structure as the production pipe.

## Geometry Preview

The GUI provides an interactive geometry plot showing the pipe profile (measured depth vs. elevation), helping visualize the pipeline trajectory.

## JSON Structure

```json
{
  "productionPipe": [
    {
      "id": 0,
      "active": true,
      "crossSectionId": 0,
      "formationId": 0,
      "inclination": -5.0,
      "externalTemperature": 4.0,
      "grouping": true,
      "discretization": [
        {"numCells": 50, "length": 20.0}
      ]
    }
  ],
  "servicePipe": [
    {
      "id": 0,
      "crossSectionId": 1,
      "formationId": 0,
      "inclination": -5.0,
      "externalTemperature": 4.0,
      "grouping": true,
      "discretization": [
        {"numCells": 50, "length": 20.0}
      ]
    }
  ]
}
```

!!! tip
    Use XY coordinate mode for complex well trajectories where it's easier to define the well path by coordinates rather than segment inclinations.
