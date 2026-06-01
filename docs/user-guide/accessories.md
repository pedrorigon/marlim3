# Accessories

Accessories are localized components that introduce active physical behavior — inflow, choking, pumping, pigging — at specific measured positions along the pipe.

## Concept

Without accessories, a simulation is a passive conduit with boundary constraints only. Accessories introduce operational control and localized physics at discrete positions within the domain.

All accessory objects share common patterns:

- They are located at a `measuredLength` position along the pipe.
- Time-dependent behavior is defined via `time` arrays paired with corresponding value arrays.
- Each object has an `id` for referencing and an `active` flag.

## Source Models

Sources introduce fluid into the production system at a given position. They represent reservoir inflow, injected fluids, or prescribed rates.

| Object | Main role | Key input |
|--------|-----------|-----------|
| `ipr` | Reservoir inflow performance relation | Static reservoir pressure, productivity index, IPR curve type |
| `liquidSource` | Prescribed liquid volumetric inflow (std conditions) | `liquidFlowRate` [sm³/d], `temperature` [°C], `prodFluidId` |
| `massSource` | Prescribed total mass inflow | `massFlowRate` [kg/s], `temperature` [°C] |
| `gasSource` | Prescribed gas inflow | Gas flow rate, temperature |
| `pressureSource` | Local pressure-coupled exchange point | Exchange pressure, direction |
| `porousRadialSource` | Radial porous-media inflow model | Permeability, reservoir pressure, radius |
| `porous2DSource` | 2D porous-media inflow model | 2D permeability field |
| `gasLiftSource` | Gas-lift valve injection device | Valve design parameters, operating pressure |

!!! note
    A `liquidSource` with given `prodFluidId` implies associated free-gas flow based on the referenced fluid's GOR and local P-T conditions.

### Example: Liquid Source

```json
{
  "liquidSource": [
    {
      "id": 0,
      "active": true,
      "prodFluidId": 0,
      "measuredLength": 2500.0,
      "time": [0, 3600],
      "liquidFlowRate": [1500, 1200],
      "temperature": [60.0, 58.0]
    }
  ]
}
```

## Flow-Control Devices

| Object | Main role | Key fields |
|--------|-----------|------------|
| `valve` | Generic two-phase restriction (Sachdeva model) | `measuredLength`, `time`, `opening`, `cd` (discharge coefficient) |
| `masterValve` | Production-tree master valve | Opening schedule, choke-like behavior |
| `masterValve2` | Secondary master valve (service tree) | Opening schedule |
| `surfaceChoke` | Surface choke at outlet | Opening schedule, choke sizing |
| `injectionChoke` | Injection-side choke | Opening schedule |

### Valve Opening Interpretation (`cvCurve`)

| `cvCurve` | Meaning of `opening` vector |
|-----------|------------------------------|
| `0` | Area ratio: valve free area / pipe area (dimensionless, 0–1) |
| `1` | Valve-stem displacement (requires `x1`/`cv1` Cv-curve) |

### Example: Valve with Time Schedule

```json
{
  "valve": [
    {
      "id": 0,
      "active": true,
      "measuredLength": 1500.0,
      "cvCurve": 0,
      "cd": 0.84,
      "time": [0, 60, 120],
      "opening": [0.0, 0.5, 1.0]
    }
  ]
}
```

## Pumping Models

| Object | Meaning | Key inputs |
|--------|---------|------------|
| `esp` | Electrical submersible pump (curve-based) | Head-flow curve, frequency, number of stages |
| `volumetricPump` | Positive-displacement pump | Displacement volume, speed |
| `pressureDrop` | Prescribed pressure increment/decrement | Δp value or schedule |

### Example: ESP

```json
{
  "esp": [
    {
      "id": 0,
      "active": true,
      "measuredLength": 3000.0,
      "time": [0],
      "frequency": [60.0],
      "stages": 100
    }
  ]
}
```

## Pigging Model

`pig` represents pig launch/receive events. The pig is tracked along its measured-length position over time, creating a moving boundary that affects pressure drop and liquid displacement.

```json
{
  "pig": [
    {
      "id": 0,
      "active": true,
      "launchTime": 100.0,
      "launchPosition": 0.0
    }
  ]
}
```

## Design Tips

- **Position alignment:** Keep accessory positions consistent with the physical well/pipeline architecture (e.g., DHSV position, ESP depth).
- **Avoid overlap:** Do not place multiple active restrictions at exactly the same measured length unless physically justified.
- **Incremental build-up:** Introduce one active device at a time during model construction to isolate its effect on the solution.
- **Opening ramps:** Avoid instantaneous valve openings (step from 0 to 1); use short ramps for numerical stability in transient mode.

!!! tip
    When combining gas-lift with an ESP, ensure the gas-lift source is upstream of the pump to avoid injecting gas directly into pump stages.
