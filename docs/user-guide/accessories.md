# Accessories

Accessories are localized components that introduce active physical behavior — inflow, choking, pumping, pigging — at specific measured-length positions along the pipe. Without accessories, a simulation is a passive conduit constrained only at its boundaries. Accessories inject operational control and localized physics at discrete interior points.

---

## Common Structure

All accessories share a consistent design:

- **Position:** Each accessory is placed at a measured-length position along the pipe.
- **Time-dependent behavior:** A `time` array paired with value arrays defines how the accessory evolves (e.g., valve opening ramp, flow-rate schedule).
- **Identification:** Each has an `id` (for referencing) and an `active` flag (for toggling without deletion).

---

## Sources — Reservoir Inflow

The most common source is an **Inflow Performance Relationship (IPR)**, representing production from a reservoir. The IPR relates well inflow to the pressure difference between the static reservoir pressure and the flowing bottomhole pressure. Production rate adjusts dynamically as system pressure changes.

Key concepts:
- **Static reservoir pressure** — far-field formation pressure driving inflow.
- **Productivity index** — the slope of the linear IPR (rate per unit drawdown).
- **IPR curve type** — selects the inflow model (linear, Vogel, Fetkovich, etc.).

> **JSON key:** `ipr` (EN) · `fonteIPR` (PT)

---

## Sources — Prescribed Flow Rates

When reservoir coupling is not needed, fluid can be injected at prescribed rates:

### Liquid Source

Injects a prescribed liquid volumetric flow rate (at standard conditions) along with its associated gas (determined by the fluid's GOR and local P-T conditions). The source fluid is identified by its production-fluid index.

> **JSON key:** `liquidSource` (EN) · `fonteLiquido` (PT)

### Mass Source

Injects a prescribed total mass flow rate. Useful when the overall mass inflow is known but not the phase split.

> **JSON key:** `massSource` (EN) · `fonteMassa` (PT)

### Gas Source

Injects prescribed gas flow into the system — used for gas-lift injection points or surface gas injection.

> **JSON key:** `gasSource` (EN) · `fonteGas` (PT)

### Pressure Source

A pressure-coupled exchange point (constant-pressure reservoir analog). Flow direction depends on local vs. source pressure difference.

> **JSON key:** `pressureSource` (EN) · `fontePressao` (PT)

### Porous Radial Source

Models radial inflow from a porous medium (Darcy-based). Used when a simplified permeability/reservoir model replaces an explicit IPR.

> **JSON key:** `porousRadialSource` (EN) · `fontePorosoRadial` (PT)

### Gas-Lift Source

Represents a gas-lift valve — injects gas when casing pressure exceeds tubing pressure plus the valve's opening threshold. Parameters include valve area, discharge coefficient, and operating-pressure settings.

> **JSON key:** `gasLiftSource` (EN) · `fonteGasLift` (PT)

---

## Flow-Control Devices — Valves and Chokes

Valves and chokes restrict flow by reducing the effective area available to the multiphase mixture. They create a local pressure drop governed by choke models (e.g., Sachdeva's critical/subcritical two-phase model).

### Generic Valve

A general-purpose two-phase restriction placed at any position. Controlled by an opening schedule and characterized by a discharge coefficient.

The **opening interpretation** depends on the Cv-curve mode:
- **Area-ratio mode (0):** Opening is the fraction of valve free area relative to pipe area (0 = closed, 1 = full bore).
- **Stem-displacement mode (1):** Opening is the physical valve stem travel; a Cv-vs-stem curve translates this to effective area.

> **JSON key:** `valve` (EN) · `valvula` (PT)
> **Opening mode:** `cvCurve` (EN) · `cvCurve` (PT) — `0` = area ratio, `1` = stem displacement

### Master Valve

The production-tree master valve (DHSV or surface tree valve). Behaves like a generic valve but is semantically distinct for operational modeling.

> **JSON key:** `masterValve` (EN) · `valvulaMestra` (PT)

### Surface Choke

An outlet restriction at the surface, typically controlling wellhead or platform arrival pressure.

> **JSON key:** `surfaceChoke` (EN) · `chokeSuperficie` (PT)

### Check Valve (Non-Return Valve)

Prevents reverse flow. When flow reversal is detected, the valve closes immediately. Commonly placed at risers or downstream of pumps.

> **JSON key:** `checkValve` (EN) · `valvulaRetencao` (PT)

---

## Pumping — Electrical Submersible Pump (ESP)

An ESP adds energy to the flow by imposing a head rise governed by pump performance curves. The key physics:

- **Head-flow curve** — at a reference frequency, maps volumetric flow to differential head.
- **Frequency** — rotational speed (time-dependent for variable-speed drives). Head scales with frequency² and flow with frequency.
- **Number of stages** — total head is per-stage head multiplied by stage count.
- **Pump degradation** — optional efficiency or head degradation factor.

> **JSON key:** `esp` (EN) · `esp` (PT)

---

## Pumping — Positive-Displacement Pump

A volumetric (positive-displacement) pump delivers a fixed volume per revolution regardless of differential pressure. It imposes a flow rate rather than a head rise.

> **JSON key:** `volumetricPump` (EN) · `bombaPD` (PT)

---

## Pumping — Prescribed Pressure Drop/Rise

A simplified model that imposes a fixed pressure increment (or decrement) at a position. Useful for representing an idealized booster or a known restriction without detailed pump curves.

> **JSON key:** `pressureDrop` (EN) · `perdaCarga` (PT)

---

## Pigging

A pig is a mechanical device launched inside the pipe that travels with the flow, sweeping liquid ahead of it. In the model, the pig is tracked as a moving boundary. It creates:

- A liquid slug ahead (displaced accumulation).
- A pressure discontinuity (pig friction against the pipe wall).
- A gas pocket behind (expanding as the pig advances).

The pig is characterized by its launch time, launch position, and friction properties.

> **JSON key:** `pig` (EN) · `pig` (PT)

---

## Practical Guidance

- **Position alignment:** Accessory positions should match the physical architecture (e.g., DHSV at 50 m below wellhead, ESP at pump-set depth).
- **Avoid overlap:** Do not place multiple active restrictions at exactly the same measured length unless physically justified.
- **Incremental build-up:** Introduce one active device at a time during model construction to isolate its effect.
- **Opening ramps:** Avoid instantaneous valve openings (step from 0 to 1); use short ramps (e.g., 10–60 s) for numerical stability during transients.
- **Gas-lift + ESP sequencing:** When combining gas-lift with an ESP, ensure the gas-lift injection point is upstream of the pump intake to avoid injecting gas directly into pump stages.

---

## Example: Liquid Source with Time Schedule

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

## Example: Valve with Opening Ramp

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

## Example: ESP at Pump-Set Depth

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

!!! tip
    When debugging a model with multiple accessories, deactivate all but one (`"active": false`) and validate in isolation before combining.
