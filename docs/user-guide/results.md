# Results

Simulation outputs in Marlim3 provide three complementary views of the solution: spatial profiles along the pipe, temporal trends at fixed locations, and radial/cross-section outputs through the pipe wall. Configuring what is written and when enables efficient post-processing without generating excessive data.

---

## Output Philosophy

Marlim3 does not dump all variables at every time step — that would produce unmanageable files for long transients. Instead, the user specifies:

1. **Which variables** to record (pressure, temperature, holdup, velocities, etc.).
2. **Where and when** to record them (time instants for profiles; measured-length positions for trends).

This design keeps output files focused and manageable while giving complete control over what is captured.

---

## Spatial Profiles

A spatial profile is a snapshot of the solution along the entire pipe length at a specific time instant. These are ideal for understanding:

- Pressure drop distribution along the system.
- Temperature gradients and insulation effectiveness.
- Liquid holdup patterns and flow-regime transitions.
- Friction and hydrostatic contributions to pressure loss.

### When are profiles written?

At each time instant listed in the profile time array. For steady-state simulations, the converged solution is written as a single profile.

> **JSON keys (production line):** `productionProfile` (EN) · `perfilProducao` (PT)
> **JSON keys (service line):** `serviceProfile` (EN) · `perfilServico` (PT)

### Selecting Variables

Each variable is toggled individually. Common choices:

| Variable concept | EN key | PT key |
|-----------------|--------|--------|
| Pressure | `pressure` | `pressao` |
| Temperature | `temperature` | `temperatura` |
| Liquid holdup | `holdup` | `holdup` |
| Flow pattern | `flowPattern` | `padraoEscoamento` |
| Mixture velocity | `mixtureVelocity` | `velocidadeMistura` |
| Liquid velocity | `liquidVelocity` | `velocidadeLiquido` |
| Gas velocity | `gasVelocity` | `velocidadeGas` |
| Liquid density | `liquidDensity` | `densidadeLiquido` |
| Gas density | `gasDensity` | `densidadeGas` |
| Friction pressure gradient | `frictionPressureGradient` | `gradientePressaoFriccao` |
| Hydrostatic pressure gradient | `hydrostaticPressureGradient` | `gradientePressaoHidrostatico` |

### Specifying Profile Times

The time array defines which time instants trigger output. For steady-state, use `[0]`.

> **JSON key (inside profile):** `time` (EN) · `tempo` (PT) — array [s]

---

## Temporal Trends

A temporal trend records the evolution of selected variables over time at fixed measured-length positions. These are ideal for understanding:

- Startup and shutdown dynamics.
- Pressure/temperature oscillations and slug frequency.
- Arrival times of thermal or pressure fronts.
- Control response to valve or pump events.

### Where are trends recorded?

At each measured-length position specified in the trend's position array. Choose positions at key locations: inlet, outlet, valve stations, ESP intake, downhole gauge positions, etc.

> **JSON keys (production line):** `productionTrend` (EN) · `tendenciaProducao` (PT)
> **JSON keys (service line):** `serviceTrend` (EN) · `tendenciaServico` (PT)

### Selecting Positions

> **JSON key (inside trend):** `measuredLength` (EN) · `comprimentoMedido` (PT) — array [m]

The same variable toggle flags apply as for profiles.

---

## Radial / Cross-Section Outputs

Cross-section outputs provide temperature (and other variables) across each concentric wall layer at selected positions and times. They are useful for:

- Verifying insulation performance (temperature drop across insulation layer).
- Identifying layer-specific thermal delays during transients.
- Validating formation coupling and cooldown behavior.
- Diagnosing hydrate or wax formation risk at inner-wall surfaces.

Two forms exist — spatial (radial profile vs. measured length at fixed times) and temporal (radial profile vs. time at fixed positions):

> **JSON keys:**
>
> - Radial along pipe: `crossProductionProfile` (EN) · `perfilCruzadoProducao` (PT)
> - Radial along pipe (service): `crossServiceProfile` (EN) · `perfilCruzadoServico` (PT)
> - Radial vs. time: `crossProductionTrend` (EN) · `tendenciaCruzadaProducao` (PT)
> - Radial vs. time (service): `crossServiceTrend` (EN) · `tendenciaCruzadaServico` (PT)

---

## Interpreting Results

### Validation Checklist

1. **Pressure monotonicity:** Pressure should decrease monotonically from inlet to outlet (for production). Abrupt drops correspond to restrictions (valves, ESPs, chokes).
2. **Thermal consistency:** Temperature profiles must be bounded between fluid inlet temperature and ambient temperature. Insulated sections should show flatter gradients.
3. **Mass balance:** At sources/sinks, flow-rate changes must match specified schedules.
4. **Steady-state baseline:** Always compare transient evolution against the initial steady-state.
5. **Holdup behavior:** Liquid holdup increases in uphill sections and decreases downhill. Jumps indicate flow-pattern transitions.

### Common Diagnostic Patterns

| Symptom | Possible cause | Suggested action |
|---------|----------------|------------------|
| Abrupt pressure spikes | Overly aggressive time step or sharp boundary changes | Reduce `maxDT` around event times |
| Unrealistic thermal oscillation | Inconsistent material/formation properties | Check material IDs, review cross-section geometry |
| Nonphysical holdup switching | Correlation mismatch or unstable control schedules | Verify flow-pattern correlations, smooth valve ramps |
| Diverging solution | Time step too large for the physical process | Reduce initial `maxDT` |
| Unexpected reverse flow | Missing check valve or incorrect separator pressure | Add check valve or verify outlet pressure |

---

## Python API Access

When using the Python package, results are accessed programmatically:

```python
branch.simulate()

# Plot spatial profiles
branch.plot_profiles()

# Access raw output data
results = branch.results
```

---

## Example: Complete Output Configuration

```json
{
  "productionProfile": {
    "time": [0, 1800, 3600],
    "pressure": true,
    "temperature": true,
    "holdup": true,
    "flowPattern": true,
    "frictionPressureGradient": true
  },
  "productionTrend": {
    "measuredLength": [0, 500, 1000, 2500],
    "pressure": true,
    "temperature": true,
    "holdup": true,
    "liquidFlowRate": true
  },
  "crossProductionProfile": {
    "time": [0, 3600],
    "temperature": true
  }
}
```

!!! tip
    Always establish a validated steady-state reference before interpreting complex transient behavior. If transient results seem unreasonable, first verify the steady-state solution is correct.
