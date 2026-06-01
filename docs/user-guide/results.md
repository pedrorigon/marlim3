# Results

This section explains how to configure, read, and interpret Marlim3 simulation outputs.

## Output Families

Marlim3 generates three complementary views of the solution:

| Family | What it shows | Use case |
|--------|--------------|----------|
| **Spatial profiles** | State variables vs. measured length at selected times | Pressure drops, temperature gradients, holdup distributions |
| **Temporal trends** | State variables vs. time at selected positions | Startup dynamics, oscillations, control response |
| **Cross/radial outputs** | State across pipe wall layers at selected positions/times | Thermal diagnostics, insulation performance |

## Spatial Profiles

Profile outputs produce snapshots of the solution along the entire pipe at specified time instants.

### Configuration Objects

| Object | Line |
|--------|------|
| `productionProfile` | Production line |
| `serviceProfile` | Service/gas-injection line |

### Common Output Variables

| Variable | Physical meaning |
|----------|-----------------|
| `pressure` | Mixture pressure along measured length |
| `temperature` | Mixture temperature |
| `holdup` | Liquid holdup (volume fraction) |
| `flowPattern` | Flow regime indicator |
| `mixtureVelocity` | Mixture superficial velocity |
| `liquidVelocity` | Liquid phase velocity |
| `gasVelocity` | Gas phase velocity |
| `liquidDensity` | In-situ liquid density |
| `gasDensity` | In-situ gas density |
| `frictionPressureGradient` | Frictional pressure-loss gradient |
| `hydrostaticPressureGradient` | Hydrostatic pressure gradient |

### Configuration Pattern

```json
{
  "productionProfile": {
    "time": [0, 1800, 3600],
    "pressure": true,
    "temperature": true,
    "holdup": true,
    "flowPattern": true,
    "frictionPressureGradient": true,
    "hydrostaticPressureGradient": true
  }
}
```

The `time` array specifies instants at which profiles are written. Set variable flags to `true` to include them in output.

## Temporal Trends

Trend outputs record solution evolution over time at fixed measured-length positions.

### Configuration Objects

| Object | Line |
|--------|------|
| `productionTrend` | Production line |
| `serviceTrend` | Service line |

### Configuration Pattern

```json
{
  "productionTrend": {
    "measuredLength": [0, 500, 1000, 2500],
    "pressure": true,
    "temperature": true,
    "holdup": true,
    "liquidFlowRate": true
  }
}
```

## Cross/Radial Outputs

Used for thermal diagnostics across the pipe wall and formation coupling:

| Object | Spatial/Temporal | Line |
|--------|-----------------|------|
| `crossProductionProfile` | Radial state vs. measured length at fixed times | Production |
| `crossServiceProfile` | Radial state vs. measured length at fixed times | Service |
| `crossProductionTrend` | Radial state vs. time at fixed positions | Production |
| `crossServiceTrend` | Radial state vs. time at fixed positions | Service |

These outputs show temperature (and potentially other variables) across each concentric layer, useful for:

- Verifying insulation performance
- Identifying layer-specific thermal delays
- Validating formation coupling

## Interpreting Results in Practice

### Validation Checklist

1. **Pressure monotonicity:** Verify pressure decreases monotonically from inlet to outlet (for production). Abrupt drops should correspond to known restrictions (valves, ESPs, chokes).
2. **Thermal consistency:** Temperature profiles should be bounded between fluid inlet temperature and ambient temperature. Check that insulation sections show flatter gradients.
3. **Mass balance:** At sources/sinks, verify that flow-rate changes match the specified source schedules.
4. **Steady-state baseline:** Always compare transient evolution against the initial steady-state to ensure changes are physically plausible.
5. **Holdup behavior:** Liquid holdup should increase in uphill sections and decrease in downhill sections. Jumps indicate flow-pattern transitions.

### Common Diagnostic Patterns

| Symptom | Possible cause | Suggested action |
|---------|----------------|------------------|
| Abrupt pressure spikes | Overly aggressive time step or sharp boundary changes | Reduce `maxDT` around event times |
| Unrealistic thermal oscillation | Inconsistent material/formation properties or too few radial nodes | Check material IDs, increase `discretization` |
| Nonphysical holdup switching | Correlation mismatch or unstable control schedules | Verify `correlationsByPattern`, smooth valve ramps |
| Diverging solution | Time step too large for the physical process | Reduce initial `maxDT`, check `simplePressureFrontier` |
| Unexpected reverse flow | Missing check valve or incorrect separator pressure | Add `checkValve = 1` or verify outlet pressure |

## Python API Access

When using the Python package, results are accessed directly:

```python
branch.simulate()

# Plot spatial profiles
branch.plot_profiles()

# Access raw output data
results = branch.results
```

!!! tip
    Always establish a validated steady-state reference before interpreting complex transient behavior. If transient results seem unreasonable, first verify the steady-state solution is correct.
