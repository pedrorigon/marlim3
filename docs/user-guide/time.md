# Time

Time configuration governs three fundamental aspects of a transient simulation: how the initial state is constructed, how the solver advances through time, and when intermediate results are stored for restart or analysis.

---

## Steady-State vs. Transient Mode

The simulation mode determines whether the solver computes a single equilibrium state or evolves the system through time:

- **Steady-state:** Solves for the equilibrium condition given fixed boundary values. No time evolution — produces spatial profiles only.
- **Transient:** Advances the conservation equations in time, capturing dynamic phenomena (startup, shutdown, slugging, valve events).

> **JSON key:** `transient` (EN) · `transiente` (PT) — inside global config, default `false`

---

## Initialization Strategy

Every transient simulation starts from an initial condition. The strategy determines where this initial state comes from:

### From Steady-State Solution (Strategy 1 — Recommended)

The solver first converges a steady-state solution using the boundary conditions at time zero, then uses that result as the starting point for transient evolution. This ensures physical consistency and is the most common approach.

### From User-Defined State (Strategy 0)

The user provides explicit initial profiles (pressure, temperature, holdup) along the pipe. Useful for academic benchmarks, synthetic test cases, and controlled verification scenarios.

### From Snapshot File (Strategy 2)

Loads a binary restart file (`.snp`) saved from a previous simulation. This enables continuation studies: run a shutdown for 1 hour, save a snapshot, then later resume from that state with different conditions.

A snapshot file name must be provided.

> **JSON key (snapshot file):** `snapshotFile` (EN) · `SnapShotArq` (PT)

### Gas-Lift Unloading (Strategy 3)

A specialized mode for gas-lift unloading analysis. It sets up a system with gas/completion-fluid interfaces at specified measured-length positions in both production and service lines, then evolves the unloading process.

Additional settings become relevant:

- **Gas-line interface position** — where the gas/fluid interface sits in the service line [m].
- **Production-line interface position** — where the gas/fluid interface sits in the production line [m].
- **Completion-fluid salinity** — for property estimation of the completion brine [g/(kg water)].
- **Discharge control** — when enabled, automatically regulates gas injection to avoid erosional velocities at gas-lift valves during unloading.

> **JSON keys (unloading):**
>
> - Gas-line interface: `gasLineInterfaceLength` (EN) · `comprimentoMedidoInterfaceLinhaGas` (PT)
> - Production interface: `prodLineInterfaceLength` (EN) · `comprimentoMedidoInterfaceLinhaProd` (PT)
> - Salinity: `fluidSalinity` (EN) · `SalinidadeFluido` (PT)
> - Control: `dischargeControl` (EN) · `controleDescarga` (PT)

> **JSON key (strategy selector):** `initialCondition` (EN) · `condicaoInicial` (PT) — default `1`
> Values: `0` = user-defined, `1` = steady-state, `2` = snapshot, `3` = unloading

---

## Initial Fluid Assignment

When initializing from user-defined conditions (strategy 0), the production line is filled with a specific fluid from the production-fluid array.

> **JSON key:** `initialFluidId` (EN) · `iniFluidoP` (PT) — default `0`

---

## Simulation Horizon

The total duration of the transient simulation, measured from time zero.

> **JSON key:** `finalTime` (EN) · `tempoFinal` (PT) — unit: seconds, inside the time object

---

## Time-Step Scheduling

The transient solver uses adaptive time stepping bounded by a user-defined maximum. A **piecewise schedule** defines how this maximum evolves over time:

- Between breakpoint `times[i]` and `times[i+1]`, the maximum allowed time step is `maxDT[i]`.
- The actual step may be smaller due to internal stability criteria (CFL condition, holdup oscillation penalization, valve events, sonic constraints, etc.).

**Design principles:**

- Use small `maxDT` during early transients (startup, initial valve movements) to capture fast dynamics.
- Progressively relax `maxDT` as the system approaches quasi-steady behavior.
- Tighten `maxDT` again around planned events (valve operations, pump trips, choke changes).

> **JSON keys (inside time object):**
>
> - Breakpoints: `times` (EN) · `tempos` (PT) — array [s]
> - Maximum step: `maxDT` (EN) · `maxDT` (PT) — array [s]

---

## Snapshot Storage (Restart Files)

At specified time instants, the simulator can write binary snapshot files (`.snp`) containing the complete system state. These files enable:

- Restarting a simulation from an intermediate point (via strategy 2).
- Branching: running multiple scenarios from the same intermediate state.
- Checkpointing long simulations for robustness.

> **JSON key:** `saveSnapshot` (EN) · `salvarSnapshot` (PT) — array of time instants [s], inside time object

---

## Segregation Mode

During extended shutdowns, gas and liquid phases separate vertically (segregation). The physical behavior during segregation differs from normal flowing conditions — the segregation mode adjusts internal model behavior to handle phase-separation dynamics more robustly.

Segregation is controlled by time windows: at each breakpoint, the mode can be activated or deactivated.

> **JSON keys (inside time object):**
>
> - Breakpoints: `segregationTime` (EN) · `tempoSegregacao` (PT) — array [s]
> - Mode flags: `segregation` (EN) · `segregacao` (PT) — array (0 = normal, 1 = segregation active)

---

## Example: Standard Transient (Initialize from Steady-State)

```json
{
  "initialConfig": {
    "transient": true,
    "initialCondition": 1,
    "initialFluidId": 0
  },
  "time": {
    "finalTime": 3600,
    "times": [0, 10, 100],
    "maxDT": [0.1, 1.0, 5.0],
    "saveSnapshot": [1800, 3600]
  }
}
```

## Example: Extended Shutdown with Segregation

```json
{
  "initialConfig": {
    "transient": true,
    "initialCondition": 1,
    "initialFluidId": 0
  },
  "time": {
    "finalTime": 86400,
    "times": [0, 60, 600, 3600],
    "maxDT": [0.5, 2.0, 10.0, 30.0],
    "saveSnapshot": [3600, 14400, 43200, 86400],
    "segregationTime": [0, 600],
    "segregation": [0, 1]
  }
}
```

## Example: Restart from Snapshot

```json
{
  "initialConfig": {
    "transient": true,
    "initialCondition": 2,
    "snapshotFile": "shutdown_3600s.snp"
  },
  "time": {
    "finalTime": 7200,
    "times": [0, 100],
    "maxDT": [1.0, 5.0]
  }
}
```

!!! tip
    A robust default workflow: initialize from steady-state (`initialCondition = 1`), use conservative early `maxDT` (0.1–0.5 s), then progressively relax after the strongest transients (valve closures, pump trips) have passed.
