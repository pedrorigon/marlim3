# Time

Time configuration controls how the simulator initializes state, advances in time, and stores intermediate results.

## Concept

Transient simulations require two fundamental decisions:

1. **Initialization:** How to define the starting state (from scratch, from steady-state, or from a previous snapshot).
2. **Time-stepping:** How large each time step is allowed to be during evolution, balancing accuracy against runtime.

These decisions strongly influence stability, runtime, and physical fidelity of results.

## Core Fields

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `initialConfig.transient` | bool | `false` | `true` = transient simulation; `false` = steady-state only |
| `initialConfig.initialCondition` | int | `1` | Initialization strategy selector (see below) |
| `initialConfig.initialFluidId` | int | `0` | Fluid ID that initially fills the production line |
| `time.finalTime` | number | — | End of simulation horizon [s] |

## Initialization Strategies

| Code | Strategy | When to use |
|------|----------|-------------|
| `0` | User-defined initial state | Controlled academic/synthetic scenarios with known initial profiles |
| `1` | Start from steady-state solution | **Most common.** Natural starting point for transient field studies |
| `2` | Load from snapshot file (`.snp`) | Restart and continuation studies from a previous simulation |
| `3` | Gas-lift unloading setup | Specialized unloading analyses with gas/liquid interface positions |

### Strategy 0: User-Defined

Requires specifying initial pressure, temperature, and holdup profiles along the pipe. Useful for benchmarking and verification.

### Strategy 1: Steady-State Initialization (Recommended)

The solver first converges a steady-state solution with the initial boundary conditions, then uses that as the transient starting point. This ensures physical consistency.

### Strategy 2: Snapshot Restart

Loads a binary `.snp` file (saved from a previous run) to continue a simulation. Requires `initialConfig.snapshotFile` to specify the file name.

### Strategy 3: Gas-Lift Unloading

Specialized mode that sets up a gas/completion-fluid interface at specified positions for unloading simulation. Additional fields become relevant:

| Field | Unit | Meaning |
|-------|------|---------|
| `initialConfig.gasLineInterfaceLength` | m | Interface position in service line (from platform) |
| `initialConfig.prodLineInterfaceLength` | m | Interface position in production line (from platform) |
| `initialConfig.fluidSalinity` | g/(kg water) | Completion-fluid salinity |
| `initialConfig.dischargeControl` | bool | Enable automatic erosion-velocity-limited unloading |

## Time-Step Scheduling

The `time` object defines a piecewise schedule for the maximum allowed time step:

| Field | Unit | Description |
|-------|------|-------------|
| `time.times` | s | Breakpoints defining time-step schedule |
| `time.maxDT` | s | Maximum allowed Δt in each interval |
| `time.finalTime` | s | Simulation end time |

The schedule works as follows: between `times[i]` and `times[i+1]`, the maximum time step is `maxDT[i]`. The actual step may be smaller due to internal stability criteria (CFL, holdup oscillation, valve events, etc.).

**Design principles:**

- Use small `maxDT` early to capture fast startup transients.
- Relax `maxDT` at later times when the system approaches quasi-steady behavior.
- Use small `maxDT` around planned events (valve operations, pump start/stop).

## Snapshot Storage

`time.saveSnapshot` defines time instants at which the simulator writes `.snp` restart files:

```json
"saveSnapshot": [1800, 3600, 7200]
```

These files can be loaded later via `initialCondition = 2` for continuation runs.

## Segregation Controls

For shutdown/segregation studies, specialized time windows control model behavior:

| Field | Unit | Description |
|-------|------|-------------|
| `time.segregationTime` | s | Time breakpoints for segregation control |
| `time.segregation` | int | `0` = normal mode; `1` = segregation mode active |

Segregation mode adjusts internal model behavior for phase-separation dynamics during extended shutdowns.

## Example JSON

### Standard Transient (Initialize from Steady-State)

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

### Extended Shutdown with Segregation

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

### Restart from Snapshot

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
