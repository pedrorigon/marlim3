# Time

Time configuration controls three things:

- How the initial transient state is built.
- How the maximum time-step cap evolves.
- When snapshots are saved for restart.

---

## Simulation Mode

`transient` switches between steady-state and transient solve.

> **JSON key:** `initialConfig.transient` (EN) · `configuracaoInicial.transiente` (PT)
> Default: `false`

---

## Initial Condition Strategy

`initialCondition` chooses the startup method:

- `0`: user-defined initial fields
- `1`: initialize from steady-state solution
- `2`: initialize from snapshot file
- `3`: gas-lift unloading initialization

> **JSON key:** `initialConfig.initialCondition` (EN) · `configuracaoInicial.condicaoInicial` (PT)
> Default: `1`

### Strategy 0: User-Defined Initial State

Use explicit initial fields in pipe cells. The initial production fluid is selected with:

> **JSON key:** `initialConfig.initialFluidId` (EN) · `configuracaoInicial.iniFluidoP` (PT)

### Strategy 2: Snapshot Restart

Loads a `.snp` state file.

> **JSON key:** `initialConfig.snapshotFile` (EN) · `configuracaoInicial.SnapShotArq` (PT)

Parser behavior: if mode 2 is selected and file does not exist, parsing fails.

### Strategy 3: Gas-Lift Unloading

Requires unloading interface information:

- `gasLineInterfaceLength` · `comprimentoMedidoInterfaceLinhaGas`
- `prodLineInterfaceLength` · `comprimentoMedidoInterfaceLinhaProd`
- `fluidSalinity` · `SalinidadeFluido`

Optional control:

- `dischargeControl` · `controleDescarga`
- `dischargeParameters` · `parametrosDescarga`

---

## Time Object

> **JSON key:** `time` (EN/PT)

Core fields:

- `finalTime` · `tempoFinal` [s]
- `times` · `tempos` [s]
- `maxDT` · `dtmax` [s]

`times` and `maxDT` define a piecewise max-step schedule. The solver may still take smaller steps for stability.

Parser checks:

- `times` must be increasing.
- `times` and `maxDT` must have the same size.

Parser defaults when omitted:

- `times = [0]`
- `maxDT = [5]`

---

## Segregation Scheduling

For shutdown/segregation scenarios, mode windows can be scheduled in the time object.

> **JSON keys:**
>
> - `segregationTime` (EN) · `tempoSegrega` (PT)
> - `segregation` (EN) · `segrega` (PT)

Values:

- `0` = normal mode
- `1` = segregation mode

Parser behavior: if these arrays are omitted, it creates a default segmentation schedule at `t=0` with mode `1`.

---

## Snapshot Save Times

Snapshot-write times are defined in `time`:

> **JSON key:** `saveSnapshot` (EN) · `gravaMomento` (PT)

Each entry is a simulation time [s] to write a snapshot file for future restart.

---

## Example: Standard Transient

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

## Example: Segregation Window

```json
{
  "time": {
    "finalTime": 86400,
    "times": [0, 60, 600, 3600],
    "maxDT": [0.5, 2.0, 10.0, 30.0],
    "segregationTime": [0, 600],
    "segregation": [0, 1]
  }
}
```
