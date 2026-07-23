# Rearchitect `SisProd` into a modular, testable and maintainable production solver

## Problem

`src/core/SisProd.cpp` (the class `SProd`) is the core of the Marlim3 flow
simulator, yet it has grown into a monolith that is expensive to maintain,
impossible to unit test and risky to evolve.

Concrete, measurable facts about the current state:

- The implementation file has roughly **26,000 lines** and the header roughly
  **1,500 lines**; it is by far the largest translation unit in `src/core`.
- The class `SProd` concentrates around **150 methods** and about **twelve
  distinct responsibilities**: input parsing and assembly, steady-state solving
  (production, gas and injection), transient evolution, mass/composition
  marching, thermal/energy marching, two-phase flow closure correlations, fluid
  property updates, equipment control (chokes, gas-lift, Master1 valves, pigs),
  generic root finding, and output/observability.
- **Every member is public** (the header has a single `public:` section). The
  orchestrator `src/core/Num4Main.cpp` reaches into the internal state in
  **700+ places**, including direct writes of physical state in the middle of a
  time step. This makes any change to the internal representation a breaking
  change.
- Several methods are enormous. The mass march (`RenovaMassPerm`) alone is about
  **1,600 lines** and dispatches on the accessory type (`acsr.tipo`) with nearly
  identical ~300-line blocks per case, replicated across `Rev`, `Comp` and
  `CompRev` variants (roughly a 4x duplication).
- Generic numerical root finders (`zbrent`, `zriddr`, `falsacorda`, `SIGN`) are
  **duplicated** inside `SProd` and in `FerramentasNumericas.cpp`.
- Error handling terminates the process: `NumError()` logs and calls `exit()`
  in 150+ places, and convergence failures are only signalled through `cout`.
  The engine is therefore not usable as a testable library.
- There are **no C++ unit tests** for the solver; the only safety net is the
  Python end-to-end regression suite that runs the compiled executable and
  compares numeric outputs.
- Resource ownership is manual (`new`/`delete`, a ~120-line destructor, a
  ~300-line deep-copy `operator=`), which is error prone.
- Shared mutable global state (`varGlob1D* vg1dSP`, the `lixo5` clock, a global
  `logger`) hinders reentrancy, deterministic testing and any future
  parallelism.

Expected behaviour after the change: the production results must remain
identical (within an agreed regression tolerance) while the code becomes
modular, encapsulated, unit-testable, observable and easier to evolve, and while
performance is preserved or improved.

## Proposed solution

Introduce a new, modular `SisProd` architecture and migrate to it incrementally
(Strangler Fig), keeping the legacy solver available as the reference until the
new one is a verified drop-in replacement.

### File organisation

- Preserve the legacy implementation intact by moving it to
  `src/core/SisProd_old.cpp` (it keeps defining `class SProd` declared in the
  unchanged `src/include/SisProd.h`, and remains the production solver and the
  integration-test reference).
- Place the new architecture in `src/core/SisProd.cpp` with its interface in
  `src/include/SisProd2.h`.
- Keep the header limited to declarations, small value types, interfaces and
  code that genuinely must be inline or template; move the actual logic into the
  `.cpp` for better encapsulation, readability and build time.
- Do not guard the new implementation with a standalone `#ifdef`/`main`. Put any
  test entry point in a separate, non-production translation unit
  (`tests/sisprod2_selftest.cpp`).

### Implementation-selection flag

- Add an explicit build option `MARLIM_USE_NEW_SISPROD` (CMake, default `OFF`)
  that defines a matching compile macro and lets the system be built with the
  legacy solver (`SisProd_old.cpp`) or the new one (`SisProd.cpp`).
- Expose `marlim::sisprod2::usingNewSisProd()` /
  `activeSisProdImplementation()` so the active implementation is observable and
  testable. Keep the default `OFF` so all existing flows and regression tests
  keep passing until the new solver is complete.

### Target components (new architecture)

Separate the twelve legacy responsibilities into small, cohesive units with
explicit ownership and dependency injection:

- Typed results (`SolveResult`, `StepResult`) that replace `cout`/`exit`.
- `SimContext` encapsulating the shared clock and configuration (replaces
  `vg1dSP`/`lixo5`).
- `Diagnostics` (null and stream implementations) for structured logging.
- `RootFinder` (Brent) injected by `std::function`, plus an independent
  bisection reference, unifying the duplicated `zbrent`/`zriddr`.
- `StandardStream` + `StreamMixing` (Strategy) that deduplicate the accessory
  switch of `RenovaMassPerm`.
- `FluidModel` and closure correlations as pure, testable kernels.
- `ProductionColumn` and a `TramoEngine` facade whose steady-state solve mirrors
  `buscaProdPfundoPerm` calling `marchaProdPerm`.
- `TrendBuffer` (RAII) replacing the `double***` buffers and manual deletes.
- A batch solver using OpenMP across independent columns (the per-cell march is
  inherently sequential and must not be parallelised).

### Dual-run comparison

- Provide a harness that runs the same flow through two executables (legacy vs
  new) and compares the numeric outputs within a tolerance
  (`scripts/compare_sisprod_impls.py`), so parity can be checked automatically as
  the new solver is wired in.

### Behavioural parity requirement

The new `SisProd.cpp` must reproduce the full behaviour of `SisProd_old.cpp`, not
merely resemble its structure. Because the new implementation is significantly
smaller, do not assume the migration is complete. Systematically compare the
legacy and the new implementation and ensure the new one covers **all** of the
following before it can replace the legacy solver behind the flag:

- Every solving flow: steady-state (production, gas, injection), transient time
  stepping, reverse flow (`*Rev`), compositional mode (`arq.flashCompleto == 2`,
  `*Comp`), water injection (`injPoc`), and every network topology (series,
  parallel, gas-lift ring, injection).
- Every business rule and special case currently guarded in the legacy code
  (accessory types, boundary conditions, choke operating modes, Master1 valve
  logic, pig tracking, gas-lift unloading, rollback via `reinicia`, moving
  averages, dynamic property tables, parkaffin/hydrate coupling).
- Error handling: the legacy `NumError`/`exit` semantics at the system boundary
  must be preserved (or converted through an adapter that keeps observable
  behaviour), while the core surfaces typed, propagable errors.
- Integrations: the Fortran black-oil/PVT routines, the OpenMP parallel regions,
  the `.mr3` snapshot format, and the trend/log outputs consumed by the GUI and
  the regression tooling.

Any missing functionality, incomplete integration, incompatibility or regression
found during this comparison must be implemented and fixed.

### Validation

- The full project must build with the new module integrated (both flag states
  must compile).
- The existing test suite must pass unchanged with the default flag.
- Add C++ unit tests for the new pure kernels (root finder vs bisection,
  stream mixing conservation, steady-state vs the closed-form hydrostatic
  reference, determinism, batch parallel == sequential) and integrate them into
  the Python test suite.
- Use the dual-run comparison to assert numeric parity between the legacy and the
  new solver on representative `.mr3` cases, within a **relative tolerance of
  1e-6** (absolute tolerance 1e-9), before flipping the default flag.

## Alternatives

- **Full rewrite in a single step.** Rejected: with ~26,000 lines of coupled
  physics, Fortran coupling and no unit tests, a big-bang rewrite cannot be
  verified and would very likely introduce numeric regressions. An incremental
  migration behind regression and dual-run checks is the only sustainable path.
- **Leave the monolith as is.** Rejected: the maintenance cost, the absence of
  unit tests and the tight coupling to `Num4Main` keep making changes slow and
  risky.
- **Refactor in place inside `SisProd.cpp`.** Rejected as the primary approach:
  it does not give a clean, testable module boundary and makes it hard to run and
  compare the two behaviours side by side. Extracting a new module while keeping
  the legacy file intact preserves a reference and enables the dual-run
  comparison.

## Additional context

- The project targets **C++11** (`CMAKE_CXX_STANDARD 11`, extensions off), so the
  new code cannot rely on `std::expected`, `std::optional` or
  `std::filesystem`; use a small result type and `std::function`.
- The build compiles a single executable via `file(GLOB ... *.cpp)` and links
  C++ with Fortran, so any new `.cpp` under `src/core` is compiled into the
  binary and must not introduce symbol clashes or a second `main`.
- The existing tests are Python end-to-end/regression tests that run the
  compiled `Marlim3` executable (`-s TRANSIENTE|INJETOR|REDE -i ... -p ... -d ...
  -o ...`) and compare outputs; they are the primary safety net for parity.
- Detailed background, the responsibility map, the phased plan and the
  comparison of the current and proposed architectures are documented under
  `docs-sisprod/` and `docs/dev-guide/SisProd-analise-arquitetural.md`.
- Open decisions that must be settled with the team and the domain experts:
  representative `ncel` and network sizes (to decide where parallelism pays off),
  which output contracts are fixed (`.mr3`, trends, logs), the meaning of magic
  constants (`1e-15`, `RGOMax`, choke thresholds), whether `exit()` on error is
  required by external scripts, and whether raising the language standard beyond
  C++11 is allowed. The regression parity tolerance is set to **1e-6 relative**.
- A living, checkpointed execution plan is maintained in
  `issues/sisprod-migration-plan.md`; it tracks per-region progress so the work
  can be resumed at any point.
