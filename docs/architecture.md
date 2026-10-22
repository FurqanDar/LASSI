# Dependency direction (Phase 3)

Allowed include / call dependency (arrows point toward lower-level modules):

```
lattice / structure  ←  energy  ←  moves  ←  analysis / I/O  ←  main / parse
```

In practice:

| Layer | Modules | May include |
|-------|---------|-------------|
| Lattice / structure | `structure*.c`, `structure.h` | `global.h` only (plus peer structure units via `structure.h`) |
| Energy | `energy.c`, `energy.h` | structure, `global.h` |
| Moves | `mcmove*.c`, `mcmove.h` | energy, structure, `cluster_mc.h` (move-facing cluster queries only), `global.h` |
| Analysis / I/O | `cluster.*`, `print*.c`, `print.h` | moves (if needed), energy, structure, `global.h` |
| Main / parse | `main.c`, `parsekey.*`, `initialize.*` | any layer above |

**Must not:** a lower layer `#include` a higher-layer header (e.g. structure must not include `print.h` / `mcmove.h` / `energy.h`; energy must not include `mcmove.h` / `print.h`; moves should not include `print.h` / `initialize.h`).

`global.h` / `globals.c` are process-wide state shared by all layers (not a dependency edge in the sense above).

## Move-facing cluster API

Cluster MC moves need a few cluster-finding helpers without depending on the full
analysis-layer `cluster.h`. Those decls live in **`cluster_mc.h`** (implemented
in `cluster.c`). `cluster.h` includes `cluster_mc.h` so analysis callers are
unchanged. `mcmove_moves.c` includes only `cluster_mc.h`.

## Known exceptions (deferred)

1. **`structure_analysis.c` → `cluster.h`** — RDF / radial helpers use cluster utilities. The file is analysis-flavored despite the `structure_` prefix; a later move/rename may place it with analysis.

Sanity checks in `structure_sanity.c` call I/O printers via **forward declarations** (not `#include "print.h"`) so structure does not take a hard include edge into I/O.

`print_fileio.c` forward-declares `Create2DInt` from initialize rather than `#include "initialize.h"`, avoiding an I/O → main/parse include edge.

## How to extend

When touching a `.c` file, prefer removing unused upward includes and avoid adding new ones. Characterization (`make characterization`) stays the behavior gate.
