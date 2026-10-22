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
| Moves | `mcmove*.c`, `mcmove.h` | energy, structure, `global.h` |
| Analysis / I/O | `cluster.*`, `print*.c`, `print.h` | moves (if needed), energy, structure, `global.h` |
| Main / parse | `main.c`, `parsekey.*`, `initialize.*` | any layer above |

**Must not:** a lower layer `#include` a higher-layer header (e.g. structure must not include `print.h` / `mcmove.h` / `energy.h`; energy must not include `mcmove.h` / `print.h`; moves should not include `print.h` / `initialize.h`).

`global.h` / `globals.c` are process-wide state shared by all layers (not a dependency edge in the sense above).

## Known exceptions (deferred to 3.2+)

These remain after Phase 3.1 because fixing them needs more than an include trim:

1. **`mcmove_moves.c` → `cluster.h`** — cluster MC moves call cluster-finding helpers. Clustering is shared by moves and analysis; untangling ownership is follow-up work.
2. **`structure_analysis.c` → `cluster.h`** — RDF / radial helpers use cluster utilities. The file is analysis-flavored despite the `structure_` prefix; a later move/rename may place it with analysis.

Sanity checks in `structure_sanity.c` call I/O printers via **forward declarations** (not `#include "print.h"`) so structure does not take a hard include edge into I/O.

`print_fileio.c` forward-declares `Create2DInt` from initialize rather than `#include "initialize.h"`, avoiding an I/O → main/parse include edge.

## How to extend

When touching a `.c` file, prefer removing unused upward includes and avoid adding new ones. Characterization (`make characterization`) stays the behavior gate.
