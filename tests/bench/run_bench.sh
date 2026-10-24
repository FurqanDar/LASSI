#!/usr/bin/env bash
# Performance bench harness for LASSI.
#
# Builds via src/compile.sh (-O3 Release-ish), runs the pinned bench deck
# N times, prints median wall time, and optionally compares against the
# checked-in baseline (fail if regression exceeds threshold).
#
# One command (preferred):
#   make bench
#   make bench-compare
#   make bench-update
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DECK_NAME="bench_diblock"
DECK_DIR="${SCRIPT_DIR}/decks/${DECK_NAME}"
BASELINE_FILE="${SCRIPT_DIR}/baseline.json"
WORK_ROOT="${BENCH_WORK_DIR:-${SCRIPT_DIR}/.work}"

# Pin locale / timezone for stable log text (timings are wall-clock).
export LC_ALL=C
export LANG=C
export TZ=UTC
unset LC_NUMERIC || true

RUNS="${BENCH_RUNS:-3}"
THRESHOLD_PCT="${BENCH_THRESHOLD_PCT:-15}"
MODE="run" # run | compare | update-baseline

usage() {
  cat >&2 <<'EOF'
Usage: run_bench.sh [--compare | --update-baseline]

  (default)           Build Release-ish binary, run bench deck N times, print median.
  --compare           Same, then fail (exit 1) if median regresses beyond threshold
                      vs tests/bench/baseline.json (default 15%; BENCH_THRESHOLD_PCT).
  --update-baseline   Same, then write tests/bench/baseline.json with this median.

Environment:
  LASSI_BIN             Prebuilt binary (skip src/compile.sh)
  BENCH_RUNS            Number of timed runs (default 3)
  BENCH_THRESHOLD_PCT   Allowed regression percent for --compare (default 15)
  BENCH_WORK_DIR        Scratch directory (default tests/bench/.work)

Note: cycle counts via `perf` are local-optional and not recorded by this harness.
EOF
}

if [[ $# -gt 1 ]]; then
  usage
  exit 2
fi
if [[ $# -eq 1 ]]; then
  case "$1" in
    --compare) MODE="compare" ;;
    --update-baseline) MODE="update-baseline" ;;
    -h|--help) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
fi

if ! [[ "${RUNS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: BENCH_RUNS must be a positive integer (got '${RUNS}')" >&2
  exit 2
fi

build_lassi() {
  if [[ -n "${LASSI_BIN:-}" ]]; then
    if [[ ! -x "${LASSI_BIN}" ]]; then
      echo "ERROR: LASSI_BIN=${LASSI_BIN} is not executable" >&2
      exit 1
    fi
    echo "Using LASSI_BIN=${LASSI_BIN}"
    return
  fi

  echo "Building lassi via src/compile.sh (Release-ish: -march=x86-64 -O3) ..."
  (cd "${REPO_ROOT}/src" && bash compile.sh)
  LASSI_BIN="${REPO_ROOT}/src/lassi"
}

require_pinned_seed() {
  local keyfile="$1"
  local seed
  seed="$(
    awk '
      /^[[:space:]]*#/ { next }
      $1 == "RANDOM_SEED" {
        print $2
        found=1
        exit
      }
      END { if (!found) exit 1 }
    ' "${keyfile}"
  )" || {
    echo "ERROR: ${keyfile} has no RANDOM_SEED keyword" >&2
    exit 1
  }
  if [[ -z "${seed}" || "${seed}" == "0" ]]; then
    echo "ERROR: ${keyfile}: bench deck must pin a non-zero RANDOM_SEED" >&2
    exit 1
  fi
  printf '%s' "${seed}"
}

# Wall time in seconds with fractional part; prefers date +%s.%N, else python3.
now_seconds() {
  if date +%s.%N 2>/dev/null | grep -Eq '^[0-9]+\.[0-9]+$'; then
    date +%s.%N
    return
  fi
  python3 -c 'import time; print(f"{time.time():.9f}")'
}

median_of() {
  # args: float times; stdout: median (average of two middle values if even)
  python3 -c '
import sys
vals = sorted(float(x) for x in sys.argv[1:])
if not vals:
    raise SystemExit("no values for median")
n = len(vals)
if n % 2:
    print(f"{vals[n // 2]:.6f}")
else:
    print(f"{(vals[n // 2 - 1] + vals[n // 2]) / 2.0:.6f}")
' "$@"
}

run_once() {
  local work_dir="$1"
  local seed="$2"
  local start end elapsed rc

  rm -rf "${work_dir}"
  mkdir -p "${work_dir}"
  cp "${DECK_DIR}/"* "${work_dir}/"

  start="$(now_seconds)"
  set +e
  (
    cd "${work_dir}"
    "${LASSI_BIN}" -k param.key >run.log 2>&1
  )
  rc=$?
  set -e
  end="$(now_seconds)"
  elapsed="$(python3 -c "print(f'{float('${end}') - float('${start}'):.6f}')")"

  if [[ "${rc}" -ne 0 ]]; then
    echo "ERROR: lassi exited ${rc}" >&2
    tail -n 40 "${work_dir}/run.log" >&2 || true
    exit 1
  fi
  if ! grep -Eq "RNG Seed[[:space:]]*=[[:space:]]*${seed}[[:space:]]*$" "${work_dir}/run.log"; then
    echo "ERROR: run log does not report RNG Seed = ${seed}" >&2
    grep -n "RNG Seed" "${work_dir}/run.log" >&2 || true
    exit 1
  fi
  printf '%s' "${elapsed}"
}

write_baseline() {
  local median="$1"
  shift
  local -a times=("$@")
  local times_json host uname_s compiler build_flags bin_sha recorded_at seed
  local n_steps preeq

  times_json="$(python3 -c 'import json,sys; print(json.dumps([float(x) for x in sys.argv[1:]]))' "${times[@]}")"
  host="$(hostname -f 2>/dev/null || hostname || echo unknown)"
  uname_s="$(uname -srm 2>/dev/null || echo unknown)"
  compiler="$(gcc --version 2>/dev/null | head -n1 || echo unknown)"
  build_flags='gcc -march=x86-64 -O3 -std=gnu11 ... -lm -fcommon (src/compile.sh)'
  bin_sha="$(sha256sum "${LASSI_BIN}" | awk '{print $1}')"
  recorded_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  seed="$(require_pinned_seed "${DECK_DIR}/param.key")"
  n_steps="$(awk '/^[[:space:]]*N_STEPS/{print $2; exit}' "${DECK_DIR}/param.key")"
  preeq="$(awk '/^[[:space:]]*PREEQ_STEPS/{print $2; exit}' "${DECK_DIR}/param.key")"

  python3 - "${BASELINE_FILE}" <<PY
import json, sys
path = sys.argv[1]
doc = {
    "deck": "${DECK_NAME}",
    "median_seconds": float("${median}"),
    "runs": int("${RUNS}"),
    "run_seconds": json.loads('${times_json}'),
    "threshold_pct_default": float("${THRESHOLD_PCT}"),
    "seed": int("${seed}"),
    "n_steps": int("${n_steps}"),
    "preeq_steps": int("${preeq}"),
    "mc_cycle_num": 1,
    "recorded_at_utc": "${recorded_at}",
    "host": "${host}",
    "uname": "${uname_s}",
    "compiler": "${compiler}",
    "build_flags": "${build_flags}",
    "lassi_sha256": "${bin_sha}",
    "metric": "wall_clock_median_seconds",
    "notes": (
        "Wall-clock median of N timed runs of decks/bench_diblock with "
        "src/compile.sh (-O3). Trajectory/energy writers disabled. "
        "Hardware-bound; compare locally or via optional CI. "
        "Cycle counts via perf are local-optional and not recorded here."
    ),
}
with open(path, "w", encoding="utf-8") as f:
    json.dump(doc, f, indent=2)
    f.write("\n")
print(f"Wrote baseline: {path}")
PY
}

compare_baseline() {
  local median="$1"
  if [[ ! -f "${BASELINE_FILE}" ]]; then
    echo "ERROR: missing baseline ${BASELINE_FILE}" >&2
    echo "Record one with: make bench-update" >&2
    exit 1
  fi

  python3 - "${BASELINE_FILE}" "${median}" "${THRESHOLD_PCT}" <<'PY'
import json, sys
path, median_s, thr_s = sys.argv[1], sys.argv[2], sys.argv[3]
median = float(median_s)
thr = float(thr_s)
with open(path, encoding="utf-8") as f:
    base = json.load(f)
ref = float(base["median_seconds"])
limit = ref * (1.0 + thr / 100.0)
delta_pct = 100.0 * (median - ref) / ref if ref else float("inf")
print(f"Baseline median: {ref:.6f} s  ({path})")
print(f"Current  median: {median:.6f} s")
print(f"Delta:           {delta_pct:+.2f}%  (fail if slower than +{thr:.1f}%)")
print(f"Fail limit:      {limit:.6f} s")
if median > limit:
    print(
        f"PERF REGRESSION: median {median:.6f}s exceeds "
        f"{ref:.6f}s + {thr:.1f}% (= {limit:.6f}s)",
        file=sys.stderr,
    )
    sys.exit(1)
print("Perf compare OK (within threshold).")
PY
}

main() {
  if [[ ! -d "${DECK_DIR}" || ! -f "${DECK_DIR}/param.key" ]]; then
    echo "ERROR: missing bench deck at ${DECK_DIR}" >&2
    exit 1
  fi

  build_lassi
  local seed
  seed="$(require_pinned_seed "${DECK_DIR}/param.key")"

  echo "Bench deck '${DECK_NAME}'  seed=${seed}  runs=${RUNS}  mode=${MODE}"
  echo "Binary: ${LASSI_BIN}"
  mkdir -p "${WORK_ROOT}"

  local -a times=()
  local i elapsed
  for ((i = 1; i <= RUNS; i++)); do
    echo "--- run ${i}/${RUNS} ---"
    elapsed="$(run_once "${WORK_ROOT}/run_${i}" "${seed}")"
    times+=("${elapsed}")
    echo "wall_seconds=${elapsed}"
  done

  local median
  median="$(median_of "${times[@]}")"
  echo "----------------------------------------"
  echo "run_seconds: ${times[*]}"
  echo "median_seconds: ${median}"
  echo "----------------------------------------"

  case "${MODE}" in
    update-baseline)
      write_baseline "${median}" "${times[@]}"
      ;;
    compare)
      compare_baseline "${median}"
      ;;
    run)
      echo "Tip: make bench-compare  # fail if >${THRESHOLD_PCT}% slower than baseline"
      echo "Tip: make bench-update   # refresh tests/bench/baseline.json"
      ;;
  esac
}

main
