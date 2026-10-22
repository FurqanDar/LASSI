#!/usr/bin/env bash
# Characterization / regression harness for LASSI.
#
# Builds via src/compile.sh (unless LASSI_BIN is set), runs pinned-seed smoke
# decks, and bit-compares selected outputs against checked-in baselines.
# Exits non-zero (CI-style) on any mismatch or missing artifact.
#
# One command (preferred):
#   make characterization
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DECKS_DIR="${SCRIPT_DIR}/decks"
BASELINES_DIR="${SCRIPT_DIR}/baselines"
WORK_ROOT="${CHAR_WORK_DIR:-${SCRIPT_DIR}/.work}"

# Pin locale / timezone so printf float text cannot drift with ambient env.
export LC_ALL=C
export LANG=C
export TZ=UTC
unset LC_NUMERIC || true

UPDATE_BASELINES=0
if [[ "${1:-}" == "--update-baselines" ]]; then
  UPDATE_BASELINES=1
elif [[ $# -gt 0 ]]; then
  echo "Usage: $0 [--update-baselines]" >&2
  exit 2
fi

DECKS=(smoke_single smoke_multi)

# Artifacts compared bit-for-bit against baselines (no wall-clock logs).
ARTIFACT_SUFFIXES=(
  trj.lammpstrj
  topo.lammpstrj
  0_restart.lammpstrj
  EQ_restart.lammpstrj
  0_energy.dat
  energy.dat
  0_mcmove.dat
  mcmove.dat
)

build_lassi() {
  if [[ -n "${LASSI_BIN:-}" ]]; then
    if [[ ! -x "${LASSI_BIN}" ]]; then
      echo "ERROR: LASSI_BIN=${LASSI_BIN} is not executable" >&2
      exit 1
    fi
    echo "Using LASSI_BIN=${LASSI_BIN}"
    return
  fi

  echo "Building lassi via src/compile.sh (pinned flags in that script) ..."
  (cd "${REPO_ROOT}/src" && bash compile.sh)
  LASSI_BIN="${REPO_ROOT}/src/lassi"
}

# Extract the RANDOM_SEED value from a keyfile; fail if missing or zero.
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

  if [[ -z "${seed}" ]]; then
    echo "ERROR: ${keyfile}: RANDOM_SEED value missing" >&2
    exit 1
  fi
  if [[ "${seed}" == "0" ]]; then
    echo "ERROR: ${keyfile}: RANDOM_SEED 0 uses time(NULL); smoke decks must pin a non-zero seed" >&2
    exit 1
  fi
  printf '%s' "${seed}"
}

run_deck() {
  local deck="$1"
  local deck_dir="${DECKS_DIR}/${deck}"
  local work_dir="${WORK_ROOT}/${deck}"
  local baseline_dir="${BASELINES_DIR}/${deck}"
  local seed
  local rc

  if [[ ! -d "${deck_dir}" ]]; then
    echo "ERROR: missing deck directory ${deck_dir}" >&2
    exit 1
  fi
  if [[ ! -f "${deck_dir}/param.key" ]]; then
    echo "ERROR: missing ${deck_dir}/param.key" >&2
    exit 1
  fi

  seed="$(require_pinned_seed "${deck_dir}/param.key")"

  rm -rf "${work_dir}"
  mkdir -p "${work_dir}"
  cp "${deck_dir}/"* "${work_dir}/"

  echo "Running deck '${deck}' with RANDOM_SEED=${seed} ..."
  set +e
  (
    cd "${work_dir}"
    "${LASSI_BIN}" -k param.key >run.log 2>&1
  )
  rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    echo "ERROR: lassi exited ${rc} for deck '${deck}'" >&2
    tail -n 40 "${work_dir}/run.log" >&2 || true
    exit 1
  fi

  # Confirm the binary used the pinned seed (guards against silent wall-clock seeding).
  if ! grep -Eq "RNG Seed[[:space:]]*=[[:space:]]*${seed}[[:space:]]*$" "${work_dir}/run.log"; then
    echo "ERROR: run log for '${deck}' does not report RNG Seed = ${seed}" >&2
    grep -n "RNG Seed" "${work_dir}/run.log" >&2 || true
    exit 1
  fi

  if [[ "${UPDATE_BASELINES}" -eq 1 ]]; then
    mkdir -p "${baseline_dir}"
    for suffix in "${ARTIFACT_SUFFIXES[@]}"; do
      local artifact="${deck}_${suffix}"
      if [[ -f "${work_dir}/${artifact}" ]]; then
        cp -a "${work_dir}/${artifact}" "${baseline_dir}/${artifact}"
      else
        echo "ERROR: expected artifact missing while updating baselines: ${artifact}" >&2
        exit 1
      fi
    done
    (
      cd "${baseline_dir}"
      sha256sum "${deck}_"* >SHA256SUMS
    )
    echo "Updated baselines for '${deck}' in ${baseline_dir}"
    return
  fi

  if [[ ! -d "${baseline_dir}" ]]; then
    echo "ERROR: missing baselines for ${deck} at ${baseline_dir}" >&2
    echo "Re-run with --update-baselines after intentional reference changes." >&2
    exit 1
  fi

  local failed=0
  for suffix in "${ARTIFACT_SUFFIXES[@]}"; do
    local artifact="${deck}_${suffix}"
    local got="${work_dir}/${artifact}"
    local want="${baseline_dir}/${artifact}"

    if [[ ! -f "${want}" ]]; then
      echo "ERROR: baseline missing: ${want}" >&2
      failed=1
      continue
    fi
    if [[ ! -f "${got}" ]]; then
      echo "ERROR: expected output missing after run: ${got}" >&2
      failed=1
      continue
    fi
    if ! cmp -s "${got}" "${want}"; then
      echo "MISMATCH: ${artifact}" >&2
      echo "  got:  ${got}" >&2
      echo "  want: ${want}" >&2
      failed=1
    else
      echo "OK: ${artifact}"
    fi
  done

  if [[ "${failed}" -ne 0 ]]; then
    echo "Characterization FAILED for deck '${deck}' (CI-style non-zero exit)." >&2
    exit 1
  fi
  echo "Deck '${deck}' matches baselines."
}

main() {
  build_lassi
  mkdir -p "${WORK_ROOT}"

  for deck in "${DECKS[@]}"; do
    run_deck "${deck}"
  done

  if [[ "${UPDATE_BASELINES}" -eq 1 ]]; then
    echo "All baselines updated."
  else
    echo "All characterization decks passed (bit-identical)."
  fi
}

main
