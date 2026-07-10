#!/usr/bin/env bash
# SELEXpipe motif-processing driver
#
# Runs totalautomodelIC on SEQ files produced by preprocessing/trim_to_kmer.py.
# Filenames are generic (no experiment-specific IDs).
#
# Original recipe (genericised):
#   ./totalautomodelIC -40N BACKGROUND_101.seq SIGNAL_101.seq \
#       1 6 6 0.25 - 20 50 0.1 SIGNALv1
#
# Argument order (from totalautomodelIC header):
#   -[seq length] background_seq signal_seq multinomial
#   min_seed_length max_seed_length length_cutoff iupac_cutoff
#   local_max_cutoff_count logo_number_cutoff logo_similarity_cutoff
#   output_prefix
#
# Requires the spacek40 binary alongside totalautomodelIC (the script
# invokes ./spacek40). Place spacek40 in --bin-dir (default: this motif/).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BACKGROUND=""
SIGNAL=""
OUTPUT_PREFIX=""
BIN_DIR="${SCRIPT_DIR}"
WORK_DIR="."
SEQ_LENGTH=40
MULTINOMIAL=1
MIN_SEED=6
MAX_SEED=6
LENGTH_CUTOFF=0.25
IUPAC_CUTOFF="-"
LOCAL_MAX_CUTOFF=20
LOGO_NUMBER=50
LOGO_SIMILARITY=0.1
ALL_GAPS=0

usage() {
  cat <<'EOF'
Usage: run_motif.sh --background BG.seq --signal SIGNAL.seq [options]

Required:
  --background PATH          Background SEQ file (e.g. BACKGROUND_101.seq)
  --signal PATH              Signal / enriched SEQ file (e.g. SIGNAL_101.seq)

Optional:
  --output-prefix NAME       Output name prefix (default: <signal_basename>v1)
  --bin-dir DIR              Directory containing totalautomodelIC and spacek40
                             (default: motif/ next to this script)
  --work-dir DIR             Working directory for motif outputs (default: .)
  --seq-length N             SEQ line length flag as -NN (default: 40 -> -40N)
  --multinomial N            Multinomial setting (default: 1)
  --min-seed N               Minimum seed length (default: 6)
  --max-seed N               Maximum seed length (default: 6)
  --length-cutoff X          Length cutoff (default: 0.25)
  --iupac-cutoff X           IUPAC cutoff, or '-' to disable (default: -)
  --local-max-cutoff N       Local-max instance cutoff (default: 20)
  --logo-number N            Max logos to consider (default: 50)
  --logo-similarity X        Logo similarity cutoff (default: 0.1)
  --all-gaps                 Pass '-allgaps' with the length flag
  -h, --help                 Show this help

Generic equivalents of the original commands:

  ./run_motif.sh \
      --background BACKGROUND_101.seq \
      --signal SIGNAL_r4_101.seq \
      --output-prefix SIGNAL_r4v1

  ./run_motif.sh \
      --background BACKGROUND_101.seq \
      --signal SIGNAL_r3_101.seq \
      --output-prefix SIGNAL_r3v1

Naming map (original -> generic):
  HTS147xlSOX3DBDc1_101.seq  ->  BACKGROUND_101.seq
  HTS147xlSOX3DBDc4_101.seq  ->  SIGNAL_r4_101.seq
  HTS147xlSOX3DBDc3_101.seq  ->  SIGNAL_r3_101.seq
  HTS147xlSOX3DBDc4v1        ->  SIGNAL_r4v1
  HTS147xlSOX3DBDc3v1        ->  SIGNAL_r3v1
EOF
}

log() {
  printf '[motif] %s\n' "$*"
}

die() {
  printf '[motif] error: %s\n' "$*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --background) BACKGROUND="${2:-}"; shift 2 ;;
    --signal) SIGNAL="${2:-}"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="${2:-}"; shift 2 ;;
    --bin-dir) BIN_DIR="${2:-}"; shift 2 ;;
    --work-dir) WORK_DIR="${2:-}"; shift 2 ;;
    --seq-length) SEQ_LENGTH="${2:-}"; shift 2 ;;
    --multinomial) MULTINOMIAL="${2:-}"; shift 2 ;;
    --min-seed) MIN_SEED="${2:-}"; shift 2 ;;
    --max-seed) MAX_SEED="${2:-}"; shift 2 ;;
    --length-cutoff) LENGTH_CUTOFF="${2:-}"; shift 2 ;;
    --iupac-cutoff) IUPAC_CUTOFF="${2:-}"; shift 2 ;;
    --local-max-cutoff) LOCAL_MAX_CUTOFF="${2:-}"; shift 2 ;;
    --logo-number) LOGO_NUMBER="${2:-}"; shift 2 ;;
    --logo-similarity) LOGO_SIMILARITY="${2:-}"; shift 2 ;;
    --all-gaps) ALL_GAPS=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "$BACKGROUND" ]] || { usage; die "--background is required"; }
[[ -n "$SIGNAL" ]] || { usage; die "--signal is required"; }
[[ -f "$BACKGROUND" ]] || die "background SEQ not found: $BACKGROUND"
[[ -f "$SIGNAL" ]] || die "signal SEQ not found: $SIGNAL"
[[ -d "$BIN_DIR" ]] || die "bin directory not found: $BIN_DIR"

BACKGROUND="$(cd "$(dirname "$BACKGROUND")" && pwd)/$(basename "$BACKGROUND")"
SIGNAL="$(cd "$(dirname "$SIGNAL")" && pwd)/$(basename "$SIGNAL")"
BIN_DIR="$(cd "$BIN_DIR" && pwd)"

TOTALAUTOMODELIC="${BIN_DIR}/totalautomodelIC"
SPACEK40="${BIN_DIR}/spacek40"

[[ -x "$TOTALAUTOMODELIC" || -f "$TOTALAUTOMODELIC" ]] || die "totalautomodelIC not found in $BIN_DIR"
[[ -f "$SPACEK40" ]] || die "spacek40 binary not found in $BIN_DIR (required; totalautomodelIC calls ./spacek40)"
chmod +x "$TOTALAUTOMODELIC" "$SPACEK40" 2>/dev/null || true

if [[ -z "$OUTPUT_PREFIX" ]]; then
  signal_base="$(basename "$SIGNAL")"
  signal_base="${signal_base%.seq}"
  OUTPUT_PREFIX="${signal_base}v1"
fi

mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"

LENGTH_FLAG="-${SEQ_LENGTH}N"
if [[ "$ALL_GAPS" -eq 1 ]]; then
  LENGTH_FLAG="-${SEQ_LENGTH}N -allgaps"
fi

log "background : $BACKGROUND"
log "signal     : $SIGNAL"
log "prefix     : $OUTPUT_PREFIX"
log "work-dir   : $WORK_DIR"
log "bin-dir    : $BIN_DIR"
log "command    : ./totalautomodelIC '${LENGTH_FLAG}' <bg> <signal> ${MULTINOMIAL} ${MIN_SEED} ${MAX_SEED} ${LENGTH_CUTOFF} ${IUPAC_CUTOFF} ${LOCAL_MAX_CUTOFF} ${LOGO_NUMBER} ${LOGO_SIMILARITY} ${OUTPUT_PREFIX}"

# totalautomodelIC invokes ./spacek40 and writes outputs relative to CWD.
# Run from WORK_DIR with binaries available via PATH-local copies/symlinks.
cd "$WORK_DIR"
ln -sfn "$TOTALAUTOMODELIC" ./totalautomodelIC
ln -sfn "$SPACEK40" ./spacek40

./totalautomodelIC \
  "$LENGTH_FLAG" \
  "$BACKGROUND" \
  "$SIGNAL" \
  "$MULTINOMIAL" \
  "$MIN_SEED" \
  "$MAX_SEED" \
  "$LENGTH_CUTOFF" \
  "$IUPAC_CUTOFF" \
  "$LOCAL_MAX_CUTOFF" \
  "$LOGO_NUMBER" \
  "$LOGO_SIMILARITY" \
  "$OUTPUT_PREFIX"

log "done"
log "  accepted logos : ${WORK_DIR}/${OUTPUT_PREFIX}_logos.svg"
log "  refined seeds  : ${WORK_DIR}/${OUTPUT_PREFIX}_refined_seeds.tmp.txt"
