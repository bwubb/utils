#!/usr/bin/env bash
# Build per-chr filtered VEP CSV: only rows whose ID is mnp_id, id1, or id2 from the plan.
# Usage: make_mnp_vep_filtered.sh [plan_dir [csv_dir [output_dir]]]
# Defaults: data/mnp for all.

set -euo pipefail

PLAN_DIR="${1:-data/mnp}"
CSV_DIR="${2:-data/mnp}"
OUT_DIR="${3:-data/mnp}"
mkdir -p "$OUT_DIR"

for plan in "$PLAN_DIR"/chr*.mnp_gt.plan.txt; do
  [ -f "$plan" ] || continue
  ch=$(basename "$plan" .mnp_gt.plan.txt | sed 's/^chr//')
  vep="$CSV_DIR/chr${ch}.mnp_gt.vep.csv"
  out="$OUT_DIR/chr${ch}.mnp_gt.vep.filtered.csv"

  if [ ! -f "$vep" ]; then
    echo "Skip chr${ch}: no VEP CSV" >&2
    continue
  fi

  # Unique IDs from plan columns 1 (mnp_id), 6 (id1), 7 (id2); skip header
  ids_file=$(mktemp)
  awk -F'\t' 'NR>1 && NF>=7 { print $1; print $6; print $7 }' "$plan" | sort -u > "$ids_file"

  # Header + matching rows only
  head -1 "$vep" > "$out"
  grep -Ff "$ids_file" "$vep" >> "$out"

  rm -f "$ids_file"
  echo "Wrote $out" >&2
done
