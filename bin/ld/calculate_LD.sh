#!/bin/bash
# calculate_LD.sh
# WARNING: Computationally intensive on large LD sets.

set -uo pipefail

# ---- Load config from runroot (preferred) or repo root fallback ----
if [[ -f "./config_shell.sh" ]]; then
  # shellcheck disable=SC1091
  source "./config_shell.sh"
else
  SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  PROJECT_ROOT=$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")
  CONFIG_FILE="$PROJECT_ROOT/config_shell.sh"
  if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
  else
    echo "Config file $CONFIG_FILE not found" >&2
    exit 1
  fi
fi

: "${SINGLE_TRAIT_LOCI:?SINGLE_TRAIT_LOCI not set}"

# Resolve PLINK1_PATH: prefer configured value if executable; otherwise fall back
# to system 'plink' on PATH. This allows using the container-installed plink.
if [[ -n "${PLINK1_PATH:-}" && -x "${PLINK1_PATH}" ]]; then
  : # PLINK1_PATH already configured and executable
elif command -v plink >/dev/null 2>&1; then
  PLINK1_PATH="$(command -v plink)"
  export PLINK1_PATH
  echo "Using PLINK from PATH: ${PLINK1_PATH}" >&2
else
  echo "ERROR: PLINK not found. Set PLINK1_PATH in config_shell.sh or ensure 'plink' is on PATH." >&2
  exit 1
fi
: "${SINGLE_TRAIT_LOCI:?SINGLE_TRAIT_LOCI not set}"
: "${SINGLE_TRAIT_LD:?SINGLE_TRAIT_LD not set}"
: "${DOWNSAMPLE:?DOWNSAMPLE (ld_reference_panel) not set}"

# Optional user template (kept as fallback; safe to leave empty for pure autodetect)
LD_BFILE_TEMPLATE="${LD_BFILE_TEMPLATE:-}"

mkdir -p "$SINGLE_TRAIT_LD"

csv_file="$SINGLE_TRAIT_LOCI"
mapfile -t lines < <(tail -n +2 "$csv_file")

# Concurrency and PLINK resource settings
MAX_PARALLEL="${PARALLEL_JOBS:-2}"
PLINK_THREADS="${PLINK_THREADS:-2}"
PLINK_MEM_MB="${PLINK_MEM_MB:-8000}"

missing_bfiles_log="$SINGLE_TRAIT_LD/missing_bfiles.txt"
missing_snps_log="$SINGLE_TRAIT_LD/missing_snps.txt"
: > "$missing_bfiles_log"
: > "$missing_snps_log"

# ----------------- Helpers -----------------

# Normalize chromosome token: strip leading "chr"/"CHR"
norm_chr() {
  local tok="$1"
  tok="${tok#chr}"; tok="${tok#CHR}"
  case "$tok" in
    23) echo "X" ;;
    24) echo "Y" ;;
    25|26|MT|Mt|mt|M|m) echo "MT" ;;
    *) echo "$tok" ;;
  esac
}

# Extract chromosome ID from filename (exactly one of: 1-22, X, Y, MT)
# Matches tokens like "...chr1_", "...chr10.", "...c1_", "...cX-"
extract_chr_from_filename() {
  local f="$1"
  local base; base="$(basename "$f")"
  local m
  # Try chrN
  if [[ "$base" =~ (^|[^A-Za-z0-9])chr([0-9]{1,2}|X|Y|MT)([^A-Za-z0-9]|$) ]]; then
    m="${BASH_REMATCH[2]}"
    echo "$(norm_chr "$m")"
    return 0
  fi
  # Try cN
  if [[ "$base" =~ (^|[^A-Za-z0-9])c([0-9]{1,2}|X|Y|MT)([^A-Za-z0-9]|$) ]]; then
    m="${BASH_REMATCH[2]}"
    echo "$(norm_chr "$m")"
    return 0
  fi
  return 1
}

# Resolve a PLINK --bfile prefix for a chromosome with exact-ID matching
resolve_bfile_for_chr() {
  local chr_raw="$1"
  local chr="$(norm_chr "$chr_raw")"
  local panel_dir="$DOWNSAMPLE"

  # 1) Template, if given and valid
  if [[ -n "$LD_BFILE_TEMPLATE" ]]; then
    local prefix="${LD_BFILE_TEMPLATE}"
    prefix="${prefix//\{panel\}/$panel_dir}"
    prefix="${prefix//\{chr\}/$chr}"
    if [[ -f "${prefix}.bed" ]]; then
      # Validate exact match (avoid chr1 picking chr10 etc.)
      if extract_chr_from_filename "${prefix}.bed" >/dev/null 2>&1; then
        local fchr; fchr="$(extract_chr_from_filename "${prefix}.bed")"
        if [[ "$fchr" == "$chr" ]]; then
          echo "$prefix"
          return 0
        fi
      fi
    fi
  fi

  # 2) Autodetect: scan all .bed once, filter to exact chr
  shopt -s nullglob
  local allbeds=("$panel_dir"/*.bed)
  shopt -u nullglob

  local exact_matches=()
  local b
  for b in "${allbeds[@]}"; do
    if extract_chr_from_filename "$b" >/dev/null 2>&1; then
      local fchr; fchr="$(extract_chr_from_filename "$b")"
      if [[ "$fchr" == "$chr" ]]; then
        exact_matches+=("$b")
      fi
    fi
  done

  if (( ${#exact_matches[@]} == 0 )); then
    return 1
  fi

  # Prefer names containing _EUR_, otherwise shortest basename, then lexicographic
  local pick=""
  for b in "${exact_matches[@]}"; do
    if [[ "$(basename "$b")" == *"_EUR_"* ]]; then
      pick="$b"; break
    fi
  done
  if [[ -z "$pick" ]]; then
    # pick shortest basename
    local best_len=99999
    for b in "${exact_matches[@]}"; do
      local len=${#b}
      if (( len < best_len )); then best_len=$len; pick="$b"; fi
    done
  fi

  echo "${pick%.bed}"
  return 0
}

# ----------------- Worker -----------------

run_one() {
  local line="$1"
  IFS=, read -r Phenotype Locus_n Locus_name SNP CHROM GENPOS ALLELE1 ALLELE0 A1FREQ MAF INFO BETA SE P GENPOS_hg19 rsid_1kg _rest <<< "$line"

  local snp_str="$SNP"
  local file_lead_snp_id; file_lead_snp_id=$(echo "$snp_str" | tr ':' '_')
  local lead_snp_id;      lead_snp_id=$(echo "$snp_str" | tr '_' ':')

  # Prefer chr from SNP token (before first ':'); fallback to CSV CHROM
  local snp_chr_token=""; [[ "$snp_str" == *:* ]] && snp_chr_token="${snp_str%%:*}"
  snp_chr_token="$(norm_chr "$snp_chr_token")"
  local csv_chr="$(norm_chr "${CHROM:-}")"

  local chr=""
  if [[ -n "$snp_chr_token" && "$snp_chr_token" != rs* ]]; then
    chr="$snp_chr_token"
  elif [[ -n "$csv_chr" ]]; then
    chr="$csv_chr"
  else
    echo "ERROR: Could not determine chromosome for SNP '$snp_str' (no SNP token or CHROM column)." >&2
    echo "$snp_str" >> "$missing_snps_log"
    return 1
  fi

  if [[ -n "$csv_chr" && -n "$snp_chr_token" && "$csv_chr" != "$snp_chr_token" ]]; then
    echo "WARN: CHROM mismatch for $snp_str — CSV=$csv_chr, SNP token=$snp_chr_token. Using $chr." >&2
  fi

  local ld_file="$SINGLE_TRAIT_LD/locus_${file_lead_snp_id}.ld"
  local lock_file="$SINGLE_TRAIT_LD/locus_${file_lead_snp_id}.lock"

  # Skip if already present and non-empty
  if [[ -s "$ld_file" ]]; then
    echo "LD exists: $ld_file (skip)"
    return 0
  fi

  # Non-blocking lock
  exec 200>"$lock_file" || true
  if ! flock -n 200; then
    echo "Another worker is generating $ld_file (skip)"
    return 0
  fi

  # Resolve exact-matching panel
  local bfile_prefix
  if ! bfile_prefix="$(resolve_bfile_for_chr "$chr")"; then
    echo "ERROR: No reference panel .bed for chr${chr} in ${DOWNSAMPLE} (and no usable template)." >&2
    echo "chr${chr} :: $snp_str" >> "$missing_bfiles_log"
    flock -u 200 || true
    return 1
  fi

  # Verify trio exists
  for ext in bed bim fam; do
    if [[ ! -f "${bfile_prefix}.${ext}" ]]; then
      echo "ERROR: Missing ${bfile_prefix}.${ext} for chr${chr}" >&2
      echo "chr${chr} :: $snp_str (missing ${ext})" >> "$missing_bfiles_log"
      flock -u 200 || true
      return 1
    fi
  done

  # Defensive: ensure chosen file really is the intended chr
  if extract_chr_from_filename "${bfile_prefix}.bed" >/dev/null 2>&1; then
    chosen_chr="$(extract_chr_from_filename "${bfile_prefix}.bed")"
    if [[ "$chosen_chr" != "$(norm_chr "$chr")" ]]; then
      echo "ERROR: Resolver picked ${bfile_prefix}.bed (chr ${chosen_chr}) for target chr ${chr} — refusing." >&2
      echo "chr${chr} :: picked ${chosen_chr} :: $snp_str" >> "$missing_bfiles_log"
      flock -u 200 || true
      return 1
    fi
  fi

  echo "PLINK: trait=$Phenotype SNP=$lead_snp_id chr=$chr bfile=$(basename "$bfile_prefix")" >&2

  # Run PLINK, capture stderr
  if ! "$PLINK1_PATH" \
    --bfile "$bfile_prefix" \
    --r2 \
    --ld-snp "$lead_snp_id" \
    --keep-allele-order \
    --ld-window-kb 4000 \
    --ld-window-r2 0.1 \
    --ld-window 99999 \
    --out "$SINGLE_TRAIT_LD/locus_${file_lead_snp_id}" \
    --threads "$PLINK_THREADS" \
    --memory "$PLINK_MEM_MB" 2> >(tee "$SINGLE_TRAIT_LD/locus_${file_lead_snp_id}.plink.err" >&2)
  then
    echo "ERROR: PLINK failed for $lead_snp_id (chr$chr)" >&2
    echo "$lead_snp_id" >> "$missing_snps_log"
    flock -u 200 || true
    return 1
  fi

  if grep -qi "No valid variants specified by --ld-snp" "$SINGLE_TRAIT_LD/locus_${file_lead_snp_id}.plink.err"; then
    echo "WARN: PLINK reported no valid variants for $lead_snp_id (chr$chr)" >&2
    echo "$lead_snp_id" >> "$missing_snps_log"
  fi

  flock -u 200 || true
  return 0
}

# ----------------- Scheduler -----------------

running=0
fail_count=0

for line in "${lines[@]}"; do
  run_one "$line" &
  ((running++))
  if (( running >= MAX_PARALLEL )); then
    if ! wait -n; then
      ((fail_count++))
    fi
    ((running--))
  fi
done

while (( running > 0 )); do
  if ! wait -n; then
    ((fail_count++))
  fi
  ((running--))
done

if (( fail_count > 0 )); then
  echo "WARNING: $fail_count PLINK job(s) failed; continuing."
fi

miss_b=$(wc -l < "$missing_bfiles_log" 2>/dev/null || echo 0)
miss_s=$(wc -l < "$missing_snps_log" 2>/dev/null || echo 0)
echo "Missing bfiles: $miss_b  |  Missing SNP hits: $miss_s" >&2

exit 0
