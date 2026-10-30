#!/usr/bin/bash
set -euo pipefail

INDEX_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.low_coverage.GRCh38DH.alignment.index"

# Parse optional --sample flag
sample_name=""
count=10  # default value
if [[ "${1:-}" == "--sample" ]]; then
  if [[ -z "${2:-}" ]]; then
    echo "Error: --sample requires a sample name argument" >&2
    exit 1
  fi
  sample_name="${2}"
  shift 2
  # When --sample is used, remaining args are [outdir] [mode]
  outdir="${1:-./downloads}"
  mode="${2:-download}"
else
  # Without --sample, remaining args are [count] [outdir] [mode]
  count="${1:-10}"
  outdir="${2:-./downloads}"
  mode="${3:-download}"
fi

if [[ "${mode}" != "download" && "${mode}" != "list" ]]; then
  echo "Usage: $0 [--sample <name>] [count=10] [outdir=./downloads] [mode=download|list]" >&2
  exit 1
fi

if [[ -z "${sample_name}" ]]; then
  if ! printf '%d\n' "${count}" >/dev/null 2>&1 || [[ "${count}" -lt 1 ]]; then
    echo "count must be a positive integer" >&2
    exit 1
  fi
fi

mkdir -p "${outdir}"
manifest="${outdir}/samples.tsv"

tmp_list="$(mktemp)"
trap 'rm -f "${tmp_list}"' EXIT

# Handle specific sample lookup vs. random sample selection
if [[ -n "${sample_name}" ]]; then
  # Look up specific sample in the index - search for sample name in any column, extract CRAM (col 1) and CRAI (col 3)
  curl -fsSL "${INDEX_URL}" | SAMPLE_NAME="${sample_name}" python3 -c "
import os
import sys
sample_name = os.environ['SAMPLE_NAME']
for line in sys.stdin:
  line = line.rstrip('\n')
  if line.startswith('#'):
    continue
  fields = [f.strip() for f in line.split('\t')]
  if len(fields) >= 3 and sample_name in line:
    print(f\"{fields[0]}\t{fields[2]}\")
" > "${tmp_list}"
  
  if [[ ! -s "${tmp_list}" ]]; then
    echo "Error: Sample '${sample_name}' not found in 1000 Genomes index" >&2
    exit 1
  fi
else
  # Select random samples
  curl -fsSL "${INDEX_URL}" | COUNT="${count}" python3 -c "
import os
import sys
import random
count = int(os.environ['COUNT'])
entries = []
for line in sys.stdin:
  line = line.rstrip('\n')
  if line.startswith('#'):
    continue
  fields = line.split('\t')
  if len(fields) >= 3:
    entries.append(f\"{fields[0]}\t{fields[2]}\")
random.shuffle(entries)
for entry in entries[:count]:
  print(entry)
" > "${tmp_list}"
fi

echo -e "sample\tcram_url\tcrai_url" > "${manifest}"

while IFS=$'\t' read -r cram_url_raw crai_url_raw; do
  cram_url="$(printf '%s' "${cram_url_raw}" | sed -E 's#^ftp:/+#https://#')"
  crai_url="$(printf '%s' "${crai_url_raw}" | sed -E 's#^ftp:/+#https://#')"

  cram_file="${cram_url##*/}"
  crai_file="${crai_url##*/}"
  sample="${cram_file%%.*}"

  echo -e "${sample}\t${cram_url}\t${crai_url}" >> "${manifest}"

  if [[ "${mode}" == "download" ]]; then
    echo "Downloading ${sample}..."
    curl -fL --retry 3 --retry-delay 2 -o "${outdir}/${cram_file}" "${cram_url}"
    curl -fL --retry 3 --retry-delay 2 -o "${outdir}/${crai_file}" "${crai_url}"
  else
    echo "${sample}"
  fi
done < "${tmp_list}"

echo "Wrote manifest: ${manifest}"
if [[ "${mode}" == "download" ]]; then
  if [[ -n "${sample_name}" ]]; then
    echo "Done. Downloaded CRAM + CRAI pair for ${sample_name}"
  else
    echo "Done. Downloaded ${count} CRAM + CRAI pairs into ${outdir}"
  fi
else
  if [[ -n "${sample_name}" ]]; then
    echo "Found sample ${sample_name}"
  else
    echo "List-only mode complete for ${count} samples."
  fi
fi
