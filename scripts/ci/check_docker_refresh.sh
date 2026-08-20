#!/usr/bin/env bash

set -euo pipefail

BASE_SHA="${1:-}"
HEAD_SHA="${2:-}"

if [[ -z "${BASE_SHA}" || -z "${HEAD_SHA}" ]]; then
  echo "Usage: $0 <base_sha> <head_sha>" >&2
  exit 0
fi

DIFF_RANGE="${BASE_SHA}...${HEAD_SHA}"
GVS_UTILS_FILE="scripts/variantstore/wdl/GvsUtils.wdl"

needs_variants_refresh=false
needs_gatk_refresh=false

triggered_variants_files=()
triggered_gatk_files=()

# This parses the git-produced diff to see which files were changed.  It's how we detect whether a Docker refresh is needed.
while IFS= read -r changed_file; do
  [[ -z "${changed_file}" ]] && continue

  if [[ "${changed_file}" == scripts/variantstore/scripts/* && "${changed_file}" != scripts/variantstore/scripts/test/* ]]; then

  # Core Java/runtime changes that typically require a GATK image rebuild.
  if [[ "${changed_file}" == src/main/* ]] ||
     [[ "${changed_file}" == build.gradle ]] ||
     [[ "${changed_file}" == settings.gradle ]] ||
     [[ "${changed_file}" == gradle.properties ]] ||
     [[ "${changed_file}" == Dockerfile ]] ||
     [[ "${changed_file}" == scripts/docker/gatkbase/* ]]; then
    needs_gatk_refresh=true
    triggered_gatk_files+=("${changed_file}")
  fi
done < <(git diff --name-only "${DIFF_RANGE}")

# This is where we check whether the docker images have been updated
gvs_utils_diff="$(git diff --unified=0 "${DIFF_RANGE}" -- "${GVS_UTILS_FILE}" || true)"

updated_variants_tag=false
updated_gatk_lite_tag=false
updated_gatk_heavy_tag=false

if grep -qE '^[+-].*String variants_docker = ' <<< "${gvs_utils_diff}"; then
  updated_variants_tag=true
fi

if grep -qE '^[+-].*String gatk_docker = ' <<< "${gvs_utils_diff}"; then
  updated_gatk_lite_tag=true
fi

if grep -qE '^[+-].*String gatk_heavy_docker = ' <<< "${gvs_utils_diff}"; then
  updated_gatk_heavy_tag=true
fi

missing_items=()

# Basic check to see if we updated them if we needed to
if [[ "${needs_variants_refresh}" == "true" && "${updated_variants_tag}" != "true" ]]; then
  missing_items+=("variants_docker")
fi

if [[ "${needs_gatk_refresh}" == "true" && "${updated_gatk_lite_tag}" != "true" ]]; then
  missing_items+=("gatk_docker")
fi

if [[ "${needs_gatk_refresh}" == "true" && "${updated_gatk_heavy_tag}" != "true" ]]; then
  missing_items+=("gatk_heavy_docker")
fi

missing_count="${#missing_items[@]}"
status="ok"
if [[ "${missing_count}" -gt 0 ]]; then
  status="warn"
fi

summary="Docker refresh check: ${status}."

if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
  {
    echo "status=${status}"
    echo "missing_count=${missing_count}"
    echo "needs_variants_refresh=${needs_variants_refresh}"
    echo "needs_gatk_refresh=${needs_gatk_refresh}"
    echo "updated_variants_tag=${updated_variants_tag}"
    echo "updated_gatk_lite_tag=${updated_gatk_lite_tag}"
    echo "updated_gatk_heavy_tag=${updated_gatk_heavy_tag}"

    echo "missing_items<<EOF"
    if [[ "${missing_count}" -gt 0 ]]; then
      printf '%s\n' "${missing_items[@]}"
    fi
    echo "EOF"

    echo "triggered_variants_files<<EOF"
    if [[ "${#triggered_variants_files[@]}" -gt 0 ]]; then
      printf '%s\n' "${triggered_variants_files[@]}"
    fi
    echo "EOF"

    echo "triggered_gatk_files<<EOF"
    if [[ "${#triggered_gatk_files[@]}" -gt 0 ]]; then
      printf '%s\n' "${triggered_gatk_files[@]}"
    fi
    echo "EOF"

    echo "summary=${summary}"
  } >> "${GITHUB_OUTPUT}"
fi

echo "${summary}"
