#!/usr/bin/env bash
#
# VS-834 Generate the `git filter-repo` command for moving files and history from the ah_var_store branch of `gatk`
# to the `variantstore` repo. Find files that were/are:
#
# * Added on ah_var_store after it branched from master.
# * Not on master (exclude files cherry picked from master).
# * Were never deleted from master (exclude files cherry picked from master that were subsequently deleted).

# Debug only
# PS4='\D{+%F %T} \w $ '
# set -o errexit -o nounset -o pipefail -o xtrace

# Make the variants branch a variable to be able to experiment with the effect of moving files around on a branch
# of ah_var_store.
# variants_branch=ah_var_store
variants_branch=vs_834_disentangle

files_on_master() {
  git ls-tree -r --name-only master
}


ah_var_store_branch_point() {
  # Finds the commit at which ah_var_store branched off master.
  # https://stackoverflow.com/a/4991675/21269164
  # The more obvious `git merge-base --fork-point master ah_var_store` does not work, exits 1 with no output.
  diff -u <(git rev-list --first-parent $variants_branch) <(git rev-list --first-parent master) |
    sed -ne 's/^ //p' | head -1
}


files_added_on_ah_var_store() {
  # Look for files added to ah_var_store since the branch point. Note that these files were not necessarily *uniquely*
  # added to ah_var_store and might represent cherry picks from master (e.g. scripts and build files for the migration
  # from Travis to GitHub Actions, VQSR Lite work, etc.)
  git diff "$(ah_var_store_branch_point)" $variants_branch --name-status | grep -E '^A' | cut -f 2-
}


files_added_on_ah_var_store_not_on_master() {
  # `comm` excluding files only on master and files on master that were also added on ah_var_store
  comm -1 -3 <(files_on_master) <(files_added_on_ah_var_store)
}


files_deleted_from_master() {
  # This intentionally does not use `git diff` as is used in `files_added_on_ah_var_store` since that would only show
  # files deleted from the branch point to the head of master. There are numerous files here (mostly related to VQSR
  # Lite) where files added to master after the branch point were cherry picked onto ah_var_store and subsequently
  # deleted from master. This `git log` finds these while the `git diff` does not.
  #
  # https://waylonwalker.com/git-find-deleted-files/#git-log-diff-filter
  # That `sort` at the end is critical to avoid confusing `comm`
  git log master --diff-filter D --pretty="format:" --name-only | sed '/^$/d' | sort
}


files_to_move() {
  # `comm` excluding files deleted from master from our previous working set of files.
  comm -1 -3 <(files_deleted_from_master) <(files_added_on_ah_var_store_not_on_master)
}

paths_to_move() {
  # Call out the paths that we want to move from gatk to variantstore.
  # Partition into paths that are unique to Variants (containing '/gvs/' or '/variantstore/') and those that are not.
  # Use perl for the non-greedy `*?` quantifier since we're trying to grab as short of a 'variantstore' or 'gvs' path
  # prefix as possible.
  files_to_move | perl -ne 'if (/(.*?(\/gvs\/|\/variantstore\/)).*/) { print "$1\n" } else { print "$_\n" }' |
    sed '/^$/d' | sort -u
}

# Files that were legitimately created on ah_var_store but currently need to stay in master barring some significant
# refactoring:
KEEP_ON_MASTER=(
  # Added on ah_var_store but referenced by BigQueryUtils which was created on master before the branch point.
  src/main/java/org/broadinstitute/hellbender/utils/bigquery/BigQueryResultAndStatistics.java
  # Same
  src/main/java/org/broadinstitute/hellbender/utils/bigquery/StorageAPIAvroReaderAndBigQueryStatistics.java
)

# https://stackoverflow.com/a/17841619/21269164
function join_by { local IFS="$1"; shift; echo "$*"; }

KEEP_ON_MASTER_FILTER=$(join_by '|' "${KEEP_ON_MASTER[@]}")
paths_to_move | grep -v -E "$KEEP_ON_MASTER_FILTER" > /tmp/gvs_paths.txt

# Build the command to be executed in *a clone of* the `gatk` repo prior to pushing to `variantstore`.

# Disable the ShellCheck warning about word splitting since that's a feature here. File names are individually quoted
# within the expression.
# shellcheck disable=SC2046

# The 'Useless cat' makes this more readable left to right.
# shellcheck disable=SC2002
echo git filter-repo $(cat /tmp/gvs_paths.txt | sed 's/^/'\''/' | sed 's/$/'\''/' | sed -E 's/^/ --path /' | tr -d '\n')
