#!/usr/bin/env bash

declare -a -r SKIP=(            # Skip these directories.
    docs
    src
    gradle
    project
    build
    NA12878
    resources_for_CI
    .github
    .gradle
    .idea
)
declare -r CROMWELL=https://github.com/broadinstitute/cromwell
declare -r WOMTOOL=$CROMWELL/releases/download/78/womtool-78.jar

declare -r BW='\033[0m'         # Black or white, the default.
blue () {
    local -r head=$1; shift
    >&2 echo -e '\033[0;36m'$head$BW "$@"
}
red () {
    >&2 echo -e '\033[0;31m'"$@"$BW
}
yellow () {
    >&2 echo -e $(tput -Txterm smul)'\033[1;33m'$1$BW$(tput -Txterm rmul)
}

check_input () {
    local -r wdl=$1 json=$2; shift 2
    local -a -r check=("$@")
    local ok=true
    blue $json: Validating against $wdl
    if ! "${check[@]}" "$wdl" -i "$json" >/dev/null; then
        red $json: Validation failed
        ok=false
    fi
    $ok
}

check_options () {
    local -r json="$1"
    local ok=true
    blue $json: Checking well-formedness
    if ! jq -er . "$json" >/dev/null; then
        red $json: Not well-formed
        ok=false
    fi
    $ok
}

check_wdl () {
    local -r wdl=$1; shift
    local -a -r check=("$@")
    local -r dir=$(dirname ${wdl}) base=${wdl##*/}
    local -r prefix=${base%.*}
    local -r options=($(find "$dir" -name "$prefix"'*.options.json' -type f))
    local ok=true option json
    blue $wdl: Validating WDL
    if ! "${check[@]}" "$wdl" >/dev/null; then
        red $wdl: Validation failed
        ok=false
    fi
    for option in "${options[@]}"; do
        check_options "$option" || ok=false
    done
    if [ -d "$dir/input_files" ]; then
        for json in $(find "$dir/input_files" -name '*.json' -type f); do
            check_input "$wdl" "$json" "${check[@]}" || ok=false
        done
    fi
    $ok
}

check_dir () {
    local -r dir=$1; shift
    local -a -r check=("$@")
    local ok=true wdl
    yellow $dir
    for wdl in $(find "$dir" -name '*.wdl' -type f); do
        check_wdl "$wdl" "${check[@]}" || ok=false
    done
    $ok
}

dirs_to_check () {
    local -a find=(find . -type d -maxdepth 1 -not -path . -not -path ./.git)
    local d; for d in "${SKIP[@]}"; do find+=(-not -path "./$d"); done
    "${find[@]}" -print
}

main () {
    local -r here=$(cd $(dirname "$0") && pwd) tool=$(mktemp)
    local -r check=(java -jar $tool validate) get=(wget -q -O $tool $WOMTOOL)
    local ok=true d
    cd $(dirname "$here")
    trap "rm -f $tool; cd ~-" EXIT
    >&2 echo ${get[@]} && "${get[@]}"
    for d in $(dirs_to_check); do
        check_dir "${d#./}" "${check[@]}" || ok=false
    done
    $ok && return 0
    red Some WDL files are bad.
    false
}

main
