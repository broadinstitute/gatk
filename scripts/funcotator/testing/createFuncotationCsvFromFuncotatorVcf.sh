#!/usr/bin/env bash

################################################################################

#Setup variables for the script:
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=1

################################################################################

f=$1

head -n 100 $f  | grep 'Funcotation fields are' | sed -e 's3.*: 33g' -e 's#">##g' -e 's#|#,#g'
grep -v '^#' $1 | perl -p -e 's#.*FUNCOTATION=##g' | perl -p -e 's#;.*##g' | perl -pe 's#,#\n#g' | perl -pe 's#\[##g' | perl -pe 's#\]##g' | perl -pe 's#\|#,#g'

