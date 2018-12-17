#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script will create a liftover chain file that goes from b37 -> hg38
#
# EXAMPLE:
#   ./createLiftoverChainFileForB37ToHg38.sh
#
# AUTHOR: Jonn Smith
#
###############################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( python -c "import os;print os.path.realpath(\"${BASH_SOURCE[0]}\")" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=1
PREREQUISITES="curl gunzip"

# Determine if this shell is interactive:
ISCALLEDBYUSER=true
[[ "${BASH_SOURCE[0]}" != "${0}" ]] && ISCALLEDBYUSER=false

# Required for the aliased call to checkPipeStatus:
shopt -s expand_aliases

################################################################################

HG19_HG38_CHAIN_FILE_URL='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
outChainFile=b37ToHg38.over.chain

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME" 
  echo -e "Create a liftover chain file for b37->hg38"
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e ""
  echo -e "Pulls a chain file for hg19->hg38 and modifies it to use"
  echo -e "b37 contig names.  The chain file is read from:"
  echo -e ""
  echo -e "${HG19_HG38_CHAIN_FILE_URL}"
  echo -e ""
  echo -e "This is OK because b37 is roughly equivalent to hg19 aside"
  echo -e "from these contig name differences."
  echo -e ""
  echo -e "Arguments:"
  echo -e ""
  echo -e "--showDiffs      Print a table of the differences between the "
  echo -e "                 b37 and hg19 sequences (compared by contig"
  echo -e "                 md5sum)"
  echo -e "--help           Show this Help menu"
  echo -e ""
  # List prerequisites, if any:
  if [[ ${#PREREQUISITES} -ne 0 ]] ; then
    echo -e "Requires the following programs to be installed:"
    for prereq in ${PREREQUISITES} ; do 
      echo -e "  ${prereq}"
    done
    echo
  fi
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  UNKNOWN ARGUMENT"
  echo -e "  4  MISSING PREREQUISITE"
  echo -e ""
}

#Display a message to std error:
function error() 
{
  echo "$1" 1>&2 
}

# Make a temporary file that will be cleaned after this script finishes. 
function makeTemp()
{
  local f
  f=$( mktemp )
  TMPFILELIST="${TMPFILELIST} $f"
  echo $f
}
TMPFILELIST=''

# Clean all temporary files made in this script:
function cleanTempVars()
{
  rm -f ${TMPFILELIST}
}

# Function run in the exit trap.
function at_exit()
{
  cleanTempVars
}

# Checks the bash built-in PIPESTATUS variable for any failures
# If given strings will output the string corresponding to the failure
# position in PIPESTATUS of any failures in the chain of commands 
# that occurred.
# This function should be called as `checkPipeStatus` as per the 
# alias below it.
function _checkPipeStatus() {

  local hadFailure=false

  for (( i = 0 ; i < ${#RC[@]} ; i++ )) ; do 
    st=${RC[i]}
    echo "st = ${st}"
    if [ $st -ne 0 ] ; then
      # If we were passed a description string for this error in the pipe, then we print it:
      let argIndex=$i+1
      description=${!argIndex}
      [[ ${#description} -ne 0 ]] && error "$description"
      hadFailure=true  
    fi
  done

  if $hadFailure ; then
    return 10
  fi
  return 0
}
alias checkPipeStatus='RC=( "${PIPESTATUS[@]}" );_checkPipeStatus'

function checkPrerequisites {
  local missing=""
  local foundMissing=false
  for c in ${@} ; do
    which $c &> /dev/null
    r=$?
    [[ $r -ne 0 ]] && foundMissing=true && missing="${missing}$c "
  done

  if $foundMissing ; then
    error "Error: the following commands could not be found:"
    error "  $missing"
    error "Please install them and try again."
    if [[ $(uname) == 'Darwin' ]] ; then
      error "NOTE: You likely must use homebrew to install them."
    fi
    return 1
  fi
  return 0
}

################################################################################

function printSequenceDictionaryComparison() {
	echo 'MD5	B37SN (Homo_sapiens_assembly19.dict)	HG19SN (ucsc.hg19.dict)'
	echo '06cbf126247d89664a4faebad130fe9c	GL000202.1	chr11_gl000202_random'
	echo '0996b4475f353ca98bacb756ac479140	GL000244.1	chrUn_gl000244'
	echo '118a25ca210cfbcdfb6c2ebb249f9680	GL000235.1	chrUn_gl000235'
	echo '131b1efc3270cc838686b54e7c34b17b	GL000238.1	chrUn_gl000238'
	echo '1c1b2cd1fccbc0a99b6a447fa24d1504	GL000226.1	chrUn_gl000226'
	echo '1d708b54644c26c7e01c2dad5426d38c	GL000218.1	chrUn_gl000218'
	echo '1d78abec37c15fe29a275eb08d5af236	GL000249.1	chrUn_gl000249'
	echo '2f8694fc47576bc81b5fe9e7de0ba49e	GL000242.1	chrUn_gl000242'
	echo '3238fb74ea87ae857f9c7508d315babb	GL000221.1	chrUn_gl000221'
	echo '325ba9e808f669dfeee210fdd7b470ac	GL000192.1	chr1_gl000192_random'
	echo '399dfa03bf32022ab52a846f7ca35b30	GL000223.1	chrUn_gl000223'
	echo '3e06b6741061ad93a8587531307057d8	GL000232.1	chrUn_gl000232'
	echo '43f69e423533e948bfae5ce1d45bd3f1	GL000206.1	chr17_gl000206_random'
	echo '445a86173da9f237d7bcf41c6cb8cc62	GL000240.1	chrUn_gl000240'
	echo '46c2032c37f2ed899eb41c0473319a69	GL000214.1	chrUn_gl000214'
	echo '563531689f3dbd691331fd6c5730a88b	GL000212.1	chrUn_gl000212'
	echo '569af3b73522fab4b40995ae4944e78e	GL000199.1	chr9_gl000199_random'
	echo '5a8e43bec9be36c7b49c84d585107776	GL000248.1	chrUn_gl000248'
	echo '5d9ec007868d517e73543b005ba48535	GL000195.1	chr7_gl000195_random'
	echo '5eb3b418480ae67a997957c909375a73	GL000215.1	chrUn_gl000215'
	echo '63945c3e6962f28ffd469719a747e73c	GL000225.1	chrUn_gl000225'
	echo '642a232d91c486ac339263820aef7fe0	GL000216.1	chrUn_gl000216'
	echo '6ac8f815bf8e845bb3031b73f812c012	GL000194.1	chr4_gl000194_random'
	echo '6d243e18dea1945fb7f2517615b8f52e	GL000217.1	chrUn_gl000217'
	echo '6f5efdd36643a9b8c8ccad6f2f1edc7b	GL000197.1	chr8_gl000197_random'
	echo '6fe9abac455169f50470f5a6b01d0f59	GL000222.1	chrUn_gl000222'
	echo '75e4c8d17cd4addf3917d1703cacaf25	GL000200.1	chr9_gl000200_random'
	echo '7daaa45c66b288847b9b32b964e623d3	GL000211.1	chrUn_gl000211'
	echo '7de00226bb7df1c57276ca6baabafd15	GL000247.1	chrUn_gl000247'
	echo '7e0e2e580297b7764e31dbc80c2540dd	X						chrX'
	echo '7fed60298a8d62ff808b74b6ce820001	GL000233.1	chrUn_gl000233'
	echo '851106a74238044126131ce2a8e5847c	GL000210.1	chr21_gl000210_random'
	echo '868e7784040da90d900d2d1b667a1383	GL000198.1	chr9_gl000198_random'
	echo '89bc61960f37d94abf0df2d481ada0ec	GL000245.1	chrUn_gl000245'
	echo '93f998536b61a56fd0ff47322a911d4b	GL000234.1	chrUn_gl000234'
	echo '96358c325fe0e70bee73436e8bb14dbd	GL000203.1	chr17_gl000203_random'
	echo '99795f15702caec4fa1c4e15f8a29c07	GL000239.1	chrUn_gl000239'
	echo '9d424fdcc98866650b58f004080a992a	GL000213.1	chrUn_gl000213'
	echo 'a4aead23f8053f2655e468bcc6ecdceb	GL000227.1	chrUn_gl000227'
	echo 'aa81be49bf3fe63a79bdc6a6f279abf6	GL000208.1	chr19_gl000208_random'
	echo 'b4eb71ee878d3706246b7c1dbef69299	GL000230.1	chrUn_gl000230'
	echo 'ba8882ce3a1efa2080e5d29b956568a4	GL000231.1	chrUn_gl000231'
	echo 'c5a17c97e2c1a0b6a9cc5a6b064b714f	GL000228.1	chrUn_gl000228'
	echo 'cc34279a7e353136741c9fce79bc4396	GL000243.1	chrUn_gl000243'
	echo 'd0f40ec87de311d8e715b52e4c7062e1	GL000229.1	chrUn_gl000229'
	echo 'd22441398d99caf673e9afb9a1908ec5	GL000205.1	chr17_gl000205_random'
	echo 'd5b2fc04f6b41b212a4198a07f450e20	GL000224.1	chrUn_gl000224'
	echo 'd75b436f50a8214ee9c2a51d30b2c2cc	GL000191.1	chr1_gl000191_random'
	echo 'd92206d1bb4c3b4019c43c0875c06dc0	GL000196.1	chr8_gl000196_random'
	echo 'dbb6e8ece0b5de29da56601613007c2a	GL000193.1	chr4_gl000193_random'
	echo 'dfb7e7ec60ffdcb85cb359ea28454ee9	GL000201.1	chr9_gl000201_random'
	echo 'e0c82e7751df73f4f6d0ed30cdc853c0	GL000237.1	chrUn_gl000237'
	echo 'e4afcd31912af9d9c2546acf1cb23af2	GL000246.1	chrUn_gl000246'
	echo 'ef4258cdc5a45c206cea8fc3e1d858cf	GL000241.1	chrUn_gl000241'
	echo 'efc49c871536fa8d79cb0a06fa739722	GL000204.1	chr17_gl000204_random'
	echo 'f3814841f1939d3ca19072d9e89f3fd7	GL000207.1	chr18_gl000207_random'
	echo 'f40598e2a5a6b26e84a3775e0d1e2c81	GL000209.1	chr19_gl000209_random'
	echo 'f977edd13bac459cb2ed4a5457dba1b3	GL000219.1	chrUn_gl000219'
	echo 'fc35de963c57bf7648429e6454f1c9db	GL000220.1	chrUn_gl000220'
	echo 'fdcd739913efa1fdc64b6c0cd7016779	GL000236.1	chrUn_gl000236'
	echo '1b22b98cdeb4a9304cb5d48026a85128	1						chr1'
	echo 'a0d9851da00400dec1098a9255ac712e	2						chr2'
	echo '23dccd106897542ad87d2765d28a19a1	4						chr4'
	echo '0740173db9ffd264d728f32784845cd7	5						chr5'
	echo '1d3a93a248d92a729ee764823acbbc6b	6						chr6'
	echo '618366e953d6aaad97dbe4777c29375e	7						chr7'
	echo '96f514a9929e410c6651697bded59aec	8						chr8'
	echo '3e273117f15e0a400f01055d9f393768	9						chr9'
	echo '988c28e000e84c26d552359af1ea2e1d	10					chr10'
	echo '98c59049a2df285c76ffb1c6db8f8b96	11					chr11'
	echo '51851ac0e1a115847ad36449b0015864	12					chr12'
	echo '283f8d7892baa81b510a015719ca7b0b	13					chr13'
	echo '98f3cae32b2a2e9524bc19813927542e	14					chr14'
	echo 'e5645a794a8238215b2cd77acb95a078	15					chr15'
	echo 'fc9b1a7b42b97a864f56b348b06095e6	16					chr16'
	echo '351f64d4f4f9ddd45b35336ad97aa6de	17					chr17'
	echo 'b15d4b2d29dde9d3e4f93d1d0f2cbc9c	18					chr18'
	echo '1aacd71f30db8e561810913e0b72636d	19					chr19'
	echo '0dec9660ec1efaaf33281c0d5ea2560f	20					chr20'
	echo '2979a6085bfe28e3ad6f552f361ed74d	21					chr21'
	echo 'a718acaa6135fdca8357d5bfe94211dd	22					chr22'
	echo '1fa3474750af0948bdf97d5a0ee52e51	Y						----'
	echo '6743bd63b3ff2b5b8985d8933c53290a	NC_007605		----'
	echo 'c68f52674c9fb33aef52dcf399755519	MT					----'
	echo 'fdfd811849cc2fadebc929bb925902e5	3						----'
	echo '094d037050cad692b57ea12c4fef790f	----				chr6_qbl_hap6'
	echo '18c17e1641ef04873b15f40f6c8659a4	----				chr6_cox_hap2'
	echo '1e86411d73e6f00a10590f976be01623	----				chrY'
	echo '2a3c677c426a10e137883ae1ffb8da3f	----				chr6_dbb_hap3'
	echo '3b6d666200e72bcc036bf88a4d7e0749	----				chr6_ssto_hap7'
	echo '641e4338fa8d52a5b781bd2a2c08d3c3	----				chr3'
	echo '9d51d4152174461cd6715c7ddc588dc8	----				chr6_mann_hap4'
	echo 'd2ed829b8a1628d16cbeee88e88e39eb	----				chrM'
	echo 'd89517b400226d3b56e753972a7cad67	----				chr17_ctg5_hap1'
	echo 'efed415dd8742349cb7aaca054675b9a	----				chr6_mcf_hap5'
	echo 'fa24f81b680df26bcfb6d69b784fbe36	----				chr4_ctg9_hap1'
	echo 'fe71bc63420d666884f37a3ad79f3317	----				chr6_apd_hap1'
}

function printDiffTable() {

  echo 'The following table summarizes the differences between'
  echo "The Broad Institute's B37 reference and UCSC's HG19 reference:"
  echo ''
  printSequenceDictionaryComparison 
  echo ''
  echo 'The table is grouped by the MD5SUM of the contents of each'
  echo 'contig.'
  echo 'Entries that are blank indicate a contig that does not have'
  echo 'a matching entry in the other reference (e.g. B37 NC_007605,'
  echo 'the Epstein-Barr virus).'
  echo 'Some major contigs contain differences (e.g. chr3).  This is'
  echo 'generally because of masked IUPAC bases (N, rather than '
  echo 'A, T, G, or C).'
  echo 'However, care should be taken to ensure that the resulting'
  echo 'liftover files are used judiciously'
  echo ''
}

################################################################################

# Do all the interactive stuff if the script is called by a user.
# The bash equivalent to `if __name__ == '__main__':
if ${ISCALLEDBYUSER} ; then
  
  #========================================

  trap at_exit EXIT 
  
  #========================================

  #Check given arguments:
  if [[ $# -gt $MAXARGS ]] ; then
    usage
    exit 1
  elif [[ $# -lt $MINARGS ]] ; then
    usage
    exit 2
  fi

  # Make sure we have all the required commands:
  checkPrerequisites $PREREQUISITES 
  r=$?
  [[ $r -ne 0 ]] && exit 4

  #----------------------------------
  #Read args:

  while [ $# -gt 0 ] ; do

    case "$1" in
      -h|--h|--help|-help|help)
        usage;
        exit 0;
        ;;
      --showDiffs)
        printDiffTable
        exit 0
        ;;
      *)
        error "${SCRIPTNAME}: invalid option: $1" 
        simpleUsage
        echo "Try \`$SCRIPTNAME --help' for more information."
        exit 3;
        ;; 
    esac

    #Get next argument in $1:
    shift
  done

  #----------------------------------
  # Do real work here.

  error "Retrieving baseline hg19->hg38 chain file..."
  hg19LiftoverFile=$( makeTemp )
  curl "${HG19_HG38_CHAIN_FILE_URL}" 2>/dev/null | gunzip > ${hg19LiftoverFile}

	error "Modifying contig names and creating chain file..."
  # Do the appropriate change here:
  awk '{                                  \
    if ( $1 == "chain" ) {                \
      if ( $3 ~ /^chrM/ ) {               \
        $3 = "MT"                         \
      }                                   \
      else if ( $3 ~ /chr11_gl000202_random/ ) {               \
        $3 = "GL000202.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000244/ ) {               \
        $3 = "GL000244.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000235/ ) {               \
        $3 = "GL000235.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000238/ ) {               \
        $3 = "GL000238.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000226/ ) {               \
        $3 = "GL000226.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000218/ ) {               \
        $3 = "GL000218.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000249/ ) {               \
        $3 = "GL000249.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000242/ ) {               \
        $3 = "GL000242.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000221/ ) {               \
        $3 = "GL000221.1"                        \
      }                                                \
      else if ( $3 ~ /chr1_gl000192_random/ ) {               \
        $3 = "GL000192.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000223/ ) {               \
        $3 = "GL000223.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000232/ ) {               \
        $3 = "GL000232.1"                        \
      }                                                \
      else if ( $3 ~ /chr17_gl000206_random/ ) {               \
        $3 = "GL000206.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000240/ ) {               \
        $3 = "GL000240.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000214/ ) {               \
        $3 = "GL000214.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000212/ ) {               \
        $3 = "GL000212.1"                        \
      }                                                \
      else if ( $3 ~ /chr9_gl000199_random/ ) {               \
        $3 = "GL000199.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000248/ ) {               \
        $3 = "GL000248.1"                        \
      }                                                \
      else if ( $3 ~ /chr7_gl000195_random/ ) {               \
        $3 = "GL000195.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000215/ ) {               \
        $3 = "GL000215.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000225/ ) {               \
        $3 = "GL000225.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000216/ ) {               \
        $3 = "GL000216.1"                        \
      }                                                \
      else if ( $3 ~ /chr4_gl000194_random/ ) {               \
        $3 = "GL000194.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000217/ ) {               \
        $3 = "GL000217.1"                        \
      }                                                \
      else if ( $3 ~ /chr8_gl000197_random/ ) {               \
        $3 = "GL000197.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000222/ ) {               \
        $3 = "GL000222.1"                        \
      }                                                \
      else if ( $3 ~ /chr9_gl000200_random/ ) {               \
        $3 = "GL000200.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000211/ ) {               \
        $3 = "GL000211.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000247/ ) {               \
        $3 = "GL000247.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000233/ ) {               \
        $3 = "GL000233.1"                        \
      }                                                \
      else if ( $3 ~ /chr21_gl000210_random/ ) {               \
        $3 = "GL000210.1"                        \
      }                                                \
      else if ( $3 ~ /chr9_gl000198_random/ ) {               \
        $3 = "GL000198.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000245/ ) {               \
        $3 = "GL000245.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000234/ ) {               \
        $3 = "GL000234.1"                        \
      }                                                \
      else if ( $3 ~ /chr17_gl000203_random/ ) {               \
        $3 = "GL000203.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000239/ ) {               \
        $3 = "GL000239.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000213/ ) {               \
        $3 = "GL000213.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000227/ ) {               \
        $3 = "GL000227.1"                        \
      }                                                \
      else if ( $3 ~ /chr19_gl000208_random/ ) {               \
        $3 = "GL000208.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000230/ ) {               \
        $3 = "GL000230.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000231/ ) {               \
        $3 = "GL000231.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000228/ ) {               \
        $3 = "GL000228.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000243/ ) {               \
        $3 = "GL000243.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000229/ ) {               \
        $3 = "GL000229.1"                        \
      }                                                \
      else if ( $3 ~ /chr17_gl000205_random/ ) {               \
        $3 = "GL000205.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000224/ ) {               \
        $3 = "GL000224.1"                        \
      }                                                \
      else if ( $3 ~ /chr1_gl000191_random/ ) {               \
        $3 = "GL000191.1"                        \
      }                                                \
      else if ( $3 ~ /chr8_gl000196_random/ ) {               \
        $3 = "GL000196.1"                        \
      }                                                \
      else if ( $3 ~ /chr4_gl000193_random/ ) {               \
        $3 = "GL000193.1"                        \
      }                                                \
      else if ( $3 ~ /chr9_gl000201_random/ ) {               \
        $3 = "GL000201.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000237/ ) {               \
        $3 = "GL000237.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000246/ ) {               \
        $3 = "GL000246.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000241/ ) {               \
        $3 = "GL000241.1"                        \
      }                                                \
      else if ( $3 ~ /chr17_gl000204_random/ ) {               \
        $3 = "GL000204.1"                        \
      }                                                \
      else if ( $3 ~ /chr18_gl000207_random/ ) {               \
        $3 = "GL000207.1"                        \
      }                                                \
      else if ( $3 ~ /chr19_gl000209_random/ ) {               \
        $3 = "GL000209.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000219/ ) {               \
        $3 = "GL000219.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000220/ ) {               \
        $3 = "GL000220.1"                        \
      }                                                \
      else if ( $3 ~ /chrUn_gl000236/ ) {               \
        $3 = "GL000236.1"                        \
      }                                               \
      else if ( $3 ~ /chr6_qbl_hap6/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr6_cox_hap2/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr6_dbb_hap3/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr6_ssto_hap7/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr6_mann_hap4/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr17_ctg5_hap1/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr6_mcf_hap5/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr4_ctg9_hap1/ ) {               \
        $3 = $3                                         \
      }                                                \
      else if ( $3 ~ /chr6_apd_hap1/ ) {               \
        $3 = $3                                         \
      }                                                \
			else if ( $3 ~ /^chr[0123456789XY]*/ ) {           \
      	$3 = substr($3,4)                                  \
			}                                   \
    }                                     \
    ; print                               \
  }' ${hg19LiftoverFile} > ${outChainFile}

	error "Created a b37 -> hg38 chain file: ${outChainFile}"

  exit 0
fi

