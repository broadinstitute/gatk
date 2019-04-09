#!/usr/bin/env bash

################################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( python -c "import os;print os.path.realpath(\"${BASH_SOURCE[0]}\")" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=4
PREREQUISITES=""
BUILD_DOCKER_SCRIPT="${SCRIPTDIR}/../../../../build_docker.sh"

# Determine if this shell is interactive:
ISCALLEDBYUSER=true
[[ "${BASH_SOURCE[0]}" != "${0}" ]] && ISCALLEDBYUSER=false

# Required for the aliased call to checkPipeStatus:
shopt -s expand_aliases

################################################################################

#DEFAULT_DOCKERHUB_REPO_NAME="broadinstitute/gatk-nightly"
DEFAULT_DOCKERHUB_REPO_NAME="us.gcr.io/broad-dsde-methods/broad-gatk-snapshots"

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME [-t DOCKERTAG] [-b BRANCH]"
  echo -e "Build a docker image and push it to docker hub."
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e ""
  echo -e "Constructs a GATK docker image then immediately uploads it to the specified location on docker hub."
  echo -e "If no BRANCH is given will use the current branch's last commit.  This will fail if the commit has not been"
  echo -e "pushed to master."
  echo -e "If no DOCKERTAG is given will commit to broadinstitute/gatk-nightly ."
  echo -e ""
  # List prerequisites, if any:
  if [[ ${#PREREQUISITES} -ne 0 ]] ; then
    echo -e "Requires the following programs to be installed:"
    for prereq in ${PREREQUISITES} ; do 
      echo -e "  ${prereq}"
    done
    echo
  fi
  echo -e "  -t dockertag the full name for where to upload the docker image"
  echo -e "  -b BRANCH the branch of gatk off of which to build an image"
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  UNKNOWN ARGUMENT"
  echo -e "  4  MISSING PREREQUISITE"
  echo -e "  5  DOCKER DAEMON NOT RUNNING"
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
function _checkPipeStatus()
{

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

function checkPrerequisites
{
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

trap at_exit EXIT 

################################################################################

# Do all the interactive stuff if the script is called by a user.
# The bash equivalent to `if __name__ == '__main__':
if ${ISCALLEDBYUSER} ; then

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
			-b|-B)
				shift
				gitbranch="${1}"
				;;
		    -t|-T)
				shift
				dockerTag="${1}"
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
    hash=$( git show ${gitbranch} | head -n 1 | grep -Eo "[A-Za-z0-9]+$" )

    # Get a docker tag if it's not set:
    if [[ ${#dockerTag} -eq 0 ]] ; then
        d=$(date +%Y-%m-%d)
        dockerTag="${DEFAULT_DOCKERHUB_REPO_NAME}:${d}-$(git describe)-SNAPSHOT"
    fi

    docker images -q > /dev/null
    r=$?
    [[ $r -ne 0 ]] && exit 5

    oldTopImage=$(docker images -q | head -n 1)

    error "============================================================================="
    error "Building docker image..."
    buildDockerTmpLocation=$( mktemp -d )
    ${BUILD_DOCKER_SCRIPT} -e ${hash} -s -u -d ${buildDockerTmpLocation}
    newTopImage=$(docker images -q | head -n 1)
    error "============================================================================="

    error "Old top docker image: ${oldTopImage}"
    error "New top docker image: ${newTopImage}"

    if [[ "${oldTopImage}" != "${newTopImage}" ]] ; then
        error "============================================================================="
        error "Pushing image to docker tag: ${dockerTag}"
        docker tag ${newTopImage} ${dockerTag}

        # Make sure we push the correct way:
        echo "${DEFAULT_DOCKERHUB_REPO_NAME}" | grep -q '^gcr.io/'
        if [[ $? -eq 0 ]] ; then
            gcloud docker -- push ${dockerTag}
        else
            docker push ${dockerTag}
        fi

        docker image rm -f ${newTopImage}
    fi
	exit 0
fi

