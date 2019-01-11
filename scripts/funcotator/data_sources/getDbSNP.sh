#!/usr/bin/env bash

################################################################################

# Creates a set of files that allows for creating annotations based off dbSNP
# Pulled directly from the ncbi FTP site.

################################################################################

#Setup variables for the script:
SCRIPTDIR="$( cd -P "$( dirname "$0" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=0

FTP_BASE_URL='ftp://ftp.ncbi.nih.gov/snp/organisms/'
BUILD_NUMBER='150'

DATA_SOURCE_NAME="dbSNP"
OUT_DIR_NAME='dbsnp'

SRC_FILE_BASE_NAME="All_"
#SRC_FILE_BASE_NAME="common_all_"

SRC_FILE_REGEX="${SRC_FILE_BASE_NAME}\\d+.vcf.gz\\s*\$"
MD5_FILE_REGEX="${SRC_FILE_BASE_NAME}\\d+.vcf.gz.md5\\s*\$"
TBI_FILE_REGEX="${SRC_FILE_BASE_NAME}\\d+.vcf.gz.tbi\\s*\$"

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME [OPTIONS] ..."
  echo -e "Creates the data sources folder for dbSnp for the GATK Funcotator tool."
}

#Define a usage function:
function usage()
{
  simpleUsage
  echo -e ""
  echo -e "Will download all data sources directly from the NCBI website:"
  echo -e "    ${FTP_BASE_URL}"
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  UNKNOWN ARGUMENT"
  echo -e "  4  BAD CHECKSUM"
  echo -e "  5  OUTPUT DIRECTORY ALREADY EXISTS"
  echo -e "  6  COULD NOT FIND BGZIP UTILITY"
  echo -e ""
}

#Display a message to std error:
function error()
{
  echo "$1" 2>&1
}

TMPFILELIST=''
function makeTemp()
{
  local f
  f=$( mktemp )
  TMPFILELIST="${TMPFILELIST} $f"
  echo $f
}

function cleanTempVars()
{
  rm -f ${TMPFILELIST}
}

function at_exit()
{
  cleanTempVars
}

################################################################################


function createConfigFile() {

    local dataSourceName=$1
    local version=$2
    local srcFile=$3
    local originLocation=$4

    echo "name = ${dataSourceName}"
    echo "version = ${version}"
    echo "src_file = ${srcFile}"
    echo "origin_location = ${originLocation}"
    echo "preprocessing_script = ${SCRIPTNAME}"
    echo ""
    echo "# Supported types:"
    echo "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID"
    echo "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location"
    echo "# gencode      -- Custom datasource class for GENCODE"
    echo "# cosmic       -- Custom datasource class for COSMIC"
    echo "# vcf          -- Custom datasource class for Variant Call Format (VCF) files"
    echo "type = vcf"
    echo ""
    echo "# Required field for GENCODE files."
    echo "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:"
    echo "gencode_fasta_path ="
    echo ""
    echo "# Required field for simpleXSV files."
    echo "# Valid values:"
    echo "#     GENE_NAME"
    echo "#     TRANSCRIPT_ID"
    echo "xsv_key ="
    echo ""
    echo "# Required field for simpleXSV files."
    echo "# The 0-based index of the column containing the key on which to match"
    echo "xsv_key_column ="
    echo ""
    echo "# Required field for simpleXSV AND locatableXSV files."
    echo "# The delimiter by which to split the XSV file into columns."
    echo "xsv_delimiter ="
    echo ""
    echo "# Required field for simpleXSV files."
    echo "# Whether to permissively match the number of columns in the header and data rows"
    echo "# Valid values:"
    echo "#     true"
    echo "#     false"
    echo "xsv_permissive_cols ="
    echo ""
    echo "# Required field for locatableXSV files."
    echo "# The 0-based index of the column containing the contig for each row"
    echo "contig_column ="
    echo ""
    echo "# Required field for locatableXSV files."
    echo "# The 0-based index of the column containing the start position for each row"
    echo "start_column ="
    echo ""
    echo "# Required field for locatableXSV files."
    echo "# The 0-based index of the column containing the end position for each row"
    echo "end_column ="
    echo ""

}

function downloadAndVerifyVcfFiles() {

    local remoteFolder=$1
    local outputFolder=$2
    local filePrefix=$3

    local listingFile=$( makeTemp )
    local indentSpace="    "
    local version=$( echo "${remoteFolder}" | sed "s#.*human_\\(.*b${BUILD_NUMBER}\\).*#\\1#g" )

    curl ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/ 2>/dev/null > ${listingFile}

    vcfFile=$( cat ${listingFile} | awk '{print $9}' | grep -E ${SRC_FILE_REGEX} )
    tbiFile=$( cat ${listingFile} | awk '{print $9}' | grep -E ${TBI_FILE_REGEX} )
    md5File=$( cat ${listingFile} | awk '{print $9}' | grep -E ${MD5_FILE_REGEX} )

    echo "${indentSpace}Retrieving MD5 sum file: ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${md5File} ... "
    wget ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${md5File}

    # Get the VCF file, then make sure that the contig names are correct for HG19 (if applicable)
    echo "${indentSpace}Retrieving VCF file: ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${vcfFile} ... "
    if [[ "${filePrefix}" == "hg19" ]] ; then
        curl ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${vcfFile} | gunzip | sed -e 's#^\([0-9][0-9]*\)#chr\1#' -e 's#^MT#chrM#' -e 's#^X#chrX#' -e 's#^Y#chrY#' | bgzip > ${vcfFile}
    else
        wget ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${vcfFile}

        echo "${indentSpace}Retrieving VCF Index file: ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${tbiFile} ... "
        wget ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${tbiFile}

        # We can only verify the checksum with hg38 because we modify the hg19 file as we stream it in:
        echo "${indentSpace}Verifying VCF checksum ..."
        if [[ "$(uname)" == "Darwin" ]] ; then
            which md5sum-lite &> /dev/null
            r=$?
            if [ $r == 0 ] ; then
                checksum=$( md5sum-lite ${vcfFile} | awk '{print $1}' | sed -e 's#^[ \t]*##g' -e 's#[ \t]*$##g' )
                expected=$( head -n1 ${md5File} | awk '{print $1}' | sed -e 's#^[ \t]*##g' -e 's#[ \t]*$##g' )

                if [[ "${checksum}" != "${expected}" ]] ; then
                    error "DOWNLOADED FILE IS CORRUPT!  (${checksum} != ${expected})"
                    error "FAILING"
                    exit 4
                fi
            else
                error "Unable to validate md5sum of file: cannot locate 'md5sum-lite'.  Use these data with caution."
            fi
        else
            which md5sum &> /dev/null
            r=$?
            if [ $r == 0 ] ; then
                checksum=$( md5sum ${vcfFile} | awk '{print $1}' | sed -e 's#^[ \t]*##g' -e 's#[ \t]*$##g' )
                expected=$( head -n1 ${md5File} | awk '{print $1}' | sed -e 's#^[ \t]*##g' -e 's#[ \t]*$##g' )

                if [[ "${checksum}" != "${expected}" ]] ; then
                    error "DOWNLOADED FILE IS CORRUPT!  (${checksum} != ${expected})"
                    error "FAILING"
                    exit 4
                fi
            else
                error "Unable to validate md5sum of file: cannot locate 'md5sum'.  Use these data with caution."
            fi
        fi
    fi

    # Now put it in the right place and clean up:
    echo "${indentSpace}Creating output directory ..."
    mkdir -p ${outputFolder}

    echo "${indentSpace}Moving files to output directory ..."
    mv ${vcfFile} ${outputFolder}/${filePrefix}_${vcfFile}
    if [[ ! "${filePrefix}" == "hg19" ]] ; then
        mv ${tbiFile} ${outputFolder}/${filePrefix}_${tbiFile}
        rm ${md5File}
    fi

    echo "${indentSpace}Creating Config File ... "
    createConfigFile "${DATA_SOURCE_NAME}" "${version}" ${filePrefix}_${vcfFile} "ftp://ftp.ncbi.nih.gov/snp/organisms/${remoteFolder}/VCF/${vcfFile}" > ${outputFolder}/${DATA_SOURCE_NAME}.config
}

################################################################################

trap at_exit EXIT

################################################################################

#Check given arguments:
if [[ $# -gt $MAXARGS ]] ; then
  usage
  exit 1
elif [[ $# -lt $MINARGS ]] ; then
  usage
  exit 2
fi

################################################################################

#Read args:
while [ $# -gt 0 ] ; do

  case "$1" in
    -h|--h|--help|-help|help)
      usage;
      exit 0;
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

################################################################################
# Do the work here:

# Make sure we don't have anything in our out folder already:
if [[ -d ${OUT_DIR_NAME} ]] ; then
    error "Output directory already exists: ${OUT_DIR_NAME} - aborting!"
    exit 5
fi

# Make sure bgzip exists:
which bgzip > /dev/null
r=$?
if [[ $r -ne 0 ]] ; then
    error "bgzip utility not found on path.  You must have bgzip installed to get the dbSNP resource.  Please install bgzip and try again. - aborting!"
    exit 6
fi

echo "Querying NCBI listing..."
tmpListing="$( makeTemp )"
curl ftp://ftp.ncbi.nih.gov/snp/organisms/ 2>/dev/null > ${tmpListing}

# Get the link for HG19:
echo "Processing HG19 ..."
hg19Version=$( cat ${tmpListing} | awk '{print $9}' | grep 'human' | grep 'GRCh37' | grep "b${BUILD_NUMBER}" )

downloadAndVerifyVcfFiles ${hg19Version} ${OUT_DIR_NAME}/hg19 "hg19"


# Get the link for HG38:
echo "Processing HG38 ..."
hg38Version=$( cat ${tmpListing} | awk '{print $9}' | grep 'human' | grep 'GRCh38' | grep "b${BUILD_NUMBER}" )

downloadAndVerifyVcfFiles ${hg38Version} ${OUT_DIR_NAME}/hg38 "hg38"

echo -e "\033[1;33;40m##################################################################################\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m# \033[1;5;37;41mWARNING\033[0;0m: You \033[4;37;40mMUST\033[0;0m index the VCF files for \033[4;37;40mHG19\033[0;0m #\033[0;0m"
echo -e "\033[1;33;40m#             before using this data source.                                     #\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m#             Use gatk IndexFeatureFile <GTF_FILE>                               #\033[0;0m"
echo -e "\033[1;33;40m#                                                                                #\033[0;0m"
echo -e "\033[1;33;40m##################################################################################\033[0;0m"

exit 0
