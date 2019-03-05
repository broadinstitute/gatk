<a name="0"></a>
## 0 - Introduction
This page explains what **Funcotator** is and how to run it.

<a name="0.1"></a>
### 0.1 - Table of Contents
0. [0.0 Introduction](#0)
    1. [0.1 Table of Contents](#0.1)
1. [1.0 Funcotator Background Information](#1)
    1. [1.1 Data Sources](#1.1) 
        1. [1.1.1 Data Source Folders](#1.1.1)
        2. [1.1.2 Pre-Packaged Data Sources](#1.1.2)
            1. [1.1.2.1 Downloading Pre-Packaged Data Sources](#1.1.2.1)
            2. [1.1.2.2 gnomAD](#1.1.2.2)
                1. [1.1.2.2.1 Enabling gnomAD](#1.1.2.2.1)
                2. [1.1.2.2.2 Included gnomAD Fields](#1.1.2.2.2)
        3. [1.1.3 Data Source Downloader Tool](#1.1.3)
        4. [1.1.4 Disabling Data Sourcesl](#1.1.4)
        5. [1.1.5 User-Defined Data Sources](#1.1.5)
            1. [1.1.5.1 Configuration File Format](#1.1.5.1)
                1. [1.1.5.1.1 Simple XSV Config File Example](#1.1.5.1.1)
                2. [1.1.5.1.2 Locatable XSV Config File Example](#1.1.5.1.2)
            2. [1.1.5.2 Cloud Data Sources](#1.1.5.2)
        6. [1.1.6 Data Source Versioning](#1.1.6)
    2. [1.2 Input Variant Data Formats](#1.2)
    3. [1.3 Output](#1.3)
        1. [1.3.1 Output Data Formats](#1.3.1)
            1. [1.3.1.1 VCF Format](#1.3.1.1)
            2. [1.3.1.2 MAF Format](#1.3.1.2)
        2. [1.3.2 Annotations for Pre-Packaged Data Sources](#1.3.2)
            1. [1.3.2.1 Gencode Annotation Specification](#1.3.2.1)
    4. [1.4 Reference Genome Versions](#1.4)
    5. [1.5 Comparisons with Oncotator](#1.5)
        1. [1.5.1 Funcotator / Oncotator Feature Comparison](#1.5.1)
        2. [1.5.2 Oncotator Bugs Compared With Funcotator](#1.5.2)
2. [2.0 Tutorial](#2)
    0. [2.0 Requirements](#2.0)
    1. [2.1 Running Funcotator in the GATK With Base Options](#2.1)
    2. [2.2 Optional Parameters](#2.2)
        1. [2.2.1 - --ignore-filtered-variants](#2.2.1)
        2. [2.2.2 - --transcript-selection-mode](#2.2.2)
        3. [2.2.1 - --transcript-list](#2.2.3)
        4. [2.2.2 - --annotation-default](#2.2.4)
        5. [2.2.1 - --annotation-override](#2.2.5)
        6. [2.2.2 - --allow-hg19-gencode-b37-contig-matching](#2.2.6)
3. [3.0 FAQ](#3)
4. [4.0 Known Issues](#4)
5. [5.0 Github](#5)
6. [6.0 Tool Documentation](#6)

----

<a name="1"></a>
## 1 - Funcotator Background Information

Funcotator (**FUNC**tional ann**OTATOR**) analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.

This tool allows a user to add their own annotations to variants based on a set of data sources.  Each data source can be customized to annotate a variant based on several matching criteria.  This allows a user to create their own custom annotations easily, without modifying any Java code.

An example Funcotator workflow based on the GATK Best Practices Somatic Pipeline is as follows:
![](https://us.v-cdn.net/5019796/uploads/editor/rt/fapszh63zvd0.png "")

<a name="1.1"></a>
### 1.1 - Data Sources 
Data sources are expected to be in folders that are specified as input arguments.  While multiple data source folders can be specified, **no two data sources can have the same name**.

<a name="1.1.1"></a>
### 1.1.1 - Data Source Folders

In each main data source folder, there should be sub-directories for each individual data source, with further sub-directories for a specific reference (e.g. _hg19_, _hg38_, etc.). In the reference-specific data source directory, there is a configuration file detailing information about the data source and how to match it to a variant. This configuration file is required.

An example of a data source directory is the following:
```
    dataSourcesFolder/
         Data_Source_1/
             hg19
                 data_source_1.config
                 data_source_1.data.file.one
                 data_source_1.data.file.two
                 data_source_1.data.file.three
                 ...
              hg38
                 data_source_1.config
                 data_source_1.data.file.one
                 data_source_1.data.file.two
                 data_source_1.data.file.three
                 ...
         Data_Source_2/
             hg19
                 data_source_2.config
                 data_source_2.data.file.one
                 data_source_2.data.file.two
                 data_source_2.data.file.three
                 ...
              hg38
                 data_source_2.config
                 data_source_2.data.file.one
                 data_source_2.data.file.two
                 data_source_2.data.file.three
                 ...
          ...
```
<a name="1.1.2"></a>
### 1.1.2 - Pre-Packaged Data Sources
The GATK includes two sets of pre-packaged data sources, allowing for [Funcotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php "Funcotator") use without (much) additional configuration.
These data source packages correspond to the **germline** and **somatic** use cases.
Broadly speaking, if you have a **germline VCF**, the **germline data sources** are what you want to use to start with.
Conversely, if you have a <strong>somatic VCF</strong>, the <strong>somatic data sources</strong> are what you want to use to start with.

<a name="1.1.2.1"></a>
### 1.1.2.1 - Downloading Pre-Packaged Data Sources
Versioned gzip archives of data source files are provided here:
● FTP: [ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/)
● Google Cloud Bucket: [gs://broad-public-datasets/funcotator/](https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator)

<a name="1.1.2.2"></a>
### 1.1.2.2 - gnomAD
The pre-packaged data sources include a subset of gnomAD, a large database of known variants.  This subset contains a greatly reduced subset of INFO fields, primarily containing allele frequency data.  gnomAD is split into two parts - one based on exome data, one based on whole genome data.  These two data sources are not equivalent and for complete coverage using gnomAD, we recommend annotating with both.
Due to the size of gnomAD, it cannot be included in the data sources package directly.  Instead, the configuration data are present and point to a Google bucket in which
the gnomAD data reside.  This will cause [Funcotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php "Funcotator") to actively connect to that bucket when it is run.  
For this reason, **gnomAD is disabled by default**.

Because [Funcotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php "Funcotator") will query the Internet **when gnomAD is enabled, performance will be impacted** by the machine's Internet connection speed.
If this degradation is significant, you can localize gnomAD to the machine running [Funcotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php "Funcotator") to improve performance (however due to the size of gnomAD this may be impractical).

<a name="1.1.2.2.1"></a>
### 1.1.2.2.1 - Enabling gnomAD
To enable gnomAD, simply change directories to your data sources directory and untar the gnomAD tar.gz files:
```
cd DATA_SOURCES_DIR
tar -zxf gnomAD_exome.tar.gz
tar -zxf gnomAD_genome.tar.gz
```

<a name="1.1.2.2.2"></a>
### 1.1.2.2.2 - Included gnomAD Fields
The fields included in the pre-packaged gnomAD subset are the following:
<table>
<tr><th>Field Name</th><th>Field Description</th></tr>
<tr><td>AF</td><td>Allele Frequency, for each ALT allele, in the same order as listed</td></tr>
<tr><td>AF_afr</td><td>Alternate allele frequency in samples of African-American ancestry</td></tr>
<tr><td>AF_afr_female</td><td>Alternate allele frequency in female samples of African-American ancestry</td></tr>
<tr><td>AF_afr_male</td><td>Alternate allele frequency in male samples of African-American ancestry</td></tr>
<tr><td>AF_amr</td><td>Alternate allele frequency in samples of Latino ancestry</td></tr>
<tr><td>AF_amr_female</td><td>Alternate allele frequency in female samples of Latino ancestry</td></tr>
<tr><td>AF_amr_male</td><td>Alternate allele frequency in male samples of Latino ancestry</td></tr>
<tr><td>AF_asj</td><td>Alternate allele frequency in samples of Ashkenazi Jewish ancestry</td></tr>
<tr><td>AF_asj_female</td><td>Alternate allele frequency in female samples of Ashkenazi Jewish ancestry</td></tr>
<tr><td>AF_asj_male</td><td>Alternate allele frequency in male samples of Ashkenazi Jewish ancestry</td></tr>
<tr><td>AF_eas</td><td>Alternate allele frequency in samples of East Asian ancestry</td></tr>
<tr><td>AF_eas_female</td><td>Alternate allele frequency in female samples of East Asian ancestry</td></tr>
<tr><td>AF_eas_jpn</td><td>Alternate allele frequency in samples of Japanese ancestry</td></tr>
<tr><td>AF_eas_kor</td><td>Alternate allele frequency in samples of Korean ancestry</td></tr>
<tr><td>AF_eas_male</td><td>Alternate allele frequency in male samples of East Asian ancestry</td></tr>
<tr><td>AF_eas_oea</td><td>Alternate allele frequency in samples of non-Korean, non-Japanese East Asian ancestry</td></tr>
<tr><td>AF_female</td><td>Alternate allele frequency in female samples</td></tr>
<tr><td>AF_fin</td><td>Alternate allele frequency in samples of Finnish ancestry</td></tr>
<tr><td>AF_fin_female</td><td>Alternate allele frequency in female samples of Finnish ancestry</td></tr>
<tr><td>AF_fin_male</td><td>Alternate allele frequency in male samples of Finnish ancestry</td></tr>
<tr><td>AF_male</td><td>Alternate allele frequency in male samples</td></tr>
<tr><td>AF_nfe</td><td>Alternate allele frequency in samples of non-Finnish European ancestry</td></tr>
<tr><td>AF_nfe_bgr</td><td>Alternate allele frequency in samples of Bulgarian ancestry</td></tr>
<tr><td>AF_nfe_est</td><td>Alternate allele frequency in samples of Estonian ancestry</td></tr>
<tr><td>AF_nfe_female</td><td>Alternate allele frequency in female samples of non-Finnish European ancestry</td></tr>
<tr><td>AF_nfe_male</td><td>Alternate allele frequency in male samples of non-Finnish European ancestry</td></tr>
<tr><td>AF_nfe_nwe</td><td>Alternate allele frequency in samples of North-Western European ancestry</td></tr>
<tr><td>AF_nfe_onf</td><td>Alternate allele frequency in samples of non-Finnish but otherwise indeterminate European ancestry</td></tr>
<tr><td>AF_nfe_seu</td><td>Alternate allele frequency in samples of Southern European ancestry</td></tr>
<tr><td>AF_nfe_swe</td><td>Alternate allele frequency in samples of Swedish ancestry</td></tr>
<tr><td>AF_oth</td><td>Alternate allele frequency in samples of uncertain ancestry</td></tr>
<tr><td>AF_oth_female</td><td>Alternate allele frequency in female samples of uncertain ancestry</td></tr>
<tr><td>AF_oth_male</td><td>Alternate allele frequency in male samples of uncertain ancestry</td></tr>
<tr><td>AF_popmax</td><td>Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry)</td></tr>
<tr><td>AF_raw</td><td>Alternate allele frequency in samples, before removing low-confidence genotypes</td></tr>
<tr><td>AF_sas</td><td>Alternate allele frequency in samples of South Asian ancestry</td></tr>
<tr><td>AF_sas_female</td><td>Alternate allele frequency in female samples of South Asian ancestry</td></tr>
<tr><td>AF_sas_male</td><td>Alternate allele frequency in male samples of South Asian ancestry</td></tr>
<tr><td>OriginalAlleles*</td><td>A list of the original alleles (including REF) of the variant prior to liftover.  If the alleles were not changed during liftover, this attribute will be omitted.</td></tr>
<tr><td>OriginalContig*</td><td>The name of the source contig/chromosome prior to liftover.</td></tr>
<tr><td>OriginalStart*</td><td>The position of the variant on the source contig prior to liftover.</td></tr>
<tr><td>ReverseComplementedAlleles*</td><td>The REF and the ALT alleles have been reverse complemented in liftover since the mapping from the previous reference to the current one was on the negative strand.</td></tr>
<tr><td>SwappedAlleles*</td><td>The REF and the ALT alleles have been swapped in liftover due to changes in the reference. It is possible that not all INFO annotations reflect this swap, and in the genotypes, only the GT, PL, and AD fields have been modified. You should check the TAGS_TO_REVERSE parameter that was used during the LiftOver to be sure.</td></tr>
</table>
\* - only available in *hg38*

<a name="1.1.3"></a>
### 1.1.3 - Data Source Downloader Tool
To improve ease-of-use of [Funcotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php "Funcotator"), there is a tool to download the pre-packaged data sources to the user's machine.
This tool is the **[FuncotatorDataSourceDownloader](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_FuncotatorDataSourceDownloader.php "FuncotatorDataSourceDownloader")** and can be run to retrieve the pre-packaged data sources from the google bucket and localize them to the machine on which it is run.
Briefly:
For **somatic** data sources:
```
./gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
```
For **germline** data sources:
```
./gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download
```

<a name="1.1.4"></a>
### 1.1.4 - Disabling Data Sources
A data source can be disabled by removing the folder containing the configuration file for that source.  This can be done on a per-reference basis.  If the entire data source should be disabled, the entire top-level data source folder can be removed.

<a name="1.1.5"></a>
### 1.1.5 - User-Defined Data Sources
Users can define their own data sources by creating a new correctly-formatted data source sub-directory in the main data sources folder. In this sub-directory, the user must create an additional folder for the reference for which the data source is valid. If the data source is valid for multiple references, then multiple reference folders should be created. Inside each reference folder, the user should place the file(s) containing the data for the data source. Additionally the user **must** create a _configuration file_ containing metadata about the data source.

There are several formats allowed for data sources:

<table>
<tr><th>Data Format Class</th><th>Data Source Description</th></tr>
<tr><td>simpleXSV</td><td>Separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID</td></tr>
<tr><td>locatableXSV</td><td>Separated value table (e.g. CSV), keyed off a genome location</td></tr>
<tr><td>gencode</td><td>Class for GENCODE data files (gtf format)</td></tr>
<tr><td>cosmic</td><td>Class for COSMIC data</td></tr>
<tr><td>vcf</td><td>Class for Variant Call Format (VCF) files</td></tr>
</table>

Two of the most useful are arbitrarily separated value (XSV) files, such as comma-separated value (CSV), tab-separated value (TSV). These files contain a table of data that can be matched to a variant by _gene name_, _transcript ID_, or _genome position_. In the case of _gene name_ and _transcript ID_, one column must contain the gene name or transcript ID for each row's data.

* For gene name, when a variant is annotated with a gene name that exactly matches an entry in the gene name column for a row, that row's other fields will be added as annotations to the variant.
* For transcript ID, when a variant is annotated with a transcript ID that exactly matches an entry in the transcript ID column for a row, that row's other fields will be added as annotations to the variant.
* For genome position, one column must contain the contig ID, another column must contain the start position (1-based, inclusive), and a column must contain the stop position (1-based, inclusive). The start and stop columns may be the same column. When a variant is annotated with a genome position that overlaps an entry in the three genome position columns for a row, that row's other fields will be added as annotations to the variant.

<a name="1.1.5.1"></a>
#### 1.1.5.1 - Configuration File Format
The configuration file is a standard Java properties-style configuration file with key-value pairs. This file name **must end in .config**.

<a name="1.1.5.1.1"></a>
##### 1.1.5.1.1 - Simple XSV

The following is an example of a Locatable XSV configuration file (for the Familial Cancer Genes data source):

```
name = Familial_Cancer_Genes
version = 20110905
src_file = Familial_Cancer_Genes.no_dupes.tsv
origin_location = oncotator_v1_ds_April052016.tar.gz
preprocessing_script = UNKNOWN

# Whether this data source is for the b37 reference.
# Required and defaults to false.
isB37DataSource = false

# Supported types:
# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
# gencode      -- Custom datasource class for GENCODE
#	cosmic       -- Custom datasource class for COSMIC
# vcf          -- Custom datasource class for Variant Call Format (VCF) files
type = simpleXSV

# Required field for GENCODE files.
# Path to the FASTA file from which to load the sequences for GENCODE transcripts:
gencode_fasta_path =

# Required field for GENCODE files.
# NCBI build version (either hg19 or hg38):
ncbi_build_version = 

# Required field for simpleXSV files.
# Valid values:
#     GENE_NAME
#     TRANSCRIPT_ID
xsv_key = GENE_NAME

# Required field for simpleXSV files.
# The 0-based index of the column containing the key on which to match
xsv_key_column = 2

# Required field for simpleXSV AND locatableXSV files.
# The delimiter by which to split the XSV file into columns.
xsv_delimiter = \t

# Required field for simpleXSV files.
# Whether to permissively match the number of columns in the header and data rows
# Valid values:
#     true
#     false
xsv_permissive_cols = true

# Required field for locatableXSV files.
# The Name or 0-based index of the column containing the contig for each row
contig_column =

# Required field for locatableXSV files.
# The Name or 0-based index of the column containing the start position for each row
start_column =

# Required field for locatableXSV files.
# The Name or 0-based index of the column containing the end position for each row
end_column =
```

<a name="1.1.5.1.2"></a>
##### 1.1.5.1.2 - Locatable XSV

The following is an example of a Locatable XSV configuration file (for the ORegAnno data source):

```
name = Oreganno
version = 20160119
src_file = oreganno.tsv
origin_location = http://www.oreganno.org/dump/ORegAnno_Combined_2016.01.19.tsv
preprocessing_script = getOreganno.py

# Whether this data source is for the b37 reference.
# Required and defaults to false.
isB37DataSource = false

# Supported types:
# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
# gencode      -- Custom datasource class for GENCODE
#	cosmic       -- Custom datasource class for COSMIC
# vcf          -- Custom datasource class for Variant Call Format (VCF) files
type = locatableXSV

# Required field for GENCODE files.
# Path to the FASTA file from which to load the sequences for GENCODE transcripts:
gencode_fasta_path =

# Required field for GENCODE files.
# NCBI build version (either hg19 or hg38):
ncbi_build_version = 

# Required field for simpleXSV files.
# Valid values:
#     GENE_NAME
#     TRANSCRIPT_ID
xsv_key =

# Required field for simpleXSV files.
# The 0-based index of the column containing the key on which to match
xsv_key_column =

# Required field for simpleXSV AND locatableXSV files.
# The delimiter by which to split the XSV file into columns.
xsv_delimiter = \t

# Required field for simpleXSV files.
# Whether to permissively match the number of columns in the header and data rows
# Valid values:
#     true
#     false
xsv_permissive_cols = true

# Required field for locatableXSV files.
# The Name or 0-based index of the column containing the contig for each row
contig_column = 1

# Required field for locatableXSV files.
# The Name or 0-based index of the column containing the start position for each row
start_column = 2

# Required field for locatableXSV files.
# The Name or 0-based index of the column containing the end position for each row
end_column = 3
```

<a name="1.1.5.2"></a>
#### 1.1.5.2 - Cloud Data Sources

Funcotator allows for data sources with source files that live on the cloud, enabling users to annotate with data sources that are not physically present on the machines running Funcotator.
To create a data source based on the cloud, create a configuration file for that data source and put the cloud URL in as the src_file property (see [Configuration File Format](#1.1.5.1) for details).
E.g.:
```
 ...
src_file = gs://broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz
...
```


<a name="1.1.6"></a>
### 1.1.6 - Data Source Versioning
Each release of the data sources contains a version number.  Newer versions of Funcotator require minimum versions of data sources in order to run.  If a new version of Funcotator is run with an older version of the data sources, an error will be thrown prompting the user to download a new release of the data sources.  

Similarly newer releases of the data source packages are not reverse compatible with older versions of Funcotator.  However, in this case Funcotator may or may not throw an error or warning.

To ensure compatibility when upgrading Funcotator, always download the latest data sources release.  Similarly, when updating data sources make sure to update Funcotator to the latest version.

<a name="1.2"></a>
### 1.2 - Input Variant Data Formats
Currently Funcotator can only accept input variants in the form of a [_VCF_](https://samtools.github.io/hts-specs/VCFv4.2.pdf "_VCF_") file.

<a name="1.3"></a>
### 1.3 - Output 

<a name="1.3.1"></a>
### 1.3.1 - Output Data Formats
Funcotator supports output in both [_VCF_](https://samtools.github.io/hts-specs/VCFv4.2.pdf "_VCF_") format and [_MAF_ ](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification "_MAF_ ") format.

<a name="1.3.1.1"></a>
#### 1.3.1.1 - VCF Output
[VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf "VCF") files will contain the annotations for each variant allele as part of a **custom _INFO_ tag** - `FUNCOTATION`.  This custom tag will contain a pipe-separated (|) list of annotations for each alternate allele on a given line of the VCF.  The VCF header will contain an INFO field comment line for the FUNCOTATION data describing the field name for each value in the pipe-separated list.  For variants with multiple alternate alleles, the INFO field will contain multiple lists of annotations (each list separated by a comma), the order of which corresponds to the alternate allele being annotated.  

For example: 

```
#fileformat=VCFv4.2
...
 #INFO=<ID=FUNCOTATION,Number=A,Type=String,Description="Functional annotation from the Funcotator tool.  Funcotation fields are: dbSNP_Val_Status|Center">
...
#CHROM  POS ID  REF ALT QUAL  FILTER  INFO
chr19 8914955 . C A 40  . FUNCOTATION=No Value|broad.mit.edu
```

In this example, the variant has one alternate allele (_A_) with two fields (_\_dbSNP_Val_Status_ and _Center_).  The values of the fields are:
<table>
<tr><th>Field Name</th><th>Field Value</th></tr>
<tr><td>dbSNP_Val_Status</td><td>No Value</td></tr>
<tr><td>Center</td><td>broad.mit.edu</td></tr>
</table>

For multiple alternate alleles: 
```
#fileformat=VCFv4.2
...
 #INFO=<ID=FUNCOTATION,Number=A,Type=String,Description="Functional annotation from the Funcotator tool.  Funcotation fields are: dbSNP_Val_Status|Center">
...
#CHROM  POS ID  REF ALT QUAL  FILTER  INFO
chr7 273846 . C A,G 40  . FUNCOTATION=No Value|broad.mit.edu,Big Value Here|brandeis.edu
```

In this example, the variant has one alternate allele (_A_) with two fields (_\_dbSNP_Val_Status_ and _Center_).  The values of the fields are:
<table>
<tr><th>Alternate Allele</th><th>Field Name</th><th>Field Value</th></tr>
<tr><td>A</td><td>dbSNP_Val_Status</td><td>No Value</td></tr>
<tr><td>A</td><td>Center</td><td>broad.mit.edu</td></tr>
<tr><td>G</td><td>dbSNP_Val_Status</td><td>Big Value Here</td></tr>
<tr><td>G</td><td>Center</td><td>brandeis.edu</td></tr>
</table>

This formatting is the result of limitations in the VCF file specification.

<a name="1.3.1.2"></a>
#### 1.3.1.2 - MAF Output
The [_MAF_ ](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification "_MAF_ ") format used in Funcotator is an extension of the standard TCGA MAF.  It is based on the MAF format specified for Oncotator [here under Output Format](https://portals.broadinstitute.org/oncotator/help/ "here under Output Format").  While the actual columns can vary (due to different data sources being used to create annotations), columns 1-67 will generally be the same.  

In the case of a variant with multiple alternate alleles, each alternate allele will be written to a separate line in the MAF file.  

<a name="1.3.2"></a>
### 1.3.2 - Annotations for Pre-Packaged Data Sources
The pre-packaged data sources will create a set of baseline, or default annotations for an input data set.
Most of these data sources copy and paste values from their source files into the output of Funcotator to create annotations.  In this sense they are trivial data sources.

<a name="1.3.2.1"></a>
#### 1.3.2.1 - Gencode Annotation Specification
Funcotator performs some processing on the input data to create the Gencode annotations.  Gencode is currently required, so Funcotator will create these annotations for all input variants.
See [this forum post](https://gatkforums.broadinstitute.org/dsde/discussion/23389/funcotator-annotation-specifications "this forum post") for the specification of Gencode annotations in Funcotator.

<a name="1.4"></a>
#### 1.4 - Reference Genome Versions
The two currently supported genomes for annotations **out of the box** are **hg19** and **hg38**.  This is due to the pre-packaged Gencode data sources being for those two references.  **Any reference genome with published Gencode data sources can be used.**

<a name="1.4.1"></a>
#### 1.4.1 - hg19 vs b37 Reference
The Broad Institute uses an alternate hg19 reference known as b37 for our sequencing.  UCSC uses the baseline hg19 reference.  These references are similar but different.  

Due to the Gencode data source being published by UCSC, the data sources all use the hg19 reference for hg19 data (as opposed to b37).  Funcotator detects when user data is from the b37 reference and forces the use of the hg19 data sources in this case.  The user is warned when this occurs.  Generally speaking this is OK, but due to the differences in the sequence data it is possible that some erroneous data will be created.  

This effect has not yet been quantified, but in most cases should not be appreciable.  For details, see [this forum post](https://gatkforums.broadinstitute.org/dsde/discussion/23390/grch37-hg19-b37-humang1kv37-human-reference-discrepancies "this forum post").

<a name="1.5"></a>
#### 1.5 - Comparisons with Oncotator
Oncotator is an older functional annotation tool developed by The Broad Institute.  Funcotator and Oncotator are fundamentally different tools with some similarities.  

While I maintain that a direct comparison should not be made, to address some inevitable questions some comparison highlights between Oncotator and Funcotator are in the following two tables:

<a name="1.5.1"></a>
#### 1.5.1 - Funcotator / Oncotator Feature Comparison
<table>
<tr><th></th><th>Funcotator</th><th>Oncotator</th><th>Notes</th></tr>
<tr><td>Override values for annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Default values for annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>VCF input</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>VCF output</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td>Annotation format b/w Funcotator and Oncotator differ.</td></tr>
<tr><td>MAF input</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>MAF output</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>TSV/maflite input</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Simple TSV output</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Removing datasources does not require developer</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>hg38 support</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Cloud datasources</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td>All data sources supported</td></tr>
<tr><td>Transcript override list</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Default config speed somatic (muts/min) (hg19)</td><td></td><td></td><td></td></tr>
<tr><td>Default config speed germline (muts/min) (hg19)</td><td></td><td bgcolor="#FF0000">A very long time....</td><td></td></tr>
<tr><td>Default config speed somatic (muts/min) (hg38)</td><td></td><td bgcolor="#FF0000">N/A</td><td></td></tr>
<tr><td>Default config speed germline (muts/min) (hg38)</td><td></td><td bgcolor="#FF0000">N/A</td><td></td></tr>
<tr><td>Documentation</td><td bgcolor="#00FF00">Tutorial; Specifications forum post; inclusion in workshop materials</td><td bgcolor="#FFFF00">Minimal support in forum</td><td></td></tr>
<tr><td>Manuscript</td><td bgcolor="#FFFF00">Planned</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>HGVS support</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>BigWig datasource support</td><td bgcolor="#FF0000">No</td><td bgcolor="#FFFF00">Linux only</td><td></td></tr>
<tr><td>Seg file input/output</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Transcript modes: canonical and most deleterious effect</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Transcript mode: ALL</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Exclude annotations/columns on CLI</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Automated datasource download tool</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Automated tool for creating datasources</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Web application</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td>Uses old version of Oncotator and datasources</td></tr>
<tr><td>Config file to specify CLI arguments</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td>GATK built-in command line arguments file</td></tr>
<tr><td>Simple MAF to VCF </td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Or </td><td></td><td></td><td></td></tr>
<tr><td>VCF to MAF conversion</td><td></td><td></td><td></td></tr>
<tr><td>Inferring ONPs</td><td bgcolor="#AAAAAA">No</td><td bgcolor="#FFFF00">Yes (Not recommended)</td><td>Mutect2 infers ONPs when calling variants.  This is not the job of a functional annotator.</td></tr>
<tr><td>Ignores filtered input variants</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Mitochondrial amino acid sequence rendering</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>gnomAD annotations</td><td bgcolor="#00FF00">Yes (cloud support)</td><td bgcolor="#FF0000">Not recommended</td><td>v2.1 support for hg19</td></tr>
<tr><td></td><td></td><td></td><td>V2.0.2 support for hg38 liftover coming soon</td></tr>
<tr><td></td><td></td><td></td><td>Must be manually enabled</td></tr>
<tr><td>UniProt ID annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Other UniProt annotations (e.g. AAxform)</td><td bgcolor="#FF0000">No</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Custom fields: t_alt_count; t_ref_count; etc</td><td bgcolor="#00FF00">MAF Output Only</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>“other_transcripts” annotation</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>Reference context annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>COSMIC annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
<tr><td>UCSC ID annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td>In Funcotator UCSC ID is part of the HGNC data source.</td></tr>
<tr><td>RefSeq ID annotations</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#00FF00">Yes</td><td></td></tr>
</table>

<a name="1.5.2"></a>
#### 1.5.2 - Oncotator Bugs Compared With Funcotator
<table>
<tr><th></th><th style=>Fixed in Funcotator</th><th>Fixed in Oncotator</th><th>Notes</th></tr>
<tr><td>Collapsing ONP counts into one number</td><td bgcolor="#AAAAAA">N/A</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Variants resulting in protein changes that do not overlap the variant codon itself are not rendered properly</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Appris ranking not properly sorted</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Using protein-coding status of gene for sorting (instead of transcript)</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>De Novo Start in UTRs not properly annotated</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Protein changes for Frame-Shift Insertions on the Negative strand incorrectly rendered</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>MNP End positions incorrectly reported</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>MNPs on the Negative strand have incorrect cDNA/codon/protein changes</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>For Negative strand indels; cDNA string is incorrect</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Negative strand splice site detection boundary check for indels is incorrect</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Inconsistent number of bases in reported reference context annotation for indels</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>5’ Flanking variants are reported with an incorrect transcript chosen for Canonical mode</td><td bgcolor="#00FF00">Yes</td><td bgcolor="#FF0000">No</td><td></td></tr>
<tr><td>Variants overlapping both introns and exons or transcript boundaries are not rendered properly</td><td bgcolor="#FFFF00">No</td><td bgcolor="#FF0000">No</td><td>Funcotator produces a ‘CANNOT_DETERMINE’ variant classification and minimal populated annotations.</td></tr>
</table>

<a name="2"></a>
## 2 - Tutorial 

 <a name="2.0"></a>
### 2.0 - Requirements
1. Java 1.8
2. A functioning GATK4 jar
3. Reference genome (fasta files) with fai and dict files. Human references can be downloaded as part of the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle "GATK resource bundle").  Other references can be used but must be provided by the user.
4. A local copy of the [Funcotator data sources](#1.1 "Funcotator data sources")
5. A VCF file containing variants to annotate.

 <a name="2.1"></a>
### 2.1 - Running Funcotator in the GATK With Base Options

Open a command line and navigate to your GATK directory.
```
cd ~/gatk
```

At this point you should choose your output format.  There are [two output format choices](#1.3 "two output format choices"), one of which must be specified. 

Additionally, you must specify a reference version.  This reference version is used verbatim to determine which data sources to use for annotations.  That is, specifying `hg19` will cause Funcotator to look in the `<data_sources_dir>/hg19` folder for data sources to use.
 
A VCF instantiation of the Funcotator tool looks like this: 
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.vcf --output-file-format VCF
```

A MAF instantiation of the Funcotator tool looks like this: 
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF
```

 <a name="2.2"></a>
### 2.2 - Optional Parameters

 <a name="2.2.1"></a>
#### 2.2.1 - --ignore-filtered-variants
This flag controls whether Funcotator will annotate filtered variants.  By **default, this flag is set to true.**
To annotate filtered variants, run Funcotator with this flag set to false:

```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF --ignore-filtered-variants false
```

 <a name="2.2.2"></a>
#### 2.2.2 - --transcript-selection-mode
This parameter determines how the primary annotated transcript is determined.  The two modes for this parameter are **BEST_EFFECT**, **CANONICAL**, and **ALL**.  By **default, Funcotator uses the CANONICAL** transcript selection mode.

The explanations and rules governing the two transcript selection modes are as follows:

_**BEST_EFFECT**_
Select a transcript to be reported with details with priority on effect according to the folowing list of selection criteria:
* Choose the transcript that is on the custom list specified by the user. If no list was specified, treat as if no transcripts were on the list (tie).
* In case of tie, choose the transcript that yields the variant classification highest on the variant classification rank list (see below).
* If still a tie, choose the transcript with highest level of curation. Note that this means lower number is better for level (see below).
* If still a tie, choose the transcript with the best appris annotation (see below).
* If still a tie, choose the transcript with the longest transcript sequence length.
* If still a tie, choose the first transcript, alphabetically.

_**CANONICAL**_

Select a transcript to be reported with details with priority on canonical order according to the folowing list of selection criteria:

* Choose the transcript that is on the custom list specified by the user. If no list was specified, treat as if all transcripts were on the list (tie).
* In case of tie, choose the transcript with highest level of curation. Note that this means lower number is better for level (see below).
* If still a tie, choose the transcript that yields the variant classification highest on the variant classification rank list (see below).
* If still a tie, choose the transcript with the best appris annotation (see below).
* If still a tie, choose the transcript with the longest transcript sequence length.
* If still a tie, choose the first transcript, alphabetically.

_**ALL**_
Same as CANONICAL, but indicates that no transcripts should be dropped.  Render all overlapping transcripts.

 <a name="2.2.3"></a>
#### 2.2.3 - --transcript-list
This parameter will restrict the reported/annotated transcripts to only include those on the given list of transcript IDs.  This list can be given as the path to a file containing one transcript ID per line OR this parameter can be given multiple times each time specifying a  transcript ID.

_When specifying transcript IDs, **transcript version numbers will be ignored**._

Using a manually specified set of transcripts for the transcript list:
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF --transcript-list TRANSCRIPT_ID1 --transcript-list TRANSCRIPT_ID2
```

Using an equivalent transcript file:
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF --transcript-list transcriptFile.txt
```

Contents of `transcriptFile.txt`:
```
TRANSCRIPT_ID1
TRANSCRIPT_ID2
```

 <a name="2.2.4"></a>
#### 2.2.4 - --annotation-default
This parameter specifies a default value for an annotation.  This default value for this annotation will be used for any annotated variant.  However if this annotation would be added by Funcotator to this variant, the Funcotator value will overwrite this default.  

To specify this annotation default, the value on the command line takes the format:
`ANNOTATION_FIELD:value`

For example, to set the _Center_ annotation to _broad.mit.edu_:
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF --annotation-default Center:broad.mit.edu
```

It is valid to provide both the`--annotation-default` and `--annotation-override` arguments to Funcotator, however the behavior of specifying an annotation-default and an annotation-overrid for the _same annotation field_ is undefined.

 <a name="2.2.5"></a>
#### 2.2.5 - --annotation-override
This parameter specifies an override value for an annotation.  If the annotation were to be added to a variant by a data source, the value for that annotation would be replaced with the value specified in the annotation override.  If the annotation would not be added by a data source it is added to the output with the given value.

To specify this annotation default, the value on the command line takes the format:
`ANNOTATION_FIELD:value`

For example, to override the _NCBI_Build_ annotation to _HG19_:
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF --annotation-override NCBI_Build:HG19
```

It is valid to provide both the `--annotation-override` and `--annotation-default` arguments to Funcotator, however the behavior of specifying an annotation-override and an annotation-default for the _same annotation field_ is undefined.

 <a name="2.2.6"></a>
#### 2.2.6 - --allow-hg19-gencode-b37-contig-matching
This flag will cause hg19 contig names to match b37 contig names, allowing a set of variants created on an hg19 reference to match a b37 reference and visa-versa.

hg19 was created by UCSC. b37 was created by the Genome Reference Consortium.  In practice these references are _very_ similar but have small differences in certain bases, as well as a different naming convention for chromosomal contigs (`chr1` in hg19 vs `1` in b37).  In 99.9% of cases the results will be identical, however for certain genomic regions the results will differ.

This flag **defaults to _true_.**

To run Funcotator without this hg19/b37 matching: 
```
./gatk Funcotator --variant variants.vcf --reference Homo_sapiens_assembly19.fasta --ref-version hg19 --data-sources-path funcotator_dataSources.v1.2.20180329 --output variants.funcotated.maf --output-file-format MAF --allow-hg19-gencode-b37-contig-matching false
```
<a name="3"></a>
## 3 - FAQ
### Why do I not get annotations from my favorite data source on my favorite variant?
This almost always happens when the data source does not overlap the variant.  Commonly a variant that is not within a gene will not be annotated by data sources because they are not in the region that the data sources cover (e.g. when the `VariantClassification` is `IGR`, `FIVE_PRIME_FLANK`, `COULD_NOT_DETERMINE`, etc.).
This can also happen if the given reference file does not match the data sources' reference (for the pre-packaged data sources either `hg19`/`b37` or `hg38`).  In this case, Funcotator will produce a large obnoxious warning:

```
  _ _ _  __        ___    ____  _   _ ___ _   _  ____   _ _ _ 
 | | | | \ \      / / \  |  _ \| \ | |_ _| \ | |/ ___| | | | |
 | | | |  \ \ /\ / / _ \ | |_) |  \| || ||  \| | |  _  | | | |
 |_|_|_|   \ V  V / ___ \|  _ <| |\  || || |\  | |_| | |_|_|_|
 (_|_|_)    \_/\_/_/   \_\_| \_\_| \_|___|_| \_|\____| (_|_|_)
--------------------------------------------------------------------------------
Only IGRs were produced for this dataset.  This STRONGLY indicates that this   
run was misconfigured.     
You MUST check your data sources to make sure they are correct for these data.
================================================================================
```

<a name="4"></a>
## 4 - Known Issues
The current list of known open issues can be found on the GATK github page [here](https://github.com/broadinstitute/gatk/labels/Funcotator).

<a name="5"></a>
## 5 - Github
Funcotator is developed as part of GATK.  The GATK github page is [here](https://github.com/broadinstitute/gatk/).

<a name="6"></a>
## 6 - Tool Documentation
Tool documentation is written in the source code for Funcotator to better explain the options for running and some details of its features.
The tool documentation for Funcotator is [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php).

----

[Back to Top](#0)

<hr>
