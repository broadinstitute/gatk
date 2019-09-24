################################################################################
# Funcotator Data Sources Package README
################################################################################

+---------------------------------------------+ 
| Data Source Version Information             |
+---------------------------------------------+ 

Version:          v1.6.20200227_ecoli
Use Case:         Ecoli
Source:           GATK Github
Alternate Source: N/A

################################################################################

+---------------------------------------------+ 
| README                                      | 
+---------------------------------------------+ 

This is a collection of data sources to be used in conjunction with Funcotator
to annotate Ecoli data (e. coli K12 mg1655 - ASM548v2).

This folder is a top-level Data Sources Folder for The Broad Institute's 
Funcotator tool.  When running Funcotator, pass the path to this directory in
as a command-line argument:

  ./gatk Funcotator --data-sources-path PATH/TO/THIS/FOLDER ...

For more information on Funcotator, see the Funcotator tool doc or tutorial:

  https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php
  https://gatkforums.broadinstitute.org/dsde/discussion/11193/funcotator-information-and-tutorial/

For more information on GATK, see the GATK development github site:

  https://github.com/broadinstitute/gatk

################################################################################

+---------------------------------------------+ 
| Data Sources                                |
+---------------------------------------------+ 

Using this Data Sources Folder will enable the following data sources:

gencode / Transcript Data Source
----------------------
This is not actually data from gencode, but is currently in the gencode folder because that is what the code looks for.
This will change very soon.
The transcript information included here is in GTF format and was found on the Ensembl website here:

Reference Genome:
  ftp://ftp.ensemblgenomes.org/pub/bacteria/release-44/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz

GTF:
	ftp://ftp.ensemblgenomes.org/pub/bacteria/release-44/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.44.gtf.gz

CDS:
  ftp://ftp.ensemblgenomes.org/pub/bacteria/release-44/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz


