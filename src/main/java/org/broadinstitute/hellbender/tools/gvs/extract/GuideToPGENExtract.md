# Extracting PGEN from GVS

## The PGEN format
PGEN is a format written for and used by version 2 of [PLINK](https://www.cog-genomics.org/plink/2.0/).  ***IT IS VERY IMPORTANT TO NOTE THAT VERSION 2 OF PLINK IS STILL IN ALPHA AND THE PGEN FORMAT IS STILL SUBJECT TO CHANGE.***  The format comprises 3 file types (or actually sometimes 4):
1. A `.pgen` file.  This is a binary file that stores a mapping of samples and sites to genotypes in a very cleverly compressed way that I can't explain super well because it's complicated.
2. A `.pvar` file. This is essentially a sites-only VCF.  It has information for each site referenced in the `.pgen` file.  PLINK also provides an option to produce/use a zstd compressed version of the file (`.pvar.zst`), and we have opted to write that for performance purposes.
3. A `.psam` file.  This is a plaintext file that contains a list of samples referenced in the `.pgen` file.  It also optionally includes some phenotype data.
4. Optionally, a `.pgi` file.  Typically, a `.pgen` file has an index at the top.  Optionally, PLINK supports using a `.pgen` file with an index in a separate `.pgi` file.

The PGEN format does not store all of the information that a VCF has.  It leaves out a lot of the fields and annotations you can store in a VCF.  As a result of this and the clever compression in the `.pgen` file, these files are typically much smaller than equivalent VCFs.

For more information on the PGEN file format, see the official spec [here](https://github.com/chrchang/plink-ng/blob/master/pgen_spec/pgen_spec.pdf).

## The code
The code for the PGEN extract can be divided into 3 parts:
1. The PGEN-JNI, a C++/JNI library that handles writing HTSJDK VariantContext objects to PGEN files,
2. ExtractCohortToPgen, a GATK tool based on ExtractCohortToVcf that processes VariantContexts and passes them to PGEN-JNI for writing, and
3. GvsExtractCallsetPgenMerged, a WDL workflow based on GvsExtractCallset that uses ExtractCohortToPgen to write a series of PGEN files and then merges them by chromosome.

### Part 1: PGEN-JNI
The PGEN-JNI library was written by Chris Norman of the GATK Engine Team and lives [here](https://github.com/broadinstitute/pgen-jni).  It is written primarily in C++ for performance purposes and also compatibility with the pgenlib library (part of the [plink repo](https://github.com/chrchang/plink-ng/tree/master)).  It builds on top of pgenlib to provide a writer for creating PGEN files and writing to them from HTSJDK VariantContext objects.

PGEN-JNI is compatible with Linux and macOS.

A build of this library is currently hosted on the Broad's artifactory repo, and that is being used as a dependency for GATK.

Ownership of the PGEN-JNI library will stay with the GATK Engine Team and we will provide support for it if y'all encounter any issues with it.

### Part 2: ExtractCohortToPgen
ExtractCohortToPgen is a GATK tool that inherits from ExtractCohort and is based very closely on ExtractCohortToVcf.  It produces 3-4 files:

1. A `.pgen` file, which contains a mapping of samples and sites to variants,
2. A `.psam` file, which contains a list of sample names,
3. A `.pvar.zst` file, which is a zstd compressed list of sites with alleles, similar to a sites-only VCF, and
4. Optionally (if specified by setting `write-mode` to `WRITE_SEPARATE_INDEX`), a `.pgi` file, which contains an index for the `.pgen` file.

It has a few arguments that are specific to it that warrant explanation.

#### pgen-chromosome-code
Plink defines a set of [chromosome codes](https://www.cog-genomics.org/plink/2.0/data#irreg_output) that correspond to different sets of contig names of chromosomes.  This tool supports two of those chromosome code options: `chrM` and `MT`, which correspond to hg38 and hg19 contig naming, respectively.  This argument is required.

#### write-mode
The PGEN writer defined in PGEN-JNI defines two different write modes.  When `WRITE_AND_COPY` is selected, a temporary .pgen file is created and written to during the running of the tool, and then once all records have been written, a new file is created with the index at the top and the contents of the temporary .pgen file appended to it.  When `WRITE_SEPARATE_INDEX` is selected, the index is instead written to a separate .pgi file.  The default is `WRITE_AND_COPY`.

#### max-alt-alleles
The PGEN format can only support up to 254 alt alleles per site.  This argument allows you to specify a limit.  The default is the max of 254.  Any sites with more alt alleles than the specified max will not be written.

#### lenient-ploidy-validation
PGEN is a bit quirky in that it requires samples to be diploid but has a special case for sex chromosomes, which are allowed to be haploid.  By default, any attempt to write a record with an unsupported ploidy will result in an exception being thrown.  If this flag is used, then ploidy failures will instead be logged and the records will be written as missing.

#### writer-log-file
The C++ code in the PGEN writer in PGEN-JNI will log sites that exceed max-alt-alleles and with unsupported ploidy (if lenient-ploidy-validation is set) to the specified log file, if this argument is set.

#### allow-empty-pgen
Empty PGEN files are not technically valid PGEN files.  However, for parallel processing purposes, it is sometimes helpful to allow the creation of empty files when there are no variants to be written.  The GvsExtractCallsetPgenMerged workflow relies on this.

If this flag is set and no variants are written, an empty .pgen, .psam, and .pvar.zst file will be written in `onShutdown()`.

By default (i.e. if this flag is not set), if there are no variants written, an exception will be thrown.  

### Part 3: GvsExtractCallsetPgenMerged
GvsExtractCallsetPgenMerged is a WDL workflow that calls ExtractCohortToPgen to extract data from GVS and write it to PGEN files, and then merges those PGEN files by chromosome.  This workflow has 3 steps:

#### Step 1: GvsExtractCallsetPgen
This is a workflow based very closely on the GvsExtractCallset workflow (which is used for extracting data from GVS and writing it to VCF files).  It has a few notable differences when compared to that workflow.

First, it features extra inputs for some of the PGEN-specific arguments described in [Part 2](#part-2-extractcohorttopgen).  These mostly have defaults which match the defaults for the ExtractCohortToPgen tool.  A default of "chrM" is also specified for chromosome code, to match the reference defaulting to hg38.

Second, the call to SplitIntervals is replaced with a call to a new task called SplitIntervalsTarred.  For very large scatter counts (I think >20k-ish?), the SplitIntervals task will fail because the shell command being run behind the scenes by the WDL `glob()` function will fail when trying to do wildcard expansion with that many files.  The best solution I was able to find for this was to just tar up all the interval files and use that as an output.

Third, as a result of this, the ExtractTask has some changes when compared to the VCF ExtractTask (beyond the fact that it now calls the PGEN extract tool instead of the VCF one).  It now accepts the interval lists tar file as an input, along with the filename for the interval list file for that particular shard.  Before it runs the extract tool, it extracts the specified interval list file from the interval lists tar.

The output of this workflow is a list of .pgen, .psam, and .pvar.zst files, along with the interval list tar and a list of the filenames within the tar (ordered to match the pgen/psam/pvar files so they can be matched up to each other more easily).

#### Step 2: SplitFilesByChromosome
This is a short task defined within GvsExtractCallsetPgenMerged.  It loops through the interval list files in the interval list tar and extracts the contig name from each.  Then, going by that contig name, it writes the names of the corresponding .pgen, .psam, and .pvar.zst files to list files named for the corresponding contigs.  

The output of this task is 3 arrays of files, containing names for .pgen, .psam, and .pvar.zst files which contain sites for each contig.  For example, the file `chr1.pgen_list` would contain a list of all the .pgen files that contain variants for sites on chromosome 1.

#### Step 3: MergePgenHierarchical
This workflow accepts a list of .pgen, .psam, and .pvar.zst files and merges them all into one file using Plink's `--pmerge-list` functionality.

For performance purposes, the merging is done in two stages.  One big, monolithic merge takes too long.

First, the file lists are sorted by index (the index in the filename).  Plink does not like (i.e. does not support) merging files with overlapping intervals, so the files need to be sorted so that when they are merged in stages, we do not create merged files with intervals that overlap.

Next, the lists are split, so we can parallelize the merge.  Each of the split lists is then sent to a task that merges the files using Plink.

Then, this is repeated to merge all the resulting files into one.  The output is a single .pgen file, .psam file, and .pvar.zst file containing all the extracted data for the corresponding chromosome.

## Testing
All of my testing (with the exception of some very small scale stuff early on) has been done in the [GVS_AoU_PGEN_Extract_Development Terra workspace](https://app.terra.bio/#workspaces/allofus-drc-prod-auxiliary/GVS_AoU_PGEN_Extract_Development) and using data from the GVS Delta callset (aou-genomics-curation-prod.aou_wgs_fullref_v2).  My test process (for the majority of tests) has been as follows:

1. Select a list of sample names from aou_wgs_fullref_v2.sample_info (excluding control samples and sample 3224672 because of the data issue mentioned [here](https://broadinstitute.slack.com/archives/CJRLP6ZSA/p1699026273329339)).
2. Use that list as an input to GvsPrepareRangesCallset to create a cohort of test data in a separate BigQuery dataset (aou-genomics-curation-prod.klydon_pgen_extract_test).
3. Run GvsExtractCallsetPgen on the newly created cohort.  (I would just run GvsExtractCallsetPgenMerged, but I like using Workflow Dashboard to monitor how the job is going and dig into it if there are any failures.  Workflow Dashboard doesn't seem to let you dig into individual tasks for workflows with sub-workflows, so it wouldn't allow me to look at individual shards running ExtractTask if I ran GvsExtractCallsetPgenMerged.  Job Manager would be an alternative for this, but it seems to be pretty much unusable for even moderately-sized jobs.)
4. Run GvsExtractCallset on the newly created cohort, making sure to use the same parameters, including scatter count.  This will generate VCF files that we can use to compare to the PGEN files created during the previous step for validation.
5. Run GvsExtractCallsetPgenMerged with the same parameters used to run GvsExtractCallsetPgen in Step 3.  This will use call-caching for the extract steps and then merge the PGEN files by chromosome.  (Running it this way is maybe not the ideal way to do this, but it's what I've been doing for reasons described in the parenthetical in Step 3).
6. Create list files, by file type, containing the gs:// URIs for the .pgen, .psam, and .pvar.zst files created in Step 3, along with the .vcf.gz files created in Step 4.  Upload them to the workspace to use for validation.
7. Run ComparePgenAndVcfScatter using the file lists as inputs.  If there are any differences, it will output files that contain those differences.  If there are no diff files generated, the files match.

ComparePgenAndVcfScatter is a workflow I wrote that converts a list of .pgen, .psam, and .pvar.zst files generated by GvsExtractCallsetPgen into .vcf.gz files and then compares those files to a list of .vcf.gz files generated by GvsExtractCallset.  

It ignores basically everything except genotypes, because PGENs do not store all the other fields and annotations that the VCFs might have.  It will also skip over any sites in the VCFs with >254 alleles because those will not be present in the PGEN files.  Any differences are written to diff files, in the form of the differing lines in the VCFs being compared.

The code for this comparison tool lives [here](https://github.com/KevinCLydon/pgen_vcf_comparator) in a repo I created under my GitHub account.  (I didn't create it under the Broad org because it's sort of half-baked and bad and not actually meant to be used by anyone other than me.)  The WDLs are copied over to the GATK repo in my working branch for now (in `scripts/variantstore/wdl/pgen_validation`), because that's an easier way to get them into Terra.  I don't know if y'all want to continue using this tool, but I'm happy to discuss it more if it would actually be useful to you.