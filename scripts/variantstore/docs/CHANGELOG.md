# Genomic Variant Store (GVS) Changelog

## 0.6.1 - 2024-11-18

### Changed

- Maximum alternate alleles for extract now defaults to 1000.
- Ingest parallelism now specified by scatter width rather than batch size.

### Fixed

- Cohort extract now chooses a more appropriate scatter width based on the number of samples loaded rather than the size of the cohort being extracted.
- Ingest now correctly handles ploidy differences in Dragen 3.7.8 VCFs.
- Corrected handling of query project vs extract project in cohort extract.

## 0.6.0 - 2024-07-25

### Added

- Blended Genome Exome (BGE) support added.
- `OUTSIDE_OF_TARGETS` genotype filter support added for Exomes and BGEs for a specified target interval list.

### Changed

- As a consequence of the `OUTSIDE_OF_TARGETS` change for Exomes and BGEs described above, the VCF header line for *all* GVS-produced VCFs (including WGS) has changed slightly. In particular the `Number` attribute of the `FORMAT` record for the `FT` genotype-level filter now has a value of `.` rather than `1`, and the `Description` attribute for it now has a value of `Genotype-level filter`.

### Fixed

- Fixed an edge case where untimely VM preemption could cause sample load status to not be recorded correctly.

## 0.5.5 - 2024-06-24

### Changed

- Inputs `vcf_files_column_name`, `vcf_index_files_column_name`, `sample_id_column_name`, and `extract_output_gcs_dir` are no longer optional.
- high_CALIBRATION_SENSITIVITY_SNP and high_CALIBRATION_SENSITIVITY_INDEL were moved from the ##FILTER entry in the header to a ##COMMENT. No change to behavior or vcf content outside of header. 
- Reduced the number of sharded vcfs coming out of beta workflow for smaller callsets. Documentation on details coming soon. 
- Bug fix to correctly handle samples with chromosomes with differing ploidy for DRAGEN 3.7.8 data.


## 0.5.4 - 2024-03-18

### Changed

VCF files generated using the Extract workflow (`GvsExtractCallset`) no longer contain the INFO field `AS_YNG`. AS_YNG contained the allele-specific YNG field (Yay/Nay/Gray). Yay meant that allele is known good (in the VQSR/VETS training set), Nay is known bad and Gray means it's not in the training set, so use the VQSR/VETS score to determine whether to filter the genotype. Although `AS_YNG` is no longer included in the VCF explicitly, it is implicitly used to filter genotypes.

## 0.5.3 - 2024-02-20

### Changed

Updated the representation of sample load states in the `sample_load_status` table from using `STARTED` and `FINISHED`
to using `HEADERS_LOADED`, `REFERENCES_LOADED` and `VARIANTS_LOADED`. The new status reading/writing code is backward
compatible with the old `STARTED` and `FINISHED` states if those have already been used within a callset. The new statuses
should allow for faster recovery if an ingest task or workflow is preempted or aborted and then restarted.

## 0.5.2 - 2024-01-17

### Changed

- Updated release process to version `GvsExtractCohortFromSampleNames` in addition to `GvsBeta`.

### Fixed

- `OutOfMemoryExceptions` conditions with default settings in `ExtractTask` addressed, memory override now passed through correctly.

## 0.5.1 - 2024-01-04

### Changed

- Modified GvsPrepareRangesCallset.wdl and GvsExtractCohortFromSampleNames.wdl to make the tables created as part of preparing the callset to *by default* have a time-to-live of 14 days.

## 0.5.0 - 2023-10-26

### Changed

- Modified loading of ref_ranges_* tables to now use GQ0 band (state = 0) instead of Missing (state = 'm') for regions where the interval list chosen for importing data extends beyond regions covered in the gVCF data. NOTE that when importing data with a drop state of 0, then those regions will not be written at all to ref_ranges_*, since we are dropping the state of 0.

## 0.4.1 - 2023-10-16

### Fixed

- Updated input VCF and VCF index localization in ingest to support Google Cloud Storage composite objects.

## 0.4.0 - 2023-10-10

### Added

- `call_PS` column added in variant data tables for phasing set support. This change makes GVS 0.4.0 incompatible with ingesting sample data into GVS schemas created by earlier versions of GVS.

### Changed

- Updated to latest version of VETS variant calling toolchain.

### Fixed

- Support for quota-restricted `GvsBeta` users was accidentally removed in a prior release and has now been restored.
- Fix bug with `GenerateImportFofnFromDataTable` task inappropriately not being marked `volatile`.

## 0.3.2 - 2023-09-11

### Changed

- Added a new Boolean `is_wgs` flag to `GvsBeta`/`GvsJointVariantCalling` that determines whether to use the WGS interval list (and weights bed) if true, or the Exome interval list if false.
  - The (optional) inputs `interval_list` and `interval_weights_bed` can be defined to override the default behavior.

## 0.3.1 - 2023-09-06

### Changed

- Added an optional input `billing_project_id` to allow for ingest from a GCS bucket with requester pays turned on where egress will be charged to the Google billing project ID passed.

## 0.3.0 - 2023-09-01

### Changed

- Made VETS tool chain the default method for filtering (used to be VQSR).
- Synced GVS development branch with GATK master, GATK-based Java tools now run on Java 17 rather than Java 8.
- If no `extract_output_gcs_dir` is specified for `GvsBeta`/`GvsJointVariantCalling`, default value will be chosen so output VCFs are collected in one convenient place. Please see the [GVS documentation for further details](../beta_docs/gvs-overview.md#input-descriptions) for further details.

## 0.2.1 - 2023-08-30

### Changed

- Updates to Beta documentation.

## 0.2.0 - 2023-08-29

### Added

- [Bulk Ingest Support](./gvs-bulk-ingest-details.md)
- [VETS variant calling support](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/release_notes/VETS_Release.pdf)
- [Exome support](./RUNNING_EXOMES_ON_GVS.md)
- Enhanced provenance for GVS workflow versions (symbolic and hashes)
- Enhanced standardization of GVS Docker image usage


## 0.1.0 - 2023-08-03

### Added

- First semantically versioned Genomic Variant Store release.
