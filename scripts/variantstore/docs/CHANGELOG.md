# Genomic Variant Store (GVS) Changelog

## 0.4.0 - 2023-??-??

### Added

- `call_PS` column added in variant data tables for phasing set support. This change makes GVS 0.4.0 incompatible with ingesting sample data into GVS schemas created by earlier versions of GVS.

### Changed

- Updated to latest version of VETS variant calling toolchain.

### Fixed

- Support for quota-restricted `GvsBeta` users was accidentally removed in a prior release and has now been restored.

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
