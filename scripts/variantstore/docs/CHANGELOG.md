# Genomic Variant Store (GVS) Changelog

## 0.3.0 - 2023-09-01

### Changed

- Made VETS tool chain the default method for filtering (used to be VQSR).
- Synced GVS development branch with GATK master, GATK-based Java tools now run on Java 17 rather than Java 8.

## 0.2.2 - 2023-08-31

### Changed

- If no `extract_output_gcs_dir` is specified for `GvsBeta`/`GvsJointVariantCalling`, choose a sensible default value so output VCFs are collected in one convenient place.

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
