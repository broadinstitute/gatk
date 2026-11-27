# Overview of the Genomic Variant Store (GVS)

The Genomic Variant Store (GVS) is a scalable cloud-based platform designed to
efficiently store, manage, and process large-scale genomic variant data. GVS
enables joint calling across hundreds of thousands of samples and provides fast
cohort extraction capabilities for population genomics studies. For a high level
overview see `gvs-product-sheet.pdf`.

GVS was created as a branch of the Genome Analysis Toolkit (GATK). The principal
GVS branch in git is called `ah_var_store`. While this GVS branch has been
updated from GATK `master` periodically, it is not clear if GVS will ever be
merged back into the main GATK branch. While GVS does build upon some of the
GATK codebase, the majority of the code in the GATK repository (even on the
`ah_var_store` branch) is not related to GVS.

## Key Features

- **Scalable Storage**: Leverages Google BigQuery for petabyte-scale genomic
  data storage with high availability
- **Joint Calling**: Performs population-scale joint genotyping across large
  cohorts
- **Cohort Extraction**: Rapid extraction of variant subsets for specific sample
  cohorts
- **Cloud-Native**: Designed for Terra platform with WDL-based workflows

## Usage Scenarios

GVS supports two primary usage scenarios as described below:

### GVS Beta Callsets

GVS Beta is intended for end users to be able to run their own joint calling
using the `GvsJointVariantCalling.wdl` as a top-level workflow. See the
documentation in the `beta_docs` directory for more information on GVS Beta
usage.

### AoU Callsets

GVS is used to produce joint callsets for the All of Us Research Program at the
scale of 400,000+ samples. AoU callsets push the scale limits of GVS and are not
suited for the `GvsJointVariantCalling.wdl` workflow. For AoU, more narrowly
scoped GVS workflows are run by the engineers of the Variants team following an
AoU-specific protocol as described in `AOU_DELIVERABLES.md`. Rather than the
collection of VCF files produced from GVS Beta, the core deliverable of an AoU
callset is a Hail VDS (Variant Dataset). More AoU-specific documentation can be
found under `scripts/variantstore/docs/aou`.

Conceptually an AoU callset is similar to a GVS Beta callset up to the point of
extraction. While a Beta callset would run `GvsPrepareCallset.wdl` and then
`GvsExtractCallset.wdl` to extract a set of VCF files, an AoU callset will run
`GvsExtractAvroFilesForHail.wdl` and then `GvsCreateVDS.wdl` to generate a Hail
VDS. The data within the VDS is used to generate an ancestry file for creating
the VAT (Variant Annotations Table) and interval lists for the "small callsets"
described in `AOU_DELIVERABLES.md`.

### AoU artifact version numbers

AoU artifact names carry a three-level version, as in
`aou_srwgs_short_variants_v9_r2_p1.vds`: `v9` is the major version, which is
AoU's call; `r2` is the release, which advances when the participant set
changes; and `p1` is the patch. Within a release the patch counter is **one
sequence shared across every artifact type**, not a separate counter per
artifact. In v9 r2, p1 is the VDS, p2 the VAT, and p3 and p4 the participant
mapping tables (VS-2023).

So artifacts of one release routinely carry different patch numbers, and that
is not evidence that they are mismatched: the v9 r2 p4 mapping table was built
from the v9 r2 p1 VDS, and pairing the two is correct. Before reporting a
version mismatch between artifacts, check the release level first, then
confirm lineage from the ticket that produced the artifact, or from the inputs
recorded in the Cromwell metadata of the workflow that built it. A differing
patch number alone should not be flagged.

The same reasoning runs the other way: an unchanged release level means an
unchanged participant set, whatever the steps in between did. Rebuilding a VDS
within a release can remove samples and later add them back. In VS-2023,
`GvsRemoveSamplesFromVDS` took two superpartitions' samples out of the r2 p1
VDS, and `GvsMergeAndRescoreVDSes` merged them back in from a rebuilt VDS, so
the r2 p5 VDS holds the same samples as r2 p1. So do not infer a changed sample
set from one step's name. Check the release level, or the net effect of the
whole procedure.

### Set provenance-bearing inputs explicitly on AoU runs

For an AoU deliverable, set explicitly in the Terra configuration any input
whose default decides what the artifact contains, even when the default is the
value you want. The configuration and the Cromwell metadata then record what
was used. Left to a default, the answer depends on which WDL version happened
to run, and a later change to the default silently changes what a rerun
produces. For example, `GvsCreateVATfromVDS`'s `manual_clinvar_path_prefix`
decides which ClinVar build annotates the VAT (VS-2052).

`git_branch_or_tag` belongs in the same category, but it does not do what its
name suggests. It does not choose the code that runs. Terra runs whichever
Dockstore version is selected in the configuration. The input is passed to
`GetToolVersions`, which clones that branch and records its HEAD as the run's
`git_hash`. So set it to the release branch *and* select that branch's version
in Terra. Set it to one while running the other, and the metadata names a commit
that is not the code that ran. Selecting the version requires the workflow to be
registered for the branch in `.dockstore.yml`. Release branches are
usually registered only for the workflows they were cut to run, so check before
the run is due.

## Key Concepts

### Hard filtered versus soft filtered variants

During the GVS joint calling process variants are evaluated for quality at both
the site and allele level. By default GVS evaluates variant quality scores using
a machine learning model called VETS (Variant Extract-Train-Score). GVS formerly
used a machine learning model called VQSR (Variant Quality Score Recalibration)
for this purpose. While still supported in GVS, VQSR is rarely used in making
GVS callsets today. See
`scripts/variantstore/docs/release_notes/VETS_Release.pdf` for more information
on VETS and the transition away from VQSR.

GVS stores all variants and their associated quality scores in its BigQuery
schema, regardless of whether these variants pass or fail quality filters. GVS
artifacts (VCF, PGEN, VDS, etc), are classified as being either "hard filtered"
or "soft filtered" depending on how they treat variants with failing quality
scores. "Hard filtered" artifacts completely exclude variants with failing
quality scores, while "soft filtered" artifacts include all variants along with
their quality scores, leaving it to the user's discretion how to handle these
variants in downstream analysis.

### Split versus unsplit multi-allelics

A multi-allelic site is a genomic location where more than one alternate allele
is observed. For example, at a given position the reference allele might be "A",
and there could be two alternate alleles "C" and "G". In an unsplit
multi-allelic representation, this would be represented as a single conceptual
row with the alternate alleles listed together (e.g., A -> C,G). In a split
multi-allelic representation, this would be represented as two separate rows
(e.g., A -> C and A -> G).

GVS contains examples of both split and unsplit multi-allelic representations.
Knowing which representation is used where is critical for a correct
understanding of the data and for writing code that interacts with GVS. For
example, the `alt_allele` table uses a split multi-allelic representation to
support efficient querying and filtering, while the `vet_%` tables use an
unsplit multi-allelic representation like the input gVCFs from which they are
derived. While Hail VDSes in general can use either split or unsplit
multi-allelic representations, the VDSes created, read, and modified by GVS are
exclusively unsplit.

Some of the GVS code in this repo does interact with split multi-allelic non-VDS
Hail artifacts; see `AOU_MARKETING_STATISTICS.md` for an example that takes a
split Exome MatrixTable as input.

### Allele representations across GVS artifacts

The same variant is written differently in different GVS artifacts. The
differences are systematic rather than incidental, and comparing alleles across
artifacts without accounting for them is a recurring source of bugs.

- **`vet_%` and `alt_allele`** hold minimized (trimmed) alleles that are *not*
  left-aligned. An indel stays wherever the caller placed it.
- **The VDS** holds unsplit multi-allelic rows in a common per-site reference
  context. `alleles[0]` must be long enough to span the longest deletion at the
  site, and every ALT is expressed relative to it. A VDS allele string is
  therefore a property of the site's allele set, not of the variant: the same
  variant is written differently at a site with different neighbours.
- **The VAT, and therefore vids,** hold alleles that are minimized *and*
  left-aligned, via `bcftools norm -f` in `GvsCreateVATFromVDS.wdl:781`.

Five consequences worth knowing before writing anything that joins across these:

1. Minimize before comparing. `hl.min_rep`, or trimming by hand, reduces a VDS
   allele to the intrinsic form that `alt_allele` stores.
2. Compare the REF/ALT **pair**, never the ALT alone. One string can be the raw
   ALT of one variant and the minimized ALT of an unrelated one. At
   chr13:32368001 in Foxtrot, `CTT` is both: the raw ALT of a 24bp deletion
   against a 27bp REF, and the minimized ALT of a 2bp insertion.
3. Minimizing is not left-aligning. `hl.min_rep` trims shared prefixes and
   suffixes but will not walk an indel through a repeat, so it cannot unify a
   left-aligned vid with a non-left-aligned `alt_allele` row. That needs real
   normalization against the reference. `GvsCreateParticipantMappingFromVDS`
   does it by running `bcftools norm` over the VDS keys that are not
   left-aligned (the `NormalizeKeys` task).
4. Indel length is invariant under normalization, so it is a safe filter when
   hunting for equivalent representations; position and allele strings are not.
   The retired synonym search keyed on it as bcftools `ILEN` over a 200bp
   window rightward from the vid's position, left-aligned being the leftmost
   equivalent form.
5. An indel can also be stored without its anchor base: a deletion whose ALT
   is REF's *last* base, not its first. Minimizing leaves such a key alone,
   and its position is not the variant's, because `bcftools norm` re-anchors
   it to the left. The AnVIL 3K VDS has 279 of them, all deletions of more
   than 200 bases, which moved between 1 and 442 bases (VS-2054). `key_class`
   in `hail_create_participant_mapping_inputs.py` calls the class
   `unanchored`, and that script fails on it by design, so
   `GvsCreateParticipantMappingFromVDS` would fail on AnVIL 3K. It ran
   genome-wide on Foxtrot, which therefore has none. Anything that maps VDS
   keys to vids has to normalize this class along with the shiftable ones.
   Treating them as ordinary keys put 179 vids at the wrong locus in the
   VS-2054 prototype, and only a comparison against the VAT's AN showed it.


### A gap in `ref_ranges` is not a no-call

AoU datasets drop GQ0 reference blocks at ingest (`drop_state` ZERO), so it is
tempting to read "no reference block at this location" as "no-call". A sample
lacks a reference block at a location for three reasons, and only the first
fits that reading:

1. a GQ0 block was dropped there;
2. the sample has a `vet` row there, which is a call unless it has GQ 0 or
   fails the filter;
3. the sample's own upstream deletion spans the location. No reference block
   is written under the deletion, so the sample is a no-call there in the VDS
   and in the VCF extract.

AN in the VAT is `hl.agg.call_stats` over the densified VDS
(`hail_create_vat_inputs.py`), after sites are hard filtered and FT-failing
genotypes are set to no-call. Per sample at a site, the first match applies:

| Sample at the site                           | Contribution to AN                  |
|----------------------------------------------|-------------------------------------|
| Variant call with GQ 0 (`import_gvs.py:262`) | 0                                   |
| Variant call failing FT                      | 0                                   |
| Any other variant call                       | ploidy of the stored genotype       |
| No variant call, a reference block covers it | ploidy from the ploidy table (2, 1) |
| Neither                                      | 0                                   |

So SQL over `ref_ranges_%` alone cannot reproduce the VAT's AN, and SQL over
`ref_ranges_%` and `vet_%` would also need the per-genotype FT rule. The
tables and a given VDS can differ as well: a VDS with dropouts has no-calls
that the tables do not. Derive per-sample call status from the VDS the VAT was
built from (VS-2054).


# High-Level Architecture (WDLs, BigQuery, Terra)

## Java Artifacts

GVS is written in a mixture of Java, Python, WDL (Workflow Description
Language), and shell scripts. The Java code can be built from the repository
root using the `./gradlew clean shadowJar` command. This compilation requires a
JDK of at least version 17. For local compilation, if
[SDKMAN!](https://sdkman.io/) is available, one of the Java 17 versions can be
used (most GATK builds seem to use Temurin).

## Docker Versioning

The `GvsUtils.wdl` workflow's `GetToolVersions` task is used to retrieve the
version of various Docker images used by GVS. This task is usually run at the
start of a GVS workflow to ensure that the correct versions of tools are being
used.

## Reference Versioning

The `GvsUtils.wdl` workflow's `GetReference` task is used to retrieve the
canonical paths to references used by GVS. This task is usually run at the start
of a GVS workflow to ensure that the correct references are being used.

## Docker artifacts

GVS uses multiple Docker images to run its workflows. Some of these are specific
to GVS as documented below.

### GATK Docker Image

This builds a GATK Docker image based on the GVS branch of the GATK codebase.
Instructions for building this image can be found in the
`scripts/variantstore/docs/Build Docker from VM/building_gatk_docker_on_a_vm.md`
file. An x86 cloud-based VM is used to build this x86 Docker image, which is
then pushed to the Broad's GAR (Google Artifact Registry) repository.

### Variants Docker Image

This is a Google Cloud SDK Alpine-based Docker image to which many useful tools,
Python packages, and scripts have been added. This can be built on a developer's
local machine using the `scripts/variantstore/scripts/build_docker.sh` script.

The Python scripts under `scripts/variantstore/scripts/` are baked into this image at `/app/`, and
WDL tasks that set `docker: variants_docker` invoke them as `python3 /app/<script>.py`. **Editing or
adding such a script has no effect on a running workflow until the image is rebuilt and the new tag
is wired in** — a common trap (the WDL keeps using the old image). So whenever you change anything
baked into the image (a script, an added tool, or a Python package):

1. Run `scripts/variantstore/scripts/build_docker.sh`. It builds the alpine-based image, runs the
   unit tests, and pushes it to GAR with a tag of the form `<ISO-8601 date>-alpine-<12-char image
   id>`, e.g. `us-central1-docker.pkg.dev/broad-dsde-methods/gvs/variants:2026-06-30-alpine-766e78c05289`.
2. Update the `variants_docker` value in the `GetToolVersions` task of `GvsUtils.wdl` to that new
   tag. Workflows resolve the image through `GetToolVersions`, so this step is what makes the rebuilt
   image actually take effect.

To run one of these scripts' unit tests locally before building, without pytest
(it is not installed on developer machines), use unittest discovery from the
scripts directory:

    cd scripts/variantstore/scripts
    python3 -m unittest discover -s test -p test_run_in_hail_cluster.py

The dotted-module form, `python3 -m unittest test.test_run_in_hail_cluster`,
fails with `ModuleNotFoundError`, because `test/` is not a package and the
tests import their script as a top-level module.

Run across the whole `test/` directory, discovery reports import errors for
`test_create_*_bqloadjson_*` (no `ijson`) and `test_parquet_loading` (no
`google-cloud-bigquery`) on a machine without those packages. CI installs them, so
these errors locally do not indicate a problem with the change.

### Nirvana Docker Image

This is an Illumina Nirvana-based Docker image build from
`scripts/variantstore/variant-annotations-table/build_docker.sh`. Since Nirvana
is no longer being updated this image hopefully will also not need to be updated
in the future.

## "location" versus "position"

GVS leverages BigQuery partitioning and clustering to optimize query performance
and cost. The core tables are either partitioned by sample ID or location,
depending on the access pattern for which they are optimized. GVS uses the term
"location" in a very specific way: the chromosome number is multiplied by
1,000,000,000,000 and then the position on the chromosome is added. For example,
chromosome 1 position 1000 becomes `1000000000000 + 1000 = 1000000001000`. GVS
locations use a chromosome number of 23 for the X chromosome and 24 for the Y
chromosome. This location encoding allows for efficient storage and querying in
BigQuery using a single partitionable value. GVS references to "position"
usually refer to position on a given chromosome and do not mean the same thing
as "location".

The canonical GVS location encoding and decoding logic lives in
`SchemaUtils.java`. Of particular interest are the following methods:

- `encodeLocation` encodes a chromosome and position into a single location
  value.
- `decodeContig`, `decodeChromosomeIndex` and `decodePosition` decode a location
  value to obtain contig, chromosome index or position respectively.

## Core BigQuery Tables

- **`sample_info`**: Holds sample name, control / non-control status, and a
  withdrawal date if applicable.
- **`sample_load_status`**: (deprecated, used for Storage Write API-based ingest
  only). This table tracks which data has been loaded for each sample. As
  described in the LoadStatusValue enum in LoadStatus.java, the values currently
  used include REFERENCES_LOADED, VARIANTS_LOADED, HEADERS_LOADED. For
  historical reasons some legacy values are also recognized.
- **`sample_chromosome_ploidy`**: This table contains a per-sample,
  per-chromosome of the ploidy observed in its input GVCF file.
- **`vet_%`**: Primary variant information (sample id, location, alleles,
  annotations). These tables are partitioned by sample id. Due to BigQuery's
  4000 partition limit there can be many `vet_%` tables in a GVS dataset, named
  `vet_001`, `vet_002`, etc. Schema information for vet tables can be found in
  the `GvsCreateTables.wdl` workflow in the variable `vet_schema_json`.
- **`ref_ranges_%`**: Reference genome locations, lengths, and qualities.
  Similar to vet, this table is partitioned by sample id and can hold at most
  4000 samples, so there will be multiple `ref_ranges_%` for larger datasets.

  Note that `ref_ranges_%` tables have two alternate schemas, compressed and
  uncompressed, to accommodate different data storage needs. Both schemas are
  described in the `GvsCreateTables.wdl` workflow. The canonical implementation
  of reference compression can be found in
  `SchemaUtils#encodeCompressedRefBlock`. Implementations of reference
  decompression can be found in `create_ranges_cohort_extract_data_table.py` and
  `GvsExtractAvroFilesForHail.wdl`, both relying on a BigQuery temporary
  function called `UnpackRefRangeInfo`.

- **`vcf_header_lines`**: Contains all unique VCF header seen while ingesting
  all samples.
- **`sample_vcf_header`**: Contains header information on a per-sample basis,
  referencing entries in the `vcf_header_lines` table.
- **`vcf_header_lines_scratch`**: Scratch table for temporary storage of VCF
  header lines during ingestion. This table is used to avoid duplicate entries
  in the `vcf_header_lines` table.
- **`alt_allele`**: This table contains allele information similar to that in
  the `vet_%` tables, but is partitioned by location to support creation of a
  filter model in joint calling. The `alt_allele` table also splits
  multi-allelic (aka "hetvar" or "1/2") genotypes into separate rows for each
  allele.
- **`filter_set_sites`**: Contains filtering information at a site level only,
  no alleles.
- **`filter_set_info`**: Contains filtering information at a position and allele
  level, including the VQSR model and other metadata.
- **`cost_observability`**: Operation cost tracking, written to by the various
  tools that comprise the GVS pipeline.

### Filter table gotchas

**`yng_status` is only ever `G` or `Y`.** `N` is recognized by the extract code
but is never written by the filter model; on Foxtrot the distribution across
`filter_set_info` was `G` 1,724,533,270 and `Y` 29,333,911. Any check keyed on
`yng_status = 'N'` is therefore vacuous, returning zero rows whether or not the
condition it tests for occurs. `Y` comes from `POSITIVE_TRAIN_SITE`
(`CreateFilteringFiles.java`), assigned by matching against left-aligned truth
resources, so input data that is not left aligned cannot be rescued by a `Y` and
is judged on calibration sensitivity alone.

**`filter_set_sites` records only sites that have filters, and does not cover
the filter model.** There is no PASS or null row; absence from the table means
the site was not filtered. On Foxtrot it held 71,032,776 rows across
`EXCESS_ALLELES`, `ExcessHet`, `LowQual`, `NO_HQ_GENOTYPES` and combinations.
Calibration sensitivity failures are *not* among them: VETS and VQSR site
failures are computed at extract time from `filter_set_info` scores, never
stored here. So this table records which sites carry site-level filter values.
Whether carrying one actually excludes the variant depends on the artifact: per
"Hard filtered versus soft filtered variants" above, the VAT is hard filtered and
drops them, while the VCF, VDS and PGEN artifacts are soft filtered and emit
every variant with its filter information attached.
`filter_set_info` answers a different question again: keyed by
`(location, ref, alt)`, it gives the calibration sensitivity score and
`yng_status` of an individual *allele*. Whether a *site* passes the filter model
is derived from its alleles by the max-over-class rule below and is not recorded
in either table.

**The filter model is applied per site and per genotype, with opposite senses.**
A site is judged by its *best* allele of a class, so a poor allele at a site
whose best allele passes is not excluded; a genotype is judged by its *worst*
non-ref allele, so a `1/2` call is filtered if either allele fails. Note that
which arithmetic expresses "best" depends on the model, and the two are
inverted. For VQSR, `vqslod` is higher-is-better, so `ExtractCohortEngine`
fails a site when `max(score) < threshold` and a genotype when
`min(score) < threshold`. For VETS — the default — `calibration_sensitivity` is
lower-is-better, so `ExtractCohortVETSEngine` overrides both and fails a site
when `min(score) > threshold` and a genotype when `max(score) > threshold`.

### `bq --apilog=false` writes a file named `false`

`--apilog` takes a file path, not a boolean, so `--apilog=false` writes the
API debug log to a file literally named `false` in the current directory.
The idiom appears throughout the WDLs, where it is harmless because the file
lands in the task's working directory and is discarded. Do not copy it into
commands run locally: there it leaves a log of several hundred KB in the
repo. Omit the flag instead.

## Key Workflows

### GvsJointVariantCalling.wdl

This is the top-level workflow for GVS Beta joint calling. It orchestrates the
joint calling process, including sample ingestion, variant calling, and
filtering. The workflow is designed to be run on the Terra platform and can be
customized for specific datasets. This workflow is _not_ used for AoU callsets,
which use a different set of workflows as described in `AOU_DELIVERABLES.md`.

### GvsBulkIngest.wdl

This workflow is used to ingest large sets of samples into GVS and is
extensively documented in `gvs-bulk-ingest-details.md`. This workflow is called
by `GvsJointVariantCalling.wdl` for Beta callsets and is called directly for AoU
callsets. One significant change for the Foxtrot callset is that this workflow
will now be directly supplied with a "FOFN" (File of File Names) TSV rather than
reading from Terra data tables. This change was motivated by hitting the scale
limits of Terra data tables and being unable to read AoU-sized sample sets, but
providing this workflow with a FOFN TSV directly is also much simpler than the
previous TSV to Terra data table to TSV approach.

### GvsPopulateAltAllele.wdl

This workflow populates the `alt_allele` table with allele information from the
`vet_%` tables. This table will enable the creation of a filter model in joint
calling in subsequent steps. The `alt_allele` table is partitioned by location
with split multi-allelics to support efficient querying and filtering.

### GvsCreateFilterSet.wdl

This workflow reads the `alt_allele` tables, scattering widely over locations to
generate a filter model and populate the `filter_set_sites` and
`filter_set_info` tables.

### GvsPrepareCallset.wdl

This workflow prepares a GVS callset for VCF or PGEN extraction. Sample,
reference, and variant data is copied to a set of "prepare tables" that are read
by GVS extract tools.

### GvsExtractCallset.wdl / GvsExtractCallsetPgenMerged.wdl

These workflows extract a GVS callset into VCF or PGEN files. They read the
prepare tables created by `GvsPrepareCallset.wdl` and generate the requested
output files. Output files are partitioned at a minimum by chromosome, but for
"wide" extracts (callsets with many large numbers of samples), there may be
hundreds of output files per chromosome.

## Nested inputs, and why a WDL 1.1 upgrade would change them

The GVS WDLs are `version 1.0`. Under 1.0, a call input that the caller leaves
unset becomes an input of the enclosing workflow, addressed by its
fully-qualified name. This is true whether the call is to a task or a
subworkflow. Terra lists these nested inputs in the workflow configuration, so
you can set them there with no plumbing. For example, `GvsQuickstartIntegration`
does not pass `use_manual_clinvar_update` to its `GvsQuickstartVATIntegration`
call, but in VS-2029 the flag could still be set from Terra when running the
top-level workflow. Terra only exposes one level of this, though: the inputs of
calls made directly by the workflow being run, not those of calls nested inside
a subworkflow. In VS-2029, setting `GvsCreateVATfromVDS`'s `use_reference_disk`
from `GvsQuickstartIntegration` (two levels down) required temporarily adding a
pass-through input to `GvsQuickstartVATIntegration`. So before adding
pass-through inputs, check whether the input belongs to a direct child call and
is therefore already exposed.

WDL 1.1 reverses the default. By default, only the top-level workflow's own
`input` section is settable. A workflow opts back in to nested inputs with
`allowNestedInputs: true` in its `meta` section. WDL 1.2 moves that key to the
workflow `hints` section, and restricts it to inputs that are optional or have
a default. The only WDL in this repo that sets it today is
`scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl`, which is not a GVS
workflow.

An upgrade to 1.1 therefore has to either add `allowNestedInputs` to every
workflow whose nested inputs are set in Terra, or plumb each of those inputs
through explicitly. Otherwise, any existing Terra method configuration that sets
a nested input should be rejected for supplying an unexpected input. This has
not been tried here. The quickstart integration workflows and the AoU workflows'
method configurations are the places to audit.

## Rerunning a quickstart integration test

Each quickstart integration subworkflow starts by creating its own BigQuery
dataset with `CreateDatasetForTest` in `GvsUtils.wdl`. The dataset name is built
from the UTC date, the branch, the short HEAD hash and a per-subworkflow suffix,
for example `quickit_2026_10_02_vs_2029_clinvar_clean_258052a_vat`. The task is
`volatile`, so it is never call-cached, and its `bq mk` fails if the dataset
already exists. A rerun of the same commit on the same UTC day therefore fails
at that step, before any of the work under test runs.

So when a quickstart run fails in `CreateDatasetForTest` with an "already
exists" error, the cause is a same-day rerun of the same commit, not the change
under test. When suggesting a rerun, say it needs a new commit or a later UTC
date. The usual way to get a new commit is an empty one, pushed so the task's
clone from GitHub picks it up (e.g. `git commit --allow-empty -m "hash bump"`,
as in VS-2029's `d179c1de6`). The empty commit takes effect as soon as it is
pushed, because `CreateDatasetForTest` gets its hash from its own `git clone` of
the branch when the task starts.

A WDL change has a separate timing concern. Terra reads the WDL from GitHub, at
`raw.githubusercontent.com/<org>/<repo>/<branch>/...`, when the run is
submitted, not from Dockstore's stored copy. GitHub serves these files with
`cache-control: max-age=300`, so allow up to 5 minutes after pushing a WDL
change before submitting. To confirm afterwards what a run executed, compare
`submittedFiles.workflow` in its metadata against the commit; its `workflowUrl`
shows where it came from. Its recorded `git_hash` is from the clone, so if the
two disagree, something was pushed between submission and the clone.

Getting past that step does not mean the run did fresh work. With call caching
on, later tasks whose inputs have not changed still reuse earlier results.
Before treating a run as evidence for a change, check the tasks that change
affects for cache hits.

## Variant Annotation Table (VAT)

### Overview

The Variant Annotation Table (VAT) is an artifact that is currently only
delivered as part of AoU callsets. See the documentation in
`scripts/variantstore/variant-annotations-table/README.md` for more information
on the VAT. The core workflow for generating the VAT is
`GvsCreateVATFromVDS.wdl`, which is run after the AoU callset has been created.

### VIDs

VIDs, or “Variant IDs” are the primary key of the VAT table. Here are three
example VIDs:

1. 1-668638-G-GA
2. 16-88461103-AGGCCCTGGTGAGGACCCAGGGATGGGGCGGGAGACCTGGTGAGGACCGAGG-A
3. Y-56694638-G-A

A VID consists of four components separated by dashes. From left to right these
are chromosome, position, reference allele, and variant allele.

In the first example above, the variant described by the VID is on chromosome 1
(chr1 in hg38 terms), position 668638, with a reference allele of G and a
variant allele of GA. Because the variant allele is longer than the reference
allele, this variant represents an insertion.

In the second example the variant is on chromosome 16 (chr16) at position
88461103 with a reference allele of
AGGCCCTGGTGAGGACCCAGGGATGGGGCGGGAGACCTGGTGAGGACCGAGG and a variant allele of A.
Because the variant allele is shorter than the reference allele, this represents
a deletion.

In the third example, the variant is on chromosome Y (chrY) at position 56694638
with a reference allele of G and a variant allele of A. This VID represents a
SNP (single nucleotide polymorphism).

### VID to Participant ID Mapping Table.

See the documentation in `AOU_DELIVERABLES.md` for information on how to
generate the VID to Participant ID Mapping Table.

### The VAT's gnomAD version depends on reference disks

`AnnotateVCF` in `GvsCreateVATfromVDS.wdl` reads its Nirvana data sources
from one of two places, and they do not hold the same data. With
`use_reference_disk = true` (the default) it reads the Broad bundle
`gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/`,
which has gnomAD 3.1.2. With it false, `Nirvana/Downloader.dll` fetches
whatever Illumina currently serves; in September 2026 that was gnomAD 4.0
(plus a newer bundled ClinVar). So the two paths produce different gnomAD
columns in the VAT, and the download path is not pinned to any version.

gnomAD 4.0 is not just newer numbers. It renames the `oth` subpopulation to
`remaining` and adds `ami` and `mid`, none of which
`create_vt_bqloadjson_from_annotations.py` knows about. On the download path
the `gnomad_oth_*` columns therefore come out empty, and `ami`/`mid` are
silently dropped. The integration truth data comes from the reference-disk
path, so nothing tests the other one.

Note that `use_reference_disk` is a WDL input, distinct from Terra's "Use
reference disks" submission option. Unchecking the option does not change the
input, and the option is not part of the call-cache hash. If a rerun with the
two mismatched gets past dataset creation (see "Rerunning a quickstart
integration test" above), what happens depends on call caching. With the option
unchecked and `use_reference_disk` still `true`, `AnnotateVCF` fails with "Could
not find reference file 'Homo_sapiens.GRCh38.Nirvana.dat'" (VS-2029,
`b3226270`). With call caching on, an `AnnotateVCF` whose inputs have not
changed instead silently reuses the previous run's shards, because the dataset
name is not one of its inputs (VS-2029, `ab40f7c3`). So when verifying a VAT
run that was meant to test a reference-disk change, check that the `AnnotateVCF`
shards were not cache hits. When suggesting such a run, recommend turning call
caching off for that submission.

### A Jasix JSON error in `AnnotateVCF` is a quoted ClinVar trait name

Nirvana writes ClinVar phenotype strings into its JSON without escaping
double quotes. A trait name that contains them, such as `"See Cases"` on
VCV001344805 at chr2:20447683, comes out as `""See Cases""`. Nirvana itself
finishes, and so does Jasix `--section genes`. Jasix `--section positions`
then fails:

    ERROR: After parsing a value an unexpected character was encountered: S.
    Path 'variants[0].clinvar[0].phenotypes[0]', line 1, position 601.

The character it names is the first letter of the trait. The failure is
deterministic, so only the shards holding such a position fail, identically
on every attempt, and `maxRetries` does not help. It does not depend on
`use_manual_clinvar_update`: the same position failed in June 2025 with the
bundled ClinVar and in VS-2054 with the 2025-07 build.

The workaround is the `sites_to_exclude` input (VS-1676), a file with one
bcftools-style region per line. It drops every allele at each listed
position before annotation, so the VAT has no row there. AnVIL 3K needs
`chr2:20447683` excluded on every VAT run. The underlying bug is unfixed, so
a callset with a variant at any other ClinVar position with a quoted trait
name fails the same way.

To find the position, download the failed shard's partial positions JSON and
take its last `"position"`. The culprit is the next row of the shard's input
VCF. Before treating the failure as new, grep the WDL for that position:
the `sites_to_exclude` help text names the one already known.

### Comparing VATs after the BigQuery dataset has expired

Integration-test VAT datasets expire after 14 days, but the workflow's TSV
outlives them in the workspace bucket. It is the `vat_complete.bgz.tsv.gz`
output of `MergeVatTSVs`, inside `GvsCreateVATFilesFromBigQuery`; find it by
walking the run's Cromwell metadata with `expand_sub_workflows=True`. Row order
is not guaranteed, so `gunzip | LC_ALL=C sort` both TSVs before comparing. If
`cmp` finds them identical, you are done. Otherwise, diff the columns keyed on
`(vid, transcript)`, taking column names from the TSV's own header line. Do not
use `scripts/variantstore/scripts/variant_annotation_table/schema/vat_schema.json`
for this: its column order is not the TSV's, so it silently mislabels every
column. In VS-2029 it reported ClinVar differences as gnomAD `fin` and `nfe`
ones. Sorting moves the header out of line 1, so read it from the unsorted
file. A per-column diff separates ClinVar changes from gnomAD changes. In
VS-2029 it showed that the flag-on VAT matched the VS-1994 VAT in every ClinVar
column, with every remaining difference in gnomAD.

### A VAT run's AC and AN outlive an `AnnotateVCF` failure

`StripCustomAnnotationsFromSitesOnlyVCF` runs before `AnnotateVCF` and
writes one `*.custom_annotations.tsv` per shard, after normalization and
deduplication. Its `AC`, `AN` and `AF` columns are what the VAT would report
as `gvs_all_ac`, `gvs_all_an` and `gvs_all_af`, followed by the same for each
subpopulation. So a run that fails in annotation still supports any
comparison that needs only those values, with no rerun. In VS-2054, 523
shards gave 125 million vids this way.

One `bq load` takes all the shards, with a wildcard over
`call-StripCustomAnnotationsFromSitesOnlyVCF/shard-*/*.custom_annotations.tsv`,
`--skip_leading_rows=7` for the header block and `--quote=''`. Load the
subpopulation columns as STRING and `SAFE_CAST` them: a subpopulation's AF
is `.` where its AN is 0, which fails a typed load with "Unable to parse".
The files can also hold a few exact duplicate rows (20 in VS-2054), so use
`SELECT DISTINCT` before joining on the vid.

### A vid from several VDS rows reports one row's AC and AN

`vet_%` alleles are not left-aligned, so the same indel can sit at different
positions in different samples and occupy several VDS rows. `bcftools norm`
brings those rows to one vid. `RemoveDuplicatesFromSitesOnlyVCF` in
`GvsCreateVATfromVDS.wdl` then keeps the row with the highest AC and drops the
others, so that the vid is reported at all. For such a vid, `gvs_all_ac`,
`gvs_all_an` and `gvs_all_af` describe one representation, not the variant.
The participant mapping takes carriers from every row, so its carrier count
can exceed what the VAT's AC implies. Do not treat the VAT's value as ground
truth when validating anything computed across all of a vid's rows; restrict
such comparisons to vids that come from a single VDS row.

A VDS key can also have no carriers at all, when every genotype carrying it
has GQ 0 or fails FT. On the AnVIL 3K VDS that was 38 of the 113 vids that
non-left-aligned keys normalize to (VS-2054). An inner join from keys to
carriers drops such a key, and with it a vid whose only keys have none. So
when building anything per key or per vid, start from the keys and LEFT JOIN
the carriers. The participant mapping's `_realigned_inputs` table lists only
keys with carriers, so it undercounts keys for the same reason.

Only the VAT's Hail pass densifies the VDS (`hl.vds.to_dense_mt` in
`hail_create_vat_inputs.py`). The participant mapping pass reads variant data
only, so anything that needs reference calls has to go in the VAT's pass or
densify again.

# Data Handling

## Redact participant IDs from anything committed to this repo

The GATK repository is public, and AoU is highly sensitive about participant
data appearing in public places. Any AoU participant ID must be redacted before
a file containing it is committed, replaced with a placeholder such as
`<PERSON_A>`. This applies everywhere, not just to data files: SQL scripts with
real IDs in an `IN` list, Python fixtures, pasted query output in a Markdown
write-up, an example invocation in a docstring, and a commit message all count.

Attach the unredacted data to the JIRA ticket for the work instead, and refer to
it from the repo by ticket number.

Redact before committing rather than afterwards. A commit is not undone by a
later one — the ID stays reachable in the history, so a slip has to be fixed by
rewriting history, which is disruptive and easy to do incompletely. Prefer
placeholders in the working file from the outset, so there is no version of it
that can be committed by accident.

Note that a redacted file and its working counterpart tend to drift apart. Where
both are needed, keep the redacted one as the file that is edited, and treat any
real-ID version as a throwaway.

# WDL Conventions

## Forwarding an optional input must preserve its undefinedness

When a wrapper workflow forwards an optional input to a sub-workflow, write the
fallback so that "nobody set this" survives the hop:

    Int? effective_width = if defined(specific_width) then specific_width
                           else general_width

Do not reach for `select_first([specific_width, general_width, <sentinel>])`. A
sentinel tail turns an undefined input into a concrete value, and callees in
this repo branch on `defined()` to mean something. `GvsImportGenomes` chooses
its own scatter width when `load_data_scatter_width` is not supplied
(`GvsImportGenomes.wdl:144-167`), so a `-1` tail would have silently taken
automatic selection away from every caller that sets neither input — while
type-checking, validating and reading perfectly well.

Two cheap checks. `womtool inputs <wrapper>.wdl` must still report the
forwarded input as `Int? (optional)` rather than `Int`. And the result is
observable at runtime: with the widths left unset, both `GvsImportGenomes`
passes scatter `GenerateParquetFilesFromInputGVCFs` to the sample count,
because automatic selection clamps to `num_samples`.

Note that a small cohort hides an explicitly set width — on three samples every
width of 3 or more yields the same three shards — so to demonstrate that an
explicit value is honoured rather than ignored, set it *below* the sample
count.
## Heredocs in command blocks

Cromwell dedents a `command <<< >>>` block by the leading whitespace common to every
non-blank line in it. So a single line at column zero -- a wrapped string, a
pasted comment -- drops that common prefix to nothing, nothing is stripped, and
the indentation the block was written with reaches bash intact.

Two things have to survive that dedent:

1. the terminator must land at column zero, or bash reads to end of file. This
   holds for every heredoc, whatever it feeds;
2. where the heredoc feeds an interpreter that rejects leading whitespace on its
   input -- in this repo that means inline Python -- the body must land at
   column zero too, or that interpreter is handed input indented by whatever the
   dedent failed to strip. A heredoc feeding `cat`, SQL, R or `yq` gets its body
   as data and does not care what column it arrives at, so only the first rule
   applies to those.

The first case fails loudly at task execution: `unexpected end of file`
reported at the generated script's last line, nowhere near the cause, and only
once the task runs in the cloud. `womtool validate` passes throughout, because a
command block is just a string to it.

The second case is quieter, and it is what the obvious repair for the first produces
-- outdent the terminator so bash is satisfied, leave the body alone. `bash -n`
then passes, since the terminator really is where bash wants it, so the script
runs and the inline Python fails on its own input instead: `IndentationError:
unexpected indent` at `"<stdin>", line 1`. The interpreter is reading the
heredoc on stdin and counts from the first line of it, so that report names
neither the WDL, nor the task, nor the real line.

Two conventions reach column zero and both are fine: indent the terminator with
the rest of the block and let the dedent strip it, or write it at column zero and
accept that the block is then never dedented at all. What fails is mixing them.
The choice is not free for a whitespace-sensitive body, though, because that body
has to agree with its own terminator -- which is why the second convention forces
inline Python against the left margin wherever this repo uses it.

Check before handing over a WDL:

```shell
scripts/variantstore/scripts/check-wdl-heredocs              # the variantstore tree
scripts/variantstore/scripts/check-wdl-heredocs FILE...      # just these
scripts/variantstore/scripts/check-wdl-heredocs --strict     # fail on FRAGILE too
```

The script does two things, and the second is not implied by its name. It flags
heredocs that will not survive the dedent, and it compiles the body of any
heredoc feeding a Python interpreter. Nothing else in the toolchain parses that
code: CI runs the Python unit tests and no linter, and a linter would not reach
this body in any case, since it is a string inside a `.wdl`. `womtool validate`
treats the command block as exactly that. So a missing colon in an embedded
script is otherwise found by running it on a cluster. WDL interpolations are stubbed out before
compiling, which makes the check structural -- it will not catch a `~{...}` that
interpolates to something Python cannot parse. It exits non-zero on a broken
block. `test/test_check_wdl_heredocs.py` runs it over every WDL under
`scripts/variantstore` as part of the Python unit tests, so CI enforces this.

"FRAGILE" means only that a block is never dedented, so it works today but breaks
the moment a heredoc whose terminator is indented with the block is added next to
it. Two blocks in the tree are in that state deliberately and need no action.

Like `reflow-md` and `check-us-spelling`, this script is deliberately
extensionless so that the Dockerfile's `COPY *.py /app/` leaves it out of the
Variants image and out of the rebuild-and-bump obligation. Do not rename it to
`check_wdl_heredocs.py`.

## `bq query` traps that pass every local check

Two failures in the `bq` CLI are invisible to `womtool validate` and
`check-wdl-heredocs`, because to both the query is only a string. Each surfaces
only when the task runs, and neither error message names its cause.

1. **A query string that begins with `-` is parsed as a flag.** That includes a
   leading `--` SQL comment. absl then tries to suggest a flag name close to the
   whole query, and its recursive edit-distance function exhausts the stack: a
   `RecursionError` traceback from `_damerau_levenshtein`, with no mention of
   the query. Start every query with a statement and put comments after it, or
   pass the query on stdin rather than as an argument.
2. **BigQuery will not create a view in a script that declares a temporary
   function**, whether or not the view uses it: `Creating views with temporary
   user-defined functions is not supported`. Because the failure comes at the
   `CREATE VIEW`, every table the script built before it already exists. A
   script using `CREATE TABLE` rather than `CREATE OR REPLACE` then fails on
   rerun until they are dropped. Create the view in a separate `bq query` call.

## `ARRAY_LENGTH` bills the whole array

BigQuery bills an array column in full for any expression that touches it,
`ARRAY_LENGTH` included. So counting carriers in the participant mapping's
`_base` table (`het_ids`, `hom_ids`, `hemi_ids`) reads nearly all of that
table's bytes, however few rows the rest of the query keeps. In VS-2054 a
query that looked like a scan of two vid columns was estimated at 18.75 TB
on Foxtrot for this reason.

A filter only reduces the bill if it lets BigQuery skip blocks. `_base` is
clustered by `vid`, so that means a filter on `vid` itself against
constants:

- `WHERE vid IN (SELECT ...)` and a join on `vid` are not reliably pruned.
  Collect the vids first and pass them as a script variable:
  `SET vids = (SELECT ARRAY_AGG(vid) FROM ...)`, then
  `WHERE vid IN UNNEST(vids)`.
- A filter on an expression over the column, such as
  `SPLIT(vid, '-')[OFFSET(0)] = '21'`, prunes nothing. The chr21 run of the
  VS-2054 query scanned the VAT's whole `vid` column.

So do the cheap narrowing first, on columns that are small, and go to the
array columns only for the rows that survive. Rewritten that way, the
VS-2054 query billed 283 GB genome-wide. The same pattern probed
`alt_allele`, which is clustered by `location`, at 1,053 locations for about
57 GB.

A dry run cannot see clustering pruning: for a clustered table it reports
an upper bound, and for a multi-statement script the console gives no
estimate at all. When suggesting such a query, tell the user to set
"Maximum bytes billed" (`--maximum_bytes_billed` in `bq`), so a query that
fails to prune stops without charge. Then ask for the bytes billed
afterward, since that is the only evidence the pruning worked.

## A `bq load` wildcard over Spark output also reads `_temporary/`

A `*` in a BigQuery source URI matches across `/`, so `gs://.../out.parquet/*.parquet`
matches every Parquet file anywhere under the directory, not just its top level.
For Spark output that includes `_temporary/`. Spark normally deletes that
directory when the job commits, but an executor lost to preemption can leave
partial files there afterward, even when the write itself succeeded. The load
then fails with `Input file is not in Parquet format`, naming a file under
`_temporary/0/_temporary/attempt_.../`. That file is from an attempt that never
committed, so the committed output is not at fault. In VS-2041 this broke the
chrY load after a run that had otherwise succeeded, though the same glob had
loaded chr21 without trouble.

So never glob a Spark output directory from its root. Name the directories that
hold committed output, such as each `partitionBy` directory, as a
comma-separated list of URIs, and leave out any that do not exist, since a URI
matching nothing fails the load. `LoadMappingInputs` in
`GvsCreateParticipantMappingFromVDS.wdl` does this. Local readers are exposed in
the same way: DuckDB pointed at `out.parquet/**/*.parquet` reads the partial file
too.

Output that is not partitioned has its part files directly in the directory,
so there are no subdirectories to name. List the files with gcloud instead,
whose `*` does not cross `/`:

    gcloud storage objects list "<dir>/part-*.parquet" --format='value(name)'

This skips `_temporary/`, and it exits 0 with empty output when nothing
matches, so an empty result can be told from an error. The names come back
without `gs://<bucket>/`, so prepend it. Do not use `value(storage_url)`,
which appends `#<generation>`. A wide table has too many parts to pass one
URI each, so collapse them to one wildcard per part-number prefix
(`part-00*.parquet`, `part-01*.parquet`, and so on). `LoadNoCallSets` in
`GvsCreateNoCallSetsFromVDS.wdl` does this.

## A Hail export from a wide table finishes on the master alone

A non-parallel `hl.export_vcf` (or `Table.export`) writes one part file per
partition, then merges them into a single file on the master. The merge is
serial and its cost tracks the part count, not the data size, so a tiny
output from a wide table can add hours of wall-clock time. An autoscaling
cluster sheds its secondary workers meanwhile, so the cost is mostly time,
plus scaling back up for the next step. In VS-2041 the genome-wide
participant mapping exported 14,459 keys from 119,189 partitions. Writing
them took a full pass over the VDS (7.4 hours at 200 workers). Merging
119,190 parts totaling 3.7 MB took another 37 minutes on the master, which
then ran for 1 hour 50 minutes with no log output while the cluster scaled
down to no secondary workers. On chr20, with 2,576 partitions, none of this
was visible.

So when the result is small, compute it in an aggregation you already run,
for example `hl.agg.filter(..., hl.agg.collect(...))` alongside the counts.
Then write it from a one-partition table built with `hl.Table.parallelize`.
`write_keys_to_normalize` in `hail_create_participant_mapping_inputs.py`
does this. On chr21 it halved the classify step (1,136 to 515 seconds) and
the merge dropped to 2 files, with identical outputs. Two actions on the
same table also mean two passes over the input, because Hail does not keep
results between actions. Before adding an export or a `take` after an
`aggregate`, check whether the aggregate can return what you need.

To tell a driver-bound phase from a worker-bound one, compare the Hail log
timestamps with the cluster's autoscaler log in Cloud Logging
(`logName="projects/<project>/logs/dataproc.googleapis.com%2Fautoscaler"`,
filtered on `resource.labels.cluster_name`). The log survives the cluster's
deletion for the 30-day retention period.

The autoscaler log does not reliably show Spot preemptions. Dataproc replaces
a preempted secondary worker, so the cluster size it records can stay at its
maximum while workers are being lost. Count preemptions from the Compute
Engine audit log instead:

    gcloud logging read 'protoPayload.methodName="compute.instances.preempted"
      AND protoPayload.resourceName:"<cluster name>"' \
      --project <project> --freshness=2d --format='value(timestamp)'

A preemption burst while the cluster is fully allocated looks like stalls:
allocation dips with nothing pending, then the backlog comes back as lost work
is recomputed.

## Why Hail is pinned at 0.2.134

Nothing in GVS itself needs 0.2.134. Each workflow installs whatever
`hail_version` names and starts its own cluster with `hailctl`, so it can run
any version. The pin matches Terra instead. As of November 2026, Terra's
Leonardo runtime ships Hail 0.2.134 on Dataproc 2.2.52
(`leo-dataproc-image-2-2-52-debian12-...`). A VDS written with that version
opens in a Terra notebook as Terra ships it, and an ad-hoc check there runs
the same Hail as the workflows. So the limit falls on what readers can open,
not on what GVS can run.

Each Hail release fixes its Dataproc image in `IMAGE_VERSION` in
`hailtop/hailctl/dataproc/start.py`:

| Hail            | Dataproc image  | Spark |
|-----------------|-----------------|-------|
| 0.2.131–0.2.135 | 2.2.5-debian12  | 3.5.0 |
| 0.2.136         | 2.2.60-debian12 | 3.5.0 |
| 0.2.137–0.2.139 | 2.3.17-debian12 | 3.5.3 |
| 0.2.140         | 3.0.1-debian13  | 4.1.2 |

The change log does not mention the image changes in 0.2.136 and 0.2.137, so
read the source at the release tag rather than relying on the log.

Ad-hoc work in a Terra runtime therefore stops at the 2.2 line. That means
0.2.135, or 0.2.136 if it is tested first. Beyond that, the newer version
has to run on a cluster you start with `hailctl` yourself. Start it from the
Terra terminal, which runs as the workspace's pet service account, with the
`--service-account` and `--subnet` flags that `run_in_hail_cluster.py`
passes. From your own account it cannot work, because you can neither act as
the pet service account nor impersonate it. In VS-2067, a cluster create in a
Terra workspace project as the pet account failed with
`iam.serviceAccounts.getAccessToken` denied. Once the cluster is up, your own
account can submit jobs to it, because that needs only `dataproc.jobs.create`.

### Hail's version numbers do not signal breaking changes

Every Hail release is numbered `0.2.N`, so the version number never says
whether a release breaks anything. Breaking changes do arrive under these
increments:

- **0.2.139 changed the on-disk index format, in one direction only.**
  hail#15445 writes each index as a single flat `partname.idx` file instead of
  a `partname.idx/` directory. Newer Hail still reads the old layout, but older
  Hail cannot read the new one, and there is no option to keep writing the old
  format. A VDS contains row indexes, so a VDS written by 0.2.139 or later
  cannot be read by anyone on an older Hail. For an AoU deliverable that means
  any researcher whose Workbench environment is older.
- **0.2.140 changed the whole runtime.** It moves to Spark 4.1, Scala 2.13,
  Java 17 or later, Dataproc 3.0 images, and Python 3.12 or later.

So before proposing a `hail_version` bump, read the changelog
(`hail/python/hail/docs/change_log.md` at the target tag) for every release
between the pinned version and the target, not only the target's entry, and
read `start.py` at the target tag for its Dataproc image, which the changelog
may not mention. Then check two things. First, that researchers' Hail can
read a VDS written by the new version. Because of the 0.2.139 index change, no
release from 0.2.139 on passes this for readers on an older Hail, so raise it
with AoU before a GVS VDS is written with one. For any other candidate, write
a small VDS with it and read it with the readers' version. Second, if the new
version needs a newer Dataproc image than Terra's, ad-hoc work moves to
clusters launched with `hailctl`. VS-2066 tracks the upgrade.

## Changing a Dataproc machine type is more than a name swap

`run_in_hail_cluster.py` starts clusters through `hailctl dataproc start`, and
three constraints apply to any machine type passed to it. All pass `womtool
validate` and fail only when the cluster is created.

1. **hailctl only knows the families in its `MACHINE_MEM` table.** It looks the
   type up there to size Spark memory, so an unlisted type fails with a
   `KeyError`. Through Hail 0.2.139 the table covers n1, n2, n2d, e2, c2, m1
   and m2, but not C3, C3D, N4 or any C4 family. The pinned version
   (`hail_version` in `GvsUtils.wdl`) is 0.2.134. Hail 0.2.140 removes the
   table and asks `gcloud compute machine-types list` for the vCPUs and memory
   instead, so it accepts any type that gcloud knows. Upgrading to get that is
   not a small change, though: see "Why Hail is pinned at 0.2.134" above. To
   try a newer family on the pinned version, add entries to the table before
   hailctl starts, since memory sizing is the only thing it uses the table
   for. For secondary workers alone there is a second route that needs no
   table entry: `--secondary-worker-machine-types`, described in "Newer
   families as Spot secondaries" below. Read the pinned version's `start.py`
   before proposing a type.
2. **Valid local SSD counts depend on the machine family and vCPU count.** n1
   accepts one SSD at any size, but N2 and N2D raise the minimum as vCPUs
   grow. N2D machines with 32 or 48 vCPUs reject 1 and accept 2, 4, 8, 16 or
   24; an N2 machine with 32 vCPUs needs at least 4. So an n2d-highmem-32
   master fails to create alongside n2d-highmem-8 workers, which accept 1:
   `Number of local SSDs for an instance of type n2d-highmem-32 should be one
   of [0, 2, 4, 8, 16, 24], while [1] is requested`. The master therefore has
   its own count, `--num-master-local-ssds` (`num_master_local_ssds` on every
   GVS workflow that takes machine types). Left unset it follows the workers'
   count, except that 1 becomes 2, so the default n2d-highmem-32 master is
   valid (VS-2055). That rule covers no other N2 or N2D master: an
   n2-highmem-32 or an n2d-highmem-64 still fails with it. Look up the allowed
   counts for both shapes in the Compute Engine local SSD documentation before
   changing either one, and set the master's count explicitly.
3. **Newer families fix the local SSD count in the machine type.** C3, C3D, C4
   and C4D get local SSD only through `-lssd` shapes, with a fixed number of
   disks: one each for `c4-highmem-8-lssd` and `c4d-highmem-8-lssd`, five for
   `c4-highmem-32-lssd`. C3 has no highmem `-lssd` shape at all, and N4, N4D
   and C4A support no local SSD. For the master and primary workers,
   `run_in_hail_cluster.py` always passes `--num-worker-local-ssds` and
   `--num-master-local-ssds`, and its "1 becomes 2" master rule is wrong for
   these shapes. Secondary workers are handled: `--secondary-worker-machine-types`
   requests 0 for an `-lssd` shape, and Dataproc still mounts the built-in SSD
   and gives it to YARN (see "Newer families as Spot secondaries" below).
   Dataproc's supported machine types page lists C3, C4, N4 and N4D, but not
   C3D or C4D (checked November 2026). C4D nonetheless ran as Spot secondaries
   in VS-2067.

An error message describes one machine shape, not its family or its
neighbors. In VS-2055 the n2d-highmem-32 error above was taken to mean that N2
and N2D with 32 vCPUs both accept 2. They do not: N2 at 32 vCPUs needs 4. That
claim went into seven WDL help strings, the script and the PR description
before review caught it. So before writing down a Compute Engine limit, look
it up in Google's table for every family and size the text names, not just
the one that failed.

n2d-highmem-8 workers and an n2d-highmem-32 master have been the default
since VS-2055, replacing n1, in the script and in every WDL that runs
`run_in_hail_cluster.py`, each of which takes its own machine types. A WDL's
own default overrides the script's, so changing the default family means
editing both.

The memory per vCPU differs too: an n2d-highmem-8 has 64 GB against an
n1-highmem-8's 52 GB. So a run on a new family is not a like-for-like
comparison of Spark settings tuned on n1.

### Newer families as Spot secondaries, and the per-family local SSD quota

Families that hailctl does not know, such as C4 and C4D, can still run as
secondary workers. `--secondary-worker-machine-types` in
`run_in_hail_cluster.py` (`secondary_worker_machine_types` in
`GvsCreateVDS.wdl`) passes them to Dataproc as ranked flexible VMs. hailctl
sizes Spark from the primary worker type, so a secondary type needs at least
that type's YARN memory. Three constraints apply:

1. C4, C4D, N4 and N4D boot only from `hyperdisk-balanced`. Their `-lssd`
   shapes carry a built-in 375 GB local SSD, so the requested local SSD count
   must be 0.
2. Dataproc gives every secondary worker the same disk configuration, so the
   listed types must agree on it. N2D cannot be the fallback for C4D, but C4
   can.
3. The `-lssd` disks count against
   `LOCAL-SSD-TOTAL-GB-PER-VM-FAMILY-per-project-region`, which defaults to
   30,000 GB per family per region. That allows 80 VMs per family. N2D's
   attached SSDs did not count against it in VS-2067: N2D clusters held 200.

So a cluster of newer families that stops short of its maximum size is not
necessarily stocked out. Before concluding that it is, check the creation
errors for `QUOTA_EXCEEDED`. In VS-2067, a C4/C4D cluster peaked at 160
secondaries, 80 of each family. 360 creates failed on the quota and only 58
on stockout. Spot CPUs come out of `PREEMPTIBLE_CPUS`, not the per-family CPU
quotas, so the CPU quotas are not the limit.

N4 and N4D are in `HYPERDISK_ONLY_FAMILIES` only so that they get the right
boot disk. They have no local SSD in any form, neither attached nor as an
`-lssd` shape, so the script gives them none, and Spark's shuffle and spill
data then lands on the Hyperdisk boot disk. VS-2067 benchmarked C4 and C4D
rather than these cheaper families for that reason: a comparison against n2d
would have changed the I/O path as well as the CPU. So do not propose N4 or
N4D as a drop-in replacement for an SSD-backed type. Treat them as untested
until a benchmark has sized the boot disk and its Hyperdisk throughput for
shuffle.

C4D's advantage over n2d has been measured on two different Hail jobs, with
similar results. On AnVIL 3K VDS creation, C4D Spot secondaries ran about 1.7
times as fast per task as n2d-highmem-8 (VS-2067). On the genome-wide no-call
sets pass, which densifies the same VDS, they used 17.0 secondary worker-hours
against n2d's 28.4, a factor of 1.67 (VS-2054). A C4D worker-hour costs about
the same as an n2d one, so the saving in cost is close to the saving in
worker-hours. Each figure is from one run per configuration, so treat it as an
estimate (see "Compare cluster configurations over repeated runs" below). C4
has only the VS-2067 measurement, about 1.5 times n2d.

A ranked list is a preference, not a guarantee: Dataproc falls back to the
next type on a stockout or the quota above. So before reading a run as a
benchmark of one family, check which types its secondaries actually were, from
the Compute Engine audit log:

    gcloud logging read 'protoPayload.methodName:"compute.instances.insert"
      AND protoPayload.resourceName:"<cluster name>-sw"' \
      --project <project> --freshness=2d \
      --format='value(protoPayload.resourceName,protoPayload.request.machineType)'

Each create is logged twice, and only the first entry carries the machine
type, so drop the lines without one. A replacement for a preempted worker
reuses its name, so creates can outnumber instance names. In the VS-2054 run
there were 54 creates under 51 names, all `c4d-highmem-8-lssd`.

## Compare cluster configurations over repeated runs

Two identical runs of the same Hail job can differ by more than the
configuration change being measured. In VS-2055, two n1 runs of AnVIL 3K VDS
creation took 2h46m and 3h14m, and two n2d runs took 1h47m and 2h00m. Most of
the n1 spread was Spot preemption: the slower run lost 106 workers while all
200 were busy, against 2 in the faster one. The rest came from straggler
stages at the end of the job.

So a single run of each configuration does not settle a comparison. When
asked to evaluate one, ask for at least two runs per configuration and count
the preemptions in each, as described in "A Hail export from a wide table
finishes on the master alone". Compare the phase where the cluster is fully
allocated with a backlog pending, as well as the total. That phase isolates
per-core speed, while the tail mostly reflects the autoscaling policy. Treat a
difference as real when the worst run of one configuration still beats the
best run of the other.

When costing the runs, price everything a secondary worker carries at Spot
rates, including its local SSD. Spot local SSD is a separate SKU, "SSD backed
Local Storage attached to Spot Preemptible VMs", at about 60% of the on-demand
rate. Pricing it on demand overstated the SSD's share of the VS-2055 costs by
about 80%. Rates change, so look them up in the Cloud Billing Catalog API
rather than reusing figures from an earlier write-up.

## Interpolated multi-line text defeats both heredoc conventions

The dedent described above is computed on the command *after* interpolation. A
`~{...}` whose value spans several lines puts its continuation lines at column
zero, so nothing is stripped. A terminator indented with the block then stays
indented, even though the block was written to the first convention. The result
is the loud failure above, `unexpected end of file`, which here also means the
task never prints what it was meant to print. `check-wdl-heredocs` cannot see
this. It stubs interpolations out with a single-line literal, so it passes such a
block, and it did pass `TerminateWorkflow` while every multi-line abort message
was being lost this way.

So do not pass interpolated text through a heredoc at all. To get a WDL string
into a file or onto stderr, let Cromwell write it:

    cat ~{write_lines([message])} >&2

bash never parses that file. A heredoc would also need its delimiter quoted
(`<<'FIN'`), or bash runs backticks and `$` in the text as substitutions, and GVS
messages routinely backtick hyphenated table names.

When a task fails with an empty `stderr` in the execution bucket, the shell's own
errors, including `unexpected end of file`, are in the task's GCP Batch log:

    gcloud logging read 'logName="projects/<terra-project>/logs/batch_task_logs"
      AND resource.labels.job_id:"job-<first 8 of the sub-workflow id>"' \
      --project <terra-project> --freshness=2d --format='value(textPayload)'

The Terra project and full job id are in the call's `jobId` in the workflow
metadata.

## Transient infrastructure failures: use `maxRetries`, not `preemptible`

`preemptible` only covers VM preemption. A task that exits non-zero for any
other reason, such as a GCS `503 Service Unavailable` surfacing from a GATK tool
after the client library has used up its own retries, fails on the first attempt
unless its `runtime` block sets `maxRetries`. The repo's convention is
`maxRetries: 3` (see `GvsExtractCallset.wdl` and `GvsUtils.wdl`).

`maxRetries` is not specific to 503s: Cromwell reruns the task on any non-zero
exit, so a genuine failure is also attempted four times before it fails. Retries
run on fresh VMs, so each one costs the full task.

A task that deletes GCS objects should judge success by what remains, not by
`gcloud storage rm`'s exit status, which can be non-zero after a successful
delete or on a retry (VS-2053). List what remains with `gcloud storage objects
list --format='value(name)'`: it exits 0 with empty output when nothing matches,
while `gcloud storage ls` exits 1, the same as for a real error.

## A Hail script's own imports must be shipped to the cluster

`run_in_hail_cluster.py` submits one script to Dataproc. Any GVS module that script
imports, such as `hail_gvs_util` or `vds_validation`, reaches the cluster only if the
WDL also passes it with `--secondary-script-path-list`, which becomes `--py-files`.
Leave one out and the job fails on the cluster with `ModuleNotFoundError: No module
named '<module>'`, after the cluster has been created. Every local check passes: the
module is in the Variants image and in `GetHailScripts`'s outputs, the unit tests
import it from the scripts directory, and `womtool validate` sees only a command
string.

So moving code into a shared module means updating every WDL whose script now
imports it, not just the scripts. In VS-1949 (#9381), `remove_samples_from_vds.py`
started importing `hail_gvs_util`, but `GvsRemoveSamplesFromVDS.wdl` was not changed.
No integration test runs that workflow, so it stayed broken for five months, until
VS-2023 ran it (fixed in VS-2062). To check, list each Hail script's imports of GVS
modules and compare them with the `--secondary-script-path-list` arguments in the
WDL that runs it:

    cd scripts/variantstore/scripts
    grep -nE '^\s*(from|import) (hail_gvs_util|vds_validation|import_gvs|import_gvs_ploidy|hail_gvs_import)\b' *.py

Allow for leading whitespace, as above: `hail_gvs_import.py` imports `import_gvs`
inside a function. The list is transitive, so `GvsCreateVDS` ships `import_gvs_ploidy`
because `import_gvs` imports it. `GvsCreateVDS` and `GvsMergeAndRescoreVDSes` ship
their modules correctly and are the examples to copy.

## Running a one-off Hail script needs no WDL and no image rebuild

`GvsRunHailScript.wdl` takes its `script` and `secondary_scripts` as GCS
paths and runs them on a Dataproc cluster through `run_in_hail_cluster.py`.
A prototype or one-off analysis therefore needs only its files copied to the
workspace bucket, including any GVS module it imports (see the section
above). Nothing comes from the Variants image, so there is no rebuild and no
tag bump. `script_arguments` is a `Map[String, String]` rendered as
`--key value`. So every argument needs a value, a key cannot repeat, and a
list has to be passed as one comma-separated value. In VS-2054 this ran the
no-call set prototype on the AnVIL 3K VDS without a new workflow.

The workflow reads the script when the task starts, so re-upload after
editing. A stale copy in the bucket runs without complaint.

## A dense MatrixTable has missing entries, and aggregations skip them

`hl.vds.to_dense_mt` gives a sample a missing *entry*, not an entry with a
missing `GT`, wherever it has neither a variant call nor a covering
reference block. Hail leaves missing entries out of entry aggregations, as
it does filtered ones. So `hl.agg.count_where(hl.is_missing(mt.GT))` counts
only the samples whose entry exists with a missing `GT`, such as a call
failing FT, and silently omits most no-calls. `call_stats` is unaffected,
because a missing entry contributes nothing to AN either way, so AN still
ties out while the no-call count is far too low. In VS-2054 a site reported
15 no-calls before the fix, against a true count several times that.

Call `mt.unfilter_entries()` before aggregating over anything that must see
every sample. To test such an aggregation, compare it against a model that
enumerates samples, not against another Hail aggregation over the same
MatrixTable.

# Documentation Conventions

## Markdown tables must be rectangular

IntelliJ wants Markdown tables to be rectangular (cell padded to its column's
widest content). A table left ragged therefore reappears as a spurious
whitespace-only diff the next time anyone opens the file in IntelliJ, and the IDE
flags it. A single cell one character wider than its column is enough to cause
this, which is nearly invisible in review — and editing one cell of an
already-aligned table is the usual way it happens.

After adding or editing a table, run:

```shell
python3 scripts/variantstore/utils/format_markdown_tables.py <file.md>
```

Add `--check` to report without rewriting (exit status 1 if anything would
change), and pass a directory to search it recursively. The script is idempotent
and never touches content inside fenced code blocks.

"Editing a table" includes edits that were not aimed at the table. A
find-and-replace that changes a word's length — a spelling pass turning
`behaviour` into `behavior`, a rename — breaks the padding of every table cell
it touches. In VS-2029 a US-spelling pass over `validation.md` left two
previously aligned tables ragged. After any bulk text edit, run the formatter
with `--check` on the files you touched.

Some existing docs predate this and are still ragged. Only format files you are
already modifying; reformatting untouched docs adds diff noise that obscures the
actual change.

## Cite a run's results by something that will still resolve

Integration-test BigQuery datasets expire 14 days after creation
(`table_ttl_seconds` in `GvsUtils.wdl`), so a doc that cites one as evidence
goes stale within weeks. In VS-2029, `validation.md` named two VAT datasets
that had expired before the doc was reviewed. The run's outputs last longer.
Every task's files stay in the workspace bucket under
`gs://<bucket>/submissions/<submission id>/` for as long as the workspace
exists, whatever the workflow.

So cite the run itself rather than a dataset it wrote to. A workspace and
submission ID may be enough. A full `gs://` path to the output file saves the
reader from walking Cromwell metadata, but it can run to several hundred
characters. Ask the user which they want.

## US English spelling

This codebase is US English. British forms slip in easily and are invisible in
review, since nothing is functionally wrong. Check before handing over prose:

```shell
scripts/variantstore/scripts/check-us-spelling <file-or-dir>...   # report; exit 1 if found
scripts/variantstore/scripts/check-us-spelling --fix <file>...    # rewrite in place
```

It uses an explicit word list rather than an `-ise` suffix rule, because "analysis",
"optimistic", "premise" and "exercise" are all correct US English and a checker that flags
them gets ignored. Only the verb shifts: `analyse` becomes `analyze` while `analysis` stays
put. Adding a word is a one-line edit to the tables in the script.

Files that legitimately contain British spellings — the word list itself, its test
fixtures, a document quoting them — opt out with `check-us-spelling: disable` anywhere in
the file, or `check-us-spelling: ok` on a single line. Both are honored by `--fix`, so an
exempt file is never rewritten.

Identifiers are deliberately out of scope: word boundaries treat underscore as a word
character, so `n_behaviour_flags` is not reported. Renaming an identifier here but not at
its definition elsewhere would turn a spelling nit into broken code.

Like `reflow-md`, this script is deliberately extensionless so that the Dockerfile's
`COPY *.py /app/` leaves it out of the Variants image and out of the rebuild-and-bump
obligation. Do not rename it to `check_us_spelling.py`.

## Editing a JIRA description through the Atlassian MCP

VS tickets live on the `broadworkbench` Atlassian site; get its cloud ID from
`getAccessibleAtlassianResources`. The MCP can only replace a description
wholesale, and it can read and write it either as markdown or as ADF (Atlassian
Document Format, JIRA's native JSON). To edit a description that already exists,
use ADF. Markdown looks like the easy route but does not round-trip.

**Markdown write-back corrupts smartlinks.** The MCP renders an inline issue chip
(an ADF `inlineCard`) as `<custom data-type="smartlink" data-id="id-N">URL</custom>`.
Sending that markdown back does not recreate the chip. JIRA stores the tag and
displays it literally. On VS-2010 this left four broken links to be repaired by
hand. The value echoed back by the edit matches what was sent, so it is not
evidence that the round-trip worked.

**ADF round-trips exactly, once it is small enough to send.** Every ADF node
carries a `localId` attribute. JIRA regenerates these, and stripping them shrank
the VS-2023 description from 90K characters to 69K, small enough to pass inline in
one `editJiraIssue` call. `scripts/variantstore/scripts/jira-adf` does the
mechanical parts:

1. Fetch with `getJiraIssue`, `responseContentFormat: "adf"`,
   `fields: ["description"]`. A long description exceeds the tool-result limit and
   is saved to a `tool-results/...txt` file. That is fine, since it never needs to
   come into context.
2. `jira-adf extract <saved result> -o desc.json`, then
   `jira-adf outline desc.json` to number the blocks. `--path 83` lists the
   children of block 83, such as the items of a bullet list.
3. Edit with guards that check the target before changing anything:
   `jira-adf delete desc.json 24:"Then two things" 83.0:"No cluster" -o desc.json`
   and `jira-adf replace desc.json "old" "new" -o desc.json`. A guard matches
   visible text only. A bare number such as `9419` would otherwise match inside a
   `localId`. For anything the subcommands do not cover, edit the JSON in Python,
   asserting on each node's text before you change it.
4. `jira-adf strip desc.json -o desc.min.json`. If the result is still too large
   to Read in one go, `jira-adf split desc.min.json -o desc.part` writes parts
   that concatenate back to it exactly. Read each part and pass the whole document
   as `fields: {"description": ...}` with `contentFormat: "adf"`.
5. Set `responseContentFormat: "markdown"` on the edit and check the returned text.
   Removed passages should be gone, and untouched landmarks should still be there:
   smartlinks, tables, the acceptance criteria.

Writing new text from scratch, such as a fresh comment, is fine as markdown. Leave
out any smartlink tags and use bare URLs or plain issue keys, which JIRA links by
itself. JIRA also hard-breaks newlines inside a paragraph, so write one line per
paragraph, as described under "One line per paragraph" below.

Do not go looking for JIRA API credentials in environment variables or `~/.netrc`
to get around the MCP. Auto mode denies it as credential exploration, and the
procedure above makes it unnecessary.

Like `reflow-md`, `jira-adf` is deliberately extensionless so that the
Dockerfile's `COPY *.py /app/` leaves it out of the Variants image. Do not rename
it to `jira_adf.py`.

# Code Review Conventions

## Do the review inline — do not use the `code-review` skill

When asked to review a branch, PR, or diff, read the diff and write the review
directly in the conversation, following the conventions in this section. Do not
invoke the `code-review` skill, and do not otherwise fork the review into a
subagent.

The skill is packaged to run forked, and what comes back does not follow the
conventions below: it leads with findings instead of an `## Action items`
section, and it does not produce the durable write-up under
`~/claude-artifacts/<project>/`. Restructuring its output afterwards costs more
than doing the review inline from the start.

A forked reviewer also cannot see the context already gathered in the
conversation — the established base ref, the diff, the surrounding files — so it
re-derives them, and its supporting claims still have to be checked against the
repo before they can be reported. A forked reviewer's findings can be sound
individually yet rest on a partial survey of the repo — correctly citing one
call site that does X while missing another that does the opposite and changes
the conclusion.

Whoever writes the review — inline or otherwise — verify each finding's
supporting claims against the repo, and reproduce a finding directly where that
is cheap, before reporting it.

## Lead with action items

A code review write-up opens with an `## Action items` section, before the scope
summary and before the findings themselves. It lists what should actually be
done, each entry pointing at the numbered finding that justifies it.

Order that list by what the PR's author can act on directly, not by strict
execution order. Their own work comes first, with the discretionary line closing
that group. Anything the branch merely runs into — pre-existing debt, a problem
in files the PR does not touch, something needing a team decision or
coordination with another branch — goes in a separate group afterward, under a
heading that says so, introduced by a sentence establishing provenance: that it
was not introduced by this branch and has not come up before, so it is something
for the team to discuss. Frame that as a shared next step, not as a limit on the
author's authority — "so we should discuss these as a team" rather than "this is
not the author's to decide".

Write those cross-cutting entries to be comprehensible cold. The author has not
been part of the conversation that produced them, so lead with the symptom and
arrive at the proposed remedy; an entry that opens with the remedy reads as an
unexplained directive. "The validation only runs on the non-Parquet ingest path,
which is not what production uses" works where "relocate the non-Parquet
coverage" does not.

When one action genuinely gates another, say so inline in the entries affected
rather than reordering the list to express it — ordering for the reader and
recording a dependency are separate jobs, and the dependency note costs a
clause.

Every finding must be accounted for in that list, including the ones needing no
work — group those into a short "discretionary" line rather than dropping them,
so nothing looks accidentally omitted.

The rationale, evidence, and reproduction detail stay in the numbered findings
below. Do not put a duplicate ordered list at the end of the document.

Both lists are numbered, so keep the two vocabularies distinct: the entries in
`## Action items` are "action N" and the write-ups below are "finding N". Never
call either one "item N" — a reference like "(item 4)" inside action item 1 is
ambiguous between the two sequences, and the reader has no way to tell which is
meant. Label the write-ups with a `## Findings` heading so the second sequence is
visibly its own list.

## One line per paragraph — do not hard-wrap

Code review write-ups get pasted into GitHub PR comments and descriptions, and
GitHub renders a newline inside a paragraph as a hard line break regardless of
how much horizontal space is available. A source document wrapped at 100 columns
therefore renders as a ragged column on GitHub rather than as flowing text.

So in a generated code review document, write each paragraph and each list item
as a single unwrapped line, however long it ends up. Headings, table rows, and
fenced code blocks are unaffected — those are already one line each by
construction, and their line breaks are significant.

Author these documents unwrapped from the start where you can. To convert a
draft that is already wrapped, use the `reflow-md` helper checked into this
repo:

```
scripts/variantstore/scripts/reflow-md FILE...            # rewrite in place
scripts/variantstore/scripts/reflow-md --dry-run FILE...  # report only
```

It folds paragraphs and list-item continuations onto single lines while leaving
headings, table rows, blank lines, and fenced-block contents alone. It compares
the whitespace-delimited token stream before and after and refuses to write if
anything other than wrapping changed, and it is idempotent — a second run
reports `already unwrapped`. Prefer it to joining lines by hand, which is easy
to get wrong around list indentation and fences, where a mistake corrupts
content rather than merely rewrapping it.

The script is deliberately named without a `.py` extension: the Dockerfile in
that directory does `COPY *.py /app/`, so an extensionless name keeps this
dev-only helper out of the Variants image and out of the rebuild-and-bump
obligation described in the Variants Docker Image section above. It is a
dependency-free Python 3 script; do not rename it to `reflow_md.py`.

JIRA hard-breaks intra-paragraph newlines the same way, so the same rule applies
to anything destined for a JIRA description or comment.

## Development scaffolding is not a review finding

Do not report temporary "point this at my feature branch so I can test it"
changes as code review findings. The most common example is adding a feature
branch name to the `branches:` filters in `.dockstore.yml` — Dockstore will not
expose a workflow for a branch that is not registered, so this edit is a
prerequisite for testing any WDL change on a branch, not an oversight.
Developers are reliable about removing these before merge, and a leftover
registration is harmless: it makes a stale branch visible in Dockstore and
nothing more.

This exclusion covers branch- and tag-registration scaffolding only. It does
**not** extend to a reference that will actually break or mislead at runtime
once merged, which remains a legitimate (and usually serious) finding:

- a Docker image tag — `variants_docker`, `gatk_docker`, etc. in `GvsUtils.wdl`'s
  `GetToolVersions` — left pointing at a personal, unpublished, or stale image
  (see the Variants Docker Image section above: a new or edited script under
  `scripts/variantstore/scripts/` does not reach a running workflow until the
  image is rebuilt and the tag is updated);
- a hardcoded project, dataset, or GCS path pointing at a developer's scratch
  resources in a code path that production runs;
- a `git_branch_or_tag` default left set to the feature branch.

The distinction is whether the leftover changes what a merged workflow does. If
it only affects what is visible or runnable on the developer's own branch, leave
it alone.

## Re-reviewing a branch starts with the previous rounds

A second or third review of the same branch is mostly triage, not fresh
reading. Before looking at the diff, read the earlier write-ups in
`~/claude-artifacts/<project>/` and classify every action item they raised as
resolved, partially resolved, or open — checking each against the current tree
rather than against what the PR body or the commit messages claim.

Most of a later round's value is in that classification. Fixes between rounds
are often partial: one of two related comments corrected, one recovery path
closed but not the other, a doc updated in one place but not where operators
actually read it. "Partially resolved" is frequently the accurate verdict.
Starting from the diff alone re-derives the open items as if they were new and
silently drops the half-done ones.

So say which earlier action each finding descends from, and add a
`## Resolved since round N` section crediting what landed. The author needs to
see their fixes were noticed, and the next round needs the running tally.

Re-check the findings an earlier round left standing, too: one may have become
moot rather than fixed — for example, a cost concern about code that turns out
to be gated off by a default nobody sets. Record a moot finding as moot, or the
next round raises it again.

# Session Retrospectives

## Offer to capture what a session taught you, every time

When a piece of work reaches a natural end — the review is delivered, the tests
pass, the image is built and pushed — say plainly that it is done, then tell the
user you would like to record what the session taught you and ask whether to go
ahead. Ask every time. Do not skip it because the session felt routine, and do
not write to the repo or to memory without the user agreeing first.

The standing preference is that this conversation happens, so the question put
to the user is *what* to write and *where*, not whether retrospection is
wanted. Keep the ask to a couple of sentences naming the specific lessons you
have in mind; a vague "shall I record some learnings?" wastes a turn because
the user cannot answer it without asking you what you mean.

If nothing from the session is worth keeping, say that instead of inventing
something. A thin retrospective is worse than none, because it dilutes the
files that later sessions actually read.

## What is worth keeping

Keep a lesson when knowing it at the start would have saved real time, and when
it will still be true next session. The strongest candidates:

- **Environment facts not derivable from the code.** A required VM shape, a
  network or IAM constraint, an auth path that works where the documented one
  does not.
- **Mistakes with a generalisable root cause.** Not "I mistyped a flag", but
  "a VM with no external IP looks healthy until the first non-Google download,
  because Private Google Access keeps gcloud working and masks the problem".
- **A procedure that took several attempts to get right**, where the working
  sequence is worth replaying rather than rediscovering.
- **Corrections the user made to how you work**, together with the reason, so
  the next session applies the principle rather than just the instance.
- **Verification techniques that caught something**, such as mutation-testing a
  new test to prove it fails without the fix.

Skip anything the code, the git history, or this file already records; anything
that mattered only inside the conversation; and narrative retellings. A lesson
is a fact plus what to do about it, not a diary entry.

## Where each kind of lesson goes

| Kind of lesson                                                   | Home                          |
|------------------------------------------------------------------|-------------------------------|
| Convention any agent working in this repo should follow          | this file                     |
| Fact about the user, their preferences, or ongoing project state | the agent's memory directory  |
| Repeatable multi-step procedure worth replaying verbatim         | a skill                       |
| Detail relevant only to one open ticket                          | the JIRA ticket, not the repo |

Propose the actual edit — a diff, or the exact text and the file it lands in —
rather than describing it in the abstract. The user is approving specific
wording in a specific place, and a description is not reviewable.

Remember that anything written to this repo is public: the redaction rules in
the Data Handling section apply to retrospectives exactly as they do to code.

# Verification Practices

## A new test is not evidence until you have seen it fail

A test that passes against fixed code tells you almost nothing on its own: it
passes whether it exercises the fix, exercises nothing, or asserts something
vacuously true. Before reporting that new coverage protects a change, break the
change and confirm the test fails — and confirm it fails for the stated reason,
not incidentally.

Temporarily invert or disable the specific line the test is meant to guard,
re-run, check the failure message names the behaviour in question, then restore
the source and re-run to confirm you are back to green. Keep the mutation in
one place at a time; mutating two guards at once cannot tell you whether both
are covered.

This is worth the extra minutes because it routinely finds hollow tests. A
typical one has an assertion that cannot distinguish the two code paths it is
named for — both paths produce the same output for the chosen input, so the
assertion holds either way. The fix is to choose an input where the paths
genuinely diverge. Nothing but a mutation run surfaces that, and such a test
reads perfectly well.

Tests guarding a *negative* ("X should not take path Y") deserve particular
suspicion, since the cheapest way to write one is also the vacuous way.

Where two call sites must agree on a predicate, a single mutation should break
tests covering both. If it breaks only one, either the coverage is incomplete
or the predicate is still duplicated.

## Run a task's `bq` commands locally before its first Terra run

`womtool validate` and `check-wdl-heredocs` treat a command block as a
string, so a wrong query in a task that only runs `bq` and `gcloud` is found
on its first run in Terra. Such a task needs no cluster, so its command can
be run from a laptop first. Extract the block, dedent it, substitute each
`~{...}` with a real value, and run it with bash against a scratch table
prefix in a dataset that expires. A dozen lines of Python do the extraction:
slice from `task <name> {` to `command <<<` and `>>>`, apply
`textwrap.dedent`, then `re.sub` over a dictionary of values, and assert
that no `~{` is left.

Strip `--apilog=false` from the rendered script, for the reason given in
"`bq --apilog=false` writes a file named `false`". Then compare the tables
it builds with ones built another way, since a query that runs is not
thereby correct. In VS-2054 this found a table function that fell back to a
vid's own location where it should have returned NULL, before the workflow
had run anywhere.

This exercises the command only. Input localization, call wiring,
conditionals and the runtime block are still untested until Cromwell runs
the workflow, so say that when reporting.

## `bcftools norm` can read the reference over https

A task that normalizes against the reference, such as `NormalizeKeys` in
`GvsCreateParticipantMappingFromVDS.wdl`, can be reproduced on a laptop
without downloading the 3 GB fasta. bcftools reads it and its `.fai` by
range request from the public bucket:

    REF=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    bcftools norm -m- --check-ref w -f "${REF}" keys.vcf.bgz -O v -o normalized.vcf

Use the `https://storage.googleapis.com/` form, not `gs://`, which needs
credentials and an htslib built with GCS support. In VS-2054 this normalized
395 keys from the AnVIL 3K VDS in seconds (bcftools 1.23.1 from Homebrew).
Check the `Lines total/split/joined/realigned` summary on stderr: the total
should equal the input's row count, with nothing removed or skipped.

## Do not trust a completion signal without checking it is fresh

When polling for a marker that means "the long job finished" — a sentinel file,
a status field, a log line — confirm the marker belongs to the run you are
watching before believing it. A leftover marker from a previous attempt is
indistinguishable from success if you only test for existence.

Check a timestamp, or delete the marker before starting and verify the deletion
actually took effect. A stale exit-code file from a failed first attempt is
easily read as a completed second attempt within seconds of launching it; the
giveaway is a file mtime that predates the run. The same
applies to `pgrep -f <script>`, which happily matches the SSH command string
that contains the script name and reports a process that does not exist.

More generally: when a check returns the answer you were hoping for
suspiciously early, spend one command confirming it means what you think.

The same holds for a report that a Terra run has finished. Check the
submission's status, and the timestamps of its tasks' `rc` files, before
verifying its outputs. In VS-2054 a run reported as finished was still
creating its Dataproc cluster; the finished one was a different submission.

## Check repository state before reporting it

Before telling the user that work is uncommitted, unpushed, or on some branch,
run `git status -sb` and `git log` and report what they show. Do not carry the
claim forward from earlier in the session. The user commits and pushes on their
own, between turns, and a conversation summary records state as of when it was
written, so a statement that was true an hour ago reads just as confidently as
one that is true now.

In the VS-2041 work, "nothing is committed yet" was repeated across several
hand-offs after the user had committed and pushed every change. One command
would have caught it.

## Quote the revision in `git show <rev>:<path>`

In zsh, `$c:scripts/foo.wdl` is parsed as a parameter expansion with a history
modifier, not as a revision followed by a path. `git show $c:path/to/file` in a
loop therefore fails with `(eval):4: bad substitution` for every iteration and
produces no output at all. Write `git show "${c}":path/to/file`.

This is worth its own note because of how the failure reads: the error names
`eval` and a line number, not git, and a loop that greps the output prints a
plausible-looking nothing for each commit. Absence of matches is easy to accept
as a real answer. When a loop over revisions reports that none of them contain
something you expect to be there, check the command ran before believing the
result.

## A space-separated file list in a zsh variable is one word

zsh does not split an unquoted `$VAR` on whitespace, unlike bash. So after
`F="a.wdl b.wdl"`, both `for f in $F` and `tool $F` see a single argument,
`a.wdl b.wdl`. The loop runs once, and a tool that takes files fails with a
"No such file or directory" error naming the whole string. Use an array:
`F=(a.wdl b.wdl)`, then `for f in $F` or `tool $F`.

Worse than the failure, a loop body that pipes through `tail` or `grep` can
print almost nothing, which looks like a clean pass. In VS-2055 a `womtool
validate` loop over seven WDLs printed the joined path once and no verdicts.
Check that the output names every file before reading it as a pass.

## Prove a rebuilt image carries the code, do not infer it

The Variants-image rebuild rule above has a verification step that costs one
command, and whoever is reviewing or authoring should run it rather than
reasoning about dates. Pull the tag named in `GetToolVersions`, extract the
script, and diff it against the branch:

    CID=$(docker create --platform linux/amd64 <variants_docker_tag>)
    docker cp "${CID}":/app/<script>.py ./img.py && docker rm "${CID}"
    git show "<rev>":scripts/variantstore/scripts/<script>.py > ./pr.py
    diff img.py pr.py

This catches what nothing else does: a tag bumped but built from a different
tree, or built before the last edit to the script. The "Detect missing Docker
image tag updates" CI check proves only that the tag string changed when a
baked-in file changed — not that the image behind the new tag holds that file's
current contents. Comparing the image's build date against the script's last
commit date is circumstantial in the same way.

It matters most when a WDL starts passing a new argument to a baked-in script.
A stale image fails with argparse's `error: unrecognized arguments`, which
surfaces only at runtime inside a Cromwell task, long after review.

## Run the Docker refresh check locally, and read its output, not its exit status

The "Detect missing Docker image tag updates" CI check is a wrapper around
`scripts/ci/check_docker_refresh.sh`. Running it locally saves a push-and-wait
cycle:

    bash scripts/ci/check_docker_refresh.sh origin/ah_var_store HEAD

Two things about it are easy to misread.

It always exits 0. The verdict is the summary line it prints, and in CI the
`status` it writes to `GITHUB_OUTPUT`. "Possible missing Docker tag updates
detected" is a failure; "required Docker tag updates are present" and "no
image-affecting changes" are passes. A script that tests `$?` passes every
time.

It diffs two commits (`base...head`), so it cannot see the working tree. A tag
bump that is edited but not yet committed reports as missing. To check one
before committing, run it against a commit object that no branch points to,
built from a scratch index so the real index is untouched:

    IDX=$(mktemp)
    cp "$(git rev-parse --git-path index)" "${IDX}"
    GIT_INDEX_FILE="${IDX}" git add scripts/variantstore/wdl/GvsUtils.wdl
    C=$(git commit-tree "$(GIT_INDEX_FILE="${IDX}" git write-tree)" -p HEAD -m tmp)
    bash scripts/ci/check_docker_refresh.sh origin/ah_var_store "${C}"

Nothing references the commit, so git eventually garbage-collects it.

## Check a rebuilt GATK image pair's shape without pulling it

`building_gatk_docker_on_a_vm.md` gives the layer counts and compressed sizes
of a correct lite/heavy pair. Its first "Gotcha" describes the failure these
numbers catch: the lite image lost or overwritten because both builds share a
tag. The registry reports these figures in seconds, with no need to pull about
3.8 GB:

    docker manifest inspect -v us-central1-docker.pkg.dev/broad-dsde-methods/gvs/gatk:<tag> |
      python3 -c 'import json,sys; d=json.load(sys.stdin); d=d if isinstance(d,list) else [d]
    for e in d:
        m=e.get("OCIManifest") or e.get("SchemaV2Manifest")
        if m and e["Descriptor"].get("platform",{}).get("architecture")!="unknown":
            print(len(m["layers"]),"layers",round(sum(l["size"] for l in m["layers"])/1e9,2),"GB")'

The `unknown` filter skips the BuildKit attestation entry that images built
under the containerd store carry. Compare the result with the tags currently
in `GetToolVersions` as well as with the doc. In VS-2057 both pairs came out at
12 layers / 1.29 GB for lite and 13 layers / 2.49 GB for heavy. A matching
shape shows the right image was pushed under each tag. It does not show that
the image holds the branch's code.
