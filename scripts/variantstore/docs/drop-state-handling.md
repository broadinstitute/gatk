# How `drop_state` is handled across ingest and extract

Observations collected while working on VS-1961 (GVS launch metadata). None of this is a defect in a
delivered callset as far as we can tell, but the handling is surprising in several places and the
correctness of at least one deliverable rests on a coincidence. Collected here to seed follow-on
tickets; see the end of the document for candidates.

## What the value means

`drop_state` names a GQ band that ingest discards. The bands come from `GQStateEnum`
(`src/main/java/org/broadinstitute/hellbender/tools/gvs/common/GQStateEnum.java`): `ZERO` through
`SIXTY` in tens, plus `VARIANT`, `STAR`, `MISSING`, `UNKNOWN`, and `NONE("")`. `NONE` is special —
it carries a null `referenceGQ` and matches no real band, so it discards nothing.

At ingest the value becomes GATK's `--ref-block-gq-to-ignore`. Two details matter:

- **GATK's own default is `SIXTY`** (`CreateVariantIngestFiles.java:63`), and the WDL task input is
  `String? drop_state` passed conditionally, so a task-level caller that omits it drops the GQ60 band.
- **Missing intervals are materialized as GQ0 unless `ZERO` is the sole discarded state**
  (`CreateVariantIngestFiles.java:469-472`). Regions of the interval list with no gVCF data get
  explicit GQ0 reference blocks written for them.

Those two rules together mean the three plausible values produce quite different `ref_ranges`
footprints:

| Value   | Bands discarded | Missing intervals | Result                                         |
|---------|-----------------|-------------------|------------------------------------------------|
| `NONE`  | none            | written as GQ0    | Largest possible `ref_ranges`                  |
| `ZERO`  | GQ0             | not written       | Smallest; absence of a row means GQ0           |
| `FORTY` | GQ40            | written as GQ0    | Gaps mean GQ40, to be inferred back at extract |

## The defaults, and why they differ

| Location                                                                                                                                   | Default   |
|--------------------------------------------------------------------------------------------------------------------------------------------|-----------|
| `GvsJointVariantCalling.wdl:28`                                                                                                            | `"FORTY"` |
| `GvsBulkIngestGenomes.wdl:44`, `GvsImportGenomes.wdl:29`                                                                                   | `"NONE"`  |
| `GvsExtractCallset.wdl:38`, `GvsExtractCallsetPgen.wdl:51`, `GvsExtractCallsetPgenMerged.wdl:43`, `GvsExtractCohortFromSampleNames.wdl:35` | `"NONE"`  |
| GATK `--ref-block-gq-to-ignore` (`CreateVariantIngestFiles.java:63`)                                                                       | `SIXTY`   |

The divergence between entry points is deliberate rather than a wiring defect:
`GvsJointVariantCalling.wdl` threads its own value down to the subworkflows that need it (`drop_state`
at `:150` and `:252`, `maximum_alternate_alleles` at `:256`), so a subworkflow default never applies
when the top-level workflow is the entry point, and each default suits how its own workflow is normally
invoked.

What is harder to justify is the `"NONE"` ingest default specifically. Its stated rationale is the
comment directly above it — `# set to "NONE" to ingest all the reference data into GVS for VDS (instead
of VCF) output` — which dates from before AoU established that GQ0 blocks can be dropped safely. No
caller in this repo relies on it:

| Caller                              | Value                                               |
|-------------------------------------|-----------------------------------------------------|
| `GvsJointVariantCalling.wdl`        | `FORTY`                                             |
| Hail/VDS quickstart integration     | `ZERO`                                              |
| VCF quickstart, benchmark workflows | `FORTY`                                             |
| AoU                                 | `ZERO`, by instruction at `AOU_DELIVERABLES.md:108` |

So the least appropriate of the three values — largest storage footprint, no consumer — is what an
operator gets by omission on exactly the invocation path AoU uses.

## At extract, `drop_state` is an assertion rather than a choice

Extract consumes the value as `--inferred-reference-state`, and
`ExtractCohortEngine.java:687` stamps that state's GQ onto every reference genotype extract synthesizes
for a sample with no row at a site. The value is therefore a claim about what ingest discarded, and
delivered genotype qualities depend on the claim being true:

| Ingested with | Extract told     | Result                                                                                  |
|---------------|------------------|-----------------------------------------------------------------------------------------|
| `ZERO`        | `FORTY`          | GQ40 confident reference calls invented where no gVCF supported any confidence          |
| `FORTY`       | `ZERO`           | GQ0 where GQ40 was intended; conservative, still wrong                                  |
| anything      | argument omitted | GATK's default `SIXTY` (`ExtractCohort.java:208`), the most confident value in the enum |

Nothing checks this. Under `GvsJointVariantCalling.wdl` the two ends cannot diverge because one input
feeds both, but they can for AoU, for sub-cohort extracts, and for any direct `GvsExtractCallset` run —
the cases where the value is supplied by hand, possibly months after ingest.
`SUB_COHORT_WORKFLOW.md:27` states the rule ("This should correspond to the same value set in
`GvsImportGenomes` (or `GvsJointVariantCalling`)"), so the invariant is known; it simply has nowhere to
live and nothing to enforce it.

`"NONE"` means the same thing at both ends — no band was discarded — and extract remaps it to `"ZERO"`
(`GvsExtractCallset.wdl:440`, `GvsExtractCallsetPgen.wdl:394`) because `GQStateEnum.NONE` has a null
`referenceGQ` and extract needs a concrete GQ to stamp. The remap is consistent with GVS's convention
that absent reference data reads as GQ0, which `GQStateEnum.java:16-18` documents directly. Note also
that the argument is interpolated from a non-optional, never-empty String, so these WDLs always pass it
and GATK's `SIXTY` default cannot apply through them.

## How the current arrangement came about

- `0ff86b12b` (2022-08-24, VS-607, "Change drop_state to NONE for Ingest/Extract") set the `"NONE"`
  default on both the ingest and extract sides at once. This is the origin of both, and it was a change
  about ingesting full reference data for VDS output rather than a decision about extract inference.
- `9cf80212d` (2023-01-23, VS-772, a sub-cohort prepare fix) added the `"NONE"` → `"ZERO"` remap.
  Before it, the WDL passed `--inferred-reference-state NONE` directly.

In today's code that pre-remap combination would fail loudly rather than quietly: `NONE` has a null
`referenceGQ` and `ExtractCohortEngine.java:687` passes it straight to `GenotypeBuilder.GQ(int)`, so it
would NPE on unboxing. That has not been verified against the engine code as it stood in January 2023,
which is the difference between "that window could not produce output" and "that window needs
auditing".

## The AoU small callsets

The "Smaller Interval Lists" and "Smaller Sample Lists" sections of `AOU_DELIVERABLES.md` do not mention
`drop_state`, so those extracts ran with the workflow default. The outcome is correct — extract logs
confirm `--inferred-reference-state ZERO` — but only through a three-part coincidence:

1. AoU ingest uses `"ZERO"` because `AOU_DELIVERABLES.md:108` tells an operator to type it in.
2. The extract workflows default to `"NONE"`, set in VS-607 for reasons unrelated to extract inference.
3. VS-772 remapped `"NONE"` to `"ZERO"`, which is what ZERO-ingested data happens to require.

Nothing in the system knows those three facts agree, and no artifact records what ingest used, so the
agreement cannot be checked — only re-derived from git history and the data.

The data-side check is cheap if it is ever wanted: query the callset's `ref_ranges_%` tables for
GQ0-state rows, whose absence confirms `ZERO` ingest. Where `use_compressed_references` is set the state
is packed into `packed_ref_data`, so this needs the unpack logic
(`SchemaUtils#encodeCompressedRefBlock` / the `UnpackRefRangeInfo` BigQuery function).

Read such a query carefully, because absence is weaker evidence than it looks. GVS bins reference GQ into
states in tens (`RefRangesCreator.java:227-243`), but the reblocking GVS expects (`-GQB 20 -GQB 30 -GQB 40
--floor-blocks`, `ReblockGVCF.java:73`) can only produce states `0`, `2`, `3` and `4`. States `1`, `5` and
`6` are therefore always absent regardless of `drop_state`, so the diagnostic is "which one of `0`, `2`,
`3`, `4` is missing", not "which state is missing". A `FORTY`-ingested callset looks like this — states
`0`, `2` and `3` populated, `4` empty:

```
gvs-internal.quickit_2026_07_21_..._vets_vcf.ref_ranges_001
  state 0: 2764630
  state 2: 2495387
  state 3: 4773472
  state 4: absent   <- drop_state = FORTY
```

Two further caveats: a callset small enough that a band is genuinely unpopulated will read as a false
positive, and `--ignore-above-gq-threshold` (`dropAboveGqThreshold`, default `false`) drops bands above the
named one as well, so more than one absent band does not imply more than one `drop_state`.

## Candidate follow-on work

1. **Remove the `drop_state` defaults from `GvsBulkIngestGenomes.wdl` and `GvsImportGenomes.wdl`**, making
   the input required so that an operator must state it. No in-repo caller relies on the default, and this
   converts a silent, expensive-to-repair mistake into an unsatisfied-input error at submit time. Leave the
   extract-side defaults alone in this change; see item 3.
2. **Document `drop_state` in the "Smaller Interval Lists" and "Smaller Sample Lists" sections of
   `AOU_DELIVERABLES.md`**, matching what `SUB_COHORT_WORKFLOW.md:27` already says.
3. **Have extract look the value up rather than ask for it.** Depends on launch metadata being recorded
   (VS-1961); see `gvs-launch-metadata-proposal.md`, where this is the worked example for the enforcement
   half of that work. The input would remain as an override that fails when it disagrees with the record,
   with band-absence inference as the fallback for datasets predating the metadata table.
4. **Decide whether the pre-remap window (2022-08-24 to 2023-01-23) needs auditing**, which reduces to
   confirming whether `NONE` could reach `GenotypeBuilder.GQ` in the engine code of that period.
