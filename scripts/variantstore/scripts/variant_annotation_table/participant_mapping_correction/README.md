# Introduction

This directory holds the artifacts from a one-off correction of the Foxtrot VID-to-participant mapping table, which
over-reported carriers. Nothing here runs as part of a workflow. It is checked in as a record of what was measured, what
was changed, and why — the correction was applied directly to a delivered BigQuery table, so the SQL is the only
description of what happened to it.

The production fix, so that a future callset does not need this treatment, is tracked separately (VS-2013) and belongs in
`GvsCreateParticipantMappingTable.wdl`, `GvsMapUnmappedVIDs.wdl` and `GvsMapDroppedDuplicateVIDs.wdl`, not here.

# What was wrong

The mapping table ships alongside the VAT, and the AoU Researcher Workbench uses it to decide which samples to pull into
an extracted VCF. A surplus person ID therefore pulls that participant's sample into extracts for a variant the callset
does not count them as carrying.

The base join in `GvsCreateParticipantMappingTable` reads `alt_allele`, which holds every call ever ingested, while the
VAT counts only genotypes that survive VDS import and the filter model. Two classes of genotype fall into that gap:

  - **GQ 0 no-calls.** `import_gvs.py` nulls the genotype when `GQ` is 0, so the person has no call in the VDS and
    contributes nothing to `AC`. The call still exists in `alt_allele` and was still mapped.

  - **Genotypes failing `FT`.** A genotype is judged by its *worst* non-ref allele and a filtered genotype stops
    contributing to `AC`. This covers `1/2` het-var calls filtered because *either* allele failed — including when the
    VID's own allele passes and the partner allele is the one that failed — and `0/1` or `1/1` calls whose one allele
    failed the filter model.

A third defect, withdrawn and control samples appearing in the mappings, was fixed in the WDLs by `11488f657` and
scrubbed from the delivered table separately. That scrub is a prerequisite for this correction, and its artifacts are
here too.

The target was not a matter of judgment: `person_ids` should list exactly the participants whose genotypes contribute to
that VID's `gvs_all_ac`.

Measured genome-wide on the delivered table: 2,341,606,369 over-reported entries out of 2,559,382,818,198 (0.0915%),
across 5,437,323 of 1,601,242,198 VIDs (0.34%), averaging 431.2 entries per affected VID. The damage is concentrated
rather than diffuse, which is why the table was patched in place rather than regenerated.

# The files, in the order they were run

| File                                           | What it does                                                                                                                                                            |
|------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `scrub_participant_mapping_table.sql`          | The earlier withdrawn/control scrub of the delivered table. Also documents the four-case taxonomy of representation mismatch that the rest of this work refers back to. |
| `verify_participant_mapping_scrub.sql`         | Post-hoc validation of that scrub.                                                                                                                                      |
| `participant-mapping-withdrawn-fix-summary.md` | Write-up of the scrub, including the two VIDs it emptied and the hand-written `INSERT`s that restored them.                                                             |
| `step4_compare_vds_ft.py`                      | Validates the SQL reconstruction of Hail's `FT` against the VDS itself rather than against arithmetic. Run on a Hail cluster.                                           |
| `step5_build_exclusions.sql`                   | The single full pass over `alt_allele` that builds the exclusion set: 2,344,756,787 `(vid, person_id)` pairs to remove.                                                 |
| `step6_patch_mapping_table.sql`                | Builds the corrected copy of the mapping table, plus the rehearsals, cost modeling and validation queries around it.                                                    |
| `participant-mapping-over-reporting-plan.md`   | The plan of record. Every measurement above is derived and cross-checked here, with results recorded inline as each step ran.                                           |

Supporting scripts, used for spot checks against the VDS rather than as pipeline stages:

| File                                | What it does                                                                             |
|-------------------------------------|------------------------------------------------------------------------------------------|
| `vds_carriers_for_vid.py`           | Lists the VDS carriers of one VID.                                                       |
| `vds_carriers_for_vid_exact.py`     | As above, but matching the allele representation exactly rather than after minimization. |
| `vds_carriers_report.py`            | Per-VID carrier report used to reconcile the mapping table against the VDS.              |
| `explain_vid_for_collaborators.sql` | Explains a single VID's carrier set end to end, for answering collaborator questions.    |

# Caveats

**Participant IDs are redacted.** Three AoU person IDs appeared in comments and prose; in these committed copies they are
`<PERSON_A>`, `<PERSON_B>` and `<PERSON_C>`. This repository is public. The substitution is confined to comments — no
executable statement referenced a participant ID — but a couple of the commented-out diagnostic queries will need the
real IDs put back before they can be re-run.

**The plan document records the reasoning as it developed, including corrections.** Several figures were revised as
measurements replaced estimates, most notably a 2.77× revision to the prevalence estimate after four defects were found in
the `FT` reconstruction. Superseded values are marked where they appear rather than deleted, so that the reasoning stays
followable — read the surrounding note before quoting a number out of it.

**Two residuals are known and were accepted for this delivery**, both recorded in the plan document and on the tickets:

  - 5,776 VIDs repaired by `GvsMapUnmappedVIDs` or `GvsMapDroppedDuplicateVIDs` are carved out of the correction, because
    their `person_ids` were built at a different allele representation than the exclusions were computed at. Applying one
    to the other could remove genuine carriers. 121 of those VIDs have exclusions at all, totaling 37,702 unadjudicated
    pairs (1.6e-05 of the correction). See VS-2012.

  - 172 VIDs lost every named participant and were dropped, which is the whole of the 172-row difference between the
    delivered 1,601,242,198 and the corrected 1,601,242,026. They had been mapped entirely to non-carriers, and the fixup
    workflows did not detect them. See VS-2014.

**Upstream artifacts for the earliest steps are not checked in here** — `build_allele_status.sql` (Step 1),
`investigate_missing_filter_rows.sql` (Step 2) and `step3_rescore_chr21.sql` (Step 3), along with the Step 4 and Step 5
reconciliation queries. Their results are recorded in the plan document. Add them if the record needs to be complete
enough to re-run from scratch.
