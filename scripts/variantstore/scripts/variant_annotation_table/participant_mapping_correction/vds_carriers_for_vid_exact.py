#!/usr/bin/env python3
"""Carriers of a VID in the delivered VDS, at the VID's own coordinate only.

Built to answer the question collaborators actually ask, which is not about left alignment:
"the mapping table says 937 for 13-32368001-C-CTT but we count 1,063 in the VDS, why is the
mapping table LOWER?"

Note the direction. Collaborators are not currently applying FT when they count the VDS, so
they get every LGT-defined genotype calling the allele. The 2026-09-05 mapping table applies
FT and drops GQ 0, so it is now SMALLER than their count -- 126 hetvar genotypes that FT-fail
on their OTHER allele, all correct exclusions. This reverses the direction of the older
complaint, where the mapping table was the larger number (1,067 against their 1,063, the
difference being 4 GQ 0 records the VDS cannot return because LGT is missing for them). A
report that explains the gap in the old direction will not match what they now see.

Separately the VAT's gvs_all_ac for this VID is 1,289, which is neither wrong nor comparable
-- it counts alleles, not people. See below.

This is the simplified sibling of `vds_carriers_for_vid.py`. That script searches a 200bp
rightward window for non-left-aligned synonyms; this one looks ONLY at the locus the VID
names and reports nothing if the variant is not there. That is the right trade for the
common case: non-left-aligned VIDs are 5,776 of 1,601,242,198 in Foxtrot (0.00036%), so a
count discrepancy is overwhelmingly likely to be genotype filtering rather than a
representation mismatch. If this script finds no matching allele at the VID's locus, that is
the signal to reach for the windowed one -- see NO MATCH below.

WHAT "CARRIER" MEANS HERE, and why the number is lower than the mapping table's

  1. NOT A GQ 0 NO-CALL. `import_gvs.py:262` nulls `LGT` when GQ is 0, so `hl.is_defined(LGT)`
     already excludes them and no `GQ > 0` test is needed. `GQ` itself stays defined for these
     -- it is set unconditionally from the vet record at `import_gvs.py:264` -- so testing `GQ`
     rather than `LGT` gets this backwards. That same asymmetry is what lets the GQ 0 count
     below be computed at all: LGT missing AND GQ defined means "had a vet record, called
     nothing", while both missing means "no record at this site".

     GQ 0 IS STILL ATTRIBUTABLE TO AN ALLELE, which is worth stating because the natural
     assumption is that it is not. `LGT` is the only field nulled at GQ 0. `LA` is built from
     the vet row's own `ref`/`alt` columns (`import_gvs.py:259`, `:286`) and never consults
     `GQ`, so it survives and still names which alleles that sample's record was about. A site
     with many alt alleles will have far more GQ 0 no-calls than belong to any one of them --
     `LA.contains(target_idx)` is what separates them, and it is reported below. That is also
     precisely the criterion by which the mapping table counted these participants: the table
     is built from `alt_allele`, which splits a hetvar vet row into one row per allele, so a
     GQ 0 record listing two alts is counted under both. LA membership and mapping-table
     inclusion are therefore the same test, which is why this count reconciles.

     The one thing that genuinely is unrecoverable is the GENOTYPE: `LA` says the record
     concerned this allele, not how many copies a call would have named. So these participants
     can be attributed but not counted as carriers, which is why they are reported separately
     rather than folded into the carrier total.

  2. NOT FT-FAILING. `FT` is a boolean entry field, True = pass, written by
     `merge_and_rescore_vdses.py:173` (and `import_gvs.py:380`). A VDS built before that step
     does not carry it, hence the guard below.

  3. ACTUALLY CARRYING THIS ALLELE. `LGT` is LOCALLY indexed and `LA` maps local -> global, so
     the called non-ref global indices are the fold at `import_gvs.py:369-371`. Comparing
     `LGT` against a global index without going through `LA` is silently wrong.

WHY FT REMOVES HETVAR CARRIERS WHOSE OWN ALLELE IS FINE

`merge_and_rescore_vdses.py:155-173` folds over ALL called non-ref alleles of the genotype:

    FT = ~any_no & (any_yes | all_ok)

so FT is a property of the GENOTYPE, not of an allele. For a `1/2` hetvar that means the call
is filtered if EITHER allele is NO, or if neither is YES and not both are OK. A participant
who carries this VID's allele -- perfectly good on its own -- is therefore dropped when their
OTHER allele fails. This is the single most likely explanation for a mapping-table-vs-VDS gap
larger than the handful of GQ 0 samples, and the `ft_failing_hetvar` count below isolates it.
The converse also occurs and is worth knowing before anyone argues the rule is too strict:
`any_yes` rescues a genotype whose bad allele is partnered with a YES one.

WHY min_rep IS STILL NEEDED even though we only look at one locus

The VDS is unsplit multi-allelic and pads `alleles[0]` to span the longest deletion at the
site, so the VID's ALT string usually is NOT the string sitting in `alleles`. This is about
padding, not left alignment, and it does not go away by pinning the locus. Compare the
minimized (REF, ALT) PAIR, never the ALT alone -- 13-32368001 is itself the cautionary case:
at that site `CTT` is both the raw ALT of a 24bp deletion against a 27bp REF and the minimized
ALT of a 2bp insertion, so matching on `CTT` alone silently picks the wrong variant.

gvs_all_ac IS NOT A PARTICIPANT COUNT

The VAT's `gvs_all_ac` counts ALLELES. A `1/1` carrier contributes 2, a `0/1` or `1/2` carrier
contributes 1, so at any site with homvar calls the AC exceeds the number of people and always
will. That is not a discrepancy and no amount of filtering work will close it. The mapping
table counts PEOPLE. Anyone reconciling one against the other without converting units is
comparing two correct numbers and concluding something is broken. This script reports both,
from the same carrier set, so the ratio is visible as the homvar burden it is.

NO WITHDRAWN/CONTROL FILTERING is applied: the delivered VDS already excludes those samples,
so its columns are the active non-control cohort and match the mapping table's population.
"""

import argparse

import hail as hl


def parse_vid(vid):
    """'13-32368001-C-CTT' -> ('chr13', 32368001, 'C', 'CTT'). VAT chromosomes are bare, VDS
    contigs are chr-prefixed, and X/Y stay alphabetic on both sides."""
    chrom, pos, ref, alt = vid.split('-')
    return f'chr{chrom}', int(pos), ref, alt


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vds-path', required=True)
    p.add_argument('--vid', required=True, help="e.g. 13-32368001-C-CTT")
    p.add_argument('--output', help="optional TSV, one row per carrier")
    p.add_argument('--reference-genome', default='GRCh38')
    args = p.parse_args()

    rg = args.reference_genome
    hl.init()

    contig, pos, ref, alt = parse_vid(args.vid)
    delta = len(alt) - len(ref)
    target_locus = hl.eval(hl.locus(contig, pos, reference_genome=rg))
    # min_rep the TARGET as well, not just the VDS side. A vid is minimized in principle, but
    # normalizing both sides costs one `hl.eval` and removes the assumption.
    target_mr = hl.eval(hl.min_rep(hl.locus(contig, pos, reference_genome=rg), [ref, alt]))

    vds = hl.vds.read_vds(args.vds_path)
    vd = vds.variant_data

    if 'FT' not in vd.entry:
        raise SystemExit(
            'This VDS has no FT entry field, so "not filtered" is unanswerable. FT is '
            'written by merge_and_rescore_vdses.py:173; a VDS produced before that step '
            'will not carry it.')

    # ---------------------------------------------------------------------------------------
    # The VID's locus, and nothing else.
    # ---------------------------------------------------------------------------------------
    # Built as an explicit closed point interval rather than parsed from a string: the string
    # form is end-EXCLUSIVE by default, so `chr13:32368001-32368001` is the empty interval and
    # would match nothing. The parent script never hit this because its 200bp window made the
    # interval non-degenerate.
    vd = hl.filter_intervals(
        vd, [hl.Interval(target_locus, target_locus,
                         includes_start=True, includes_end=True)])

    # ---------------------------------------------------------------------------------------
    # Which global alt index at this site IS this VID?
    # ---------------------------------------------------------------------------------------
    # Indel length is invariant under trimming, so this prefilter discards nothing real and
    # runs before the expensive min_rep rather than after.
    vd = vd.annotate_rows(
        target_idx=(hl.range(1, hl.len(vd.alleles))
                    .filter(lambda i: (hl.len(vd.alleles[i]) - hl.len(vd.alleles[0])) == delta)
                    .filter(lambda i: hl.bind(
                        lambda mr: (mr.locus == target_mr.locus)
                        & (mr.alleles[0] == target_mr.alleles[0])
                        & (mr.alleles[1] == target_mr.alleles[1]),
                        hl.min_rep(vd.locus, [vd.alleles[0], vd.alleles[i]])))
                    .first()))

    # `alleles` is part of the rows-table key, so it comes back from collect() on its own and
    # naming it here would be rejected as an attempt to overwrite a key field.
    hits = vd.rows().select('target_idx').collect()
    matched = [h for h in hits if h.target_idx is not None]

    if not matched:
        print(f'\nNO MATCH for {args.vid} at {contig}:{pos}.')
        for h in hits:
            print(f'  site present with alleles {list(h.alleles)}, none equal to '
                  f'{target_mr.alleles} after min_rep')
        if not hits:
            print('  no VDS row at this locus at all')
        print('\nThis VID may be one of the 5,776 whose stored representation is not '
              'left-aligned. Re-run vds_carriers_for_vid.py --window 200 before concluding '
              'the variant is absent.')
        return

    row = matched[0]
    print(f'\n{args.vid} -> {contig}:{pos}  {ref}/{alt}')
    print(f'  site alleles: {list(row.alleles)}')
    print(f'  matched global alt index: {row.target_idx}')

    vd = vd.filter_rows(hl.is_defined(vd.target_idx))

    # ---------------------------------------------------------------------------------------
    # The carriers, and the decomposition that explains the mapping table's larger number.
    # ---------------------------------------------------------------------------------------
    # Built on the entries Table rather than with `filter_entries`, which only sets entries to
    # missing. One row, so this is n_samples rows.
    ent = vd.entries()
    called = (hl.range(ent.LGT.ploidy)
              .map(lambda i: ent.LA[ent.LGT[i]])
              .filter(lambda x: x != 0))
    # Materialized as fields rather than kept as loose expressions, so that the filtered table
    # built for --output below can refer to them by name. Expressions built against `ent`
    # cannot be used on the Table that `ent.filter(...)` returns.
    ent = ent.annotate(
        has_allele=hl.is_defined(ent.LGT) & called.contains(ent.target_idx),
        # Two DISTINCT non-ref alleles called, i.e. a 1/2-style genotype. `1/1` is not hetvar.
        hetvar=hl.len(hl.set(called)) > 1,
        n_copies=called.filter(lambda i: i == ent.target_idx).length(),
        called_alt_indices=hl.delimit(called.map(hl.str), ','),
        # A vet record that called nothing. Both fields come from the same vet row, but only
        # LGT is nulled at GQ 0 (import_gvs.py:262) while GQ is not (:264), so the pair
        # separates "had a record, called nothing" from "no record at this site".
        gq0=hl.is_missing(ent.LGT) & hl.is_defined(ent.GQ))

    s = ent.aggregate(hl.struct(
        with_allele=hl.agg.count_where(ent.has_allele),
        carriers=hl.agg.count_where(ent.has_allele & ent.FT),
        ft_failing=hl.agg.count_where(ent.has_allele & ~ent.FT),
        ft_failing_hetvar=hl.agg.count_where(ent.has_allele & ~ent.FT & ent.hetvar),
        ft_missing=hl.agg.count_where(ent.has_allele & hl.is_missing(ent.FT)),
        # The VAT's gvs_all_ac counts ALLELES, not people: a 1/1 carrier contributes 2. Summed
        # over exactly the carrier set, so `ac` and `carriers` are the same population in
        # different units and their ratio is the homvar burden, not a discrepancy.
        ac=hl.agg.filter(ent.has_allele & ent.FT, hl.agg.sum(ent.n_copies)),
        carriers_1_copy=hl.agg.count_where(ent.has_allele & ent.FT & (ent.n_copies == 1)),
        carriers_2_copies=hl.agg.count_where(ent.has_allele & ent.FT & (ent.n_copies == 2)),
        carriers_hetvar=hl.agg.count_where(ent.has_allele & ent.FT & ent.hetvar),
        # Counted independently of FT rather than as pass+fail, so that an FT-missing entry
        # shows up as the two not summing rather than being quietly absorbed.
        hetvar_total=hl.agg.count_where(ent.has_allele & ent.hetvar),
        gq0_at_site=hl.agg.count_where(ent.gq0),
        gq0_this_allele=hl.agg.count_where(ent.gq0 & ent.LA.contains(ent.target_idx)),
    ))

    print(f'\n  genotypes calling this allele (LGT defined)   {s.with_allele:>8,}')
    print(f'    FT pass  -> CARRIERS                        {s.carriers:>8,}')
    print(f'    FT fail  -> excluded                        {s.ft_failing:>8,}')
    print(f'      of which hetvar (two distinct alts)       {s.ft_failing_hetvar:>8,}')
    if s.ft_missing:
        print(f'    FT MISSING -> unexpected, investigate       {s.ft_missing:>8,}')
    print(f'      1 copy  (0/1 or 1/2)                    {s.carriers_1_copy:>8,}')
    print(f'      2 copies (1/1)                          {s.carriers_2_copies:>8,}')
    print(f'    ALLELE COUNT over those carriers          {s.ac:>8,}')
    print("      compare to the VAT's gvs_all_ac; alleles, not people")

    # Broken out because "we dropped 126 people over an allele you did not ask about" invites
    # "so is hetvar filtering just broken?", and the passing count is the answer to that.
    if s.hetvar_total:
        pass_pct = 100.0 * s.carriers_hetvar / s.hetvar_total
        fail_pct = 100.0 * s.ft_failing_hetvar / s.hetvar_total
        print(f'\n  HETVAR (1/2) genotypes calling this allele    {s.hetvar_total:>8,}')
        print(f'    FT pass                                   {s.carriers_hetvar:>8,}'
              f'   ({pass_pct:.1f}%)')
        print(f'    FT fail                                   {s.ft_failing_hetvar:>8,}'
              f'   ({fail_pct:.1f}%)')
        if s.carriers_hetvar + s.ft_failing_hetvar != s.hetvar_total:
            print('    NB: pass + fail != total, so some entries have FT missing')

    print(f'\n  GQ 0 no-calls at this site                    {s.gq0_at_site:>8,}')
    print(f'    whose vet record listed THIS allele         {s.gq0_this_allele:>8,}')
    print(f'    listing only other alleles at the site      '
          f'{s.gq0_at_site - s.gq0_this_allele:>8,}')
    print('    (LGT is nulled at GQ 0, but LA is not, so the allele is still attributable;')
    print('     LA membership is the same criterion the mapping table counted them by)')

    # Reconciliation against the pre-fix mapping table. That table was built from `alt_allele`,
    # which holds a row per called allele per sample regardless of GQ or FT, so its count for
    # this VID should be the LGT-defined callers plus the GQ 0 records that listed the allele.
    # Printed rather than asserted: if the collaborator's number is something else, the gap is
    # itself the finding. The one known way this identity can slip is a GQ 0 record whose vet
    # GT names no allele -- LA would still count it here while alt_allele would not.
    implied_old = s.with_allele + s.gq0_this_allele
    print(f'\n  RECONCILIATION (people, not alleles)')
    print(f'    a VDS count that ignores FT, as collaborators')
    print(f'      currently do                             {s.with_allele:>8,}')
    print(f'    pre-fix mapping table should have said     {implied_old:>8,}')
    print(f'    2026-09-05 mapping table / VDS carriers    {s.carriers:>8,}')
    print(f'    removed                                    '
          f'{implied_old - s.carriers:>8,}'
          f'   ({s.ft_failing:,} FT-failing + {s.gq0_this_allele:,} GQ 0)')
    # The fix moved the mapping table below the collaborators' own number, so the gap they see
    # now has the opposite sign from the one they originally reported. Said explicitly because
    # an explanation written for the old direction will not match what they are looking at.
    if s.carriers < s.with_allele:
        print(f'    -> the mapping table is now {s.with_allele - s.carriers:,} LOWER than an '
              f'FT-ignoring VDS count,')
        print(f'       having been {implied_old - s.with_allele:,} higher before the fix')

    print(f'\n{args.vid}: {s.carriers} carrier(s) with an unfiltered, non-GQ-0 genotype.')
    if s.ft_failing_hetvar:
        print(f'Note for collaborators: {s.ft_failing_hetvar} of the {s.ft_failing} FT-failing '
              'genotypes are hetvar. FT is per-genotype, not per-allele '
              '(merge_and_rescore_vdses.py:155-173), so those participants are excluded on '
              "account of their OTHER allele, not this one.")

    if args.output:
        out = ent.filter(ent.has_allele)
        out = out.select('GQ', 'FT', 'n_copies', 'hetvar', 'called_alt_indices')
        out.export(args.output)
        print(f'\nWrote {args.output}')


if __name__ == '__main__':
    main()
