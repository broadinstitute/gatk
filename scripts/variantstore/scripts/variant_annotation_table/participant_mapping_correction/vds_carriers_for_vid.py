#!/usr/bin/env python3
"""Participants who genuinely carry a VID in the delivered VDS.

The Hail-side counterpart to `explain_vid_for_collaborators.sql`. That query decomposes the
mapping table's `delivered` into `corrected + gq0_only + ft_only + both` from BigQuery; this
one recovers `corrected` independently, from the VDS, so the two can be shown to agree
without either being derived from the other.

"Genuinely carries" means three things at once, and the first two are cheaper than they look:

  1. NOT A GQ 0 NO-CALL. `import_gvs.py:262` nulls `LGT` when GQ is 0, so `hl.is_defined(LGT)`
     already excludes them and no `GQ > 0` test is needed. Note that `GQ` itself stays defined
     for these -- it is set unconditionally from the vet record at `import_gvs.py:264` -- so
     testing `GQ` rather than `LGT` gets this backwards.
  2. NOT FT-FAILING. `FT` is a boolean entry field, True = pass, written by
     `merge_and_rescore_vdses.py:173` (and `import_gvs.py:381`). A VDS built before that step
     does not carry it, hence the guard below.
  3. ACTUALLY CARRYING THIS ALLELE. `LGT` is LOCALLY indexed and `LA` maps local -> global, so
     the called non-ref global indices are the fold at `import_gvs.py:369-371`. Comparing
     `LGT` against a global index without going through `LA` is silently wrong.

No withdrawn/control filtering is applied: the delivered VDS already excludes those samples
(see the shortfall note in `step4_compare_vds_ft.py`), so its columns are the active
non-control cohort of 535,662 and match the mapping table's population as delivered.
"""

import argparse

import hail as hl


def parse_vid(vid, rg):
    """'4-76917565-A-AG' -> (contig, pos, ref, alt). VAT chromosomes are bare, VDS contigs
    are chr-prefixed, and X/Y stay alphabetic on both sides."""
    chrom, pos, ref, alt = vid.split('-')
    return f'chr{chrom}', int(pos), ref, alt


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vds-path', required=True)
    p.add_argument('--vid', required=True, help="e.g. 4-76917565-A-AG")
    p.add_argument('--output', help="optional TSV of carrier sample names")
    p.add_argument('--window', type=int, default=0,
                   help="bp to search RIGHTWARD from the vid position for an equivalent "
                        "representation. 0 (default) means exact locus only. Use ~200 to "
                        "hunt a case-4 mismatch; see WHY THE WINDOW IS RIGHTWARD below.")
    p.add_argument('--reference-genome', default='GRCh38')
    args = p.parse_args()

    rg = args.reference_genome
    hl.init()

    contig, pos, ref, alt = parse_vid(args.vid, rg)
    delta = len(alt) - len(ref)
    # min_rep the TARGET as well, not just the VDS side. A vid is minimized in principle, but
    # normalizing both sides costs one `hl.eval` and removes the assumption. Same as
    # check_vds_for_alleles.py:128.
    target_mr = hl.eval(hl.min_rep(hl.locus(contig, pos, reference_genome=rg), [ref, alt]))

    vds = hl.vds.read_vds(args.vds_path)
    vd = vds.variant_data

    if 'FT' not in vd.entry:
        raise SystemExit(
            'This VDS has no FT entry field, so "not filtered" is unanswerable. FT is '
            'written by merge_and_rescore_vdses.py:173; a VDS produced before that step '
            'will not carry it.')

    # -------------------------------------------------------------------------------------
    # Narrow to the site.
    # -------------------------------------------------------------------------------------
    # WHY THE WINDOW IS RIGHTWARD. A vid is minimized AND left-aligned; `alt_allele` and the
    # VDS are not left-aligned. Left-aligned is the LEFTMOST equivalent form, so an equivalent
    # non-left-aligned representation can only sit at or to the RIGHT of the vid's position.
    # The production synonym search keys on the same asymmetry -- bcftools ILEN over a 200bp
    # window rightward (`generate_bcftools_searches_for_variant_synonyms.py:38-45`).
    interval = hl.parse_locus_interval(
        f'{contig}:{pos}-{pos + max(args.window, 0)}', reference_genome=rg)
    vd = hl.filter_intervals(vd, [interval])

    # -------------------------------------------------------------------------------------
    # Which global alt index at each surviving row IS this vid?
    # -------------------------------------------------------------------------------------
    # `hl.min_rep` strips the padding the VDS adds so that alleles[0] spans the longest
    # deletion at the site, reducing each (ref, alt_i) pair to its intrinsic form. Compare the
    # PAIR, never the alt alone: one string can be the raw ALT of one variant and the
    # minimized ALT of an unrelated one at the same site.
    #
    # WHAT THIS WILL NOT CATCH. min_rep trims shared prefixes and suffixes; it does not
    # left-align. It cannot walk an indel through a homopolymer or an Alu repeat, so it will
    # not unify a left-aligned vid with a non-left-aligned VDS allele. That is the case-4
    # mismatch behind the 172 emptied VIDs, and a `--window 0` run returning zero carriers for
    # a VID whose VAT `gvs_all_ac` is nonzero is its signature. Re-run with `--window 200` and
    # read `equivalent_by_ilen` below.
    # Indel length is invariant under trimming -- min_rep removes the same prefix and suffix
    # from both alleles -- so this prefilter discards nothing real, and it runs before the
    # expensive min_rep rather than after. Position and allele strings are NOT invariant; do
    # not prefilter on those. Ordering and idiom follow check_vds_for_alleles.py:141-152.
    vd = vd.annotate_rows(
        equivalent_by_ilen=(hl.range(1, hl.len(vd.alleles))
                            .filter(lambda i: (hl.len(vd.alleles[i])
                                               - hl.len(vd.alleles[0])) == delta)))
    vd = vd.annotate_rows(
        target_idx=(vd.equivalent_by_ilen
                    .filter(lambda i: hl.bind(
                        lambda mr: (mr.locus == target_mr.locus)
                        & (mr.alleles[0] == target_mr.alleles[0])
                        & (mr.alleles[1] == target_mr.alleles[1]),
                        hl.min_rep(vd.locus, [vd.alleles[0], vd.alleles[i]])))
                    .first()))

    hits = vd.rows().select('target_idx', 'equivalent_by_ilen').collect()
    exact = [h for h in hits if h.target_idx is not None]
    print(f'\n{len(hits)} row(s) in the window; {len(exact)} match {args.vid} exactly '
          f'after min_rep.')
    if not exact:
        near = [(h.locus, h.equivalent_by_ilen) for h in hits if h.equivalent_by_ilen]
        print('NO EXACT MATCH. Candidates by indel length (verify by hand, these are not '
              f'matches): {near if near else "none"}')
        if args.window == 0:
            print('Window was 0. Re-run with --window 200 before concluding anything.')
        return

    vd = vd.filter_rows(hl.is_defined(vd.target_idx))

    # -------------------------------------------------------------------------------------
    # The carriers.
    # -------------------------------------------------------------------------------------
    # Built on the entries Table rather than with `filter_entries`, which only sets entries to
    # missing -- `entries()` would still emit one row per (row, col) pair. Streamed, so the
    # product is never materialized. Same reasoning as step4_compare_vds_ft.py.
    ent = vd.entries()
    called = (hl.range(ent.LGT.ploidy)
              .map(lambda i: ent.LA[ent.LGT[i]])
              .filter(lambda x: x != 0))
    ent = ent.filter(hl.is_defined(ent.LGT)        # excludes GQ 0 no-calls
                     & ent.FT                      # excludes FT-failing genotypes
                     & called.contains(ent.target_idx))

    ent = ent.select('GQ', 'FT', n_copies=called.filter(
        lambda i: i == ent.target_idx).length())

    n = ent.count()
    print(f'\n{args.vid}: {n} carrier(s) with an unfiltered, non-GQ-0 genotype.')
    print('Compare against `corrected` in foxtrot.mapping_correction_audit for this vid. '
          'Expect agreement; `delivered` there is the FT-unaware number a collaborator sees.')

    if args.output:
        ent.export(args.output)
        print(f'Wrote {args.output}')
    else:
        ent.show(20)


if __name__ == '__main__':
    main()
