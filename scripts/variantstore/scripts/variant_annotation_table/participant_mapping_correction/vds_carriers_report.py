#!/usr/bin/env python3
"""Batch carrier report for a list of VIDs: one TSV row per VID, one pass over the VDS.

Answers, for each VID you give it, "who does the VDS say carries this, and does that match the
participant mapping table you received?" -- reading the VDS directly rather than trusting either
artifact. The batch sibling of `vds_carriers_for_vid_exact.py`: same semantics and definitions,
but it takes a VID list, makes ONE pass over the VDS instead of one per VID, and emits a TSV.

TWO THINGS THAT MAKE A NAIVE VDS COUNT DISAGREE WITH THE MAPPING TABLE

  1. FT. A count of every LGT-defined genotype calling the allele ignores genotype filtering.
     The mapping table from 2026-09-05 onward applies `FT`, so it is SMALLER than such a count
     -- the reverse of the pre-fix situation, where the mapping table was LARGER because it also
     carried GQ 0 no-calls. Both directions are reported (`vds_count_no_ft` vs `cb_after_fix`)
     so the sign of the gap is never in doubt.

  2. HETVARS ARE INVISIBLE IN A SPLIT REPRESENTATION. Splitting a multi-allelic site turns a
     `1/2` into `0/1` at each allele's row, so a split view counts hetvar carriers as ordinary
     hets. At 13-32368001-C-CTT a split count gives 711 where the unsplit truth is 283 `0/1`
     plus 428 `1/2`. This matters because FT is a property of the GENOTYPE, not the allele
     (`merge_and_rescore_vdses.py:155-173` folds `FT = ~any_no & (any_yes | all_ok)` over every
     called non-ref allele), so a hetvar is filtered when its OTHER allele fails. At that VID
     all 126 FT removals are hetvars. In a split view those participants look like ordinary het
     carriers vanishing for no reason, which is the question this report exists to pre-empt.

COLUMN GROUPS

  the answer          cb_after_fix (= carriers, what the corrected mapping table should hold)
  what changed        cb_before_fix, removed_total, and the reasons: removed_gq0,
                      removed_ft_fail, removed_ft_fail_hetvar
  unfiltered counts   vds_count_no_ft, vds_one_copy_no_ft, vds_two_copies_no_ft -- what a
                      count that ignores FT would give, so a pre-existing tally can be joined
                      and confirmed to be looking at the same data before anything else is read
  context             hetvar_total/pass/fail, gvs_all_ac, carriers_one_copy/two_copies

CHECKING. Nothing here is validated by an identity among its own columns -- `cb_before_fix` is
derived from `vds_count_no_ft` and `removed_gq0`, so subtracting them back out is arithmetic,
not evidence. The checks that mean something compare against a source outside this script:

  agrees_mapping       our carrier count vs the array length in the mapping table you were
                       sent, via --mapping-tsv. This is the check that answers "is the file I
                       received right", and the one to run first.
  agrees_cb_before_fix our reconstructed pre-fix count vs a pre-fix cb_participants tally
  agrees_vds_count     our FT-ignoring count vs an existing vds_participant_count tally
  agrees_gvs_all_ac    our AC vs the VAT's, requires --vat-ac-tsv. The strongest of these:
                       gvs_all_ac is produced by an entirely separate pipeline, so agreement
                       says the FT and GQ 0 handling reproduces what the VAT actually did,
                       rather than merely being self-consistent.
  ft_all_defined       internal, and the only one that is: no entry has a missing FT.

`all_checks_pass` ANDs whichever of those are available. It is absent entirely when no
comparison input is supplied, because with nothing external to check against there is nothing
honest to put in it.

Each VID is looked up at its own coordinate only. A VID whose stored representation differs
comes back `vds_found=false` rather than being silently counted as absent; rerun those few
through `vds_carriers_for_vid.py --window 200`.

Usage:
  vds_carriers_report.py --vds-path gs://... --vids-file my_vids.tsv --output report.tsv \\
      [--mapping-tsv mapping_lengths.tsv] [--vat-ac-tsv vat_ac.tsv] \\
      [--compare-tsv an_existing_tally.tsv]

`--vids-file` accepts a bare list of VIDs, one per line, or any TSV with a `vid` column.
`--mapping-tsv` wants `vid` plus the number of person_ids for that VID; see the delivery notes
for how to produce it from the Parquet. `--compare-tsv` merges an existing tally in beside ours
and flags disagreement.
"""

import argparse
import csv
import sys

import hail as hl


def min_rep_py(pos, ref, alt):
    """Pure-Python `hl.min_rep`: trim shared suffix, then shared prefix, never below length 1
    on either side. Done here rather than with `hl.eval` to avoid one Hail round trip per VID."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def parse_vid(vid):
    """'13-32368001-C-CTT' -> ('chr13', 32368001, 'C', 'CTT')."""
    chrom, pos, ref, alt = vid.rsplit('-', 3)
    return f'chr{chrom}', int(pos), ref, alt


def hopen(path, mode='r'):
    """Open a local path or a URI. Builtin open() cannot touch gs://, and every path here can
    reasonably be one -- the VDS lives in a bucket, so the report usually wants to as well.
    Requires hl.init() to have run already, which is why main() initializes Hail first."""
    return hl.hadoop_open(path, mode) if '://' in path else open(path, mode)


def read_vids(path):
    """One VID per line, or a delimited file with a `vid` column."""
    with hopen(path) as f:
        rows = [ln.rstrip('\n').split('\t') for ln in f if ln.strip()]
    if not rows:
        raise SystemExit(f'{path} is empty')
    header, idx = rows[0], 0
    if 'vid' in header:
        idx = header.index('vid')
        rows = rows[1:]
    elif header[0].lower() == 'vid':
        rows = rows[1:]
    vids, seen = [], set()
    for r in rows:
        v = r[idx].strip()
        if v and v not in seen:
            seen.add(v)
            vids.append(v)
    return vids


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vds-path', required=True)
    p.add_argument('--vids-file', required=True)
    p.add_argument('--output', required=True, help='TSV to write')
    p.add_argument('--compare-tsv', help="an existing tally, to merge in beside ours")
    p.add_argument('--vat-ac-tsv',
                   help="TSV of vid + gvs_all_ac from the VAT, to validate the carrier "
                        "definition against an independently built artifact")
    p.add_argument('--mapping-tsv',
                   help="TSV of vid + the number of person_ids that VID has in the mapping "
                        "table you received. Checks the delivered file against the VDS.")
    p.add_argument('--reference-genome', default='GRCh38')
    args = p.parse_args()

    rg = args.reference_genome
    # Before reading anything: the input paths may themselves be gs:// URIs, and hopen() goes
    # through Hail's filesystem, which needs a live context.
    hl.init()

    vids = read_vids(args.vids_file)
    print(f'{len(vids)} distinct VIDs from {args.vids_file}', file=sys.stderr)

    # -------------------------------------------------------------------------------------
    # Targets, grouped by locus. Several VIDs commonly share one -- in the chr13 poly-T tract
    # 100 VIDs sit on 36 loci -- so this is a locus -> [target] dict, not a flat list.
    # -------------------------------------------------------------------------------------
    by_locus = {}
    for vid in vids:
        contig, pos, ref, alt = parse_vid(vid)
        mpos, mref, malt = min_rep_py(pos, ref, alt)
        by_locus.setdefault(f'{contig}:{pos}', []).append(
            hl.Struct(vid=vid, min_pos=mpos, min_ref=mref, min_alt=malt,
                      delta=len(alt) - len(ref)))
    print(f'{len(by_locus)} distinct loci', file=sys.stderr)

    target_type = hl.tarray(hl.tstruct(vid=hl.tstr, min_pos=hl.tint32, min_ref=hl.tstr,
                                       min_alt=hl.tstr, delta=hl.tint32))
    targets_lit = hl.literal(by_locus, hl.tdict(hl.tstr, target_type))

    vds = hl.vds.read_vds(args.vds_path)
    vd = vds.variant_data
    if 'FT' not in vd.entry:
        raise SystemExit('This VDS has no FT entry field; see merge_and_rescore_vdses.py:173.')

    # Closed point intervals. The string form `chr13:X-X` is end-exclusive and therefore empty.
    # Built with hl.Locus rather than hl.eval(hl.locus(...)): the latter is a Hail round trip
    # each, which is one thing when there is a single VID and 36 of them when there is a list.
    intervals = []
    for key in by_locus:
        contig, pos = key.split(':')
        loc = hl.Locus(contig, int(pos), reference_genome=rg)
        intervals.append(hl.Interval(loc, loc, includes_start=True, includes_end=True))
    vd = hl.filter_intervals(vd, intervals)

    # -------------------------------------------------------------------------------------
    # Match each target to a global alt index at its locus.
    # -------------------------------------------------------------------------------------
    # Padding, not left alignment, is why this needs min_rep at all: the VDS pads alleles[0]
    # to span the longest deletion at the site, so a VID's ALT string is usually not the string
    # in `alleles`. Compare the minimized PAIR -- at chr13:32368001 the string `CTT` is both the
    # raw ALT of a 24bp deletion and the minimized ALT of a 2bp insertion.
    row_targets = targets_lit.get(hl.str(vd.locus), hl.empty_array(target_type.element_type))
    vd = vd.annotate_rows(
        matches=row_targets.map(lambda t: hl.struct(
            vid=t.vid,
            # Indel length is invariant under trimming, so this prefilter is free and exact.
            target_idx=(hl.range(1, hl.len(vd.alleles))
                        .filter(lambda i: (hl.len(vd.alleles[i]) - hl.len(vd.alleles[0]))
                                == t.delta)
                        .filter(lambda i: hl.bind(
                            lambda mr: (mr.locus.position == t.min_pos)
                            & (mr.alleles[0] == t.min_ref)
                            & (mr.alleles[1] == t.min_alt),
                            hl.min_rep(vd.locus, [vd.alleles[0], vd.alleles[i]])))
                        .first()))))
    vd = vd.filter_rows(hl.len(vd.matches) > 0)

    # -------------------------------------------------------------------------------------
    # Aggregate per VID.
    # -------------------------------------------------------------------------------------
    ent = vd.entries()
    # Before exploding, drop samples with no vet record at the site. GQ is set unconditionally
    # from the vet row (import_gvs.py:264) while LGT is nulled at GQ 0 (:262), so `GQ` defined
    # is exactly "has a record" and this is both the right filter and the one that keeps the
    # explode from multiplying 535k mostly-empty entries by the VIDs at each locus.
    ent = ent.filter(hl.is_defined(ent.GQ))
    ent = ent.explode('matches')
    ent = ent.filter(hl.is_defined(ent.matches.target_idx))
    ent = ent.annotate(vid=ent.matches.vid, ti=ent.matches.target_idx)

    called = (hl.range(ent.LGT.ploidy)
              .map(lambda i: ent.LA[ent.LGT[i]])
              .filter(lambda x: x != 0))
    ent = ent.annotate(
        has_allele=hl.is_defined(ent.LGT) & called.contains(ent.ti),
        hetvar=hl.len(hl.set(called)) > 1,
        n_copies=called.filter(lambda i: i == ent.ti).length(),
        # LGT missing but GQ present: a record that called nothing. LA survives -- it is built
        # from the vet row's own ref/alt (import_gvs.py:259, :286) and never consults GQ -- so
        # the allele is still attributable, and LA membership is the same criterion by which
        # the pre-fix mapping table counted these participants.
        gq0=hl.is_missing(ent.LGT))

    res = ent.group_by(ent.vid).aggregate(
        vds_count_no_ft=hl.agg.count_where(ent.has_allele),
        vds_one_copy_no_ft=hl.agg.count_where(ent.has_allele & (ent.n_copies == 1)),
        vds_two_copies_no_ft=hl.agg.count_where(ent.has_allele & (ent.n_copies == 2)),
        cb_after_fix=hl.agg.count_where(ent.has_allele & ent.FT),
        removed_ft_fail=hl.agg.count_where(ent.has_allele & ~ent.FT),
        removed_ft_fail_hetvar=hl.agg.count_where(ent.has_allele & ~ent.FT & ent.hetvar),
        removed_gq0=hl.agg.count_where(ent.gq0 & ent.LA.contains(ent.ti)),
        hetvar_total=hl.agg.count_where(ent.has_allele & ent.hetvar),
        hetvar_pass=hl.agg.count_where(ent.has_allele & ent.FT & ent.hetvar),
        carriers_one_copy=hl.agg.count_where(ent.has_allele & ent.FT & (ent.n_copies == 1)),
        carriers_two_copies=hl.agg.count_where(ent.has_allele & ent.FT & (ent.n_copies == 2)),
        # gvs_all_ac counts ALLELES, not people: a 1/1 carrier contributes 2. Summed over
        # exactly the carrier set, so it is comparable to the VAT and not to any people column.
        gvs_all_ac=hl.agg.filter(ent.has_allele & ent.FT, hl.agg.sum(ent.n_copies)),
        ft_missing=hl.agg.count_where(ent.has_allele & hl.is_missing(ent.FT)),
        gq0_at_site=hl.agg.count_where(ent.gq0),
    )

    # 100 rows. Collect and write in Python so unmatched VIDs and their file can be merged in
    # without a second Hail join.
    found = {r['vid']: r for r in res.collect()}
    print(f'{len(found)} of {len(vids)} VIDs matched at their own coordinate', file=sys.stderr)

    theirs = {}
    if args.compare_tsv:
        with hopen(args.compare_tsv) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                theirs[row['vid']] = row

    # VAT allele counts, if supplied. This is the one check that validates the CARRIER
    # DEFINITION rather than the arithmetic: gvs_all_ac is computed by an entirely separate
    # pipeline, so agreement means our FT/GQ-0 handling reproduces what the VAT did.
    vat_ac = {}
    if args.vat_ac_tsv:
        with hopen(args.vat_ac_tsv) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                ac = row.get('gvs_all_ac')
                if ac is None:  # tolerate a headerless-ish two-column export
                    ac = list(row.values())[1]
                vat_ac[row['vid']] = int(ac)

    if not vat_ac and theirs and 'gvs_all_ac' in next(iter(theirs.values())):
        # Their comparison file may already carry the VAT's AC, in which case it is the same
        # external reference --vat-ac-tsv would supply and there is no reason to make anyone
        # run a separate export for it.
        vat_ac = {v: int(t['gvs_all_ac']) for v, t in theirs.items() if t.get('gvs_all_ac')}
        print(f'using gvs_all_ac from {args.compare_tsv} ({len(vat_ac)} VIDs)', file=sys.stderr)

    # What the delivered mapping table actually holds. Checking against this is the point of
    # the exercise for anyone who received the file: every other comparison here validates our
    # reasoning, and this one validates the artifact.
    mapping_n = {}
    if args.mapping_tsv:
        with hopen(args.mapping_tsv) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                n = next((row[k] for k in ('n_person_ids', 'n', 'count', 'person_ids_length')
                          if row.get(k) not in (None, '')), None)
                if n is None:  # tolerate any two-column vid/count export
                    n = list(row.values())[1]
                mapping_n[row['vid']] = int(n)
        print(f'{len(mapping_n)} VIDs from {args.mapping_tsv}', file=sys.stderr)

    cols = ['vid', 'vds_found',
            'cb_before_fix', 'cb_after_fix', 'vds_count_no_ft',
            'removed_total', 'removed_gq0', 'removed_ft_fail', 'removed_ft_fail_hetvar',
            'vds_one_copy_no_ft', 'vds_two_copies_no_ft',
            'hetvar_total', 'hetvar_pass', 'hetvar_fail',
            'carriers_one_copy', 'carriers_two_copies', 'gvs_all_ac',
            'gq0_at_site', 'ft_missing', 'ft_all_defined']
    if mapping_n:
        cols += ['mapping_n_person_ids', 'agrees_mapping']
    if theirs:
        cols += ['their_cb_participants', 'their_vds_participant_count',
                 'their_vds_n_ref_alt', 'their_vds_n_alt_alt',
                 'agrees_cb_before_fix', 'agrees_vds_count']
    if vat_ac:
        cols += ['vat_gvs_all_ac', 'agrees_gvs_all_ac']
    if mapping_n or theirs or vat_ac:
        cols += ['all_checks_pass']

    n_bad, n_missing = 0, 0
    with hopen(args.output, 'w') as f:
        # lineterminator is explicit because csv's default dialect writes \r\n, and the usual
        # `newline=''` guard is not available on a Hail filesystem handle.
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore',
                           lineterminator='\n')
        w.writeheader()
        for vid in vids:
            r = found.get(vid)
            if r is None:
                n_missing += 1
                out = {c: '' for c in cols}
                out['vid'], out['vds_found'] = vid, 'false'
            else:
                cb_before = r['vds_count_no_ft'] + r['removed_gq0']
                out = {c: r[c] for c in cols if c in r}
                out.update(vid=vid, vds_found='true', cb_before_fix=cb_before,
                           removed_total=cb_before - r['cb_after_fix'],
                           hetvar_fail=r['removed_ft_fail_hetvar'],
                           ft_all_defined=str(r['ft_missing'] == 0).lower())

            # Every check below compares against a source OUTSIDE this script. An internal
            # identity among these columns would be arithmetic, not evidence: cb_before_fix is
            # derived from vds_count_no_ft and removed_gq0, so subtracting them back out proves
            # nothing.
            checks = [] if r is None else [r['ft_missing'] == 0]
            if mapping_n:
                m = mapping_n.get(vid)
                if m is not None:
                    out['mapping_n_person_ids'] = m
                    if r is not None:
                        agrees = m == r['cb_after_fix']
                        out['agrees_mapping'] = str(agrees).lower()
                        checks.append(agrees)
            if theirs:
                t = theirs.get(vid)
                if t:
                    out['their_cb_participants'] = t['cb_participants']
                    out['their_vds_participant_count'] = t['vds_participant_count']
                    out['their_vds_n_ref_alt'] = t['vds_n_ref_alt']
                    out['their_vds_n_alt_alt'] = t['vds_n_alt_alt']
                    if r is not None:
                        a = int(t['cb_participants']) == out['cb_before_fix']
                        b = int(t['vds_participant_count']) == r['vds_count_no_ft']
                        out['agrees_cb_before_fix'] = str(a).lower()
                        out['agrees_vds_count'] = str(b).lower()
                        checks += [a, b]
            if vat_ac:
                v = vat_ac.get(vid)
                if v is not None:
                    out['vat_gvs_all_ac'] = v
                    if r is not None:
                        c = v == r['gvs_all_ac']
                        out['agrees_gvs_all_ac'] = str(c).lower()
                        checks.append(c)
            if (mapping_n or theirs or vat_ac) and r is not None:
                ok = all(checks)
                out['all_checks_pass'] = str(ok).lower()
                n_bad += 0 if ok else 1
            w.writerow(out)

    print(f'\nWrote {args.output}', file=sys.stderr)
    if n_missing:
        print(f'{n_missing} VID(s) did not match at their own coordinate (vds_found=false). '
              'Rerun those through vds_carriers_for_vid.py --window 200.', file=sys.stderr)
    if n_bad:
        print(f'WARNING: {n_bad} row(s) have all_checks_pass=false.', file=sys.stderr)
    elif not (args.mapping_tsv or args.compare_tsv or args.vat_ac_tsv):
        print('NOTE: none of --mapping-tsv, --compare-tsv or --vat-ac-tsv given, so nothing in '
              'this report was checked against an external source.', file=sys.stderr)


if __name__ == '__main__':
    main()
