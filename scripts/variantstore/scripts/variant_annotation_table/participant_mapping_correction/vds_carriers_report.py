#!/usr/bin/env python3
"""Batch carrier report for a list of VIDs: one TSV row per VID, one pass over the VDS.

For each VID you give it, the script reads the VDS and computes `cb_after_fix`: the participants
the VDS calls as carrying that allele with a passing FT, which is what the v9_r2_p4 mappings
should contain. Every other column exists to explain that one -- what the count was before the
correction, what came out of it, and for what reason. The batch sibling of
`vds_carriers_for_vid_exact.py`: same semantics and definitions, but it takes a VID list, makes
ONE pass over the VDS instead of one per VID, and emits a TSV.

VIDs in, counts out. It does not read your mapping file and does not join anything to it -- join
the output to your own data on `vid` in whatever you already use. `--vids-file` takes an existing
spreadsheet directly as long as it has a `vid` column, so the usual workflow is to hand it that
file and merge the result back afterwards.

TWO THINGS THAT MAKE A NAIVE VDS COUNT DISAGREE WITH THE MAPPINGS

  1. FT. A count of every LGT-defined genotype calling the allele ignores genotype filtering.
     The v9_r2_p4 mappings apply `FT`, so they are SMALLER than such a count -- the reverse of
     v9_r2_p3, which was LARGER because it also carried GQ 0 no-calls. Both counts are reported
     (`vds_count_no_ft` and `cb_after_fix`) so the sign of the gap is never in doubt.

  2. HETVARS ARE INVISIBLE IN A SPLIT REPRESENTATION. Splitting a multi-allelic site turns a
     `1/2` into `0/1` at each allele's row, so a split view counts hetvar carriers as ordinary
     hets. At 13-32368001-C-CTT a split count gives 711 where the unsplit truth is 283 `0/1`
     plus 428 `1/2`. This matters because FT is a property of the GENOTYPE, not the allele
     (`merge_and_rescore_vdses.py:155-173` folds `FT = ~any_no & (any_yes | all_ok)` over every
     called non-ref allele), so a hetvar is filtered when its OTHER allele fails. At that VID
     all 126 FT removals are hetvars. In a split view those participants look like ordinary het
     carriers vanishing for no reason, which is the question this report exists to pre-empt.

COLUMNS

  vid                     as given, echoed so the output can be joined straight back
  vds_found               false when the VID was not found at its own coordinate; every count
                          is then blank. See the note on representation below.
  cb_after_fix            THE ANSWER: participants the VDS calls as carrying this allele with a
                          passing FT. This is what the v9_r2_p4 mappings should contain.
  cb_before_fix           what the v9_r2_p3 mappings contained, reconstructed
  removed_total           cb_before_fix - cb_after_fix, split by reason into the next two
  removed_gq0             ...GQ 0 no-calls: a record exists but the VDS calls no genotype
  removed_ft_fail         ...genotypes calling this allele that fail FT
  removed_ft_fail_hetvar  of removed_ft_fail, how many were 1/2 het-vars. A large share here
                          means most of the removals are the partner-allele case in (2) above.
  vds_count_no_ft         carriers counted while ignoring FT, i.e. what a straightforward VDS
                          count gives. Equals cb_after_fix + removed_ft_fail.
  hetvar_total            of the carriers, how many are 1/2 -- how exposed this VID is to (2)
  hetvar_pass             of those, how many pass FT
  carriers_one_copy       cb_after_fix split by 0/1 and 1/1
  carriers_two_copies
  gvs_all_ac              allele count over the carriers, so a 1/1 counts twice. Computed here
                          from the VDS.
  vat_gvs_all_ac          only with --vat-ac-file: the VAT's own column of that name
  gvs_all_ac_diff         gvs_all_ac - vat_gvs_all_ac, so 0 means they agree. This is the one
                          number here checked against something built independently -- the VAT
                          comes from a separate pipeline, so 0 means the FT and GQ 0 handling
                          reproduces what the VAT did rather than merely being self-consistent.
  gq0_at_site             GQ 0 records at the site, whether or not they call this allele.
                          Context for removed_gq0, which is the subset that does.

Each VID is looked up at its own coordinate only. A VID whose stored representation differs
comes back `vds_found=false` rather than being silently counted as absent; rerun those few
through `vds_carriers_for_vid.py --window 200`.

Usage:
  vds_carriers_report.py --vds-path gs://... --vids-file my_vids.tsv --output report.tsv \\
      [--vat-ac-file vat_ac.csv]

`--vids-file` accepts a bare list of VIDs, one per line, or any file with a `vid` column. Both
input options read TSV or CSV -- the delimiter is detected, so a spreadsheet or a BigQuery
export can be handed over as-is. Run it on a Hail cluster; it reads the VDS.

PRODUCING --vat-ac-file. Worth the trouble: it is the only column here checked against something
this script did not compute. Query the VAT for the same VIDs and save the result:

    SELECT DISTINCT vid, gvs_all_ac
    FROM `<project>.<dataset>.<vat_table>`
    WHERE vid IN UNNEST(['13-32368001-C-CTT', '1-668638-G-GA']);

DISTINCT is not optional. The VAT holds one row per (vid, transcript), so a variant overlapping
several transcripts repeats `gvs_all_ac` once per transcript. It is a per-variant value, so
after DISTINCT there should be exactly one row per vid; if there is more than one, stop and
work out why before trusting anything downstream of it.

For more VIDs than are comfortable in an IN list, load them as a table and join:

    SELECT DISTINCT v.vid, v.gvs_all_ac
    FROM `<project>.<dataset>.<vat_table>` v
    JOIN `<project>.<dataset>.my_vids` m USING (vid);
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


def read_lines(path):
    """Whole file, blank lines dropped, plus its delimiter. BigQuery hands you CSV and our own
    exports are TSV; guessing wrong silently collapses a row into a single column, so detect it
    from the header rather than making the caller convert. Reads it all -- these are VID lists,
    not the VDS."""
    with hopen(path) as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip()]
    if not lines:
        raise SystemExit(f'{path} is empty')
    return lines, ('\t' if '\t' in lines[0] else ',')


def read_vids(path):
    """One VID per line, or a delimited file with a `vid` column."""
    lines, delim = read_lines(path)
    rows = [ln.split(delim) for ln in lines]
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
    p.add_argument('--vat-ac-file',
                   help="vid + gvs_all_ac from the VAT, TSV or CSV. Adds a vat_gvs_all_ac column "
                        "beside the gvs_all_ac this script computes from the VDS. See the SQL "
                        "in the module docstring for how to produce it.")
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
        # Structurally zero, and not reported as a column for that reason: FT is a fold over
        # hl.range(LGT.ploidy) (merge_and_rescore_vdses.py:173, import_gvs.py:380), so it is
        # missing exactly when LGT is, and has_allele already requires LGT defined. Kept as a
        # tripwire in case this is ever pointed at a VDS whose FT came from somewhere else.
        ft_missing=hl.agg.count_where(ent.has_allele & hl.is_missing(ent.FT)),
        gq0_at_site=hl.agg.count_where(ent.gq0),
    )

    # A row per VID, so a handful. Collect and write in Python: the VIDs that matched nothing
    # have to be emitted too, and reintroducing them is an outer join Hail would charge for.
    found = {r['vid']: r for r in res.collect()}
    print(f'{len(found)} of {len(vids)} VIDs matched at their own coordinate', file=sys.stderr)

    # The VAT's own allele count, if supplied. This is worth having beside ours because the VAT
    # is built by an entirely separate pipeline: agreement means the FT and GQ 0 handling here
    # reproduces what the VAT actually did, rather than merely being self-consistent.
    vat_ac = {}
    if args.vat_ac_file:
        lines, delim = read_lines(args.vat_ac_file)
        for row in csv.DictReader(lines, delimiter=delim):
            ac = row.get('gvs_all_ac')
            if ac is None:  # tolerate a two-column export that named the columns differently
                ac = list(row.values())[1]
            vat_ac[row['vid']] = int(ac)
        # One row per vid after DISTINCT; more means the query kept the transcript explosion,
        # and the dict silently kept whichever row came last.
        if len(vat_ac) != len(lines) - 1:
            raise SystemExit(f'{args.vat_ac_file} has {len(lines) - 1} rows for {len(vat_ac)} '
                             'distinct VIDs. Add DISTINCT to the query -- the VAT holds one row '
                             'per (vid, transcript).')
        print(f'{len(vat_ac)} VIDs from {args.vat_ac_file}', file=sys.stderr)

    cols = ['vid', 'vds_found',
            'cb_before_fix', 'cb_after_fix',
            'removed_total', 'removed_gq0', 'removed_ft_fail', 'removed_ft_fail_hetvar',
            'vds_count_no_ft',
            'hetvar_total', 'hetvar_pass',
            'carriers_one_copy', 'carriers_two_copies', 'gvs_all_ac']
    if vat_ac:
        cols += ['vat_gvs_all_ac', 'gvs_all_ac_diff']
    cols += ['gq0_at_site']

    n_ac_differs, n_ft_missing, n_missing = 0, 0, 0
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
                           removed_total=cb_before - r['cb_after_fix'])
                n_ft_missing += 1 if r['ft_missing'] else 0

            # Reported as a signed difference rather than a boolean: it reads as 0 down the
            # column when all is well, and on the day it does not, the size and direction of
            # the gap are the first thing anyone would want.
            if vat_ac:
                v = vat_ac.get(vid)
                if v is not None:
                    out['vat_gvs_all_ac'] = v
                    if r is not None:
                        out['gvs_all_ac_diff'] = r['gvs_all_ac'] - v
                        n_ac_differs += 1 if r['gvs_all_ac'] != v else 0
            w.writerow(out)

    print(f'\nWrote {args.output}', file=sys.stderr)
    if n_missing:
        print(f'{n_missing} VID(s) did not match at their own coordinate (vds_found=false). '
              'Rerun those through vds_carriers_for_vid.py --window 200.', file=sys.stderr)
    if n_ft_missing:
        print(f'WARNING: {n_ft_missing} VID(s) have entries with an undefined FT, which should '
              'not be possible. Their counts are not trustworthy; check how this VDS was built.',
              file=sys.stderr)
    if vat_ac:
        print(f'{n_ac_differs} VID(s) where gvs_all_ac differs from the VAT'
              f"'s{' -- investigate' if n_ac_differs else ''}", file=sys.stderr)


if __name__ == '__main__':
    main()
