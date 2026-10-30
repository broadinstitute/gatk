#!/usr/bin/env python3
"""Batch carrier report for a list of VIDs: one TSV row per VID, one pass over the VDS.

For each input VID, the script reads the VDS and computes `cb_after_fix`: the participants
the VDS calls as carrying that allele with a passing FT, which is what the v9_r2_p4 mappings
should contain. Every other column exists to explain that one -- what the count was before the
correction, what came out of it, and for what reason. The batch sibling of
`vds_carriers_for_vid_exact.py`: same semantics and definitions, but it takes a VID list, makes
ONE pass over the VDS instead of one per VID, and emits a TSV.

TWO THINGS THAT MAKE A NAIVE VDS COUNT DISAGREE WITH THE MAPPINGS

  1. FT. A count of every LGT-defined genotype calling the allele ignores genotype filtering.
     The v9_r2_p4 mappings apply `FT`, so they are SMALLER than such a count -- the reverse of
     v9_r2_p3, which was LARGER because it also carried GQ 0 no-calls. Both counts are reported
     (`vds_count_no_ft` and `cb_after_fix`) so the sign of the gap is never in doubt.

  2. A VID IS ONE ALLELE; FT IS ONE GENOTYPE. The mappings and the VAT are keyed per allele, so
     nothing in a VID-keyed view tells can show that the participant's genotype at that site was `1/2`.
     FT is a property of the genotype (`merge_and_rescore_vdses.py:155-173` folds
     `FT = ~any_no & (any_yes | all_ok)` over every called non-ref allele), so a hetvar is
     filtered when its OTHER allele fails -- a different allele that has its own VID.
     At 13-32368001-C-CTT, 428 of the 711 carriers are `1/2`, and all 126
     FT removals are among them. Working per-VID those participants look like ordinary het
     carriers vanishing for no reason, which is the question this report exists to preempt.
     `removed_ft_fail_hetvar` counts these removals.

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
  gvs_all_ac_no_ft        the same allele count with FT ignored. The AC-scale twin of
                          vds_count_no_ft, and the one to use when comparing a naive VDS count
                          against the VAT -- vds_count_no_ft is people, gvs_all_ac is alleles,
                          and differencing across that boundary mixes in hom-alt doubling.
  vat_gvs_all_ac          only with --vat-file: the VAT's own column of that name
  gvs_all_ac_diff         gvs_all_ac - vat_gvs_all_ac, so 0 means they agree. Read this as an
                          IMPLEMENTATION check, not as independent corroboration: the VAT's
                          gvs_all_ac is hl.agg.call_stats(GT, alleles).AC computed on the same
                          VDS after the same FT-nulling (hail_create_vat_inputs.py:60, :130),
                          so it is the same quantity under the same definition, arrived at by
                          different code. A 0 says the LA/LGT local-allele decoding here, the
                          min_rep allele matching, and the reading of FT all reproduce the dense
                          GT path the VAT took -- any of which could be wrong, so it is worth
                          having. What it CANNOT do is confirm that FT is the right definition
                          of a carrier, because gvs_all_ac is itself defined by FT.
  gvs_all_ac_no_ft_diff   gvs_all_ac_no_ft - vat_gvs_all_ac: the error a naive VDS count would
                          have made against the delivered AC. Expected to be nonzero wherever
                          the VID has any FT-failing genotype, so unlike gvs_all_ac_diff it is
                          not a warning sign and is not counted in the summary line at the end.
  vat_gvs_all_sc          only when the --vat-file export also has a gvs_all_sc column: the
                          VAT's own PEOPLE-scale carrier count, AC - Hom. Same provenance as
                          vat_gvs_all_ac and so not a second opinion, but it is on the scale the
                          mapping table counts on, which makes it the right thing to line the
                          mapping up against.
  gvs_all_sc_diff         cb_after_fix - vat_gvs_all_sc. Expected 0. A structured nonzero on
                          chrX/chrY non-PAR would be the hemizygous question flagged at
                          create_vt_bqloadjson_from_annotations.py:304, not an error here.
  mapping_n_participants  only with --mapping-count-file: how many DISTINCT participants the
                          delivered mapping table named for that release. Distinct rather than
                          ARRAY_LENGTH on purpose -- person_ids is built by a bare ARRAY_AGG and
                          can repeat a participant; see the note under PRODUCING below.
  mapping_vs_vds_ft       mapping_n_participants - cb_after_fix. THE HEADLINE. How many more
                          participants the delivered mapping named than the VDS calls with a
                          passing FT. This is the only comparison in the report between two
                          artifacts that do not share a derivation: the mapping table is built
                          from alt_allele joined to sample_info in BigQuery and never reads the
                          VDS. A positive number in a pre-v9 release is the demonstration that
                          over-reporting is not new.
  mapping_vs_vds_no_ft    mapping_n_participants - vds_count_no_ft. Near zero is the diagnosis:
                          it says the mapping behaved as an FT-ignoring carrier count. What is
                          left over is withdrawn and control samples -- excluded from the build
                          only as of VS-2000 (2026-09-09), so present in every earlier release --
                          plus any allele-representation mismatch between alt_allele and the vid.
  gq0_at_site             GQ 0 records at the site, whether or not they call this allele.
                          Context for removed_gq0, which is the subset that does.

Each VID is looked up at its own coordinate only. A VID whose stored representation differs
comes back `vds_found=false` rather than being silently counted as absent; rerun those few
through `vds_carriers_for_vid.py --window 200`.

Usage:
  vds_carriers_report.py --vds-path gs://... --vids-file my_vids.tsv --output report.tsv \\
      [--vat-file vat.csv] [--mapping-count-file mapping_counts.csv]

`--vids-file` accepts a bare list of VIDs, one per line, or any file with a `vid` column. Both
input options read TSV or CSV -- the delimiter is detected, so a spreadsheet or a BigQuery
export can be handed over as-is. Run it on a Hail cluster; it reads the VDS.

PRODUCING --vat-file. Query the VAT for the same VIDs and save the result:

    SELECT DISTINCT vid, gvs_all_ac, gvs_all_sc
    FROM `<project>.<dataset>.<vat_table>`
    WHERE vid IN UNNEST(['13-32368001-C-CTT', '1-668638-G-GA']);

Take both columns every time. They are two scales of the same thing -- gvs_all_ac counts alleles,
gvs_all_sc counts people (AC - Hom) -- and which one you need depends on what you are comparing
against, so there is no reason to make a second trip to the VAT to find out. Only gvs_all_ac is
required; gvs_all_sc simply adds the people-scale columns when present.

DISTINCT is not optional. The VAT holds one row per (vid, transcript), so a variant overlapping
several transcripts repeats both values once per transcript. They are per-variant values, so
after DISTINCT there should be exactly one row per vid; the script stops and names the vid if
there is more than one, rather than silently keeping whichever transcript sorted last.

For more VIDs than are comfortable in an IN list, load them as a table and join:

    SELECT DISTINCT v.vid, v.gvs_all_ac, v.gvs_all_sc
    FROM `<project>.<dataset>.<vat_table>` v
    JOIN `<project>.<dataset>.my_vids` m USING (vid);

PRODUCING --mapping-count-file. The one input not derived from the VDS, so worth the trouble:

    SELECT vid,
           (SELECT COUNT(DISTINCT p) FROM UNNEST(person_ids) p) AS n_participants,
           ARRAY_LENGTH(person_ids) AS n_entries
    FROM `<project>.<dataset>.<participant_mapping_table>`
    WHERE vid IN UNNEST(['13-32368001-C-CTT', '1-668638-G-GA']);

Same join-to-a-table form as above if the VID list is long.

COUNT DISTINCT, not ARRAY_LENGTH. None of the three builds de-duplicates -- a bare ARRAY_AGG at
GvsCreateParticipantMappingTable.wdl:60, GvsMapUnmappedVIDs.wdl:137 and
GvsMapDroppedDuplicateVIDs.wdl:149 -- so a participant whose data was ingested twice appears
twice in person_ids. ARRAY_LENGTH would charge that inflation to FT, which is a different defect
with a different ticket (VS-2011 duplicate batches; 24,477 VIDs on chr21 alone). Selecting
n_entries beside it keeps the distinction visible: where the two disagree, the gap is
duplication, not over-reporting. The script reads n_participants and ignores n_entries.

One row per vid in the mapping table, so no row-level DISTINCT is needed and a repeated vid is an
error worth stopping on.

Point this at the mapping table AS DELIVERED for the release in question -- for v9 that means the
pre-correction table if the question is what shipped, and the corrected one if the question is
whether the correction worked. They are different tables and the whole point of the column is
lost if the wrong one is used.
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


def read_vid_keyed(path, column, fallback_to_second=True):
    """{vid: int} from a two-column-ish export. `column` is looked up by name; with
    fallback_to_second, a file that named it something else still works as long as the value is
    the second column -- BigQuery exports get renamed by hand often enough to be worth tolerating.
    Returns {} when the column is absent and there is no fallback.

    Duplicate vids are an error, not a last-writer-wins: the VAT holds one row per
    (vid, transcript), so a query missing DISTINCT silently keeps whichever transcript sorted
    last. The counts are per-variant, so the collapse is invisible in the output."""
    lines, delim = read_lines(path)
    rows = list(csv.DictReader(lines, delimiter=delim))
    if rows and column not in rows[0]:
        if not fallback_to_second:
            return {}
        if len(rows[0]) < 2:
            raise SystemExit(f'{path} has no {column!r} column and only one column to fall back to')
    out = {}
    for row in rows:
        v = row.get(column)
        if v is None:
            v = list(row.values())[1]
        if v is None or v == '':
            continue
        if row['vid'] in out:
            raise SystemExit(f'{path} has more than one row for {row["vid"]}. Add DISTINCT to the '
                             'query -- the VAT holds one row per (vid, transcript).')
        out[row['vid']] = int(v)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vds-path', required=True)
    p.add_argument('--vids-file', required=True)
    p.add_argument('--output', required=True, help='TSV to write')
    p.add_argument('--vat-file', '--vat-ac-file', dest='vat_file',
                   help="vid + gvs_all_ac + gvs_all_sc from the VAT, TSV or CSV. Adds the VAT's "
                        "own values beside the ones this script computes from the VDS. Only "
                        "gvs_all_ac is required; gvs_all_sc adds the people-scale columns when "
                        "present. --vat-ac-file is accepted as the old spelling. See the SQL in "
                        "the module docstring for how to produce it.")
    p.add_argument('--mapping-count-file',
                   help="vid + n_participants from the delivered participant mapping table, TSV "
                        "or CSV. This is the one input here NOT derived from the VDS -- the "
                        "mapping table is built from alt_allele joined to sample_info in "
                        "BigQuery -- so it is the only genuine cross-artifact comparison the "
                        "report can make. See the SQL in the module docstring.")
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
        # The same allele count with FT ignored -- what counting the VDS naively arrives at.
        # On the AC scale specifically, because that is the only scale on which it can be
        # differenced against the VAT's gvs_all_ac: vds_count_no_ft is this same population
        # counted as PEOPLE, where a 1/1 carrier is 1 rather than 2, so differencing that
        # against the VAT would mix the FT effect with hom-alt doubling.
        gvs_all_ac_no_ft=hl.agg.filter(ent.has_allele, hl.agg.sum(ent.n_copies)),
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
    # computes its AC by a different route -- call_stats over a densified MT rather than the
    # LA/LGT arithmetic here -- so agreement means this script's local-allele decoding and
    # allele matching reproduce what the VAT actually did. Note the limit of that: it is the
    # same VDS and the same FT-based definition, so this is a check on the reimplementation,
    # not independent evidence that the definition is right. See gvs_all_ac_diff in the
    # docstring; the genuinely cross-artifact check is BigQuery's filter model against the
    # VDS's FT, which is what step4_compare_vds_ft.py does.
    vat_ac, vat_sc, mapping_n = {}, {}, {}
    if args.vat_file:
        vat_ac = read_vid_keyed(args.vat_file, 'gvs_all_ac')
        # Same file, second column: the two are scales of one quantity and are always queried
        # together, so there is no second flag for it. Absent gvs_all_sc just drops those columns.
        vat_sc = read_vid_keyed(args.vat_file, 'gvs_all_sc', fallback_to_second=False)
        print(f'{len(vat_ac)} VIDs from {args.vat_file}'
              f'{" (with gvs_all_sc)" if vat_sc else " (no gvs_all_sc column)"}', file=sys.stderr)
    if args.mapping_count_file:
        mapping_n = read_vid_keyed(args.mapping_count_file, 'n_participants')
        print(f'{len(mapping_n)} VIDs from {args.mapping_count_file}', file=sys.stderr)

    cols = ['vid', 'vds_found',
            'cb_before_fix', 'cb_after_fix',
            'removed_total', 'removed_gq0', 'removed_ft_fail', 'removed_ft_fail_hetvar',
            'vds_count_no_ft',
            'hetvar_total', 'hetvar_pass',
            'carriers_one_copy', 'carriers_two_copies', 'gvs_all_ac', 'gvs_all_ac_no_ft']
    if vat_ac:
        cols += ['vat_gvs_all_ac', 'gvs_all_ac_diff', 'gvs_all_ac_no_ft_diff']
    if vat_sc:
        cols += ['vat_gvs_all_sc', 'gvs_all_sc_diff']
    if mapping_n:
        cols += ['mapping_n_participants', 'mapping_vs_vds_ft', 'mapping_vs_vds_no_ft']
    cols += ['gq0_at_site']

    n_ac_differs, n_ft_missing, n_missing, n_mapping_over = 0, 0, 0, 0
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
                        # Deliberately NOT counted into n_ac_differs: this one is EXPECTED to be
                        # nonzero wherever the VID has any FT-failing genotype. It is the size of
                        # the error a naive VDS count would have made, not a discrepancy to chase.
                        out['gvs_all_ac_no_ft_diff'] = r['gvs_all_ac_no_ft'] - v
                        n_ac_differs += 1 if r['gvs_all_ac'] != v else 0

            # The VAT's own people-scale carrier count, AC - Hom. Same provenance as
            # vat_gvs_all_ac, so this is another implementation check, not a second opinion --
            # but it is on the people scale, which is what the mapping table counts, so it is
            # the right thing to line the mapping up against.
            if vat_sc:
                v = vat_sc.get(vid)
                if v is not None:
                    out['vat_gvs_all_sc'] = v
                    if r is not None:
                        out['gvs_all_sc_diff'] = r['cb_after_fix'] - v

            # THE INDEPENDENT NUMBER. Everything else in this row comes from the VDS, directly
            # or via the VAT. The mapping table is built from alt_allele joined to sample_info
            # in BigQuery and never reads the VDS, so these two diffs are the only cross-artifact
            # measurements here.
            if mapping_n:
                v = mapping_n.get(vid)
                if v is not None:
                    out['mapping_n_participants'] = v
                    if r is not None:
                        # What the delivered mapping over-reported against the carriers the VDS
                        # actually calls with a passing FT. This is the VS-2012 quantity.
                        out['mapping_vs_vds_ft'] = v - r['cb_after_fix']
                        # Near zero is the diagnosis: it says the mapping behaved as an
                        # FT-ignoring carrier count. Residual is withdrawn/control samples
                        # (excluded from the build only as of VS-2000, 2026-09-09) plus any
                        # allele-representation mismatch between alt_allele and the vid.
                        out['mapping_vs_vds_no_ft'] = v - r['vds_count_no_ft']
                        n_mapping_over += 1 if v > r['cb_after_fix'] else 0
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
    if mapping_n:
        print(f'{n_mapping_over} VID(s) where the delivered mapping lists MORE participants than '
              'the VDS calls with a passing FT (mapping_vs_vds_ft > 0)', file=sys.stderr)


if __name__ == '__main__':
    main()
