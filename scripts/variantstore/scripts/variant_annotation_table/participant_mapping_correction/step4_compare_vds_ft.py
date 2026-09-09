import argparse
import hail as hl

###
# Step 4 -- adjudicate the SQL FT reconstruction against the delivered VDS.
#
# Reads the two TSV exports produced by step4_extract_genotypes.sql and compares them,
# genotype by genotype, against the VDS's own FT flag over the same (region x sample
# subset). Nothing here recomputes FT: the whole point is that the VDS already carries the
# authoritative value, written by merge_and_rescore_vdses.py:173, and the SQL is what is
# on trial.
#
# WHAT IT SETTLES
#
#   1. GQ 0 -> no call. import_gvs.py:262 nulls LGT when GQ is 0, so FT is missing too.
#      Every BigQuery row with call_GQ = 0 should appear as "BigQuery row, no VDS call",
#      and nothing else should.
#   2. The per-allele predicates. is_yes / is_ok in the SQL mirror import_gvs.py:356-363
#      including its coalesce-to-passing defaults. Disagreements on single-alt genotypes
#      isolate this, since correction 3 cannot reach them.
#   3. Correction 3 -- FT over every called non-ref allele, not only the VAT ones. Where
#      ft_sql and ft_sql_vat_only differ, exactly one can match the VDS. This is the claim
#      the 2.77x revision in Step 3 rests on and the reason this script exists.
#
# KEYING
#
# On (sample, locus). Safe for two reasons that are worth stating rather than assuming.
# VDS padding extends alleles[0] rightward from the same start position, so a variant's
# locus is identical in both artifacts even though its allele strings are not. And QUERY E
# established that no (sample_id, location) on chr21 is two separate genotypes, so the
# BigQuery side is one row per key.
#
# COST, AND WHERE TO RUN IT
#
# Hail cannot column-prune a VDS read. Entries are stored row-major, so filter_cols to a
# couple of thousand samples still decodes all ~500,000 of them; it happens after the read
# and only shrinks what is streamed downstream. The read is therefore set by the INTERVAL
# alone, which is why step4_extract_genotypes.sql buys its statistics with samples and
# keeps the window narrow. At the 2 Mb default that is about 2.5 billion entry decodes.
# filter_intervals prunes on the row index so only the overlapping partitions are read,
# and the reference data is never touched.
#
# Run it in a Terra notebook backed by a Dataproc cluster -- the configuration in
# scripts/variantstore/docs/aou/vds/cluster/AoU VDS Cluster Configuration.md -- rather than
# on a single-node runtime or wrapped in a WDL.
#
#   Not single-node. Peak memory is set by the widest row, not by the interval: a row here
#   carries ~500,000 entries of LGT/LA/LAD/GQ/RGQ/FT, tens of MB decoded, and every core
#   running a task holds one. A many-core single node is exactly the shape that turns that
#   into an OOM, and a notebook OOM costs the session.
#
#   Not a WDL either. This is a one-shot investigative comparison whose output is meant to
#   be read and poked at, not a pipeline stage. Wrapping it in run_in_hail_cluster.py buys
#   reproducible cluster provisioning at the cost of Dockstore registration and a round
#   trip per iteration, and there is nothing here to re-run on a schedule.
#
# Terra does not offer autoscaling Dataproc, so the cluster is fixed and has to be sized up
# front. At the 2 Mb / 1-in-200 default: driver n1-highmem-8, 8 workers of n1-highmem-8,
# 200 GB standard disk each, no preemptibles. Highmem specifically -- 6.5 GB per core is
# what makes the wide-row decode safe, and a standard-memory worker of the same core count
# is the configuration that OOMs. Worker count is a wall-clock dial rather than a cost one,
# since a fixed cluster is billed per VM-hour and halving the time by doubling the workers
# comes out roughly even; 8 should land in the 15-25 minute range. Scale workers with the
# WINDOW, not with the sample fraction, for the row-major reason above.
#
# Note that cluster doc's startup script hardcodes spark.driver.memory to match the driver
# machine, so it has to be adjusted if the machine type changes.
#
# EXPECTED NON-FINDINGS
#
# vds_only genotypes are not necessarily a defect. populate_alt_allele_table.py:34,40
# ingests two closed lists of GTs; 1/3, 2/3 and anything higher-arity are simply absent
# from alt_allele. A VDS call with such a GT has no BigQuery row by construction. That is
# an UNDER-reporting direction -- a carrier the mapping table misses -- which is outside
# this investigation's scope but worth counting, so vds_only is broken down by LGT.
###

TRUNC = 24


def main():
    parser = argparse.ArgumentParser(
        description='Compare the SQL FT reconstruction against the delivered VDS.')
    parser.add_argument('--vds-path', type=str, required=True,
                        help='gs:// path to the delivered VDS.')
    parser.add_argument('--genotypes', type=str, required=True, metavar='GLOB',
                        help='gs:// glob for the G3 export, e.g. '
                             'gs://BUCKET/PREFIX/step4_genotypes-*.tsv')
    parser.add_argument('--samples', type=str, required=True, metavar='GLOB',
                        help='gs:// glob for the G2 export.')
    parser.add_argument('--contig', type=str, default='chr21')
    parser.add_argument('--start', type=int, default=20000000,
                        help='Window start, inclusive. Must match G1/G3.')
    parser.add_argument('--end', type=int, default=22000000,
                        help='Window end, inclusive. Must match G1/G3.')
    parser.add_argument('--out-prefix', type=str, default=None, metavar='GS_PREFIX',
                        help='gs:// prefix under which to write the disagreement and '
                             'adjudication TSVs. Omit to report counts only.')
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--temp-path', type=str, default=None)
    parser.add_argument('--max-show', type=int, default=40)
    parser.add_argument('--trunc', type=int, default=TRUNC)
    args = parser.parse_args()

    if args.out_prefix and not args.out_prefix.startswith('gs://'):
        # Table.export() writes through the Hadoop FileSystem layer, not Python's, so a
        # local path resolves against whatever machine runs the driver -- on Dataproc a
        # cluster VM, not your workstation. It does not raise; it writes somewhere else.
        raise SystemExit('--out-prefix must be a gs:// path.')

    if args.temp_path:
        hl.init(tmp_dir=args.temp_path)
    else:
        hl.init()

    rg = args.reference

    # -----------------------------------------------------------------------------
    # The BigQuery side.
    # -----------------------------------------------------------------------------
    bq = hl.import_table(
        args.genotypes,
        delimiter='\t',
        types={'position': hl.tint32, 'call_GQ': hl.tint32, 'n_alleles': hl.tint32,
               'n_in_vat': hl.tint32, 'n_non_vat': hl.tint32,
               'n_no_filter_row': hl.tint32, 'any_yes': hl.tbool, 'all_ok': hl.tbool,
               'ft_sql': hl.tbool, 'ft_sql_vat_only': hl.tbool})
    bq = bq.key_by(locus=hl.locus(bq.contig, bq.position, reference_genome=rg),
                   s=bq.sample_name)
    bq = bq.annotate(bq_present=True)
    n_bq = bq.count()

    sample_ht = hl.import_table(args.samples, delimiter='\t')
    sample_list = sample_ht.sample_name.collect()
    sample_set = hl.literal(hl.set(sample_list))
    print(f'\nBigQuery side: {n_bq} genotypes over {len(sample_list)} selected samples.')

    # -----------------------------------------------------------------------------
    # The VDS side.
    # -----------------------------------------------------------------------------
    vds = hl.vds.read_vds(args.vds_path)
    vd = vds.variant_data

    interval = hl.parse_locus_interval(
        f'{args.contig}:{args.start}-{args.end}', reference_genome=rg)
    vd = hl.filter_intervals(vd, [interval])
    vd = vd.filter_cols(sample_set.contains(vd.s))
    n_cols = vd.count_cols()
    if n_cols != len(sample_list):
        print(f'\nNOTE: {len(sample_list)} samples selected in BigQuery but {n_cols} '
              f'present as VDS columns. The extract already excludes withdrawn and '
              f'control samples, so a shortfall here is a finding, not expected noise.')

    if 'FT' not in vd.entry:
        raise SystemExit(
            'This VDS has no FT entry field, so there is nothing to validate against. '
            'FT is written by merge_and_rescore_vdses.py:173 (and import_gvs.py:381); a '
            'VDS produced before that step will not carry it.')

    # Global indices of the called non-ref alleles. LGT is locally indexed and LA maps
    # local -> global, exactly as import_gvs.py:369-371 does it before folding FT.
    called = (hl.range(vd.LGT.ploidy)
              .map(lambda i: vd.LA[vd.LGT[i]])
              .filter(lambda x: x != 0))

    # Hail's own SNP verdict on the PADDED site alleles, and the calibration sensitivity
    # it actually looked up. as_vets survives into the delivered VDS with yng_status
    # dropped (merge_and_rescore_vdses.py:145), so this is the exact input to allele_OK
    # rather than a reconstruction of it. Guarded because a VDS built by another path may
    # not carry it.
    has_as_vets = 'as_vets' in vd.row
    if not has_as_vets:
        print('\nNOTE: this VDS has no as_vets row field, so per-allele calibration '
              'sensitivity cannot be reported. FT comparison is unaffected.')

    ent_fields = {
        'vds_present': True,
        'vds_FT': vd.FT,
        'vds_GQ': vd.GQ,
        'vds_has_lgt': hl.is_defined(vd.LGT),
        'vds_lgt': hl.str(vd.LGT),
        'vds_ref_len': hl.len(vd.alleles[0]),
        'vds_n_alt': hl.len(vd.alleles) - 1,
        'vds_called_alts': hl.delimit(
            called.map(lambda i: vd.alleles[i][:args.trunc]), ';'),
        # The padded-site question, measured rather than guessed. A true SNP at a site
        # padded by another sample's deletion is written ACT -> GCT; whether hl.is_snp
        # still calls that a SNP decides whether Hail applied 0.997 or 0.990, where the
        # SQL -- keyed on the minimized alleles -- always applies 0.997.
        'vds_called_is_snp': hl.delimit(
            called.map(lambda i: hl.str(hl.is_snp(vd.alleles[0], vd.alleles[i]))), ';'),
    }
    if has_as_vets:
        ent_fields['vds_called_cal_sens'] = hl.delimit(
            called.map(lambda i: hl.str(
                vd.as_vets.get(vd.alleles[i]).calibration_sensitivity)), ';')

    vd = vd.annotate_entries(**ent_fields)

    # Keep entries where this sample had a vet record here. LGT alone is not enough: a
    # GQ 0 call has LGT missing but is exactly the population claim 1 is about. GQ is set
    # unconditionally from the vet record (import_gvs.py:264), and the LGT expression's
    # own is_missing(GQ) guard shows GQ can itself be absent, so test both.
    #
    # Filtered on the Table rather than with filter_entries, which only sets entries to
    # missing -- entries() would still emit one row per (row, col) pair. Streamed, so the
    # product is never materialized.
    ent = vd.entries()
    ent = ent.filter(hl.is_defined(ent.GQ) | hl.is_defined(ent.LGT))
    # Re-key before selecting so that the VDS's own `alleles` field, part of the entries()
    # key, becomes droppable. Left in place it would collide with the SQL side's `alleles`
    # column in the join and be silently renamed.
    ent = ent.key_by('locus', 's').select(*ent_fields.keys())
    # Checkpoint before anything asks for a count. Without this the VDS read is repeated
    # for every downstream action -- ent.count(), the aggregate, and each export -- which
    # on a fixed cluster is the difference between one pass over the window and four or
    # five. Hail is lazy and will not cache this for you.
    ent = ent.checkpoint(hl.utils.new_temp_file('step4_vds_entries', 'ht'))
    n_vds = ent.count()
    print(f'VDS side: {n_vds} genotypes in {args.contig}:{args.start}-{args.end} '
          f'across {n_cols} columns.')

    # -----------------------------------------------------------------------------
    # Compare.
    # -----------------------------------------------------------------------------
    res = bq.join(ent, how='outer')
    res = res.annotate(bq_present=hl.is_defined(res.bq_present),
                       vds_present=hl.is_defined(res.vds_present))
    # Same reasoning: the aggregate below and the two filtered exports are three separate
    # actions over this join.
    res = res.checkpoint(hl.utils.new_temp_file('step4_joined', 'ht'))

    both = res.bq_present & res.vds_present
    # The comparable population: both sides have the genotype and the VDS has a verdict.
    # FT is missing exactly when LGT is, which is the GQ 0 case handled separately.
    cmp = both & hl.is_defined(res.vds_FT)
    adjudicable = cmp & (res.ft_sql != res.ft_sql_vat_only)

    summary = res.aggregate(hl.struct(
        n_keys=hl.agg.count(),
        n_both=hl.agg.count_where(both),
        n_bq_only=hl.agg.count_where(res.bq_present & ~res.vds_present),
        n_vds_only=hl.agg.count_where(~res.bq_present & res.vds_present),

        # Claim 1.
        n_bq_gq0=hl.agg.count_where(res.bq_present & (res.call_GQ == 0)),
        n_bq_gq0_no_vds_call=hl.agg.count_where(
            res.bq_present & (res.call_GQ == 0) & (~res.vds_present | ~res.vds_has_lgt)),
        n_bq_gq0_vds_called=hl.agg.count_where(
            res.bq_present & (res.call_GQ == 0) & res.vds_present & res.vds_has_lgt),
        n_bq_nonzero_gq_no_vds_call=hl.agg.count_where(
            res.bq_present & (res.call_GQ != 0)
            & (~res.vds_present | ~res.vds_has_lgt)),

        # Claims 2 and 3.
        n_cmp=hl.agg.count_where(cmp),
        n_agree=hl.agg.count_where(cmp & (res.vds_FT == res.ft_sql)),
        n_disagree=hl.agg.count_where(cmp & (res.vds_FT != res.ft_sql)),
        n_sql_pass_vds_fail=hl.agg.count_where(cmp & res.ft_sql & ~res.vds_FT),
        n_sql_fail_vds_pass=hl.agg.count_where(cmp & ~res.ft_sql & res.vds_FT),

        n_adjudicable=hl.agg.count_where(adjudicable),
        n_adj_corrected_wins=hl.agg.count_where(
            adjudicable & (res.vds_FT == res.ft_sql)),
        n_adj_vat_only_wins=hl.agg.count_where(
            adjudicable & (res.vds_FT == res.ft_sql_vat_only)),

        by_stratum=hl.agg.filter(res.bq_present, hl.agg.group_by(res.stratum, hl.struct(
            n=hl.agg.count(),
            n_cmp=hl.agg.count_where(cmp),
            n_agree=hl.agg.count_where(cmp & (res.vds_FT == res.ft_sql)),
            n_disagree=hl.agg.count_where(cmp & (res.vds_FT != res.ft_sql)),
            n_adjudicable=hl.agg.count_where(adjudicable),
            n_adj_corrected_wins=hl.agg.count_where(
                adjudicable & (res.vds_FT == res.ft_sql))))),

        # Padded sites, where the SQL's minimized-allele SNP test and Hail's padded-allele
        # one can part company. Reported unconditionally so its size is known even when
        # nothing disagrees.
        n_cmp_padded_site=hl.agg.count_where(cmp & (res.vds_ref_len > 1)),
        n_disagree_padded_site=hl.agg.count_where(
            cmp & (res.vds_FT != res.ft_sql) & (res.vds_ref_len > 1)),

        vds_only_by_lgt=hl.agg.filter(
            ~res.bq_present & res.vds_present & res.vds_has_lgt,
            hl.agg.counter(res.vds_lgt)),
    ))

    def pct(num, den):
        return f'{100.0 * num / den:.4f}%' if den else 'n/a'

    print(f'\n{"=" * 78}')
    print(f'STEP 4 -- SQL FT vs delivered VDS FT, {args.contig}:{args.start}-{args.end}')
    print(f'{"=" * 78}')
    print(f'  join keys (sample, locus)      : {summary.n_keys}')
    print(f'    both sides                   : {summary.n_both}')
    print(f'    BigQuery only                : {summary.n_bq_only}')
    print(f'    VDS only                     : {summary.n_vds_only}')

    print(f'\n--- claim 1: GQ 0 calls have no VDS genotype')
    print(f'  BigQuery rows with call_GQ = 0 : {summary.n_bq_gq0}')
    print(f'    no VDS call (expected)       : {summary.n_bq_gq0_no_vds_call}')
    print(f'    VDS DID call (unexpected)    : {summary.n_bq_gq0_vds_called}')
    print(f'  call_GQ != 0 but no VDS call   : {summary.n_bq_nonzero_gq_no_vds_call}')
    if summary.n_bq_gq0_vds_called or summary.n_bq_nonzero_gq_no_vds_call:
        print('  ^ either of these nonzero means the GQ 0 exclusion rule is not what the '
              'plan assumes.')

    print(f'\n--- claims 2 and 3: FT agreement')
    print(f'  comparable genotypes           : {summary.n_cmp}')
    print(f'    agree                        : {summary.n_agree} '
          f'({pct(summary.n_agree, summary.n_cmp)})')
    print(f'    disagree                     : {summary.n_disagree} '
          f'({pct(summary.n_disagree, summary.n_cmp)})')
    print(f'      SQL passes, VDS fails      : {summary.n_sql_pass_vds_fail} '
          f'(SQL under-excludes -- corrected rate is too LOW)')
    print(f'      SQL fails, VDS passes      : {summary.n_sql_fail_vds_pass} '
          f'(SQL over-excludes -- corrected rate is too HIGH)')

    print(f'\n--- correction 3, adjudicated')
    print(f'  genotypes where the two formulas differ : {summary.n_adjudicable}')
    print(f'    VDS agrees with ft_sql (corrected)    : {summary.n_adj_corrected_wins}')
    print(f'    VDS agrees with ft_sql_vat_only       : {summary.n_adj_vat_only_wins}')
    if summary.n_adjudicable == 0:
        print('  ^ zero. The window adjudicates nothing about correction 3; widen it or '
              'lower sample_modulus and rerun. Do NOT read the overall agreement rate as '
              'validating the 2.77x revision.')

    print(f'\n--- padded sites (VDS alleles[0] longer than 1 base)')
    print(f'  comparable genotypes at padded sites    : '
          f'{summary.n_cmp_padded_site}')
    print(f'  of which disagree                       : '
          f'{summary.n_disagree_padded_site}')
    print('  A disagreement concentrated here points at the SNP/indel threshold: the SQL '
          'tests LENGTH(ref) = LENGTH(allele) = 1 on minimized alleles, Hail calls '
          'hl.is_snp on the padded ones. vds_called_is_snp in the export below is Hail\'s '
          'own verdict per called allele.')

    print(f'\n--- by stratum (BigQuery side)')
    print(f'  {"stratum":<24} {"n":>10} {"n_cmp":>10} {"agree":>10} {"disagree":>10} '
          f'{"adjud":>8} {"corr_wins":>10}')
    for stratum in sorted(summary.by_stratum):
        v = summary.by_stratum[stratum]
        print(f'  {stratum:<24} {v.n:>10} {v.n_cmp:>10} {v.n_agree:>10} '
              f'{v.n_disagree:>10} {v.n_adjudicable:>8} {v.n_adj_corrected_wins:>10}')

    if summary.vds_only_by_lgt:
        print(f'\n--- VDS-only genotypes by LGT')
        print('  alt_allele ingests two closed GT lists (populate_alt_allele_table.py:'
              '34,40); 1/3, 2/3 and higher-arity calls are absent by construction. These '
              'are carriers the mapping table MISSES -- the opposite direction from this '
              'investigation, and out of its scope, but they should be recorded.')
        for lgt, n in sorted(summary.vds_only_by_lgt.items(),
                             key=lambda kv: -kv[1])[:args.max_show]:
            print(f'    {lgt:<12} {n}')

    # -----------------------------------------------------------------------------
    # Itemize the disagreements.
    # -----------------------------------------------------------------------------
    if args.out_prefix:
        prefix = args.out_prefix.rstrip('/')
        for name, ht in (('disagreements', res.filter(cmp & (res.vds_FT != res.ft_sql))),
                         ('adjudicable', res.filter(adjudicable))):
            path = f'{prefix}/step4_{name}.tsv'
            n = ht.count()
            ht.export(path)
            if hl.hadoop_exists(path):
                size = hl.hadoop_stat(path)['size_bytes']
                print(f'\nVERIFIED: wrote {n} {name} rows to {path} ({size} bytes)')
            else:
                raise RuntimeError(
                    f'export() reported no error but {path} does not exist. This is '
                    f'almost always path resolution: Hail wrote via the Hadoop '
                    f'FileSystem layer to a location relative to the driver machine.')

        print('\nRead the disagreements against vds_called_is_snp, vds_called_cal_sens '
              'and vds_ref_len. A row where the SQL called an allele a SNP and Hail did '
              'not, with a calibration sensitivity between 0.990 and 0.997, is the '
              'padded-site threshold divergence and nothing more.')
    else:
        print('\nNo --out-prefix given, so the disagreeing genotypes were not written '
              'out. Rerun with one to itemize them; the counts above cannot distinguish '
              'a systematic cause from a scatter.')


if __name__ == '__main__':
    main()
