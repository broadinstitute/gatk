"""Prototype for VS-2054: per-locus no-call sample sets from a VDS, for computing AN over an arbitrary cohort.

For every variant locus that survives the VAT's filters, this records which samples are no-calls there, using the
same genotypes the VAT computes AN from: the three transforms of hail_create_vat_inputs.main, then
hl.vds.to_dense_mt. A sample is a no-call at a locus when its dense GT is missing.

A set is stored in one of two forms, whichever is smaller:

  nocall_sample_idxs  sorted sample indices, when there are at most n_samples // 32 of them
  nocall_bitmap_hex   one bit per sample as a hex string, otherwise. Sample i is bit (0x80 >> (i % 8)) of byte
                      i // 8, so BigQuery's FROM_HEX turns the string straight into the BYTES bitmap.

A sample's index is its position in the dense MatrixTable's columns, written to samples.tsv. The VDS does not carry
GVS sample_ids, so mapping indices to them is left to whatever loads the output.

The sets say who is called, not with what ploidy. A lookup assumes each called sample contributes a default: 2,
or on non-PAR chrX the ploidy most of that sample's calls there have (chrx_default_ploidy in samples.tsv).
n_nondefault_ploidy counts the called samples at a locus for which that is wrong, and should be 0 wherever the
lookup is to be exact.

Vids. The VAT splits each row with split_multi(left_aligned=True), which fails if minimizing an allele moves its
locus, so every allele's key (contig-position-ref-alt, minimized) sits at its row's own locus. That key is the vid
unless the indel can still move left (see hail_create_participant_mapping_inputs.py), in which case bcftools norm
gives the vid, at another position. So a vid is looked up at its own position unless it is the image of such a
key. Those keys go to shifted_keys.parquet with their carriers, and to keys_to_normalize.vcf.bgz for bcftools
norm. alt_alleles lists every key at a locus as ref-alt, which is enough to tell whether a vid also has a row of
its own at its position.

With --num_test_cohorts, the same pass also records the AN that hl.agg.call_stats gives for each of several random
cohorts (expected_an), so AN computed from the sets can be checked against it wherever the output ends up.
"""

import argparse
import json
import random
from datetime import datetime

import hail as hl

import hail_create_vat_inputs
from hail_create_participant_mapping_inputs import SHIFTABLE_CLASSES, key_class

BITS_PER_WORD = 32

# Keys bcftools norm moves. 'unanchored' is an indel written without its anchor base, such as REF C...G with ALT G:
# minimizing leaves it alone, but bcftools norm re-anchors it at least one base to the left. The AnVIL 3K VDS has 279
# of them, all deletions of more than 200 bases. The participant mapping script fails on this class instead.
NORMALIZED_CLASSES = SHIFTABLE_CLASSES + ['unanchored']

# GVS location encoding: chromosome index * 1e12 + position, as in SchemaUtils.encodeLocation.
CONTIG_INDEX = {**{f'chr{i}': i for i in range(1, 23)}, 'chrX': 23, 'chrY': 24}
LOCATION_MULTIPLIER = 1_000_000_000_000


def filtered_dense_mt(vds):
    """The dense MatrixTable the VAT computes AC, AN and AF from."""
    for transform in (hail_create_vat_inputs.remove_too_many_alt_allele_sites,
                      hail_create_vat_inputs.hard_filter_non_passing_sites,
                      hail_create_vat_inputs.failing_gts_to_no_call):
        vds = transform(vds)
    return hl.vds.to_dense_mt(vds)


def pick_test_cohorts(n_samples, num_cohorts, seed):
    """Random cohorts as sorted lists of sample indices, sized from a single sample up to nearly everyone."""
    rng = random.Random(seed)
    sizes = [1, 2, 10, n_samples // 100, n_samples // 10, n_samples // 2, n_samples - 1]
    sizes = [s for s in sizes if 1 <= s <= n_samples]
    return [sorted(rng.sample(range(n_samples), sizes[i % len(sizes)])) for i in range(num_cohorts)]


def sample_table(mt):
    """One dict per sample, in column order: its name and its default ploidy on non-PAR chrX."""
    x = mt.filter_rows(mt.locus.in_x_nonpar())
    x = x.annotate_cols(_ploidies=hl.agg.filter(hl.is_defined(x.GT), hl.agg.counter(x.GT.ploidy)))
    # cols() sorts by column key unless the columns are unkeyed first.
    samples = []
    for col in x.key_cols_by().cols().collect():
        ploidies = dict(col._ploidies)
        samples.append(dict(
            sample_name=col.s,
            chrx_nonpar_haploid_calls=ploidies.get(1, 0),
            chrx_nonpar_diploid_calls=ploidies.get(2, 0),
            chrx_default_ploidy=max(ploidies, key=ploidies.get) if ploidies else 2,
        ))
    return samples


def split_alt_alleles(mt):
    """Per row, one struct per alternate allele the VAT would give a vid: its index and minimized alleles."""
    ref = mt.alleles[0]

    def minimized(a):
        m = hl.min_rep(mt.locus, [ref, mt.alleles[a]])
        # split_multi(left_aligned=True) fails the VAT's run on such an allele, so this never fires on a VDS that
        # has a VAT.
        return (hl.case()
                .when(m.locus == mt.locus,
                      hl.struct(a_index=a, ref=m.alleles[0], alt=m.alleles[1],
                                shiftable=hl.literal(NORMALIZED_CLASSES).contains(key_class(m.alleles[0], m.alleles[1]))))
                .or_error('minimizing an allele moved its locus at ' + hl.str(mt.locus)))

    # The VAT drops spanning deletions, and sites with an N in REF (see carrier_entries in the mapping script).
    return (hl.range(1, hl.len(mt.alleles))
            .filter(lambda a: (mt.alleles[a] != '*') & ~ref.contains('N'))
            .map(minimized))


def nocall_sets(mt, n_samples, cohorts=None, chrx_default_ploidy=None):
    """One row per locus of the dense MatrixTable `mt`: its no-call set, and AN over all samples."""
    n_words = (n_samples + BITS_PER_WORD - 1) // BITS_PER_WORD
    max_sparse = n_samples // BITS_PER_WORD

    mt = mt.add_col_index('_idx')
    mt = mt.annotate_cols(_idx=hl.int32(mt._idx))
    if cohorts is not None:
        membership = [[False] * n_samples for _ in cohorts]
        for c, cohort in enumerate(cohorts):
            for i in cohort:
                membership[c][i] = True
        per_sample = hl.literal([[membership[c][i] for c in range(len(cohorts))] for i in range(n_samples)])
        mt = mt.annotate_cols(_cohorts=per_sample[mt._idx])

    # A sample with neither a variant call nor a covering reference block has a missing entry in the dense
    # MatrixTable, and Hail leaves missing entries out of aggregations as it does filtered ones. They are most of the
    # no-calls, so make them visible first.
    mt = mt.unfilter_entries()
    mt = mt.annotate_rows(_alts=split_alt_alleles(mt))
    mt = mt.annotate_rows(_shiftable=mt._alts.filter(lambda x: x.shiftable).map(lambda x: x.a_index))
    default_ploidy = hl.literal(2)
    if chrx_default_ploidy is not None:
        mt = mt.annotate_cols(_chrx_ploidy=hl.literal(chrx_default_ploidy)[mt._idx])
        default_ploidy = hl.if_else(mt.locus.in_x_nonpar(), mt._chrx_ploidy, 2)
    is_nocall = hl.is_missing(mt.GT)
    # Bits within a word run from the most significant down, so that printing each word as eight hex digits gives
    # the bitmap's bytes in order.
    bit = hl.bit_lshift(hl.int64(1), BITS_PER_WORD - 1 - mt._idx % BITS_PER_WORD)
    aggregations = dict(
        _alts=mt._alts,
        _shiftable=mt._shiftable,
        # Carriers of the alleles whose key is not its vid. Few rows have one, so this is empty nearly everywhere.
        _carriers=hl.agg.array_agg(
            lambda a: hl.agg.filter(hl.coalesce(mt.GT.contains_allele(a), False), hl.agg.collect(mt._idx)),
            mt._shiftable),
        an=hl.agg.call_stats(mt.GT, mt.alleles).AN,
        n_nocall=hl.int32(hl.agg.count_where(is_nocall)),
        n_nondefault_ploidy=hl.int32(hl.agg.count_where(~is_nocall & (mt.GT.ploidy != default_ploidy))),
        # Bounded: only used when the set is sparse, so there is never a need to collect a dense set's indices.
        _sparse=hl.agg.filter(is_nocall, hl.agg.take(mt._idx, max_sparse)),
        _words=hl.agg.filter(is_nocall, hl.agg.group_by(mt._idx // BITS_PER_WORD, hl.agg.sum(bit))),
    )
    if cohorts is not None:
        aggregations['expected_an'] = hl.agg.array_agg(
            lambda in_cohort: hl.agg.filter(in_cohort, hl.agg.call_stats(mt.GT, mt.alleles).AN), mt._cohorts)

    ht = mt.select_rows(**aggregations).rows()
    is_sparse = ht.n_nocall <= max_sparse
    ht = ht.annotate(
        nocall_sample_idxs=hl.or_missing(is_sparse, hl.sorted(ht._sparse)),
        nocall_bitmap_hex=hl.or_missing(
            ~is_sparse,
            hl.delimit(hl.range(n_words).map(lambda w: hl.format('%08x', ht._words.get(w, hl.int64(0)))), '')),
    )
    vid_contig = ht.locus.contig.replace('^chr', '')
    ht = ht.annotate(
        alt_alleles=ht._alts.map(lambda x: x.ref + '-' + x.alt),
        shifted_keys=ht._alts.filter(lambda x: x.shiftable).map(lambda x: hl.struct(
            vds_key=hl.format('%s-%s-%s-%s', vid_contig, ht.locus.position, x.ref, x.alt),
            ref=x.ref, alt=x.alt,
            carrier_idxs=hl.sorted(ht._carriers[ht._shiftable.index(x.a_index)]))),
    )
    ht = ht.annotate(
        contig=ht.locus.contig,
        position=ht.locus.position,
        location=hl.literal(CONTIG_INDEX).get(ht.locus.contig) * LOCATION_MULTIPLIER + ht.locus.position,
    )
    # The VDS is unsplit, so there is one row per locus and the alleles add nothing to the key.
    return ht.key_by('locus').drop('alleles', '_sparse', '_words', '_alts', '_shiftable', '_carriers')


def default_ploidy_summary(ht, n_samples):
    """Per contig, how many loci have AN equal to 2 * (called samples): every called sample diploid."""
    return ht.aggregate(hl.agg.group_by(ht.contig, hl.struct(
        loci=hl.agg.count(),
        sparse=hl.agg.count_where(hl.is_defined(ht.nocall_sample_idxs)),
        an_is_twice_called=hl.agg.count_where(ht.an == 2 * (n_samples - ht.n_nocall)),
        all_default_ploidy=hl.agg.count_where(ht.n_nondefault_ploidy == 0),
    )))


def write_shifted_keys(ht, output_path):
    """The keys bcftools norm has to turn into vids: their carriers as Parquet, and a sites-only VCF with the key in
    ID. Returns how many there are."""
    keys = ht.filter(hl.len(ht.shifted_keys) > 0).key_by()
    keys = keys.select('locus', 'location', 'shifted_keys').explode('shifted_keys')
    keys = keys.select('locus', 'location', **keys.shifted_keys)
    keys = keys.checkpoint(f'{output_path}/shifted_keys.ht', overwrite=True)
    keys.drop('locus', 'ref', 'alt').to_spark().write.mode('overwrite').parquet(f'{output_path}/shifted_keys.parquet')
    # About 10^4 keys genome-wide on Foxtrot, so collecting them is cheap; the carriers stay in the table. Exporting
    # from one partition avoids merging a part file per partition (see write_keys_to_normalize in the mapping script).
    sites = keys.select(locus=keys.locus, alleles=[keys.ref, keys.alt], rsid=keys.vds_key).collect()
    schema = hl.tstruct(locus=keys.locus.dtype, alleles=hl.tarray(hl.tstr), rsid=hl.tstr)
    hl.export_vcf(hl.Table.parallelize(sites, schema, n_partitions=1).key_by('locus', 'alleles'),
                  f'{output_path}/keys_to_normalize.vcf.bgz')
    return len(sites)


def main(vds, output_path, num_test_cohorts=0, seed=0):
    mt = filtered_dense_mt(vds)
    samples = sample_table(mt)
    n_samples = len(samples)
    cohorts = pick_test_cohorts(n_samples, num_test_cohorts, seed) if num_test_cohorts else None

    ht_path = f'{output_path}/nocall_sets.ht'
    nocall_sets(mt, n_samples, cohorts, [s['chrx_default_ploidy'] for s in samples]).write(ht_path, overwrite=True)
    ht = hl.read_table(ht_path)
    n_shifted_keys = write_shifted_keys(ht, output_path)
    ht.key_by().drop('locus', 'shifted_keys').to_spark().write.mode('overwrite').parquet(f'{output_path}/nocall_sets.parquet')

    with hl.hadoop_open(f'{output_path}/samples.tsv', 'w') as out:
        fields = list(samples[0])
        out.write('\t'.join(['sample_idx'] + fields) + '\n')
        for i, s in enumerate(samples):
            out.write('\t'.join([str(i)] + [str(s[f]) for f in fields]) + '\n')
    if cohorts is not None:
        with hl.hadoop_open(f'{output_path}/test_cohorts.json', 'w') as out:
            json.dump(cohorts, out)

    summary = {
        'n_samples': n_samples,
        'max_sparse': n_samples // BITS_PER_WORD,
        'shifted_keys': n_shifted_keys,
        'by_contig': {contig: dict(counts) for contig, counts in default_ploidy_summary(ht, n_samples).items()},
    }
    with hl.hadoop_open(f'{output_path}/summary.json', 'w') as out:
        json.dump(summary, out, indent=2, sort_keys=True)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Create per-locus no-call sets from a VDS')
    parser.add_argument('--vds_input_path', type=str, required=True, help='Input VDS path')
    parser.add_argument('--output_path', type=str, required=True, help='Directory to write outputs under')
    parser.add_argument('--temp_path', type=str, required=True, help='Path to temporary directory')
    parser.add_argument('--contigs', type=str,
                        help='Restrict to these contigs, comma-separated, e.g. chr21,chrX. Default: the whole VDS')
    parser.add_argument('--num_test_cohorts', type=int, default=0,
                        help='Also record call_stats AN for this many random cohorts')
    parser.add_argument('--seed', type=int, default=0, help='Seed for choosing the test cohorts')
    args = parser.parse_args()

    time_stamp = datetime.today().strftime('%Y-%m-%d_%H-%M-%S')
    hl.init(tmp_dir=f'{args.temp_path.rstrip("/")}/hail_tmp_create_nocall_sets_{time_stamp}')
    input_vds = hl.vds.read_vds(args.vds_input_path)
    if args.contigs:
        input_vds = hl.vds.filter_intervals(
            input_vds, [hl.parse_locus_interval(c, reference_genome='GRCh38') for c in args.contigs.split(',')])
    main(input_vds, args.output_path.rstrip('/'), args.num_test_cohorts, args.seed)
