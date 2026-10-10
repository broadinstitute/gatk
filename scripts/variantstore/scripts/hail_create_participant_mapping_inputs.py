"""
Participant mapping inputs from the VDS (VS-2013, VS-2041): every carrier genotype as a flat
(vds_key, person_id, zygosity) row, plus the few keys that need normalizing before they name a VAT VID.

Carriers are the genotypes the VAT counted. The entry transforms are imported from
hail_create_vat_inputs.py rather than copied, and applied in the same order, so a carrier here is a
carrier in gvs_all_sc. Reference blocks contribute hom-ref and nothing else, so this does not
densify; the VS-2041 chr21 and chrX probes matched the VAT's het, hom and hemi counts exactly this way.

The VDS key is contig-position-ref-alt after split_multi(left_aligned=True), with the chr prefix stripped,
which is the VAT VID for all but a few keys in 10^5. The exceptions are indels written somewhere other
than their left-most position. Every split key is minimized, so an indel carries its anchor base as the
shared first base of REF and ALT, and it can move left exactly when the last base of the inserted or
deleted sequence equals that anchor: 21-43583452-G-GG can, 21-43583450-A-AG cannot. Those keys are
written to a sites-only VCF with the key in ID, for bcftools norm to give each the VID the VAT gave it.
On chr21 this flags 248 keys, exactly the VDS keys that are not VAT VIDs. The test covers only SNVs and
anchored indels, so the job fails before the expensive write if any key is anything else.

Outputs under --output:
  keys_to_normalize.vcf.bgz   sites-only, key in ID, every flagged key whether or not it has carriers
  pairs.parquet               vds_key STRING, person_id INT64, zygosity STRING (het, hom or hemi), split by
                              Spark into needs_normalization=true/ and needs_normalization=false/ so the
                              few pairs on flagged keys can be loaded without scanning the rest
  summary.json                counts and timings

The pairs are checkpointed as a Hail Table before Spark converts them to Parquet. Spark aborts a whole
Parquet write if an executor is lost after its task is authorized to commit (OutputCommitCoordinator,
unconditionally, since a retry could duplicate that task's output), which preemptible workers make
likely over a long write. Hail's own writes are not subject to this, so the expensive pass over the VDS
completes once, and only the conversion is retried. A run that fails after the checkpoint resumes from it
when relaunched with the same arguments; the checkpoint is deleted once the Parquet is written.

The VDS must hold exactly the active samples in --sample-info, an export of sample_info: none withdrawn,
no controls, none unknown to it, and none missing. The VAT is built from the same VDS, so the mapping's
checks against it cannot see a wrong sample set. The job fails before its first pass over the entries
rather than filtering, which would leave the mapping disagreeing with the VAT.

Zygosity mirrors call_stats: a homozygote needs ploidy 2, so a haploid carrier (male chrX and chrY
outside the PARs) is hemi. Per key, n_hom == Hom and n_het + n_hom + n_hemi == AC - Hom (gvs_all_sc).

GvsCreateParticipantMappingFromVDS.wdl runs this on a Dataproc cluster, then builds the mapping tables
from its outputs. Its interval input, e.g. chr21, restricts the run to one locus interval.
"""
import argparse
import hashlib
import json
import sys
import time

import hail as hl

from hail_create_vat_inputs import (
    failing_gts_to_no_call,
    hard_filter_non_passing_sites,
    remove_too_many_alt_allele_sites,
)

SHIFTABLE_CLASSES = ['insertion_shiftable', 'deletion_shiftable']
SUPPORTED_CLASSES = ['snv', 'insertion', 'deletion'] + SHIFTABLE_CLASSES
PARQUET_WRITE_ATTEMPTS = 3


def carrier_entries(vds_path, intervals):
    vds = hl.vds.read_vds(vds_path)
    if intervals:
        vds = hl.vds.filter_intervals(
            vds, [hl.parse_locus_interval(i, reference_genome='GRCh38') for i in intervals],
            split_reference_blocks=False)
    # Same three transforms, in the same order, as hail_create_vat_inputs.py:187-191.
    for transform in [remove_too_many_alt_allele_sites, hard_filter_non_passing_sites, failing_gts_to_no_call]:
        vds = transform(vds)

    vd = vds.variant_data.select_rows().select_entries('GT')
    # Matches hail_create_vat_inputs.py:136. Entries keep the unsplit GT, so carriers are tested against a_index.
    vd = hl.split_multi(vd, left_aligned=True)
    # The VAT drops spanning deletions before it assigns VIDs (hail_create_vat_inputs.py:172).
    vd = vd.filter_rows(vd.alleles[1] != '*')
    # Nirvana cannot handle N, so RemoveDuplicatesFromSitesOnlyVCF in GvsCreateVATfromVDS.wdl drops any site with an N in
    # REF before normalizing; such a key never becomes a VAT VID.
    vd = vd.filter_rows(~vd.alleles[0].contains('N'))
    ref, alt = vd.alleles[0], vd.alleles[1]
    # VAT VID format: contig without the chr prefix, X/Y as letters.
    return vd.annotate_rows(vds_key=hl.format('%s-%s-%s-%s', vd.locus.contig.replace('^chr', ''),
                                              vd.locus.position, ref, alt),
                            key_class=key_class(ref, alt))


def key_class(ref, alt):
    def last(s):
        return s[hl.len(s) - 1]

    return (hl.case()
            .when((hl.len(ref) == 1) & (hl.len(alt) == 1), 'snv')
            .when((hl.len(ref) > 1) & (hl.len(alt) > 1), 'mnp_or_complex')
            .when(ref[0] != alt[0], 'unanchored')
            .when(hl.len(alt) > hl.len(ref), hl.if_else(last(alt) == ref[0], 'insertion_shiftable', 'insertion'))
            .default(hl.if_else(last(ref) == alt[0], 'deletion_shiftable', 'deletion')))


def write_keys_to_normalize(vd, out):
    """Classify every key, fail on any the anchor test does not cover, and export the shiftable ones."""
    rows = vd.rows()
    classes = rows.aggregate(hl.agg.counter(rows.key_class))
    unsupported = {c: n for c, n in classes.items() if c not in SUPPORTED_CLASSES}
    if unsupported:
        examples = rows.filter(~hl.literal(SUPPORTED_CLASSES).contains(rows.key_class)).vds_key.take(10)
        raise ValueError(f'keys the left-alignment test does not cover: {unsupported}, e.g. {examples}')

    to_normalize = rows.filter(hl.literal(SHIFTABLE_CLASSES).contains(rows.key_class))
    hl.export_vcf(to_normalize.select(rsid=to_normalize.vds_key), f'{out}/keys_to_normalize.vcf.bgz')
    return classes


def carrier_pairs(vd):
    a = vd.a_index
    gt = vd.GT
    # Mirrors call_stats: a homozygote needs ploidy > 1 (CallStatsAggregator.scala:201).
    zygosity = (hl.case()
                .when(gt.ploidy == 1, 'hemi')
                .when(gt.ploidy == 2, hl.if_else((gt[0] == a) & (gt[1] == a), 'hom', 'het'))
                .or_error('expected ploidy 1 or 2 for a carrier, found: ' + hl.str(gt.ploidy)))
    # Annotate before filtering: filter_entries returns a new MatrixTable, which a, gt and zygosity do not refer to.
    # or_missing evaluates zygosity only for carriers, so or_error cannot fire on anyone else.
    vd = vd.annotate_entries(zygosity=hl.or_missing(gt.contains_allele(a), zygosity))
    vd = vd.filter_entries(hl.is_defined(vd.zygosity))
    et = vd.entries().key_by()
    return et.select(vds_key=et.vds_key, person_id=hl.parse_int64(et.s), zygosity=et.zygosity,
                     needs_normalization=hl.literal(SHIFTABLE_CLASSES).contains(et.key_class))


def write_pairs_parquet(pairs, out):
    """Convert the checkpointed pairs to Parquet, retrying only the commit-coordination abort. mode='overwrite'
    clears a failed attempt's partial output first, so a retry cannot leave duplicates behind."""
    for attempt in range(1, PARQUET_WRITE_ATTEMPTS + 1):
        try:
            pairs.to_spark().write.partitionBy('needs_normalization').parquet(f'{out}/pairs.parquet', mode='overwrite')
            return attempt
        except Exception as e:
            if 'Authorized committer' not in str(e) or attempt == PARQUET_WRITE_ATTEMPTS:
                raise
            print(f'Parquet write attempt {attempt} lost an authorized committer; retrying', file=sys.stderr)


def check_sample_set(vds_samples, sample_info_path, summary):
    """Fail unless the VDS holds exactly the active samples in sample_info."""
    # The VDS column key is sample_info.sample_name (import_gvs.py), so the names compare directly.
    status, is_loaded = {}, {}
    with hl.hadoop_open(sample_info_path) as f:
        for line in f:
            if line.strip():
                name, status[name], is_loaded[name] = line.rstrip('\n').split('\t')
    summary['sample_info'] = {s: sum(1 for v in status.values() if v == s) for s in ['active', 'withdrawn', 'control']}
    active = {name for name, s in status.items() if s == 'active'}
    # is_loaded explains a missing sample but does not decide what is expected: Parquet callsets once left it unset.
    problems = {
        'withdrawn in the VDS': sorted(n for n in vds_samples if status.get(n) == 'withdrawn'),
        'controls in the VDS': sorted(n for n in vds_samples if status.get(n) == 'control'),
        'in the VDS but not in sample_info': sorted(vds_samples - status.keys()),
        'active, is_loaded, missing from the VDS': sorted(n for n in active - vds_samples if is_loaded[n] == 'true'),
        'active, not is_loaded, missing from the VDS': sorted(n for n in active - vds_samples if is_loaded[n] != 'true'),
    }
    problems = {k: v for k, v in problems.items() if v}
    if problems:
        raise ValueError('the VDS samples are not the active samples in sample_info: ' +
                         '; '.join(f'{len(v)} {k}, e.g. {v[:10]}' for k, v in problems.items()))


def main():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('--vds', required=True)
    parser.add_argument('--output', required=True, help='GCS prefix for the outputs and the summary')
    parser.add_argument('--temp-path', required=True)
    parser.add_argument('--sample-info', required=True,
                        help='TSV of sample_name, status (active, withdrawn or control) and is_loaded, without a header')
    parser.add_argument('--interval', action='append', default=None,
                        help='Locus interval, repeatable (default: the whole VDS)')
    args = parser.parse_args()
    out = args.output.rstrip('/')

    hl.init(tmp_dir=f"{args.temp_path.rstrip('/')}/hail_tmp_participant_mapping_inputs")
    vd = carrier_entries(args.vds, args.interval)

    summary = {'vds': args.vds, 'intervals': args.interval, 'hail_version': hl.version(),
               'variant_partitions': vd.n_partitions()}

    # person_id is mt.s parsed as an integer; a name that does not parse would silently become a null.
    cols = vd.cols()
    summary['samples'], unparseable = cols.aggregate(
        (hl.agg.count(), hl.agg.filter(hl.is_missing(hl.parse_int64(cols.s)), hl.agg.take(cols.s, 10))))
    if unparseable:
        raise ValueError(f'sample names that are not integer person IDs: {unparseable}')

    check_sample_set(set(cols.s.collect()), args.sample_info, summary)

    # Resume state for this run's arguments, deleted on success, so a later run never picks up another's output.
    run_id = hashlib.sha256(json.dumps([args.vds, args.interval, out]).encode()).hexdigest()[:16]
    resume_dir = f"{args.temp_path.rstrip('/')}/participant_mapping_inputs_resume/{run_id}"
    classes_path, checkpoint_path = f'{resume_dir}/key_classes.json', f'{resume_dir}/pairs.ht'

    start = time.time()
    if hl.hadoop_exists(classes_path):
        with hl.hadoop_open(classes_path) as f:
            summary['key_classes'] = json.load(f)
        summary['classify_resumed'] = True
    else:
        summary['key_classes'] = write_keys_to_normalize(vd, out)
        with hl.hadoop_open(classes_path, 'w') as f:
            json.dump(summary['key_classes'], f)
    summary['keys_to_normalize'] = sum(summary['key_classes'].get(c, 0) for c in SHIFTABLE_CLASSES)
    summary['classify_seconds'] = round(time.time() - start, 1)

    start = time.time()
    if hl.hadoop_exists(f'{checkpoint_path}/_SUCCESS'):
        pairs = hl.read_table(checkpoint_path)
        summary['pairs_resumed'] = True
    else:
        pairs = carrier_pairs(vd).checkpoint(checkpoint_path, overwrite=True)
    summary['pairs_checkpoint_seconds'] = round(time.time() - start, 1)

    start = time.time()
    summary['parquet_write_attempts'] = write_pairs_parquet(pairs, out)
    summary['parquet_seconds'] = round(time.time() - start, 1)
    hl.current_backend().fs.rmtree(resume_dir)

    # Hail's progress bar is left mid-line on stderr; end it so the JSON starts on its own line.
    sys.stderr.write('\n')
    sys.stderr.flush()
    print(json.dumps(summary, indent=2))
    with hl.hadoop_open(f'{out}/summary.json', 'w') as f:
        json.dump(summary, f, indent=2)


if __name__ == '__main__':
    main()
