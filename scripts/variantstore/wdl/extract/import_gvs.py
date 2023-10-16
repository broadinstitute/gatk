import os
from typing import List

import hail as hl
from hail.vds.combiner.combine import merge_alleles
from hail.genetics.reference_genome import reference_genome_type
from hail.typecheck import typecheck, sequenceof, numeric


@typecheck(refs=sequenceof(sequenceof(str)),
           vets=sequenceof(sequenceof(str)),
           sample_mapping=sequenceof(str),
           site_filtering_data=sequenceof(str),
           vqsr_filtering_data=sequenceof(str),
           vqsr_tranche_data=sequenceof(str),
           final_path=str,
           tmp_dir=str,
           truth_sensitivity_snp_threshold=float,
           truth_sensitivity_indel_threshold=float,
           reference_genome=reference_genome_type,
           partitions_per_sample=numeric,
           intermediate_resume_point=int,
           skip_final_merge=bool,
           unphase=bool,
           ref_block_max_length=int,
           use_classic_vqsr=bool
           )
def import_gvs(refs: 'List[List[str]]',
               vets: 'List[List[str]]',
               sample_mapping: 'List[str]',
               site_filtering_data: 'List[str]',
               vqsr_filtering_data: 'List[str]',
               vqsr_tranche_data: 'List[str]',
               final_path: 'str',
               tmp_dir: 'str',
               truth_sensitivity_snp_threshold: 'float' = 0.997,
               truth_sensitivity_indel_threshold: 'float' = 0.990,
               reference_genome='GRCh38',
               partitions_per_sample=0.35,
               intermediate_resume_point=0,
               skip_final_merge=False,
               unphase=False,
               ref_block_max_length: 'int' = 1000,
               use_classic_vqsr=False
               ):
    """Import a collection of Avro files exported from GVS.

    This function is used to import Avro files exported from BigQuery for
    the GVS database used for the All of Us project. The resulting data type is
    a :class:`.VariantDataset`, which is a modern representation of sparse cohort-level
    data in Hail with reference blocks instead of a dense VCF-like representation.

    This function accepts inputs where the reference and variant data is broken into
    sample groups (with identical blocking for reference and variant data). Each sample
    group table is represented by a group of Avro files. Data must be sorted by the
    genomic coordinate (location) field both **within and between Avro files**. Order
    of samples at a genomic coordinate is not required.

    The ``tmp_dir`` argument should refer to a path in network visible storage
    (Google Bucket, etc) preferably with a lifecycle policy to delete temporary
    data after a short duration (e.g. 5 days).

    **Data transformations**

      - `refs` -- The `state` field is transformed to an int32 `GQ` field. The reference base
        is added in from a FASTA file to ensure compatibility with Hail functionality that
        requires an allele at every locus.
      - `vets` -- The `GT` field is parsed as a Hail :class:`.tcall` and renamed to `LGT`
        to denote its indexing by local alleles. The `alt` field is split by commas and recorded
        from strings to integers that refer to an index (0 being ref) into the array of all
        alleles discovered at a locus. The `GQ` and `RGQ` fields are converted to 32-bit integers.
      - `site_filtering_data` -- The site filters are converted from a comma-delimited string
        to a `set<str>` value in Hail. This set contains all unique filters applied, and is
        an empty set for loci with no record in the input site filtering table.
      - `vqsr_filtering_data` -- The VQSR data records information specific to a
        (locus, alternate allele) pair. This data is read in as a dictionary in the row
        scope of the resulting variant table, where keys of the dictionary are an alternate
        allele string, where the dictionary values are the full records from the input VQSR
        filtering table, minus the `ref` and `alt` fields.
      - `vqsr_tranche_data` -- The VQSR tranche data is recorded as an array of records in the
        globals of the resulting variant data table.

    Execution notes
    ---------------
    Currently this method executes three queries per sample block -- one to collect sample IDs,
    one to analyze sortedness of the reference table, and one to analyze sortedness of the variant
    table. These three queries are in addition to:
      - 3 queries executed to import the filtering/tranche data
      - 1 query to merge and write a temporary reference table
      - 1 query to merge and write a temporary variant table
      - 1 query to repartition and write the final reference component of the VDS.
      - 1 query to repartition and write the final variant component of the VDS.

    The three extra queries per sample block can be eliminated with the right information.
    The necessary information is (1) the sample IDs corresponding to each sample block,
    and (2) the chromosomal coordinate start/end of each Avro file.

    Note
    ----
    The reference genome must have a sequence file loaded. This can be added with
    `:meth:`.ReferenceGenome.add_sequence`.

    Parameters
    ----------
    refs
        Paths to reference Avro files. The outer list has one entry for each sample group,
        and the inner lists contain all files in a sample group.
    vets : List[List[str]]
        Paths to variant Avro files. The outer list has one entry for each sample group,
        and the inner lists contain all files in a sample group.
    sample_mapping : List[str]
        Paths to sample mapping Avro files.
    site_filtering_data : List[str]
        Paths to site filtering files.
    vqsr_filtering_data : List[str]
        Paths to VQSR filtering files.
    vqsr_tranche_data : List[str]
        Paths to VQSR tranche files.
    final_path : :class:`str`
        Desired path to final VariantDataset on disk.
    tmp_dir : :class:`str`
        Path to network-visible temporary directory/bucket for intermediate file storage.
    truth_sensitivity_snp_threshold : :class:`float`
        VQSR sensitivity threshold for SNPs.
    truth_sensitivity_indel_threshold : :class:`float`
        VQSR sensitivity threshold for Indels.
    reference_genome : :class:`str` or :class:`.ReferenceGenome`
        Name or object referring to reference genome.
    partitions_per_sample : :class:`int` or :class:`float`
        Number of partitions per sample in the final VDS. Can be fractional (total rounds down).
    intermediate_resume_point : :class:`int`
        Index at which to resume sample-group imports. Default 0 (entire import)
    skip_final_merge : :class:`bool`
        Skip final merge if true.
    unphase : :class:`bool`
        Unphase VET genotypes on final merge.
    ref_block_max_length : :class:`int`
        Maximum reference block length.
    use_classic_vqsr : :class:`bool`
        Expect input Avro files to have been generated from VQSR 'Classic' data

    Script workflow:
    ---------------
    * Load VQSR data into table.
    * Convert each 4k sample group into a VDS.
    * Run the combiner on the VDSes (which does a hierarchical merge).
    * Annotate with VQSR / filter.
    * Compute FT genotype annotation.
    * Remove phase on LGT.
    * Add GT (computed from LGT).
    * Write final VDS.

    Returns
    -------
    :class:`.VariantDataset`

    """


    from hail.utils.java import info
    vds_paths = []
    assert len(refs) == len(vets)

    if not reference_genome.has_sequence():
        raise ValueError(f"reference genome {reference_genome.name!r} has no sequence file."
                         f"\n  Add the sequence with `<rg>.add_sequence(...)`")

    def translate_locus(location):
        """Translates an int64-encoded locus into a locus object."""
        factor = 1000000000000
        chrom = hl.literal(reference_genome.contigs[:26])[hl.int32(location / factor) - 1]
        pos = hl.int32(location % factor)
        return hl.locus(chrom, pos, reference_genome=reference_genome)

    def translate_state(state_var):
        """Translates a char-encoded GQ to int32."""
        return hl.literal({'m': 0, '0': 0, '1': 10, '2': 20, '3': 30, '4': 40, '5': 50, '6': 60})[state_var]

    def convert_array_with_id_keys_to_dense_array(arr, ids, drop=[]):
        """Converts a coordinate-represented sparse array into a dense array used in MatrixTables."""
        sdict = hl.dict(arr.map(lambda x: (x.sample_id, x.drop('sample_id', *drop))))
        return hl.rbind(sdict, lambda sdict: ids.map(lambda x: sdict.get(x)))

    info('import_gvs: Importing and collecting sample mapping lookup table')

    samp = hl.import_avro(sample_mapping)
    sample_mapping_dict = samp.aggregate(hl.dict(hl.agg.collect((samp.sample_id, samp.sample_name))))

    site_path = os.path.join(tmp_dir, 'site_filters.ht')
    vqsr_path = os.path.join(tmp_dir, 'vqsr.ht')

    if intermediate_resume_point > 0:
        info('import_gvs: skipping site and VQSR filter import')
    else:
        info('import_gvs: Importing and writing site filters to temporary storage')
        site = hl.import_avro(site_filtering_data)
        site = site.transmute(
            locus=translate_locus(site.location),
            filters=hl.set(site.filters.split(','))
        )
        site = site.key_by('locus')
        site.write(site_path, overwrite=True)

        info('import_gvs: Importing and writing VQSR filter data to temporary storage')
        vqsr = hl.import_avro(vqsr_filtering_data)
        vqsr = vqsr.transmute(
            locus=translate_locus(vqsr.location)
        )
        vqsr = vqsr.key_by('locus')
        vqsr.write(vqsr_path, overwrite=True)

    if use_classic_vqsr:
        info('vqsr_classic: Loading tranche data')
        tranche = hl.import_avro(vqsr_tranche_data)

    n_samples = 0

    with hl._with_flags(use_new_shuffle='1'):
        for idx in range(len(refs)):
            ref_group = refs[idx]
            var_group = vets[idx]
            path = os.path.join(tmp_dir, f'sample_group_{idx+1}.vds')
            vds_paths.append(path)

            if idx < intermediate_resume_point:
                n_samples += hl.vds.read_vds(path).n_samples()
                info(f'import_gvs: skipping group {idx+1}/{len(refs)}...')
                continue

            info(f'import_gvs: scanning group {idx+1}/{len(refs)}...')
            ref_ht = hl.import_avro(ref_group)

            # Note -- availability of sample IDs statically would make import more efficient
            info(f'import_gvs: collecting sample IDs...')
            sample_ids = sorted(list(ref_ht.aggregate(hl.agg.collect_as_set(ref_ht.sample_id))))
            samples = [sample_mapping_dict[s] for s in sample_ids]
            samples_lit = hl.literal(samples, hl.tarray(hl.tstr))
            sample_ids_lit = hl.literal(sample_ids, hl.tarray(hl.tint32))
            n_new_samples = len(sample_ids)
            n_samples += n_new_samples
            assert n_new_samples == len(samples), (n_new_samples, len(samples))

            new_loc = translate_locus(ref_ht.location)

            # transform fields to Hail expectations (locus object, GQ int32, end as absolute instead of relative
            ref_ht = ref_ht.transmute(locus=new_loc,
                                      GQ=translate_state(ref_ht.state),
                                      END=new_loc.position + hl.int32(ref_ht.length) - 1)
            ref_ht = ref_ht.key_by('locus')
            ref_ht = ref_ht.group_by(ref_ht.locus).aggregate(data_per_sample=hl.agg.collect(ref_ht.row.drop('locus')))
            ref_ht = ref_ht.transmute(entries=convert_array_with_id_keys_to_dense_array(ref_ht.data_per_sample, sample_ids_lit))

            # vds column keys assume string sample IDs
            ref_ht = ref_ht.annotate_globals(col_data=samples_lit.map(lambda s: hl.struct(s=s)),
                                             ref_block_max_length=ref_block_max_length)
            ref_mt = ref_ht._unlocalize_entries('entries', 'col_data', col_key=['s'])

            var_ht = hl.import_avro(var_group)
            var_ht = var_ht.transmute(locus=translate_locus(var_ht.location),
                                      local_alleles=hl.array([var_ht.ref]).extend(var_ht.alt.split(',')),
                                      LGT=hl.parse_call(var_ht.GT),
                                      LAD=var_ht.AD.split(',').map(lambda x: hl.int32(x)),
                                      GQ=hl.int32(var_ht.GQ),
                                      RGQ=hl.int32(var_ht.RGQ))
            var_ht = var_ht.key_by('locus')
            var_ht = var_ht.group_by(var_ht.locus).aggregate(data_per_sample=hl.agg.collect(var_ht.row.drop('locus')))

            alleles_list = hl.array(hl.set(var_ht.data_per_sample.map(lambda x:x.local_alleles)))

            alleles_and_translation = merge_alleles(alleles_list)
            alleles = alleles_and_translation[0]
            allele_to_index = hl.dict(hl.enumerate(alleles, index_first=False))
            local_allele_lookup = hl.dict(hl.zip(alleles_list, hl.rbind(allele_to_index,
                                                                        lambda ai: alleles_and_translation[1].map(
                                                                            lambda norm_alleles: norm_alleles.map(
                                                                                lambda a: allele_to_index[a])))))
            var_ht = var_ht.annotate(
                alleles=alleles,
                local_allele_lookup=local_allele_lookup)

            var_ht = var_ht.transmute(entries=convert_array_with_id_keys_to_dense_array(var_ht.data_per_sample, sample_ids_lit))
            var_ht = var_ht.annotate_globals(col_data=samples_lit.map(lambda s: hl.struct(s=hl.str(s))))
            var_mt = var_ht._unlocalize_entries('entries', 'col_data', col_key=['s'])

            # replace 'local_alleles' strings with indices using our global lookup table
            var_mt = var_mt.transmute_entries(LA=var_mt.local_allele_lookup[var_mt.local_alleles])
            var_mt = var_mt.drop('local_allele_lookup')
            var_mt = var_mt._key_rows_by_assert_sorted('locus', 'alleles')

            info(f'import_gvs: writing intermediate VDS for sample group {idx+1} with {n_new_samples} samples...')
            vds = hl.vds.VariantDataset(ref_mt, var_mt)
            vds.write(path, overwrite=True)

    if skip_final_merge:
        info("import_gvs: skipping final merge")
        return

    # compute partitioning for the final VDS
    total_partitions = int(partitions_per_sample * n_samples)
    first_ref_mt = hl.read_matrix_table(hl.vds.VariantDataset._reference_path(vds_paths[0]))
    target_records = first_ref_mt.count_rows() // total_partitions
    info(f'import_gvs: using target_records (records per partition) of {target_records} for VDS merge')

    target_final_intervals = first_ref_mt._calculate_new_partitions(total_partitions)

    with hl._with_flags(no_whole_stage_codegen='1'):

        merge_tmp = os.path.join(tmp_dir, 'merge_tmp.vds')
        hl.current_backend().fs.rmtree(merge_tmp)
        info(f'import_gvs: calling Hail VDS combiner for merging {len(vds_paths)} intermediates')
        combiner = hl.vds.new_combiner(output_path=merge_tmp,
                                       vds_paths=vds_paths,
                                       target_records=target_records,
                                       temp_path=tmp_dir,
                                       use_genome_default_intervals=True)
        combiner.run()
        combined = hl.vds.read_vds(merge_tmp, intervals=target_final_intervals)

        rd = combined.reference_data
        vd = combined.variant_data

        # read site and vqsr data with same intervals for efficient joins
        site = hl.read_table(site_path, _intervals=target_final_intervals)
        vqsr = hl.read_table(vqsr_path, _intervals=target_final_intervals)

        vd = vd.annotate_rows(filters=hl.coalesce(site[vd.locus].filters, hl.empty_set(hl.tstr)))

        # vqsr ref/alt come in normalized individually, so need to renormalize to the dataset ref allele
        vd = vd.annotate_rows(as_vqsr = hl.dict(vqsr.index(vd.locus, all_matches=True)
                                                .map(lambda record: (record.alt + vd.alleles[0][hl.len(record.ref):], record.drop('ref', 'alt')))))

        if use_classic_vqsr:
            vd = vd.annotate_globals(tranche_data=tranche.collect(_localize=False),
                                     truth_sensitivity_snp_threshold=truth_sensitivity_snp_threshold,
                                     truth_sensitivity_indel_threshold=truth_sensitivity_indel_threshold)

            sorted_tranche_data = hl.sorted(vd.tranche_data, key=lambda x: x.truth_sensitivity)
            vd = vd.annotate_globals(snp_vqslod_threshold=
                                     sorted_tranche_data.filter(lambda x: (x.model == 'SNP') & (
                                             x.truth_sensitivity >= truth_sensitivity_snp_threshold))
                                     .head().min_vqslod
                                     ,
                                     indel_vqslod_threshold=sorted_tranche_data.filter(lambda x: (x.model == 'INDEL') & (
                                             x.truth_sensitivity >= truth_sensitivity_indel_threshold))
                                     .head().min_vqslod
                                     )

            is_snp = vd.alleles[1:].map(lambda alt: hl.is_snp(vd.alleles[0], alt))
            vd = vd.annotate_rows(
                allele_NO=vd.alleles[1:].map(
                    lambda allele: hl.coalesce(vd.as_vqsr.get(allele).yng_status == 'N', False)),
                allele_YES=vd.alleles[1:].map(
                    lambda allele: hl.coalesce(vd.as_vqsr.get(allele).yng_status == 'Y', True)),
                allele_is_snp=is_snp,
                allele_OK=hl._zip_func(is_snp, vd.alleles[1:],
                                        f=lambda is_snp, alt:
                                        hl.coalesce(vd.as_vqsr.get(alt).vqslod >=
                                                    hl.if_else(is_snp, vd.snp_vqslod_threshold, vd.indel_vqslod_threshold),
                                                    True))
            )
        else:
            vd = vd.annotate_globals(truth_sensitivity_snp_threshold=truth_sensitivity_snp_threshold,
                                     truth_sensitivity_indel_threshold=truth_sensitivity_indel_threshold)
            is_snp = vd.alleles[1:].map(lambda alt: hl.is_snp(vd.alleles[0], alt))
            vd = vd.annotate_rows(
                allele_NO=vd.alleles[1:].map(
                    lambda allele: hl.coalesce(vd.as_vqsr.get(allele).yng_status == 'N', False)),
                allele_YES=vd.alleles[1:].map(
                    lambda allele: hl.coalesce(vd.as_vqsr.get(allele).yng_status == 'Y', True)),
                allele_is_snp=is_snp,
                allele_OK=hl._zip_func(is_snp, vd.alleles[1:],
                                      f=lambda is_snp, alt:
                                      hl.coalesce(vd.as_vqsr.get(alt).calibration_sensitivity <=
                                                  hl.if_else(is_snp, vd.truth_sensitivity_snp_threshold, vd.truth_sensitivity_indel_threshold),
                                                  True))
            )

        lgt = vd.LGT
        la = vd.LA
        allele_NO = vd.allele_NO
        allele_YES = vd.allele_YES
        allele_OK = vd.allele_OK
        allele_is_snp=vd.allele_is_snp
        ft = (hl.range(lgt.ploidy)
              .map(lambda idx: la[lgt[idx]])
              .filter(lambda x: x != 0)
              .fold(lambda acc, called_idx: hl.struct(
            any_no=acc.any_no | allele_NO[called_idx - 1],
            any_yes=acc.any_yes | allele_YES[called_idx - 1],
            any_snp=acc.any_snp | allele_is_snp[called_idx - 1],
            any_indel=acc.any_indel | ~allele_is_snp[called_idx - 1],
            any_snp_ok=acc.any_snp_ok | (allele_is_snp[called_idx - 1] & allele_OK[called_idx - 1]),
            any_indel_ok=acc.any_indel_ok | (~allele_is_snp[called_idx - 1] & allele_OK[called_idx - 1]),
        ), hl.struct(any_no=False, any_yes=False, any_snp=False, any_indel=False, any_snp_ok=False, any_indel_ok=False)))

        vd = vd.annotate_entries(FT=~ft.any_no & (ft.any_yes | ((~ft.any_snp | ft.any_snp_ok) & (~ft.any_indel | ft.any_indel_ok))))

        vd = vd.drop('allele_NO', 'allele_YES', 'allele_is_snp', 'allele_OK')

        if unphase:
            vd = vd.annotate_entries(LGT=vd.LGT.unphase())
        vd = vd.select_entries(GT=hl.vds.lgt_to_gt(vd.LGT, vd.LA), **vd.entry)

        hl.vds.VariantDataset(
            reference_data=rd,
            variant_data=vd,
        ).write(final_path, overwrite=True)
