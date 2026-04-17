import argparse

import hail as hl


def remove_samples_from_vds(
        input_vds_path: str,
        samples_to_remove_path: str,
        output_vds_path: str,
) -> None:
    """
    Remove a set of samples from a Hail VDS and write the result to a new path.

    Steps performed:
    1. Filter out the specified samples (and remove dead alleles).
    2. Drop variant-data rows that no longer carry any alternate alleles
       (monomorphic reference rows left behind after sample removal).
    3. Recalculate GT from LGT + LA so that per-genotype GT values remain
       accurate after the cohort composition changes.

    Parameters
    ----------
    input_vds_path : str
        GCS (or local) path to the input VDS.
    samples_to_remove_path : str
        Path to a comma-delimited file listing research IDs to remove.
        The file must have a single column with a header of ``research_id``,
        one sample ID per line.
    output_vds_path : str
        GCS (or local) path where the filtered VDS will be written.
    """
    # Load the input VDS.
    input_vds = hl.vds.read_vds(input_vds_path)

    # Load the samples-to-remove table and key by 's' (the Hail sample-column key).
    samples_to_remove_table = hl.import_table(samples_to_remove_path, delimiter=',')

    # Fail fast if the removal file contains duplicate research_ids.
    n_to_remove = samples_to_remove_table.count()
    n_distinct = samples_to_remove_table.distinct().count()
    if n_to_remove != n_distinct:
        raise ValueError(
            f"samples_to_remove file contains {n_to_remove - n_distinct} duplicate research_id(s). "
            "Please deduplicate before proceeding."
        )

    samples_to_remove_table = samples_to_remove_table.key_by(s=samples_to_remove_table.research_id)

    # Compute the intersection of the removal list with the VDS columns so we know
    # how many samples will actually be removed, and print informative counts.
    n_input = input_vds.variant_data.cols().count()
    n_in_vds = samples_to_remove_table.semi_join(input_vds.variant_data.cols()).count()
    n_not_in_vds = n_to_remove - n_in_vds
    print(f"Input VDS sample count:                {n_input}")
    print(f"Samples in removal file:               {n_to_remove}")
    print(f"Samples in removal file found in VDS:  {n_in_vds}")
    print(f"Samples in removal file not in VDS:    {n_not_in_vds}")

    if n_in_vds == 0:
        raise ValueError(
            "None of the samples listed in the removal file are present in the VDS. Aborting."
        )

    # Step 1: filter samples and remove dead alleles.
    filtered_vds = hl.vds.filter_samples(
        input_vds,
        samples_to_remove_table,
        keep=False,
        remove_dead_alleles=True,
    )

    # Step 2: drop rows that no longer have any non-reference calls.
    # This removes vestigial monomorphic reference rows in the variant data
    # (e.g. sites where only the removed samples carried an alt allele).
    filtered_vds = hl.vds.VariantDataset(
        filtered_vds.reference_data,
        filtered_vds.variant_data.filter_rows(
            hl.agg.any(filtered_vds.variant_data.LGT.is_non_ref())
        ),
    )

    # Step 3: recalculate GT from LGT + LA.
    # Hail stores GT at import time; it is NOT updated automatically when
    # samples are removed, so we must recompute it ourselves.
    # annotate_entries overwrites GT in-place, preserving all other entry fields.
    vd = filtered_vds.variant_data
    vd = vd.annotate_entries(GT=hl.vds.lgt_to_gt(vd.LGT, vd.LA))
    filtered_vds = hl.vds.VariantDataset(filtered_vds.reference_data, vd)

    # Verify the arithmetic adds up using the intersection count.
    n_output = filtered_vds.variant_data.cols().count()
    print(f"Output VDS sample count:               {n_output}")
    assert n_output + n_in_vds == n_input, (
        f"Sample count mismatch: {n_output} + {n_in_vds} != {n_input}"
    )
    print("Sample count check passed.")

    # Write the filtered VDS.
    filtered_vds.write(output_vds_path)
    print(f"Filtered VDS written to: {output_vds_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            'Remove a set of samples from a Hail VDS. '
            'Monomorphic reference rows are dropped and GT is recalculated from LGT.'
        ),
    )

    parser.add_argument(
        '--input-vds-path',
        type=str,
        required=True,
        help='Path to the input VDS.',
    )
    parser.add_argument(
        '--samples-to-remove-path',
        type=str,
        required=True,
        help=(
            "Path to a file listing research IDs of samples to remove. "
            "Must be a comma-delimited file with a single column whose header is 'research_id', "
            "one sample ID per line."
        ),
    )
    parser.add_argument(
        '--output-vds-path',
        type=str,
        required=True,
        help='Path to write the output (filtered) VDS.',
    )
    parser.add_argument(
        '--temp-path',
        type=str,
        required=False,
        default=None,
        help='Optional path to a temporary directory for Hail intermediate files.',
    )

    args = parser.parse_args()

    hl.init(
        idempotent=True,
        tmp_dir=args.temp_path,
    )
    hl.default_reference('GRCh38')

    remove_samples_from_vds(
        input_vds_path=args.input_vds_path,
        samples_to_remove_path=args.samples_to_remove_path,
        output_vds_path=args.output_vds_path,
    )

