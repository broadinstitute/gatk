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
        Path to a CSV/TSV file listing research IDs to remove.
        The file must have a header row whose first column is ``research_id``.
        One sample ID per line.
    output_vds_path : str
        GCS (or local) path where the filtered VDS will be written.
    """
    # Load the input VDS.
    input_vds = hl.vds.read_vds(input_vds_path)

    # Load the samples-to-remove table and key by 's' (the Hail sample-column key).
    samples_to_remove_table = hl.import_table(samples_to_remove_path, delimiter=',')
    samples_to_remove_table = samples_to_remove_table.key_by(s=samples_to_remove_table.research_id)

    # Sanity-check counts before proceeding.
    n_input = input_vds.variant_data.cols().count()
    n_to_remove = samples_to_remove_table.count()
    print(f"Input VDS sample count:       {n_input}")
    print(f"Samples to remove:            {n_to_remove}")

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
    vd = filtered_vds.variant_data
    vd = vd.drop(vd.GT)
    vd = vd.select_entries(GT=hl.vds.lgt_to_gt(vd.LGT, vd.LA), **vd.entry)
    filtered_vds = hl.vds.VariantDataset(filtered_vds.reference_data, vd)

    # Verify the arithmetic adds up.
    n_output = filtered_vds.variant_data.cols().count()
    print(f"Output VDS sample count:      {n_output}")
    assert n_output + n_to_remove == n_input, (
        f"Sample count mismatch: {n_output} + {n_to_remove} != {n_input}"
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
            "Must be a comma-delimited file with a header whose first column is 'research_id', "
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

