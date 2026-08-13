import argparse

import hail as hl

from hail_gvs_util import filter_samples_and_remove_monomorphic_rows, load_samples_to_remove


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

    Parameters
    ----------
    input_vds_path : str
        GCS (or local) path to the input VDS.
    samples_to_remove_path : str
        Path to a single-column file listing research IDs to remove.
        The file must have a header of ``research_id``, one research ID per line.
    output_vds_path : str
        GCS (or local) path where the filtered VDS will be written.
    """
    # Load and validate the samples-to-remove table (keyed by 's').
    samples_to_remove_table = load_samples_to_remove(samples_to_remove_path)
    n_to_remove = samples_to_remove_table.count()

    # Load the input VDS.
    input_vds = hl.vds.read_vds(input_vds_path)

    # Compute the intersection of the removal list with the VDS columns so we know
    # how many samples will actually be removed, and print informative counts.
    input_cols = input_vds.variant_data.cols()
    n_input = input_cols.count()
    n_in_vds = samples_to_remove_table.semi_join(input_cols).count()
    n_not_in_vds = n_to_remove - n_in_vds
    print(f"Input VDS sample count:                {n_input}")
    print(f"Samples in removal file:               {n_to_remove}")
    print(f"Samples in removal file found in VDS:  {n_in_vds}")
    print(f"Samples in removal file not in VDS:    {n_not_in_vds}")

    if n_in_vds == 0:
        raise ValueError(
            "None of the samples listed in the removal file are present in the VDS. Aborting."
        )

    # Filter samples and remove monomorphic reference rows.
    filtered_vds = filter_samples_and_remove_monomorphic_rows(input_vds, samples_to_remove_table)

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
            'Monomorphic reference rows are dropped.'
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
            "Path to a single-column file listing research IDs of samples to remove. "
            "Must have a header of 'research_id', with one research ID per line."
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
