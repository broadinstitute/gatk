from datetime import datetime
from google.cloud import storage

import csv
import hail as hl
import re
import tempfile


def parse_ancestry_file(ancestry_file):
    start = datetime.now()
    current_time = start.strftime("%H:%M:%S")
    print("Start Time =", current_time)
    # Now parse the ancestry file to get it ready for the subpopulation work
    sample_id_to_sub_population = {}

    reader = csv.reader(ancestry_file, delimiter='\t')
    next(reader) # skip header
    for row in reader:
        key = row[0]
        value = row[4]
        sample_id_to_sub_population[key] = value
    # print(sample_id_to_sub_population_map)
    # need to make it looks like: sample_id_to_sub_population =
    # {"ERS4367795":"eur","ERS4367796":"eas","ERS4367797":"eur","ERS4367798":"afr","ERS4367799":"oth","ERS4367800":"oth","ERS4367801":"oth","ERS4367802":"oth","ERS4367803":"oth","ERS4367804":"oth","ERS4367805":"oth"}
    return sample_id_to_sub_population


def download_ancestry_file(gcs_ancestry_file):
    """
    Download the specified ancestry file from GCS to a local temporary file. This temporary file should be explicitly
    deleted once we done with it.
    """
    client = storage.Client()
    gcs_re = re.compile("^gs://(?P<bucket_name>[^/]+)/(?P<blob_name>.*)$")
    match = gcs_re.match(gcs_ancestry_file)

    if not match:
        raise ValueError(f"'{ancestry_file}' does not look like a GCS path")

    bucket_name, blob_name = match.groups()
    bucket = client.get_bucket(bucket_name)
    blob = bucket.get_blob(blob_name)
    fd, temp_file = tempfile.mkstemp()
    # Close open descriptor, do not remove temporary file.
    fd.close()

    blob.download_to_filename(temp_file)
    return temp_file


def hard_filter_non_passing_sites(vds):
    """
    Hard filter out non-passing sites
    TODO ok wait, DO WE want to do this? I guess it depends on what we are handing off and when.
    note: no AC/AN and AF for filtered out positions
    """
    vd = vds.variant_data
    filtered_vd = vd.filter_rows(hl.len(vd.filters) == 0)
    # now we apply it back to the vds
    filtered_vds = hl.vds.VariantDataset(vds.reference_data, filtered_vd)
    return filtered_vds


def replace_lgt_with_gt(vds):
    """
    Replace LGT with GT for easier calculations later.
    """
    filtered_vd = vds.variant_data
    filtered_vd = filtered_vd.annotate_entries(GT=hl.vds.lgt_to_gt(filtered_vd.LGT, filtered_vd.LA))
    return hl.vds.VariantDataset(vds.reference_data, filtered_vd)


def failing_gts_to_no_call(vds):
    """
     Respect the FT flag by setting all failing GTs to a no call

     Logic for assigning non passing GTs as no-calls to ensure that AC,AN and AF respect the FT flag.

     FT is True => GT keeps its current value
     FT is False => GT assigned no-call
     FT is missing => GT keeps its current value
    """
    vd = vds.variant_data
    filtered_vd = vd.annotate_entries(GT=hl.or_missing(hl.coalesce(vd.FT, True), vd.GT))
    # TODO drop GT after it is used since it will not line up with the LGT
    return hl.vds.VariantDataset(vds.reference_data, filtered_vd)


def gq0_to_no_call(vds):
    """
    Turn the GQ0s into no calls so that ANs are correct
    """
    rd = vds.reference_data
    # would be better to drop these once its a dense mt?
    rd = rd.filter_entries(rd.GQ > 0)
    return hl.vds.VariantDataset(rd, vds.variant_data)


def mt_ac_an_af(mt, ancestry_file):
    """
    Create a DENSE MATRIX TABLE to calculate AC, AN, AF and TODO: Sample Count
    """

    sample_id_to_sub_population = parse_ancestry_file(ancestry_file)

    mt = mt.annotate_cols(pop=hl.literal(sample_id_to_sub_population)[mt.s])
    return mt.select_rows(
        ac_an_af=hl.agg.call_stats(mt.GT, mt.alleles),
        call_stats_by_pop=hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles))
    )


def vds_ac_an_af(mt, vds):
    """
    This is what we will use to create the TSV for Nirvana custom annotations
    Now we join this back to the VDS
    """
    qc_data = mt.rows()
    filtered_vd = vds.variant_data
    filtered_vd = filtered_vd.annotate_rows(ac_an_af=qc_data[filtered_vd.row_key])
    return hl.vds.VariantDataset(vds.reference_data, filtered_vd)


def write_sites_only_vcf(vds, sites_only_vcf_path):
    # TODO we will want to drop some cols because there's no reason to carry around some of this stuff
    hl.export_vcf(vds.variant_data.rows(), sites_only_vcf_path)


def write_vat_custom_annotations(mt, vat_custom_annotations_tsv_path):
    """
    Create the VAT inputs:
    1. Create a Sites only VCF
    2. Create a VAT TSV file with subpopulation data
    """
    # TODO create the desired TSV-- We're gonna pull in Tim for this
    # Do we need to track the dropped sites?
    # filter out sites with too many alt alleles and trim extraneous INFO and FORMAT fields
    # normalize, left align and split multi allelic sites to new lines, remove duplicate lines
    # filter out spanning deletions and variants with an AC of 0
    # drop GTs if we have not already since LGT doesn't match
    # figure out how to calculate the Sample Count (is it just AC_ref + AC_var?)
    #
    hl.export_table(mt.rows(), vat_custom_annotations_tsv_path)


def main(vds, ancestry_file_location, vat_custom_annotations_tsv_path):
    transforms = [
        hard_filter_non_passing_sites,
        replace_lgt_with_gt,
        failing_gts_to_no_call,
        gq0_to_no_call
    ]

    for transform in transforms:
        vds = transform(vds)

    mt = hl.vds.to_dense_mt(vds)

    with open(ancestry_file_location, 'r') as ancestry_file:
        mt = mt_ac_an_af(mt, ancestry_file)

    vds = vds_ac_an_af(mt, vds)
    write_sites_only_vcf(vds, sites_only_vcf_path)
    write_vat_custom_annotations(mt, vat_custom_annotations_tsv_path)


def annotate_entry_filter_flag(mt):
    """
    Annotate FT flag for entries
    """
    fail_case = 'FAIL'
    return mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))


def write_tie_out_vcf(vds):
    """
    This is for writing tieout VCFs for toy-sized data only. Do not use for AoU-sized data as it would write a giant
    MatrixTable!
    """
    mt = hl.vds.to_dense_mt(vds)
    mt = annotate_entry_filter_flag(mt)

    # Write VCF
    hl.export_vcf(mt, vcf_output_path)


if __name__ == '__main__':
    vds_path = "@VDS_INPUT_PATH@"
    ancestry_file = "@ANCESTRY_INPUT_PATH@"
    vcf_output_path = "@VCF_OUTPUT_PATH@"
    sites_only_vcf_path = "@SITES_ONLY_VCF_OUTPUT_PATH@"
    vat_custom_annotations_path = "@VAT_CUSTOM_ANNOTATIONS_OUTPUT_PATH@"

    vds = hl.vds.read_vds(vds_path)
    local_ancestry_file = download_ancestry_file(ancestry_file)

    main(vds, local_ancestry_file, vat_custom_annotations_path)
