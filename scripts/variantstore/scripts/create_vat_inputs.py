import csv
import re

from google.cloud import storage
import tempfile


def parse_ancestry_file(ancestry_file):
    """
    Parse the specified TSV input file to create a sample to ancestry dictionary. The result should look
    something like:

    {"ERS4367795":"eur","ERS4367796":"eas","ERS4367797":"eur","ERS4367798":"afr","ERS4367799":"oth", ... }

    """
    sample_id_to_sub_population = {}

    reader = csv.reader(ancestry_file, delimiter='\t')
    # skip header
    next(reader)
    for row in reader:
        key = row[0]
        value = row[4]
        sample_id_to_sub_population[key] = value

    return sample_id_to_sub_population


def download_ancestry_file(gcs_ancestry_file):
    """
    Download the specified ancestry file from GCS to a local temporary file. This temporary file should be explicitly
    deleted once we are done with it.
    """
    client = storage.Client()
    gcs_re = re.compile("^gs://(?P<bucket_name>[^/]+)/(?P<blob_name>.*)$")
    match = gcs_re.match(gcs_ancestry_file)

    if not match:
        raise ValueError(f"'{gcs_ancestry_file}' does not look like a GCS path")

    bucket_name, blob_name = match.groups()
    bucket = client.get_bucket(bucket_name)
    blob = bucket.get_blob(blob_name)
    fd, temp_file = tempfile.mkstemp()
    # Close open descriptor, do not remove temporary file.
    # fd.close()

    blob.download_to_filename(temp_file)
    return temp_file
