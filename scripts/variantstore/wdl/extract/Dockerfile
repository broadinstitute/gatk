FROM gcr.io/google.com/cloudsdktool/cloud-sdk:402.0.0

# Copy the application's requirements.txt and run pip to install
COPY requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt
RUN apt-get update && apt-get -y upgrade && apt-get -y install bcftools jq
RUN pip install firecloud

# Copy the application source code.
COPY create_variant_annotation_table.py /app
COPY create_ranges_cohort_extract_data_table.py /app
COPY extract_subpop.py /app
COPY populate_alt_allele_table.py /app
COPY alt_allele_positions.sql /app
COPY alt_allele_temp_function.sql /app
COPY utils.py /app
COPY add_max_as_vqslod.py /app
COPY curate_input_array_files.py /app
COPY scale_xy_bed_values.py /app
COPY workflow_compute_costs.py /app
COPY hail_gvs_import.py /app
COPY hail_create_vat_inputs.py /app

# Copy the schema files needed for VAT creation.
COPY variant_annotation_table/schema/vat_schema.json /data/variant_annotation_table/schema/
COPY variant_annotation_table/schema/vt_schema.json /data/variant_annotation_table/schema/
COPY variant_annotation_table/schema/genes_schema.json /data/variant_annotation_table/schema/

WORKDIR /app
