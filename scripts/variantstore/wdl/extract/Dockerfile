FROM gcr.io/google.com/cloudsdktool/cloud-sdk:349.0.0

# Copy the application's requirements.txt and run pip to install
COPY requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt
RUN apt-get update && apt-get -y upgrade && apt-get -y install bcftools

# Copy the application source code.
COPY create_cohort_extract_data_table.py /app
COPY create_variant_annotation_table.py /app
COPY create_ranges_cohort_extract_data_table.py /app
COPY extract_subpop.py /app
COPY populate_alt_allele_table.py /app
COPY alt_allele_positions.sql /app
COPY alt_allele_temp_function.sql /app
COPY utils.py /app

WORKDIR /app
ENTRYPOINT ["/bin/bash"]
