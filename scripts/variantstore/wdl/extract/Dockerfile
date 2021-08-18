FROM gcr.io/google.com/cloudsdktool/cloud-sdk:349.0.0

# Copy the application's requirements.txt and run pip to install
ADD requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt
RUN apt-get update && apt-get -y upgrade && apt-get -y install bcftools

# Add the application source code.
ADD create_cohort_extract_data_table.py /app
ADD create_variant_annotation_table.py /app
ADD extract_subpop.py /app
ADD populate_alt_allele_table.py /app
ADD alt_allele_positions.sql /app
ADD alt_allele_temp_function.sql /app

WORKDIR /app
ENTRYPOINT ["/bin/bash"]
