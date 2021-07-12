FROM gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0

# Copy the application's requirements.txt and run pip to install
ADD requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt

# Add the application source code.
ADD create_cohort_extract_data_table.py /app
ADD create_variant_annotation_table.py /app

WORKDIR /app
ENTRYPOINT ["/bin/bash"]
