FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-slim as build

# Copy the application's requirements.txt and run pip to install
COPY requirements.txt /app/requirements.txt
RUN apt-get update && apt-get -y install cmake bcftools jq
RUN curl -O https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py
RUN pip3 install -r /app/requirements.txt

# Copy the application source code.
COPY *.py /app
COPY *.sql /app

# Copy the schema files needed for VAT creation.
COPY *.json /data/variant_annotation_table/schema/

WORKDIR /app
