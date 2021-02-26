FROM python:3.7

# Copy the application's requirements.txt and run pip to install
ADD requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt

# Add the application source code.
ADD raw_array_cohort_extract.py /app
ADD ngs_cohort_extract.py /app

# install google SDK
RUN curl -sSL https://sdk.cloud.google.com | bash

WORKDIR /app
