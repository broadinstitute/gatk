# Run Dataflow
The following steps will run a Dataflow pipeline remotely, which in turn will, tensorize fields of a type
specified by the user (e.g. categorical, continuous) and write them onto a GCS bucket in the form of
one `hd5` file per sample id.

* Clone the repo and cd into it:
```
    git clone git@github.com:broadinstitute/ml.git
    cd ml
```

* Create and activate the right Python environment:
```
    conda env create -f env/ml4h_dataflow.yml
    conda activate ml4h_dataflow
```

* Make sure you are authenticated by Google Cloud:
```
    gcloud auth application-default login
```

* Re install ml4h if you have made any changes:
```
    pip install .
```


* Run with the help option to see the list of command line arguments.
```
    python tensorize/tensorize_main.py -h
```

* **Note** that Google requires the `id` consist of only the
characters `[-a-z0-9]`, i.e. starting with a letter and ending with a letter or number.

* Run the application to submit the pipeline to Dataflow to be executed remotely provided the
command line argument `--beam_runner` is set to `DataflowRunner`. Set it to `DirectRunner` for local execution.
For example:
```
    python ml4h/tensorize/tensorize_main.py \
    --id my-pipeline-run \
    --tensor_type continuous \
    --bigquery_dataset ukbb_dev \
    --beam_runner DataflowRunner \
    --repo_root /Users/kyuksel/github/ml \
    --gcs_output_path my-pipeline-run/tensors_ukbb_dev_continous
```

* The pipeline can be run multiple times to tensorize different types of fields. This will populate the per-sample tensors
in specified GCS buckets. In order to unify them, they can be downloaded via `gsutil` as shown below
and merged using `merge_hd5s.py` script.
```
    gsutil -m cp -r <gcs bucket with tensors> <local directory>
```
