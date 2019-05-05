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
    conda env create -f env/ml4cvd_dataflow.yml
    conda activate ml4cvd_dataflow
```

* Edit the parameters within the `USER VARIABLES` section in `tensorize/tensorize/defines.py` to update the fields as described in that file

* Run the application to submit the pipeline to Dataflow to be executed remotely provided the
`RUNNER` in `defines.py` is set to `DataflowRunner`). Set it to `DirectRunner` for local execution.
```
    python tensorize/tensorize_main.py
```

* The pipeline can be run multiple times to tensorize different types of fields. This will populate the per-sample tensors
in specified GCS buckets. In order to unify them, they can be downloaded via `gsutil` as shown below
and merged using `merge_hd5s.py` script.
```
    gsutil -m cp -r <gcs bucket with tensors> <local directory>
```
