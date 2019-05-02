# Run Dataflow
The following steps will run a Dataflow pipeline remotely, which in turn will, read the table `lubitz.coding`, 
calculate counts of `coding_file_id` and write them onto a GCS bucket.

* Clone the repo and cd into it:
```
    git clone git@github.com:broadinstitute/ml.git
    cd ml
```

* Create and activate the right Python environment:
```
    conda env create -f env/ml4cvd_osx64_dataflow.yml
    conda activate ml4cvd_osx64_dataflow
```

* Edit the parameters within the `USER VARIABLES` section in `tensorize/tensorize/defines.py` to update the fields as described in that file

* Run the application to submit the pipeline to Dataflow to be executed remotely provided the
`RUNNER` in `defines.py` is set to `DataflowRunner`). Set it to `DirectRunner` for local execution.
```
    python tensorize/tensorize_main.py
```
