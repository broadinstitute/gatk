# Run Dataflow
* Clone the repo and cd into it:
```
    git clone git@github.com:broadinstitute/ml.git
    cd ml
```

* Activate the necessary branch:
```
    git checkout tensorize_with_dataflow
```

* Create and activate the right Python environment:
```
    conda env create -f env/ml4cvd_osx64_py35.yml
    source activate ml4cvd_py35
```

* Edit `tensorize/tensorize/defines.py` to update the file paths and the run name

* Submit the application to Dataflow to run remotely:
```
    python tensorize/tensorize_main.py
```
