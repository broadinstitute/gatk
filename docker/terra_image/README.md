# Terra image

TODO(@lucidtronix, @deflaux) update these instructions to match deployment to the public Docker image location

To build and push:
```
gcloud --project uk-biobank-sek-data builds submit \
  --timeout 20m \
  --tag gcr.io/uk-biobank-sek-data/ml4h_terra:`date +"%Y%m%d_%H%M%S"` .
```

Terra notebooks list which container to use. To update them all in-place, run a command similar to the following:
```
cd notebooks
find . -name "*.ipynb" -type f -print0 | \
  xargs -0 perl -i -pe \
  's/gcr.io\/uk-biobank-sek-data\/ml4h_terra:\d{8}_\d{6}/gcr.io\/uk-biobank-sek-data\/ml4h_terra:20200623_145127/g'
```
