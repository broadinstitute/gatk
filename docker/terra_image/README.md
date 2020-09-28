# Terra image

To build and push:
```
mv ml4h ml4hBAK_$(date +"%Y%m%d_%H%M%S") \
  && mv config configBAK_$(date +"%Y%m%d_%H%M%S") \
  && cp -r ../../ml4h . \
  && cp -r ../vm_boot_images/config . \
  && gcloud --project uk-biobank-sek-data builds submit \
  --timeout 20m \
  --tag gcr.io/uk-biobank-sek-data/ml4h_terra:`date +"%Y%m%d_%H%M%S"` .
```
Notes:

1. We're running a `cp` command to make the python package code files
available to docker.
    * TODO(deflaux) instead clone from GitHub once the repository is public.
1. Terra notebooks list which container to use. To update them all in-place, run a command similar to the following:
```
cd notebooks
find . -name "*.ipynb" -type f -print0 | \
  xargs -0 perl -i -pe \
  's/gcr.io\/uk-biobank-sek-data\/ml4h_terra:\d{8}_\d{6}/gcr.io\/uk-biobank-sek-data\/ml4h_terra:20200623_145127/g'
```
