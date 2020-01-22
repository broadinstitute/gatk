# Terra image

To build and push:
```
gcloud --project uk-biobank-sek-data builds submit \
  --timeout 20m \
  --tag gcr.io/uk-biobank-sek-data/ml4cvd_terra:`date +"%Y%m%d_%H%M%S"` .
```

After successful testing of the new image, if it is backwardly compatible, tag it as `latest`:
```
gcloud --project uk-biobank-sek-data container images add-tag \
  gcr.io/uk-biobank-sek-data/ml4cvd_terra:<YYYYMMDD_HHMMSS> \
  gcr.io/uk-biobank-sek-data/ml4cvd_terra:latest

```
