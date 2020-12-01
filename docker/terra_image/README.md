# Terra image

GitHub Action [docker-publish.yml](../../.github/workflows/docker-publish.yml) is used to publish a public copy of this container to [ghcr.io/broadinstitute/ml4h/ml4h_terra](https://github.com/orgs/broadinstitute/packages/container/package/ml4h%2Fml4h_terra).

If you wish to build your own container, you can use a command similar to the following to build and push to Google Container Registry:
```
gcloud --project YOUR-PROJECT-ID builds submit \
  --timeout 20m \
  --tag gcr.io/YOUR-PROJECT-ID/ml4h_terra:`date +"%Y%m%d_%H%M%S"` .
```

