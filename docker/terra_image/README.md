# Terra image

To build and push:
```
mv ml4cvd ml4cvdBAK_$(date +"%Y%m%d_%H%M%S") \
  && mv config configBAK_$(date +"%Y%m%d_%H%M%S") \
  && cp -r ../../ml4cvd . \
  && cp -r ../vm_boot_images/config . \
  && gcloud --project uk-biobank-sek-data builds submit \
  --timeout 20m \
  --tag gcr.io/uk-biobank-sek-data/ml4cvd_terra:`date +"%Y%m%d_%H%M%S"` .
```
Note:

1. We're running a `cp` command to make the python package code files
available to docker.
    * We also place it in the home directory instead of pip installing the package to make it easier to edit the code.
    * See also [pyukbb.sh](https://github.com/broadinstitute/ml/blob/master/docker/vm_boot_images/config/pyukbb.sh)
      and [Intel MKL Dockerfile for local CPUs](https://github.com/broadinstitute/ml/pull/127/files)
      for alternate approaches to make the ml4cvd python package available within the environment.
1. We're commenting out `import vtk` to bypass a Docker build issue (see [Dockerfile](./Dockerfile)).
    * Since the code we're running on Terra does not currently make use of `vtk`, this has no impact.
    * But its messy to do this and the plan is to take another pass to reconcile the
      differences/similarties between this Docker base image and that of the ml4cvd
      vm boot images after the [Tensorflow 2.0](https://github.com/broadinstitute/ml/pull/94)
      update is in place.


After successful testing of the new image, if it is backwardly compatible, tag it as `latest`:
```
gcloud --project uk-biobank-sek-data container images add-tag \
  gcr.io/uk-biobank-sek-data/ml4cvd_terra:<YYYYMMDD_HHMMSS> \
  gcr.io/uk-biobank-sek-data/ml4cvd_terra:latest
```
