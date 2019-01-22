# How to update the GATK base docker image:

1. choose a new version number for the base image and manually update the version in `scripts/docker/gatkbase/build_docker_base.sh`
2. build the gatkbase image using that script and upload it to the [gatk-dev docker repo](https://hub.docker.com/r/broadinstitute/gatk-dev/) or [gcr-gatk-snapshots](us.gcr.io/broad-dsde-methods/broad-gatk-snapshots)
   * cd to scripts/docker/gatkbase
   * run `./build_docker_base.sh`
   * `docker tag broadinstitute/gatk-dev:your-version-rc1` or whatever the correct tag is for where you want it uploaded
   * `docker push tagname`
3. update the Dockerfile in the main gatk to use the image you pushed
4. commit the changes to the two docker files and to a new pull request
5. wait for the tests to pass and show it to a reviewer
6. push the base image to the official [gatk repo](https://hub.docker.com/r/broadinstitute/gatk) with the right name
7. update the main docker to point to the official version you just released
8. wait for tests to pass in travis and merge


