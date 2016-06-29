Docker for GATK4-protected (Developers)
--------------------------------------

*This docker script is meant for the GATK4 Broad dev team to create and push docker images.*

*Users should get gatk-protected docker images from dockerhub (https://hub.docker.com/r/broadinstitute/gatk-protected/)*

In this directory lives a simple, *unsupported* docker script.

This doc does not give usage info on Docker and assumes you are already versed in its use.

Please see the excellent docker documentation if you need info about docker:  https://docs.docker.com/

Notes:
- Image is built on ``ubuntu:14.04``
- Once image is built, the gatk-protected.jar (symlink) is found in ``/root`` (i.e. ``$HOME``).
- HDF5 jni library is in ``/usr/lib/jni``, since this is the default for Ubuntu 14.04.  This has been added to the ``JAVA_LIBRARY_PATH`` environment variable.

#### Create docker image and push the image to the gatk-protected dockerhub (Seriously, only for Broad GATK4 dev team)

*Please allow 1.5 hours for completion*

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_TAG=1.0.0.0-alpha1.2.1

sudo bash build_docker.sh -e ${GITHUB_TAG} -p
```

#### Create docker image (do not push to dockerhub)

From this directory, run:

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_TAG=1.0.0.0-alpha1.2.1

sudo bash build_docker.sh -e ${GITHUB_TAG}

```

#### Create docker image from a github hash (pushing to dockerhub prohibited)

From this directory, run:

```bash
# REPLACE VALUE OF GITHUB_TAG WITH DESIRED VERSION
export GITHUB_HASH=e454ac88c7791c0f8b385b3e82138ec52c61ef48

sudo bash build_docker.sh -e ${GITHUB_HASH} -s
```

#### Run gradle test

*Note that the unit tests are run during a build of the Dockerfile*

Currently, this is a bit manual.

```bash
# bash is the default command for the docker image.
sudo docker pull broadinstitute/gatk-protected
sudo docker run -i -t broadinstitute/gatk-protected:latest

# On the docker prompt
cd /root/gatk-protected
./gradlew test

# To leave the docker prompt:
exit
```

#### See the GATK-protected version from the docker prompt
```bash
# Start a docker image
sudo docker pull broadinstitute/gatk-protected
sudo docker run -i -t broadinstitute/gatk-protected:latest

# In the image prompt:
cat GATK_PROTECTED_VERSION
```

#### Delete all containers
```bash
# Credit:  http://stackoverflow.com/questions/17236796/how-to-remove-old-docker-containers
sudo docker rm `sudo docker ps --no-trunc -aq`
```