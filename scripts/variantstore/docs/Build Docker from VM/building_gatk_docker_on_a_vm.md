# Build the GATK Docker image on a cloud VM

The instructions here are written specifically for building an `ah_var_store` version of the GATK Docker image from an
Azure virtual machine, though much of this would likely apply to Google Cloud (or other clouds) and other Docker images.


## Allocate the VM

Through the [Azure portal](https://portal.azure.com/) allocate a VM in the Variants subscription. I used a Standard
E4-2ads v5 (2 vcpus, 32 GiB memory) for this purpose:

![Azure VM for building Docker image](./Azure%20VM%20for%20building%20Docker%20image.png)

Generate a new SSH key when prompted and save this to a safe location as you will need it soon. 

Once the machine has been created and started, go to its page in the Azure portal and select Connect -> SSH:


![Connect to Azure VM](./Azure%20VM%20Connect%20SSH.png)

SSH to the VM following the instructions on this page. Once connected, first become root as nearly every step requires
root access:

```
sudo bash
```

Now as root:

```
branch=<your branch name here>

apt-get update

# Install git and Docker dependencies
apt-get install --assume-yes git-core git-lfs

# Switch to 150 GiB data disk
cd /mnt
mkdir gitrepos && cd gitrepos
git clone https://github.com/broadinstitute/gatk.git --depth 1 --branch ${branch} --single-branch
cd gatk

# Run a helper script with lots more commands for building a GATK Docker image:
./scripts/variantstore/azure/gatk_docker_setup.sh

# Log in to Google Cloud
gcloud init

FULL_IMAGE_ID=$(cat /tmp/idfile.txt)

# Take the slice of this full Docker image ID that corresponds with the output of `docker images`:
IMAGE_ID=${FULL_IMAGE_ID:7:12}

# The GATK Docker image is based on gatkbase.
IMAGE_TYPE="gatkbase"
TAG=$(python3 ./scripts/variantstore/wdl/extract/build_docker_tag.py --image-id "${IMAGE_ID}" --image-type "${IMAGE_TYPE}")

BASE_REPO="broad-dsde-methods/gvs"
REPO_WITH_TAG="${BASE_REPO}/gatk:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"

# Configure the credential helper for GAR
gcloud auth configure-docker us-central1-docker.pkg.dev

# Tag and push
GAR_TAG="us-central1-docker.pkg.dev/${REPO_WITH_TAG}"
docker tag "${REPO_WITH_TAG}" "${GAR_TAG}"

docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""
```

Don't forget to shut down (and possibly delete) your VM once you're done!
