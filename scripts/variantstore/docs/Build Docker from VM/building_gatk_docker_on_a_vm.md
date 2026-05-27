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

# Configure the credential helper for GAR (only needs to be done once)
gcloud auth configure-docker us-central1-docker.pkg.dev

BASE_REPO="broad-dsde-methods/gvs"

# --- Lite image (gatk_docker: no Conda/ML stack, used for most GVS tasks) ---
FULL_IMAGE_ID_LITE=$(cat /tmp/idfile_lite.txt)
IMAGE_ID_LITE=${FULL_IMAGE_ID_LITE:7:12}
IMAGE_TYPE_LITE="gatkbase-lite"
TAG_LITE=$(python3 ./scripts/variantstore/scripts/build_docker_tag.py --image-id "${IMAGE_ID_LITE}" --image-type "${IMAGE_TYPE_LITE}")
REPO_WITH_TAG_LITE="${BASE_REPO}/gatk:${TAG_LITE}"
docker tag "${IMAGE_ID_LITE}" "${REPO_WITH_TAG_LITE}"
GAR_TAG_LITE="us-central1-docker.pkg.dev/${REPO_WITH_TAG_LITE}"
docker tag "${REPO_WITH_TAG_LITE}" "${GAR_TAG_LITE}"
docker push "${GAR_TAG_LITE}"
echo "Lite image pushed to \"${GAR_TAG_LITE}\""

# --- Heavy image (gatk_heavy_docker: full gatkbase with Conda/ML stack, used for VETS/VQSR) ---
FULL_IMAGE_ID_HEAVY=$(cat /tmp/idfile_heavy.txt)
IMAGE_ID_HEAVY=${FULL_IMAGE_ID_HEAVY:7:12}
IMAGE_TYPE_HEAVY="gatkbase"
TAG_HEAVY=$(python3 ./scripts/variantstore/scripts/build_docker_tag.py --image-id "${IMAGE_ID_HEAVY}" --image-type "${IMAGE_TYPE_HEAVY}")
REPO_WITH_TAG_HEAVY="${BASE_REPO}/gatk:${TAG_HEAVY}"
docker tag "${IMAGE_ID_HEAVY}" "${REPO_WITH_TAG_HEAVY}"
GAR_TAG_HEAVY="us-central1-docker.pkg.dev/${REPO_WITH_TAG_HEAVY}"
docker tag "${REPO_WITH_TAG_HEAVY}" "${GAR_TAG_HEAVY}"
docker push "${GAR_TAG_HEAVY}"
echo "Heavy image pushed to \"${GAR_TAG_HEAVY}\""
```

Don't forget to shut down (and possibly delete) your VM once you're done!
