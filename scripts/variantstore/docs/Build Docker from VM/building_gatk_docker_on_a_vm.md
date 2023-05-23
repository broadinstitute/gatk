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
apt-get update

# Install git and Docker dependencies
apt-get install --assume-yes git-core git-lfs

# Switch to 150 GiB data disk
cd /mnt
mkdir gitrepos && cd gitrepos
git clone https://github.com/broadinstitute/gatk.git --depth 1 --branch ah_var_store --single-branch
cd gatk

# Run a helper script with lots more commands for building a GATK Docker image:
./scripts/variantstore/azure/gatk_docker_setup.sh

# Log in to Google Cloud
gcloud init

# Configure the credential helper for GCR
gcloud auth configure-docker

# Tag and push
DOCKER_TAG=us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_$(date -I | sed 's/-/_/g')
DOCKER_IMAGE_ID=$(docker images | head -n 2 | tail -n 1 | awk '{print $3}')
docker tag $DOCKER_IMAGE_ID $DOCKER_TAG
docker push $DOCKER_TAG
```

Don't forget to shut down (and possibly delete) your VM once you're done!
