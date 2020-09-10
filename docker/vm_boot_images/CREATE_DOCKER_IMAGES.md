# Create and Use VM Boot Images
Everytime we want new persistent disks automatically added to our VMs, we need to update the `mount`s in
`scripts/vm_image/ml4cvd-image.sh`, and the `--disk` arguments in `scripts/vm_launch/launch_instance.sh` and
`scripts/vm_launch/launch_dl_instance.sh`, and then follow the steps listed in this section.

* Verify that `scripts/vm_image/ml4cvd-image.sh` has the
desired auto-mounting specified under the `# Mount the persistent disks` section that should look something like
    ```
    # Mount the persistent disks
    sudo mkdir -p /mnt/disks/data
    echo "UUID=3c62f761-3d8a-42ef-a029-1bfc6fd9be3f /mnt/disks/data ext4 ro,norecovery,discard,defaults,nofail" | sudo tee -a /etc/fstab
    ```

* Verify that `scripts/vm_image/dl-image.sh` has the desired disks specified to be attached, which should look something
like
    ```
    --disk=name=data,device-name=data,mode=ro,boot=no,auto-delete=no
    ```
* Make sure you have installed the
[google cloud tools (gcloud)](https://cloud.google.com/storage/docs/gsutil_install) on your laptop.
With [Homebrew](https://brew.sh/), you can use 
    ```
    brew cask install google-cloud-sdk
    ```

* Set up some environment variables to use throughout the rest of the section:
    ```
    export PROJECT=broad-ml4cvd
    export SERVICE_ACCOUNT=783282864357-compute@developer.gserviceaccount.com 
    export ZONE=us-central1-a
    export DATE=`date +%Y-%m-%d`
    export BOOT_DISK_SIZE=10GB
    export BOOT_DISK_TYPE=pd-standard
    export MACHINE_TYPE=n1-standard-1
    export BASE_IMAGE=ubuntu-1804-bionic-v20190429
    export BASE_IMAGE_PROJECT=ubuntu-os-cloud
    export CPU_IMAGE=ml4h-image
    export GPU_IMAGE=dl-image
    export CPU_VM=${USER}-create-cpu-image
    export CPU_TEST_VM=${USER}-test-cpu-image
    export GPU_VM=${USER}-create-gpu-image
    export GPU_TEST_VM=${USER}-test-gpu-image
    export ACCELERATOR=nvidia-tesla-k80
    ```

* If you want to create a **CPU-only (non-GPU)** VM image, set further environment variables as follows:
    ```
    export VM=${CPU_VM}
    export TEST_VM=${CPU_TEST_VM}
    export IMAGE=${CPU_IMAGE}
    export IMAGE_PROJECT=${BASE_IMAGE_PROJECT}
    ```
  If you want to create a **GPU** VM image, set those environment variables as follows:
    ```
    export VM=${GPU_VM}
    export TEST_VM=${GPU_TEST_VM}
    export IMAGE=${GPU_IMAGE}
    export IMAGE_PROJECT=${PROJECT}
    ```
  

* If you're creating a CPU image, create a fresh VM using an `Ubuntu` image such as `Ubuntu 18.04`. Note that not all `Ubuntu` 
images work; for example, `18.10` did not have `gcsfuse` as of 5/10/19.
    ```    
    gcloud compute instances create ${VM} \
        --project=${PROJECT} \
        --zone=${ZONE} \
        --machine-type=${MACHINE_TYPE} \
        --subnet=default \
        --maintenance-policy=MIGRATE \
        --service-account=${SERVICE_ACCOUNT} \
        --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append \
        --boot-disk-size=${BOOT_DISK_SIZE} \
        --boot-disk-type=${BOOT_DISK_TYPE} \
        --boot-disk-device-name=${VM} \
        --image-project=${IMAGE_PROJECT} \
        --image=${BASE_IMAGE}        
    ```
  If you're creating a GPU image, create a fresh VM via the command below (note that the last two lines differ from above):
    ```
    gcloud compute instances create ${VM} \
        --project=${PROJECT} \
        --zone=${ZONE} \
        --machine-type=${MACHINE_TYPE} \
        --subnet=default \
        --maintenance-policy=TERMINATE \
        --service-account=${SERVICE_ACCOUNT} \
        --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append \
        --boot-disk-size=${BOOT_DISK_SIZE} \
        --boot-disk-type=${BOOT_DISK_TYPE} \
        --boot-disk-device-name=${VM} \
        --image-project=${IMAGE_PROJECT} \
        --image-family=${CPU_IMAGE} \
        --accelerator=type=${ACCELERATOR},count=1
    ```

* Clone and go to the `ml` repo on your laptop:
    ```
    git clone git@github.com:broadinstitute/ml.git && cd ml
    ```

* Copy the scripts that will install the image content, to your VM
    ```
    gcloud compute scp scripts/vm_image/* ${VM}:/home/${USER}
    ```

* Log into your VM:
    ```
    gcloud --project ${PROJECT} compute ssh ${VM} --zone ${ZONE}
    ```

* If you're creating a **CPU** image, run the following script (**without sudo**):  
    ```
    ./ml4cvd-image.sh
    ```
  If you're creating a **GPU** image, first run:
    ```
    ./dl-image-part-1.sh
    ```
  After the script finishes, it reboots the VM so you can expect to see a message like:
    ```
    Connection to [IP_ADDRESS] closed by remote host
    ```
  It will have logged you out so first, log back in, and then run the part-2 script:
    ```
    ./dl-image-part-2.sh
    ```

* Delete the scripts from your VM:
    ```
    rm *.sh
    ```

* Exit out of the VM and stop it before attempting to create an image off of its boot disk:
    ```
    gcloud compute instances stop ${VM} \
        --project=${PROJECT} \
        --zone=${ZONE}
    ```

* Create the image:
    ```    
    gcloud compute images create ${IMAGE}-${DATE} \
        --project=${PROJECT} \
        --family=${IMAGE} \
        --source-disk=${VM} \
        --source-disk-zone=${ZONE}
    ```

* Verify that you can view your newly created image on GCP Console's [Images Page](https://console.cloud.google.com/compute/images?_ga=2.132530574.-1060415104.1522950615&project=broad-ml4cvd&folder&organizationId=548622027621&imagessize=50&imagesquery=%255B%255D).

* Delete the VM:
    ```
    gcloud -q compute instances delete ${VM} \
        --project=${PROJECT} \
        --zone=${ZONE}
    ```

* If you built a CPU base image, launch a test instance with the new image:
    ```
    scripts/vm_launch/launch_instance.sh ${TEST_VM}
    ```
  If you built a GPU image, run the following script instead:
    ```
    scripts/vm_launch/launch_dl_instance.sh ${TEST_VM}
    ``` 

* Login to `TEST_VM`:
    ```
    gcloud --project ${PROJECT} compute ssh ${TEST_VM} --zone ${ZONE}
    ```

* If you created a CPU image, that new image will be selected automatically because it is the latest
  one in the specified image family. Now, verify that the correct disk(s) are mounted and bucket(s) are `gcsfuse`d.
    ```
    ls /mnt/ml4h && ls /mnt/disks/*
    ``` 
    
* Set up your SSH keys on the VM for GitHub to be able to clone the `ml` repo.

* Clone and go to the `ml` repo on your laptop:
    ```
    git clone git@github.com:broadinstitute/ml.git && cd ml
    ```
    
* Because we don't know everyone's username, you need to run one more script to make sure
that you are added as a docker user and that you have permission to pull down our docker
instances from GCP's gcr.io. Run this while you're logged into your VM:
    ```
    scripts/vm_launch/run_once.sh
    ```

  Note that you may see warnings like below, but these are expected:
    ```
    WARNING: Unable to execute `docker version`: exit status 1
    This is expected if `docker` is not installed, or if `dockerd` cannot be reached...
    Configuring docker-credential-gcr as a registry-specific credential helper. This is only supported by Docker client versions 1.13+
    /home/username/.docker/config.json configured to use this credential helper for GCR registries
    ```
  You need to log out after that (`exit`) then ssh back in so everything takes effect.

* Run the tests (`ml/ml4h/DATA_MODELING_TESTS.md`)

* Exit out of the VM and delete it:
    ```
    gcloud -q compute instances delete ${TEST_VM} \
            --project=${PROJECT} \
            --zone=${ZONE}
    ```

