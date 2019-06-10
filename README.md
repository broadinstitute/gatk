# ml4cvd
`ml4cvd` is a set of tooling to make it easy to work on UK Biobank data on a Google Cloud machine.

## How-to:

### Set up a VM
Clone this repo and cd into it:
```
git clone git@github.com:broadinstitute/ml.git && cd ml
```
Make sure you have installed the [google cloud tools (gcloud)](https://cloud.google.com/storage/docs/gsutil_install). With [Homebrew](https://brew.sh/), you can use 
```
brew cask install google-cloud-sdk
```

If you don't have your gcloud already configed -- set the project to broad-ml4cvd
```gcloud config set project broad-ml4cvd```

**TODO**: Refactor instance creation scripts to centralized location with other useful scripts.

To create a VM without a GPU run:
```
./scripts/vm_launch/launch-instance.sh ${USER}-cpu
```
With GPU (not recommended unless you need something beefy and expensive)
```
./scripts/vm_launch/launch-dl-instance.sh ${USER}-gpu
```
This will take a few moments to run, after which you will have a VM in the cloud.  Remember to shut it off from the command line or [console](https://console.cloud.google.com/compute/instances?project=broad-ml4cvd) when you are not using it!  

Now ssh onto your instance:
```
gcloud --project broad-ml4cvd compute ssh ${USER}-gpu --zone us-central1-a
```

Because we don't know everyone's username, you need to run one more script to make sure that you are added as a docker user and that you have permission to pull down our docker instances from GCP's gcr.io. Run this while you're logged into your VM:
```
/scripts/vm_launch/run-once.sh
```

Note that you may see warnings like below, but these are expected:
```
WARNING: Unable to execute `docker version`: exit status 1
This is expected if `docker` is not installed, or if `dockerd` cannot be reached...
Configuring docker-credential-gcr as a registry-specific credential helper. This is only supported by Docker client versions 1.13+
/home/username/.docker/config.json configured to use this credential helper for GCR registries
```

You need to log out after that (`exit`) then ssh back in so everything takes effect.


Next, clone this repo onto your instance (if you have Two-Factor authentication setup you need to generate an SSH key on your VM and add it to your github settings as described [here](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/#platform-linux)):
```
git clone git@github.com:broadinstitute/ml.git
cd ml
```

### Do deep learning with TensorFlow
Once you have a virtual machine and an environment setup it is time to start learning.
The first step is to create training data by writing tensors to the disk.  

To write tensors with default categorical and continuous phenotypes, and no MRI or EKG data, `source activate` your
`conda` environment and run:
```
${HOME}/ml/scripts/tf.sh ${HOME}/ml/ml4cvd/recipes.py --mode tensorize --tensors ${HOME}/my_tensors/ --max_sample_id 1003000 --mri_field_id  --xml_field_id
```
This should take about a minute to run and will output the SQL queries as well as the counts for the phenotype categories and responses that it finds.  Now let's train a model:
```
${HOME}/ml/scripts/tf.sh ${HOME}/ml/ml4cvd/recipes.py --mode train --tensors ${HOME}/my_tensors/ --input_tensors categorical-phenotypes-94 --output_tensors coronary_artery_disease_soft --id my_first_mlp_for_cvd
```
This model should achieve about 75% validation set accuracy on predicting from the phenotypes whether this person was labelled with an ICD code corresponding to cardivascular disease.

### Run a notebook
Now let's run a Jupyter notebook.  On your VM run:
```
${HOME}/ml/scripts/dl-jupyter.sh 
```
This will start a notebook server on your VM. If you a Docker error like
```
docker: Error response from daemon: driver failed programming external connectivity on endpoint agitated_joliot (1fa914cb1fe9530f6599092c655b7036c2f9c5b362aa0438711cb2c405f3f354): Bind for 0.0.0.0:8888 failed: port is already allocated.
```
overwrite the default port (8888) like so
```
${HOME}/ml/scripts/dl-jupyter.sh 8889
```
The command also outputs two command lines in red.
Copy the line that looks like this:
```
ssh -i ~/.ssh/google_compute_engine -nNT -L 8888:localhost:8888 <YOUR VM's IP ADDRESS>
```
Open a terminal on your local machine and paste that command.  

If you get a public key error run: `gcloud compute config-ssh`

Now open a browser on your laptop and go to the URL `http://localhost:8888`

### Run tests
The command below will run all (integration and unit) tests:
```
${HOME}/ml/scripts/tf.sh -t ${HOME}/ml/ml4cvd/tests.py
```

Most of the tests need a GPU-enabled machine so it's best to run them on your VM. Note that spurious failures with
assertion errors such as
```
FAIL: test_train_mri_sax_zoom (__main__.TestTrainingModels)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/kyuksel/ml4cvd/tests.py", line 220, in test_train_mri_sax_zoom
    self.assertAlmostEqual(performances[k], expected[k], delta=delta)
AssertionError: 0.2925531914893617 != 0.7606382978723405 within 0.4 delta

----------------------------------------------------------------------
```
are expected since models are not trained much in the interest of time and they may not have learnt enough. They can
still be valuable as smoke tests since they exercise a lot of the codebase.

### Set up a local development environment
Setting up a Python environment with the right dependencies can be challenging. The first step below enables running 
the `ml4cvd` Python code on a Mac (tested on `High Sierra`) with `Conda`. The second one is for those who prefer taking 
advantage of the `PyCharm` IDE.
 
#### Conda
* Download onto your laptop the Miniconda `bash` or `.pkg` installer for `Python 3.7` and `Mac OS X` 
from [here](https://conda.io/en/latest/miniconda.html), and run it. If you installed Python via a package manager
such as `Homebrew`, you may want to uninstall that first, to avoid potential conflicts.
* On your laptop, at the root directory of your `ml4cvd` GitHub clone, load the `ml4cvd` environment via
    ```
    conda env create -f envs/ml4cvd_osx64.yml
    ``` 
    If you get an error, try updating your `Conda` via
    ```
    sudo conda update -n base -c defaults conda
    ```
    The version used at the time of this writing was `4.6.1`.
* Activate the environment:
    ```
    source activate ml4cvd
    ``` 
You may now run code on your `Terminal`, like so
```
python recipes.py --mode ...
``` 
**Note** that *recipe*s require having the right input files in place and running them without proper inputs will not
yield meaningful results.   

#### PyCharm
* Install PyCharm either directly from [here](https://www.jetbrains.com/pycharm/download/#section=mac), or download 
the [Toolbox App](https://www.jetbrains.com/toolbox/app/) and have the app install PyCharm. The latter makes 
PyCharm upgrades easier. It also allows you to manage your JetBrains IDEs from a single place if you have multiple
(e.g. IntelliJ for Java/Scala).
* Launch PyCharm.
* (Optional) Import the custom [settings](https://drive.google.com/open?id=1YvNVgVEH-rzsCJtrJ0mCi1nyAxG8Xync) as 
described [here](https://www.jetbrains.com/help/pycharm/exporting-and-importing-settings.html).
* Open the project on PyCharm from the `File` menu by pointing to where you have your GitHub repo.
* Next, configure your Python interpreter to use the Conda environment you set up previously:
    * Open `Preferences` from `PyCharm -> Preferences...`.
    * On the upcoming `Preferences` window's left-hand side, expand `Project: ml4cvd` if it isn't already.
    * Highlight `Project Interpreter`.
    * On the right-hand side of the window, where it says `Project Interpreter`, find and select your `python`
    binary installed by `Conda`. It should be a path like `~/conda/miniconda3/envs/ml4cvd/bin/python` where `conda`
    is the directory you may have selected when installing `Conda`. 
    * For a test run:
        * Open `recipes.py` (shortcut `Shift+Cmd+N` if you imported the custom settings).
        * Right-click on `if __name__=='__main__'` and select `Run recipes`.
        * You can specify input arguments by expanding the `Parameters` text box on the window
         that can be opened using the menu `Run -> Edit Configurations...`.    

### Create and Use VM Boot Images
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
    export CPU_IMAGE=ml4cvd-image
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
    ls /mnt/ml4cvd && ls /mnt/disks/*
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

* [Run the tests](#run-tests)

* Exit out of the VM and delete it:
    ```
    gcloud -q compute instances delete ${TEST_VM} \
            --project=${PROJECT} \
            --zone=${ZONE}
    ```

## Data
When you are on an instance, you can access data as follows:

### Phenotypic SQLite database
Data for 500k people containing almost everything available in the UK Biobank Showcase:

```/mnt/disks/data/raw/sql/ukbb7089.r10data.db``` 

To access the data using `sqlite`:

```sqlite3 /mnt/disks/data/raw/sql/ukbb7089.r10data.db```

The data can also be accessed through [BigQuery](https://console.cloud.google.com/bigquery?project=broad-ml4cvd&p=broad-ml4cvd&page=project).


### Cardiac MRI
212,158 individual zip files in ~20k people. Dicom-formatted files inside:

```/mnt/disks/data/raw/mris/cardiac/*.zip```

### Liver MRI
10,132 individual zip files in ~10k people. Dicom-formatted files inside:

```/mnt/disks/data/raw/mris/liver/*.zip```

### ECG: XML
119,097 ECGGs (12-lead resting and 3-lead exercise):

```/mnt/disks/data/raw/ecgs/{resting,exercise}/*.xml```

### Direct Genotypes 
~800k/person:

```/mnt/imputed_v2```

### Imputed Genotypes
90 million/person:

```/mnt/imputed_v3```

### Running an instance
There are 3 instance types envisioned:
* Standard ml4cvd: for running R, python, etc: anything that does not require a GPU
* TensorFlow ml4cvd: for running TF on a GPU instance
* PyTorch ml4cvd: for running PyTorch on a GPU instance

If you don't need a GPU, we get a *lot* more bang for our buck by sticking with a non-GPU instance.

All instances are based on a Ubuntu 18.04 base image.

### Data availability
One of the key features of using these instances is that they come with pre-mapped data, including:
* The UK Biobank genotype and imputation data (~3 TB)
* The UK Biobank cardiac MRI and liver MRI data (~3 TB)
* Pre-parsed UK Biobank tabular phenotypic data (~0.1 TB)


### Preinstalled tooling
We've preinstalled common tools that we use often. We have also built Docker images to make it easy to start playing, not getting lost in dependency hell.
* All instances come with Jupyter notebooks configured as Docker instances, with key data pre-mapped with launch scripts
* Deep learning instances have the right NVidia drivers pre-installed
* Deep learning instances have Docker images with GPU-aware configurations that are ready to go with PyTorch or TensorFlow

## Tensorization
We create tensors from the raw UK Biobank data and feed those tensors into our models for training. 
The steps below describe how to perform tensorization in Google Cloud Platform (GCP). We assume you have created a VM,
installed the command line tool `google-cloud-sdk` on your laptop, and know how to log onto your VM. If this is not
the case, please refer to [Set Up a VM](#set-up-a-vm) on how to do those. You can use the `launch-instance.sh` script
as described in that section, to also create VMs with more CPUs and/or memory by specifying the
[instance type](https://cloud.google.com/compute/docs/machine-types) as an argument, like so:

```
./scripts/vm_launch/launch-instance.sh ${USER}-cpu n1-highcpu-64
``` 

### Create a disk
Our current practice is to write the output (tensors) on a persistent GCP disk. 
If you already have a disk you can write to, you can skip this section. 

Note that GCP limits disks from being
attached to more than one VM in `read-write` mode at any given time. Furthermore, if a VM is attached to a disk in
`read-write` mode, no other VM is allowed to attach to the disk, even in `read-only` mode.

The following command line creates a solid state drive (SSD) persistent disk of size `100GB` named `my-disk` 
(disk names must contain lowercase letters, numbers, or hyphens only):
```
gcloud beta compute disks create my-disk --project broad-ml4cvd --type pd-ssd --size 100GB --zone us-central1-a
```

If you don't have `gcloud Beta Commands` installed (in addition to `google-cloud-sdk`), you will get a prompt saying

```
You do not currently have this command group installed.  Using it 
requires the installation of components: [beta]
.
.
.
Do you want to continue (Y/n)?
```

If so, enter `Y`, and it will start creating the disk after installing `gcloud Beta Commands`. 

Upon successful disk creation, you will see a note saying

```
Created [https://www.googleapis.com/compute/beta/projects/broad-ml4cvd/zones/us-central1-a/disks/my-disk].
NAME          ZONE           SIZE_GB  TYPE    STATUS
my-disk       us-central1-a  100      pd-ssd  READY

New disks are unformatted. You must format and mount a disk before it
can be used.
```

We *will* format the disk but first, we need to attach it to our VM.

### Attach the disk to a VM
Running

```
gcloud compute instances attach-disk my-vm --disk my-disk --device-name my-disk --mode rw --zone us-central1-a
```

will attach `my-disk` to `my-vm`. 

Now we can format and mount the disk.   

### Format and mount the disk
To be able to format and mount a disk, we need to know what device name it is assigned to. To find that out,
let's log onto our VM:

```
gcloud --project broad-ml4cvd compute ssh my-vm --zone us-central1-a
```

and run 

```
ls -l /dev/disk/by-id
```

That will output something like

```
total 0
lrwxrwxrwx 1 root root  9 Feb 11 19:13 google-mri-october -> ../../sdb
lrwxrwxrwx 1 root root  9 Feb 15 21:42 google-my-disk -> ../../sdd
lrwxrwxrwx 1 root root  9 Feb 11 19:13 google-persistent-disk-0 -> ../../sda
lrwxrwxrwx 1 root root 10 Feb 11 19:13 google-persistent-disk-0-part1 -> ../../sda1
lrwxrwxrwx 1 root root  9 Feb 11 19:13 scsi-0Google_PersistentDisk_mri-october -> ../../sdb
lrwxrwxrwx 1 root root  9 Feb 15 21:42 scsi-0Google_PersistentDisk_my-disk -> ../../sdd
lrwxrwxrwx 1 root root  9 Feb 11 19:13 scsi-0Google_PersistentDisk_persistent-disk-0 -> ../../sda
lrwxrwxrwx 1 root root 10 Feb 11 19:13 scsi-0Google_PersistentDisk_persistent-disk-0-part1 -> ../../sda1
``` 

The line that contains the name of our disk (`my-disk`) 

```
lrwxrwxrwx 1 root root  9 Feb 15 21:42 google-my-disk -> ../../sdd
```

indicates that our disk was assigned to the device named `sdd`. The subsequent steps will use this device name.
Make sure to replace it with yours, if different.

Let's format `my-disk` via

```
sudo mkfs.ext4 -m 0 -F -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/sdd
```   

Next, we'll create a directory that will serve as the mount point for the new disk:

```
sudo mkdir -p /mnt/disks/my-disk
```

The directory name doesn't have to match the name of the disk but it's one less name to keep track of that way.

Finally, we can do the mounting:

```
sudo mount -o norecovery,discard,defaults /dev/sdd /mnt/disks/my-disk
```

We will also add the persistent disk to the `/etc/fstab` file so that the device automatically mounts again 
if/when the VM restarts:

```
echo UUID=`sudo blkid -s UUID -o value /dev/sdd` /mnt/disks/my-disk ext4 norecovery,discard,defaults,nofail 0 2 | sudo tee -a /etc/fstab
```

If you detach this persistent disk or create a snapshot from the boot disk for this instance, edit the `/etc/fstab`
file and remove the entry for the disk. Even with the `nofail` or `nobootwait` options in place,
keep the `/etc/fstab` file in sync with the devices that are attached to your instance.

Voilà! You now have a shiny new disk where you can persist the tensors you will generate next.

### Tensorize away
It is recommended to run the tensorization in a screen multiplexer such as
[tmux](https://hackernoon.com/a-gentle-introduction-to-tmux-8d784c404340) which should come pre-installed on
your VM, especially for long running executions so your tensorization process is not halted when you disconnect from
your VM. To create a new session named, say `tensorization`, you can run `tmux new -s tensorization`. Once done going
through the steps below, pressing `ctrl+b d` will let you detach from your session. When you come back, run 
`tmux attach -t tensorization` to get back into it. See this [cheat sheet](https://tmuxcheatsheet.com) for more nifty
shortcuts.

Go to your `ml4cvd` clone's directory and run the tensorization script like below. To learn what options mean, feel
free to run `scripts/tensorize.sh` without any arguments or with `-h`.

```
scripts/tensorize.sh  -t /mnt/disks/my-disk/ -i log -n 4 -s 1000000 -e 1030000 -x 20205 -m 20209
``` 

This creates and writes tensors to `/mnt/disks/my-disk` (logs can be found at `/mnt/disks/my-disk/log`)
by running `4` jobs in parallel (recommended to match that number to your VM's number of vCPUs) 
starting with the sample ID `1000000` and ending with the sample ID `1004000` using both the `EKG` (field ID `20205`) 
and the `MRI` (field ID `20209`) data. The tensors will be `.hd5` files named after corresponding sample IDs 
(e.g. `/mnt/disks/my-disk/1002798.hd5`).

If you're happy with the results and would like to share them with your collaborators, don't forget to make your disk
`read-only` so they can attach their VMs to it as well. You can do this by

* (on VM) unmounting your disk directory:
```
    sudo umount /mnt/disks/my-disk
```

* (on laptop) detaching the disk:
```
    gcloud compute instances detach-disk my-vm --disk my-disk --zone us-central1-a
```

* (on laptop) reattaching the disk as `read-only`:
```
    gcloud compute instances attach-disk my-vm --disk my-disk --device-name my-disk --mode ro --zone us-central1-a
```

* (on VM) re-mounting the disk:
```
    sudo mount -o norecovery,discard,defaults /dev/sdd /mnt/disks/my-disk
```

It is also possible to do the disk mode changing dance on [Google Cloud Console](https://console.cloud.google.com) by
* searching for your VM's name
* clicking `EDIT` from the top bar
* finding the disk under `Additional disks`
* selecting `Read only` next to the name of your disk

You still want to have `umount`ed your disk beforehand to avoid potential issues later.

Finally, to save cost, don't forget to stop your VM when you're done using it. 

Happy model training/tuning!
