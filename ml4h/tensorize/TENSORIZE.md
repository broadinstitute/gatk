# Tensorization
We create tensors from the raw UK Biobank data and feed those tensors into our models for training. 
The steps below describe how to perform tensorization in Google Cloud Platform (GCP). We assume you have created a VM,
installed the command line tool `google-cloud-sdk` on your laptop, and know how to log onto your VM. If this is not
the case, please refer to [Set Up a VM](#set-up-a-vm) on how to do those. You can use the `launch-instance.sh` script
as described in that section, to also create VMs with more CPUs and/or memory by specifying the
[instance type](https://cloud.google.com/compute/docs/machine-types) as an argument, like so:

```
./scripts/vm_launch/launch-instance.sh ${USER}-cpu n1-highcpu-64
``` 

## Create a disk
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

## Attach the disk to a VM
Running

```
gcloud compute instances attach-disk my-vm --disk my-disk --device-name my-disk --mode rw --zone us-central1-a
```

will attach `my-disk` to `my-vm`. 

Now we can format and mount the disk.   

## Format and mount the disk
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

## Tensorize away
It is recommended to run the tensorization in a screen multiplexer such as
[tmux](https://hackernoon.com/a-gentle-introduction-to-tmux-8d784c404340) which should come pre-installed on
your VM, especially for long running executions so your tensorization process is not halted when you disconnect from
your VM. To create a new session named, say `tensorization`, you can run `tmux new -s tensorization`. Once done going
through the steps below, pressing `ctrl+b d` will let you detach from your session. When you come back, run 
`tmux attach -t tensorization` to get back into it. See this [cheat sheet](https://tmuxcheatsheet.com) for more nifty
shortcuts.

Go to your `ml4h` clone's directory and run the tensorization script like below. To learn what options mean, feel
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

