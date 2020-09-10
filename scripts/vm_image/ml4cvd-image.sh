#!/usr/bin/env bash

# server-conf-scripts are for configuration of a *fresh* VM and should not be
# treated as startup scripts. (They are not idempotent.)

GCP_BUCKET="ml4h"

# We assume we are running as a regular user, not root.

# Enable gcsfuse to allow mounting of the google storage bucket as if it were a drive
export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -

# Install frequently-used packages
# First, update apt since we have added new repos (above)
sudo apt update

sudo apt install -y r-base r-base-core unzip wget bzip2 python sqlite3 gcsfuse

# Make gcsfuse auto-mount to /mnt/${GCP_BUCKET} in the future. Modify fstab to
# do this automatically. Via
# https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/mounting.md
# and https://serverfault.com/a/830726/118452 to enable easier mount with read and
# write access by non-root users.
echo "${GCP_BUCKET} /mnt/${GCP_BUCKET} gcsfuse rw,allow_other,implicit_dirs,default_permissions,file_mode=777,dir_mode=777" | sudo tee -a /etc/fstab
echo "fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b /mnt/imputed_v2 gcsfuse ro,allow_other,implicit_dirs,default_permissions,file_mode=777,dir_mode=777" | sudo tee -a /etc/fstab
echo "fc-7d5088b4-7673-45b5-95c2-17ae00a04183 /mnt/imputed_v3 gcsfuse ro,allow_other,implicit_dirs,default_permissions,file_mode=777,dir_mode=777" | sudo tee -a /etc/fstab

# Mount the persistent disks
sudo mkdir -p /mnt/disks/survey-tensors2
sudo mkdir -p /mnt/disks/ecg-text3
sudo mkdir -p /mnt/disks/pix-size-tensors
echo "UUID=65bba926-210b-48ee-aa0e-241599fad8d5 /mnt/disks/survey-tensors2 ext4 ro,norecovery,discard,defaults,nofail 0 2" | sudo tee -a /etc/fstab
echo "UUID=f27bb394-9fcd-41a8-8374-333cae177af8 /mnt/disks/ecg-text3 ext4 ro,norecovery,discard,defaults,nofail 0 2" | sudo tee -a /etc/fstab
echo "UUID=46f2f929-44d4-4925-800e-ec08bf3a5a92 /mnt/disks/pix-size-tensors ext4 ro,norecovery,discard,defaults,nofail 0 2" | sudo tee -a /etc/fstab

# Other packages that jpp uses
sudo /usr/bin/env Rscript -<<EOF
list.of.packages <- c('ggplot2','poweRlaw', 'Hmisc', 'speedglm', 'data.table', 'CMplot')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
EOF

#
# 2018/08/16 additions
#

# Enable docker (assumes Ubuntu, of any supported version)
# See https://docs.docker.com/install/linux/docker-ce/ubuntu/#set-up-the-repository
sudo apt-get remove -y docker docker-engine docker.io
sudo apt install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt update -y
sudo apt install -y docker-ce
sudo systemctl enable docker
sudo groupadd -f docker

# Manually install gcr
# Via https://cloud.google.com/container-registry/docs/advanced-authentication#standalone_docker_credential_helper
VERSION=1.5.0
OS=linux
ARCH=amd64
curl -fsSL "https://github.com/GoogleCloudPlatform/docker-credential-gcr/releases/download/v${VERSION}/docker-credential-gcr_${OS}_${ARCH}-${VERSION}.tar.gz" \
  | tar xz --to-stdout ./docker-credential-gcr | sudo tee -a /usr/bin/docker-credential-gcr 1>/dev/null && sudo chmod +x /usr/bin/docker-credential-gcr
docker-credential-gcr configure-docker

#
# 2018/08/26 additions
#
sudo apt-get install -y python-setuptools

# Enable dsub

# Note: for now, need to use my repo because it supports min-cpu-platform
#git clone https://github.com/DataBiosphere/dsub
git clone https://github.com/carbocation/dsub
cd dsub

# Note: for now, need to use my add-min-cpu-platform because it supports min-cpu-platform
git checkout add-min-cpu-platform

sudo python setup.py install
source bash_tab_complete
cd ..

#
# Do last
#

# Cleanup apt cache
sudo apt autoremove -y
