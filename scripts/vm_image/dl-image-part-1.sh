#!/usr/bin/env bash

# dl-image-part-{1,2} setup scripts are separated from the default because it seems that there
# may be issues with GPU-enabled docker machines if you run on a non-GPU
# instance. See https://github.com/NVIDIA/nvidia-docker/issues/652.
# These scripts should not be run on a blank Ubuntu 18.04 VM, but instead on a ml4h-image.

# TODO Peg the NVIDIA drivers to specific versions that are compatible with the versions
# of related software defined within the Docker image

# 2018/09/11 additions
# Enable NVidia-docker
# Via https://askubuntu.com/a/1036265/411855
sudo add-apt-repository -y ppa:graphics-drivers/ppa
sudo apt update
sudo apt-get install -y ubuntu-drivers-common
sudo ubuntu-drivers autoinstall

# Reboot. Then return and move on to part 2
sudo reboot
