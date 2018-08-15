#!/bin/bash

# This script's purpose is for use with jitpack.io - a repository to publish snapshot automatically
# This script downloads git-lfs and pull needed sources to build GATK in the jitpack environment

GIT_LFS_VERSION="2.5.1"
GIT_LFS_LINK=https://github.com/github/git-lfs/releases/download/v${GIT_LFS_VERSION}/git-lfs-linux-amd64-v${GIT_LFS_VERSION}.tar.gz
GIT_LFS="./git-lfs"
echo "Downloading and untarring git-lfs binary"
wget -qO- $GIT_LFS_LINK | tar xvz git-lfs

echo "Fetching LFS files."
$GIT_LFS install
$GIT_LFS pull --include src/main/resources/large
