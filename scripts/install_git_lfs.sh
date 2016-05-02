#!/bin/bash

GIT_LFS_VERSION="1.1.2"
GIT_LFS_LINK=https://github.com/github/git-lfs/releases/download/v${GIT_LFS_VERSION}/git-lfs-linux-amd64-${GIT_LFS_VERSION}.tar.gz
GIT_LFS="git-lfs-${GIT_LFS_VERSION}/git-lfs"
echo "downloading and untarring git-lfs binary" 
wget -qO- $GIT_LFS_LINK | tar xvz

echo "ls"
ls

echo "resetting travis remote"
git remote set-url origin "https://github.com/broadinstitute/gatk.git"

echo "git lfs install"
GIT_TRACE=1 $GIT_LFS install

echo "fetch"
GIT_TRACE=1 $GIT_LFS fetch

echo "checkout"
GIT_TRACE=1 $GIT_LFS checkout

echo "ls-files"
GIT_TRACE=1 $GIT_LFS ls-files
ls -lah src/test/resources/large/
md5sum src/test/resources/large/*.*

echo "logs"
GIT_TRACE=1 $GIT_LFS logs last
