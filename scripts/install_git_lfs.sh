#!/bin/bash

GIT_LFS_VERSION=v1.1.1

echo "cloning"
git clone https://github.com/github/git-lfs.git

echo "cd into git-lfs"
cd git-lfs

echo "checkout commit e8a5cfb85e27c1284a73a2eb8269bd5f5c4d0955"
git checkout $GIT_LFS_VERSION

echo "scripts/bootstrap"
script/bootstrap

echo "ls bin"
ls bin

echo "cd back down"
cd ..

echo "resetting travis remote"
git remote set-url origin "git@github.com:broadinstitute/gatk.git"

echo "install"
git-lfs/bin/git-lfs install

echo "pull"
git-lfs/bin/git-lfs pull

echo "ls-files"
git-lfs/bin/git-lfs ls-files
ls -lah src/test/resources/large/
md5sum src/test/resources/large/*.*
