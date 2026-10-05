# Build the GATK Docker image on a cloud VM

The instructions here are written specifically for building an `ah_var_store` version of the GATK Docker image from an
Azure virtual machine, though much of this would likely apply to Google Cloud (or other clouds) and other Docker images.
See "Building on Google Cloud instead" at the end for a GCP recipe, and read the two "Gotcha" sections before you run
the push block -- on a current Docker they are the difference between a working build and half a working build.


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

# Run a helper script with lots more commands for building a GATK Docker image.
# This also tags each image as soon as it is built -- see the first "Gotcha" section below for why that matters.
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

## Gotcha: the two builds share a tag, and on newer Docker the lite image does not survive

`gatk_docker_setup.sh` now tags each image immediately after building it, so this trap should not bite you if you
followed the steps above. The rest of this section explains what it is guarding against, since the same failure
turns up whenever the two builds are run by hand.

Both builds pass the same `-e ${GITHUB_HASH}`, and `build_docker.sh` turns that into
`docker build -t ${REPO_PRJ}:${GITHUB_TAG}` for *both* of them. The heavy build runs second and takes that tag over,
which leaves the lite image with no reference pointing at it.

Whether that matters depends on the Docker version of the machine you are building on:

- **Classic Docker (graph driver).** The now-untagged lite image lingers as a dangling image, `docker tag
  "${IMAGE_ID_LITE}" ...` in the push block still resolves it by ID, and everything above works as written. This is
  what the Azure VMs have historically given us, which is why this procedure has been fine for years.
- **Docker 29+ with the containerd image store**, which is what a current Ubuntu 22.04 image gives you by default on
  GCP. The unreferenced lite image is not retained, `docker inspect ${IMAGE_ID_LITE}` reports `no such object`, and
  the push block fails at the first `docker tag`. You get to the end of a successful ~30 minute build with only the
  heavy image to show for it.

Check which store you are on with `docker info --format '{{.DriverStatus}}'`; `io.containerd.snapshotter.v1` means you
are exposed to this.

The fix, which `gatk_docker_setup.sh` now applies, is to give each image its unique GAR tag *immediately after it is
built* rather than relying on both surviving until the end. The tag/push block above re-derives the same tags and is
harmlessly idempotent, so it still works whether or not the images were already tagged.

If you are recovering a run that lost the lite image -- an older copy of the setup script, or builds run by hand --
tag the heavy image first. Otherwise rebuilding the lite image takes `${REPO_PRJ}:${GITHUB_TAG}` back and puts the
heavy image in exactly the same position:

```
# 1. Protect what you already have before rebuilding anything
IMAGE_ID_HEAVY=$(cut -c8-19 /tmp/idfile_heavy.txt)
TAG_HEAVY=$(python3 ./scripts/variantstore/scripts/build_docker_tag.py --image-id "${IMAGE_ID_HEAVY}" --image-type gatkbase)
docker tag "${IMAGE_ID_HEAVY}" "us-central1-docker.pkg.dev/broad-dsde-methods/gvs/gatk:${TAG_HEAVY}"

# 2. Rebuild only the lite image, with the same flags gatk_docker_setup.sh uses
export GITHUB_HASH=$(git rev-parse HEAD)
bash build_docker.sh -m -u -e ${GITHUB_HASH} -s -d /mnt/staging-tmp
cp /tmp/idfile.txt /tmp/idfile_lite.txt

# 3. Tag it straight away, then push both
IMAGE_ID_LITE=$(cut -c8-19 /tmp/idfile_lite.txt)
TAG_LITE=$(python3 ./scripts/variantstore/scripts/build_docker_tag.py --image-id "${IMAGE_ID_LITE}" --image-type gatkbase-lite)
docker tag "${IMAGE_ID_LITE}" "us-central1-docker.pkg.dev/broad-dsde-methods/gvs/gatk:${TAG_LITE}"
```


## Gotcha: what the hex suffix in a tag means depends on who built it

`build_docker_tag.py` embeds whatever Docker calls the image ID, and the two image stores mean different things by
that. Under the classic graph driver it is the image's *config* digest; under the containerd store it is the
*manifest* digest. So `2026-08-19-gatkbase-lite-087565f1a432` (Azure) and `2026-09-29-gatkbase-lite-64383fd3ade2`
(GCP) are not the same kind of hash, and neither is derivable from the other. Do not try to match them up.

Relatedly, an image built under the containerd store pushes as an **OCI image index** rather than a single Docker v2
manifest, and carries an extra BuildKit attestation manifest that shows up as an `unknown/unknown` platform entry.
Modern Docker and Cromwell pull both forms without complaint, but tooling that fetches manifests directly needs to
send an `Accept` header that includes `application/vnd.oci.image.index.v1+json` or it will get `MANIFEST_UNKNOWN`.

To verify an image is what you think it is, compare layer count and compressed size against a known-good tag rather
than reading anything into the name. At the time of writing a correct pair looks like: lite = 12 layers / 1.29 GB
compressed, heavy = 13 layers / 2.49 GB compressed, the difference being a single ~1.2 GB Conda/ML layer.

## Building on Google Cloud instead

The procedure above works on GCP with one adaptation to VM allocation and the two gotchas above. An
`n2-standard-8` finishes both images in about half an hour, comfortably faster than the 2-vCPU Azure box.

```
gcloud compute instances create gatk-docker-build-<ticket> \
  --project gvs-internal --zone us-central1-a \
  --machine-type n2-standard-8 \
  --image-family ubuntu-2204-lts --image-project ubuntu-os-cloud \
  --boot-disk-size 300GB --boot-disk-type pd-balanced \
  --network gvs-network \
  --tags gatk-docker-build \
  --scopes cloud-platform --metadata enable-oslogin=TRUE
```

Notes specific to GCP:

- **GCP has no Azure-style ephemeral `/mnt` data disk.** Size the boot disk instead (300 GB is ample -- a full run
  uses about 36 GB) and the `cd /mnt` steps above work unchanged.
- **The VM needs outbound internet, and `gvs-internal` has no Cloud NAT.** If you create it with `--no-address` the
  build fails in a thoroughly confusing way: Private Google Access keeps `gcloud` and some apt mirrors working, so
  the machine looks healthy right up until the first non-Google download, where `apt-get` cannot find `git-lfs` and
  the `git clone` hangs. Either give the VM an external IP (as above, by omitting `--no-address`) or stand up Cloud
  NAT. Sanity-check with `curl -so /dev/null -w '%{http_code}' https://download.docker.com` before starting a build.
- **SSH from outside the Broad network** needs IAP, since `gvs-network`'s SSH rule only admits Broad ranges. Add a
  tag-scoped rule and delete it with the VM:
  ```
  gcloud compute firewall-rules create gatk-docker-build-iap-ssh --project gvs-internal \
    --network gvs-network --direction INGRESS --action allow --rules tcp:22 \
    --source-ranges 35.235.240.0/20 --target-tags gatk-docker-build
  gcloud compute ssh gatk-docker-build-<ticket> --project gvs-internal --zone us-central1-a --tunnel-through-iap
  ```
- **The VM's service account cannot push to `broad-dsde-methods` GAR.** Rather than `gcloud init` on the VM, forward
  a short-lived token from a machine that already has upload rights. It expires in about an hour and leaves no
  credential behind:
  ```
  gcloud auth print-access-token | gcloud compute ssh <vm> --project gvs-internal --zone us-central1-a \
    --tunnel-through-iap --command 'read -r T; echo "$T" | sudo docker login -u oauth2accesstoken \
    --password-stdin https://us-central1-docker.pkg.dev'
  ```
  Do **not** run `gcloud auth configure-docker` if you do this -- it installs a credential helper that authenticates
  as the VM's service account and overrides the token login.
- Long builds outlive SSH sessions, and IAP tunnels drop. Launch under `setsid nohup ... </dev/null` writing to a log
  file, and poll the log rather than holding a session open.

Don't forget to shut down (and possibly delete) your VM once you're done! Delete the temporary firewall rule too.
