# Install Docker
# https://docs.docker.com/engine/install/ubuntu/
apt-get install --assume-yes ca-certificates curl gnupg lsb-release

mkdir -m 0755 -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
apt-get update
apt-get install --assume-yes docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Make sure this all worked with a Docker 'Hello World!' test.
docker run hello-world
# required cleanup, boot disk space is very tight by default
docker system prune --force

# Temurin Java 17
# https://adoptium.net/installation/linux/
apt install -y wget apt-transport-https
mkdir -p /etc/apt/keyrings
wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc
echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
apt-get -qq update
apt-get -qq install temurin-17-jdk

export GITHUB_HASH=$(git rev-parse HEAD)
export STAGING_DIR=/mnt/staging-tmp

# Both builds below pass the same -e ${GITHUB_HASH}, so build_docker.sh runs `docker build -t
# broadinstitute/gatk:${GITHUB_HASH}` for each of them and the second one takes that tag over. That leaves the
# first image with no reference pointing at it. Under the classic Docker graph driver it survives as a dangling
# image and can still be addressed by ID at push time, but under the containerd image store (the default from
# Docker 29, which is what a current Ubuntu 22.04 VM gets) it is not retained and its ID stops resolving -- so
# the push step fails after a successful build. Giving each image its own tag as soon as it is built keeps both
# alive regardless of image store. See the "Gotcha" sections in
# scripts/variantstore/docs/Build Docker from VM/building_gatk_docker_on_a_vm.md
GAR_GATK_REPO="us-central1-docker.pkg.dev/broad-dsde-methods/gvs/gatk"

# Build the lite image (no Conda/ML stack, used for most GVS tasks)
bash build_docker.sh -m -u -e ${GITHUB_HASH} -s -d ${STAGING_DIR}
cp /tmp/idfile.txt /tmp/idfile_lite.txt
IMAGE_ID_LITE=$(cut -c8-19 /tmp/idfile_lite.txt)
TAG_LITE=$(python3 ./scripts/variantstore/scripts/build_docker_tag.py --image-id "${IMAGE_ID_LITE}" --image-type "gatkbase-lite")
docker tag "${IMAGE_ID_LITE}" "${GAR_GATK_REPO}:${TAG_LITE}"
echo "Lite image ${IMAGE_ID_LITE} tagged as ${GAR_GATK_REPO}:${TAG_LITE}"

# Build the heavy image (full gatkbase with Conda/ML stack, used for VETS/VQSR)
bash build_docker.sh -u -e ${GITHUB_HASH} -s -d ${STAGING_DIR}
cp /tmp/idfile.txt /tmp/idfile_heavy.txt
IMAGE_ID_HEAVY=$(cut -c8-19 /tmp/idfile_heavy.txt)
TAG_HEAVY=$(python3 ./scripts/variantstore/scripts/build_docker_tag.py --image-id "${IMAGE_ID_HEAVY}" --image-type "gatkbase")
docker tag "${IMAGE_ID_HEAVY}" "${GAR_GATK_REPO}:${TAG_HEAVY}"
echo "Heavy image ${IMAGE_ID_HEAVY} tagged as ${GAR_GATK_REPO}:${TAG_HEAVY}"

# Install gcloud
# https://cloud.google.com/sdk/docs/install#deb
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
apt-get update && apt-get install google-cloud-cli
