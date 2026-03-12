The Dockerfile in this directory is used to build a Docker image for Ensembl VEP 115 with GRCh38 LOFTEE support.
This build process should be run on an x86 machine such as a GCE or Azure VM. Attempts to build this on an Apple Silicon
Mac resulted in frequent random segmentation faults and a final image that was mysteriously missing required components.

Steps to build:

- Clone the Ensembl VEP repo [here](https://github.com/Ensembl/ensembl-vep) and make sure the appropriate branch is checked out.
For Ensembl VEP 115 this branch is `releases/115`, which at the time of this writing also happens to be the default branch.

- Copy the Dockerfile in this directory over the `docker/Dockerfile` in the Ensembl VEP repo.
  It is very possible that adjustments will need to be made to this Dockerfile for Ensembl VEP releases after 115. `diff`ing this Dockerfile
  with the one in the `docker/Dockerfile` for Ensembl VEP 115 should illustrate what needs to be added for LOFTEE support.

- Copy the script `build_vep_loftee_docker.sh` in this directory to the root of the local Ensembl VEP repo.

- Run the script from the root of the local Ensembl VEP repo. No arguments are required.
  The script should build and push the image and print out its tag to standard output.

- Update the `vep_loftee_docker` output of `GetToolVersions` in `GvsUtils.wdl` to reference the updated tag of this image.
