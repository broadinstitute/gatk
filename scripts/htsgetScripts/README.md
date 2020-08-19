This folder contains scripts and files necessary for starting a local htsget reference server so that htsget functionality can be tested.

`start-htsget-test-server.sh` starts a docker image running the reference server, serving files from the directory `src/test/resources/htsjdk/samtools/BAMFileIndexTest/`.

`htsget_config.json` configures the mapping between path name patterns and resource sources used by the reference server.