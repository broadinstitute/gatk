# Download 10 random 1000 Genomes CRAMs

This folder includes a script that picks random samples from the official 1000 Genomes low-coverage GRCh38 alignment index and downloads CRAM + CRAI files.

## Prerequisites

- **bash** (the script uses bash-specific features like `[[` and `pipefail`)
- **curl** (for downloading files)
- **Python 3** (for robust cross-platform data parsing)

## Quick start

```bash
cd /Users/hatcher/Projects/1000_genome_crams
chmod +x ./get_1kg_crams.sh
./get_1kg_crams.sh 10 ./downloads download
```

## See sample names first (no download)

```bash
cd /Users/hatcher/Projects/1000_genome_crams
./get_1kg_crams.sh 10 ./downloads list
```

This also writes a manifest at `./downloads/samples.tsv` with sample name + CRAM/CRAI URLs.

## Search for a specific sample

To download CRAM files for a specific sample (e.g., `NA12878`):

```bash
./get_1kg_crams.sh --sample NA12878 ./downloads download
```

To look up a sample without downloading:

```bash
./get_1kg_crams.sh --sample NA12878 ./downloads list
```

If the sample is not found in the 1000 Genomes index, the script will exit with an error message.

## Notes

- Source index: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.low_coverage.GRCh38DH.alignment.index`
- Each CRAM is large (often multiple GB), so ensure enough disk/network before full download.