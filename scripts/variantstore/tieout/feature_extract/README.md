# NOTE: the current best practice for comparing feature extract is through a full tieout of filtered reslts

# ACMG35 Feature Extract Tieout

The goal is to compare data from the WARP pipeline going into VQSR with the BigQuery approach

## Data

The 35 sample gVCF used for this analysis are listed in `warp_samples.tsv`

Get a shard of pre-hard filtered extracted data from the WARP run:

```bash
WORFLOW_ID=68026724-7fdc-4a0e-bcfa-6a5e4a86cc0a
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-TotallyRadicalGatherVcfs/shard-0/*.gnarly.vcf.gz bq_validation_35.0.gnarly.vcf.gz
gunzip bq_validation_35.0.gnarly.vcf.gz
```

Extract the same region from BQ (using the `spec-ops-aou.kc_acmg_tieout_v6` which was loaded with spanning deletions in the VET for tieout purposes).  The interval below is the same interval as shard-0 from gnarly.

```bash
gatk --java-options "-Xmx4g" \
    ExtractFeatures \
        --ref-version 38  \
        -R /Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta  \
        -O acmg_feature_extract_debug.vcf.gz \
        --local-sort-max-records-in-ram 1000000 \
        --sample-table spec-ops-aou.kc_acmg_tieout_v6.metadata \
        --alt-allele-table  spec-ops-aou.kc_acmg_tieout_v6.alt_allele \
        --min-location 1000000000000 --max-location 1000035055462 \
        --project-id spec-ops-aou
```

From there WARP pipeline, there are two places we could get our source files from with different caveats:


1. The input to VQSR -- this is sites-only and the final input, so it would be perfect EXCEPT that it has had hard-filtering applied also which includes ExcessHet which is not currently part of the BQ process

2. The inputs to the scatered HardFilterAndMakeSitesOnlyVcf, which is not sites-only (which just makes it bigger) and is scattered... but has not been hard filtered yet.  Interestingly enough... this is turns out to be the same file that we used to compare the cohort extract.

We are starting with option (2) which can be found at the following path:

gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-TotallyRadicalGatherVcfs/shard-0/bq_validation_v5_35.0.gnarly.vcf.gz

## Comparison

NOTE: WARP currently has 1 row per SITE with multiple alleles (and values per allele).  BQ breaks this out to one row per allele.  Question out to Laura if this is equivalent... if it is we just need to fix up the tieout scripts.  If it's not, we need to combine these in the Java code.  Currently we just ignore multi-allelic sites.

First you have to clean each of the VCFs to split up multiallelics and normalize alleles

```bash
./clean_vcf.sh <reference-fasta> <input-vcf> <output-prefix>
```

For example

```bash
./clean_vcf.sh /Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta bq_validation_v5_35.0.gnarly.vcf.gz warp
./clean_vcf.sh /Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta acmg_feature_extract_debug.vcf.gz bq
```

The to run the comparisonscript:

```bash
python compare_feature_data.py <first-vcf> <second-vcf> [<intervals-to-exclude>]
```

So for example,

```bash
python compare_feature_data.py warp.clean.vcf.gz bq.clean.vcf.gz
```

## Debugging

The query that calculates the base data for the metrics is in `src/main/resources/org/broadinstitute/hellbender/tools/variantdb/nextgen/feature_extract.sql`.  Running that against a single position will gather the raw data used to calculate.

### Genomics DB

Download the GenomicsDB TAR for this shard

```bash
WORKFLOW_ID=68026724-7fdc-4a0e-bcfa-6a5e4a86cc0a
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-ImportGVCFs/shard-0/attempt-2/genomicsdb.tar .

tar -xf genomicsdb.tar
WORKSPACE=genomicsdb

reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

gatk --java-options "-Xms8g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  GnarlyGenotyper \
  -R $reference \
  -O debug_warp.vcf \
  -V gendb://$WORKSPACE \
  --only-output-calls-starting-in-intervals \
  -stand-call-conf 10 \
  -L chr1:602222
```




 



