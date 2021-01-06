# ACMG35 Tieout

The goal is to compare data from the WARP pipeline (from Gnarly from GenomicsDB but before hard-filtering) and compare the fields that will be included in a _cohort export_ (position, alleles, GT, GQ, RGQ).  

## Data

Get a shard of pre-hard filtered extracted data from the WARP run:

```bash
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/52245f96-6b2d-47c3-a191-4119cc31f319/call-TotallyRadicalGatherVcfs/shard-0/bq_validation_35.0.gnarly.vcf.gz .

gunzip bq_validation_35.0.gnarly.vcf.gz
```

Extract the same region from BQ (using the `spec-ops-aou.kc_acmg_tieout` which was loaded with spanning deletions in the VET for tieout purposes).  The interval below is the same interval as shard-0 from gnarly.

```bash
gatk ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT -R ~/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta \
  -O acmg_35_chr1.vcf --local-sort-max-records-in-ram 1000000 --sample-table spec-ops-aou.kc_acmg.metadata  --project-id spec-ops-aou \
    --cohort-extract-table spec-ops-aou.kc_acmg_tieout.exported_cohort_35_test \
      -L chr1:1-35055461
```

## Comparison

### Position Comparison

This shell script compares the two VCFs at the line/position level and output some statistics about that comparison.  Importantly, it also creates two files (`/tmp/only.f1.pos` and `/temp/only.f2.pos`) which contain a list of position that are only found in the first or second files respectively.  This is important so we can ignore those in the site-level comparison

```bash
./find_common_position.sh bq_validation_35.0.gnarly.vcf acmg_35_chr1.vcf
```

### Allele / Genotype Comparison

```bash
python compare_data.py <first-vcf> <second-vcf> <file-of-positions-to-exclude>
```

So for example,

```bash
python compare_data.py bq_validation_35.0.gnarly.vcf acmg_35_chr1.vcf /tmp/only.f1.pos
```



## Debugging Tips

Download the GenomicsDB TAR for this shard

```bash
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/52245f96-6b2d-47c3-a191-4119cc31f319/call-ImportGVCFs/shard-0/genomicsdb.tar .

tar -xf genomicsdb.tar
WORKSPACE=genomicsdb

reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

gatk --java-options "-Xms8g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  GnarlyGenotyper \
  -R $reference \
  -O bq_validation_35.0.0.vcf.gz \
  -V gendb://$WORKSPACE \
  --only-output-calls-starting-in-intervals \
  -stand-call-conf 10 \
  -L chr1:691395-691500
```
