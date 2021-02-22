# ACMG35 Tieout

The goal is to compare data from the WARP pipeline (from Gnarly from GenomicsDB but before hard-filtering) and compare the fields that will be included in a _cohort export_ (position, alleles, GT, GQ, RGQ).  

## Data

The 35 sample gVCF used for this analysis are listed in `warp_samples.tsv`

Get a shard of pre-hard filtered extracted data from the WARP run:

```bash
WORFLOW_ID=e2644477-e3d3-4a7a-afcb-12e31ce1f85c
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORFLOW_ID}/call-TotallyRadicalGatherVcfs/shard-0/*.gnarly.vcf.gz bq_validation_35.0.gnarly.vcf.gz
gunzip bq_validation_35.0.gnarly.vcf.gz
```

Extract the same region from BQ (using the `spec-ops-aou.kc_acmg_tieout_v2` which was loaded with spanning deletions in the VET for tieout purposes).  The interval below is the same interval as shard-0 from gnarly.

```bash
gatk ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT -R ~/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta \
  -O acmg_35_chr1.vcf --local-sort-max-records-in-ram 1000000 --sample-table spec-ops-aou.kc_acmg_tieout_v2.metadata  --project-id spec-ops-aou \
    --cohort-extract-table spec-ops-aou.kc_acmg_tieout_v2.exported_cohort_35_test_pl \
      -L chr1:1-35055461
```

## Comparison

### Allele / Genotype Comparison

The code in this python script is iterative, and a bit janky.  As we find cases we understand we can exclude them from the comparison with a TODO note about what the difference is.

For example, we are currently ignoring:

* spaning deletions (possibly incorrectly)
* phasing (e.g. 0/1 vs 0|1)


To run the script:

```bash
python compare_data.py <first-vcf> <second-vcf>
```

So for example,

```bash
python compare_data.py bq_validation_35.0.gnarly.vcf.gz acmg_35_chr1.vcf.gz
```

The output has grown organically to help with debugging discrepancies.  It will look something like this for a difference:

```text
DIFF on Genotypes for SM-GXZUI at chr1:891951 with C and C
['.', '.'] vs ['T', 'T']
{'GT': './.'}
{'GT': '0/0', 'GQ': '50'}
--------------
DIFF on Genotypes for SM-GXZUP at chr1:960406 with A and A
['.', '.'] vs ['G', 'G']
{'GT': './.'}
{'GT': '0/0', 'GQ': '50'}
--------------
```

Where each discrepancy is separated by `------------`

## Notes

* the order of samples in the VCFs is different... why would this be?  Shouldn't the be lexically sorted?
* the order of alleles in the VCFs is different... same question as above?

## Debugging Tips

Download the GenomicsDB TAR for this shard

```bash
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/feeca595-85a8-4d2b-b1ef-7f6dc64d714c/call-ImportGVCFs/shard-0/genomicsdb.tar .

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
  -L chr1:1020061
```

```
reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"
dataset="spec-ops-aou.kc_acmg_tieout_v5"

GATK_LOCAL_JAR=/Users/kcibul/projects/gatk/build/libs/gatk-package-4.1.9.0-191-g61b5dfb-SNAPSHOT-local.jar

gatk --java-options "-Xms8g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
  -R $reference \
  -O acmg_35_debug.vcf \
  --local-sort-max-records-in-ram 1000000 \
	--print-debug-information \
  --sample-table ${dataset}.metadata  \
  --project-id spec-ops-aou \
  --cohort-extract-table ${dataset}.exported_cohort_35_test \
  -L chr1:10336
```


## WGS Tieout (Cohort Extract)

`#WORFLOW_ID=e2644477-e3d3-4a7a-afcb-12e31ce1f85c`
`#WORFLOW_ID=b9b42238-6b4a-420f-ae82-586742af9f7d`
WORKFLOW_ID=68026724-7fdc-4a0e-bcfa-6a5e4a86cc0a

gsutil ls gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-TotallyRadicalGatherVcfs/shard-*/*.gnarly.vcf.gz > warp_vcfs.list

### crazy command to sort these files by 
cat warp_vcfs.list  | rev | cut -d"." -f4 | rev | paste - warp_vcfs.list | sort -n | cut -f2 > warp_vcfs.sorted.list

~/gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I warp_vcfs.sorted.list \
--output warp.vcf.gz

tabix warp.vcf.gz

### GVS
`#WORKFLOW_ID=a78e6345-d53a-40ad-ab40-80c509bb5a8d`
`#WORKFLOW_ID=ab6e11c6-9efc-4893-918f-41f537cf3803`
WORKFLOW_ID=fc109a53-bf34-4cd3-bdae-f0b0e6c7c112
gsutil ls gs://broad-dsp-spec-ops-cromwell-execution/NgsCohortExtract/${WORKFLOW_ID}/call-ExtractTask/shard-*/acmg_35_*.vcf.gz > gvs_vcfs.list

cat gvs_vcfs.list  | rev | cut -d"." -f3 | cut -d_ -f1 | rev | paste - gvs_vcfs.list | sort -n | cut -f2 > gvs_vcfs.sorted.list

~/gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I gvs_vcfs.sorted.list \
--output gvs.vcf.gz

tabix gvs.vcf.gz




V4
Compared 15 156 671 positions
    296 Genotypes
    296 GQ
  13900 position
	
v5
Compared 15 156 671 positions	
  21140 Genotypes
  21140 GQ
  13900 position
