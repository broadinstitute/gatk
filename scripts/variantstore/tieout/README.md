# ACMG35 Tieout

There are many challenges in comparing two VCFs -- at a high level, there are two complexities arising from flexibility in the VCF specification around samples and alleles.  There is no defined ordering.  For samples this isn't so bad we just need to compare FORMAT field data based on the appropriate sample (not the column number).  For alleles, since ordering is not defined it means two different GTs (e.g. 0/1 and 0/2) may actually resolve to the exact same genotype because the alleles that they reference may be in a different order.  The code has to handle this as well.

## Comparison Details

### Global comparison

- Skip all `##` header lines
- Verify both VCFs have the same set of samples
- Compare each site

### Site-level comparison**

- assert identical chrom, position, id and ref fields
- SKIP any position with a `*` allele in the ALT of the WARP VCF
- assert set of ALT alleles is the same
- Compare sample-level data

### Sample-level Comparison

- if both have a GQ, assert it is identical
- if both have an RGQ, assert it is identical
- if WARP has a PL, and BQ has an RGQ, assert PL[0] == RGQ
- compare GT alleles (after deconvoluting).  EITHER they are an exact match OR they are different but with equivalent PLs OR warp has no PL but both have GQ == 0

## Data

The tieout was performed using the data in the "ACMG_Cell_Line_Validation_specops_reblock_v2" workspace with 35 reblocked gVCFs.

## Regions

GVCFs are generated using the calling region defined in `gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list`.  The Genomic Variant Store also loads and processes data using that same interval list.

While the WARP pipeline uses those same GVCFs, it processes using a different set of intervals which are defined in `gs://gcp-public-data--broad-references/hg38/v0/hg38.even.handcurated.20k.intervals`.  This list of intervals is largely the same as the GVCF but was trimmed down and subdivided by hand for better performance.  There are two sets of differences:

Sites in WARP but not in GVCF
The WARP intervals include regions where the reference is N, whereas the GVCFs do not attempt to call in these regions.  Thus these extra regions will never contain variants in a WARP callset.

Sites in GVCF but not in WARP:
Initially the WARP intervals completely excluded chrY for performance reasons, and then had selected sites added back in due to user demand.  The GVCFs however call in the majority of the calling regions.  

Therefore we exclude these additional chrY sites from the comparison since there will be valid, valuable data there but WARP is not calling it.

To obtain the list of sites to be excluded:

```bash
gatk IntervalListTools \
--ACTION SUBTRACT \
-I wgs_calling_regions.hg38.interval_list \
-SI hg38.even.handcurated.20k.interval_list \
-O wgs_not_in_handcurated.interval_list

cat wgs_not_in_handcurated.interval_list | grep -v "@" | awk '{ print $1":"$2"-"$3 }' > chrY.excludes.intervals
```

## Running the Comparison

To run the script:

```bash
python compare_data.py <warp-vcf-gz> <gvs-vcf-gz> <file-of-excluded-intervals>
```

So for example,

```bash
python compare_data.py bq_validation_35.0.gnarly.vcf.gz acmg_35_chr1.vcf.gz chrY.excludes.intervals
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

## Gathering WGS Data

Since the data set is small enough we can gather together all the shards into a single gVCF for each pipeline to compare results.  For GVS, this is straighforward because we can just gather the final outputs.  For WARP however, the final output also has hard filtering and VQSR applied to the results, so we need to grab an earlier output specifically from the `TotallyRadicalGatherVcfs` tasks.  The instructions below assume the pipelines were run on the SpecOps Cromwell server which outputs to the `gs://broad-dsp-spec-ops-cromwell-execution` bucket.

### Gather WARP pipeline results

```bash
WORKFLOW_ID=6d57ab66-8668-495f-8947-3b0ee7387389

gsutil ls gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-TotallyRadicalGatherVcfs/shard-*/*.gnarly.vcf.gz > warp_vcfs.list

# sort by shard number
cat warp_vcfs.list  | rev | cut -d"." -f4 | rev | paste - warp_vcfs.list | sort -n | cut -f2 > warp_vcfs.sorted.list

# block gather
~/gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I warp_vcfs.sorted.list \
--output warp.vcf.gz

#index
tabix warp.vcf.gz
```

### Gather GVS pipeline results

NOTE: The CohortExtract should be run with the workflow setting `"NgsCohortExtract.emit_pls":"true"`

```bash
WORKFLOW_ID=76d21253-d0f5-42f0-a6f1-068bcd3e4e4b
gsutil ls gs://broad-dsp-spec-ops-cromwell-execution/NgsCohortExtract/${WORKFLOW_ID}/call-ExtractTask/shard-*/acmg_35_*.vcf.gz > gvs_vcfs.list

# sort by shard number
cat gvs_vcfs.list  | rev | cut -d"." -f3 | cut -d_ -f1 | rev | paste - gvs_vcfs.list | sort -n | cut -f2 > gvs_vcfs.sorted.list

# block gather
~/gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I gvs_vcfs.sorted.list \
--output gvs.vcf.gz

#index
tabix gvs.vcf.gz
```

## Debugging Tips

### Running GnarlyGenotyper against GenomicsDB

Download the GenomicsDB TAR for a shard, and run GnarlyGenotyper

```bash
WORKFLOW_ID=68026724-7fdc-4a0e-bcfa-6a5e4a86cc0a

gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-ImportGVCFs/shard-0/genomicsdb.tar .

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

### Running ExtractCohort against BQ

```bash
reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"
dataset="spec-ops-aou.kc_acmg_tieout_v6"

gatk --java-options "-Xms2g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
  -R $reference \
  -O acmg_35_debug.vcf \
  --local-sort-max-records-in-ram 1000000 \
  --print-debug-information \
  --sample-table ${dataset}.metadata  \
  --project-id spec-ops-aou \
  --cohort-extract-table ${dataset}.exported_cohort_35_test \
  -L chr1:55398671
```
