# NOTES
 - excess het sites are PRESENT but FILTERED (ExcessHet) in WARP pipeline
 - excess het sites are NOT MARKED OR FILTERED IN BQ PIPELINE
 - python strips out INFO field annotations with no `=` part (no effect)
 - Previously there was an attempt to tie out the feature extract input to VQSR between WARP and GVS.  However, that was abandoned in favor of the truth-sample based comparison here due to non-determinism and bugs in GenomicsDB that were hard to compensate for.  However, the docs and scripts for doing this are in the branch `gvs_vqsr_input_tieout` under `scripts/variantstore/tieout/feature_extract`
 
## One time tasks
```
REFERENCE="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

gsutil -m cp gs://dsde-palantir/NISTTruthData/* .
gatk IntervalListToBed -I NIST_highconfidenceregions.interval_list -O NIST_highconfidenceregions.bed
rtg format --output human_REF_SDF $REFERENCE
```
## WARP -- take full results, extract two samples, and strip down to chr20
```
WORKFLOW_ID=36e7547e-3253-4888-9b1c-2a7437401aee
gsutil -m cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-FinalGatherVcf/warp_tieout_acmg_cohort_v1.vcf.gz .

~/gatk SelectVariants -V warp_tieout_acmg_cohort_v1.vcf.gz --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.warp.vcf.gz
bcftools view -O z NA12878.warp.vcf.gz chr20 > NA12878.warp.chr20.vcf.gz
tabix NA12878.warp.chr20.vcf.gz

~/gatk SelectVariants -V warp_tieout_acmg_cohort_v1.vcf.gz --sample-name CHMI_CHMI3_WGS1 --select-type-to-exclude NO_VARIATION -O SYNDIP.warp.vcf.gz
bcftools view -O z SYNDIP.warp.vcf.gz chr20 > SYNDIP.warp.chr20.vcf.gz
tabix SYNDIP.warp.chr20.vcf.gz

gsutil -m cp NA12878* gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/
gsutil -m cp SYNDIP* gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/

```

## BIGQUERY 

First, create a full cohort extract (as described in README.md) using the `gvs_tieout_acmg_v3` (baseline), or otherwise desired, filtering model.  

```
~/gatk SelectVariants -V gvs.bq.all.vcf.gz --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.bq.all.vcf.gz
bcftools view -O z NA12878.bq.all.vcf.gz chr20 > NA12878.bq.all.chr20.vcf.gz
tabix NA12878.bq.all.chr20.vcf.gz

~/gatk SelectVariants -V gvs.bq.all.vcf.gz --sample-name CHMI_CHMI3_WGS1 --select-type-to-exclude NO_VARIATION -O SYNDIP.bq.all.vcf.gz
bcftools view -O z SYNDIP.bq.all.vcf.gz chr20 > SYNDIP.bq.all.chr20.vcf.gz
tabix SYNDIP.bq.all.chr20.vcf.gz

gsutil -m cp gvs.bq.all.vcf.gz* gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/
gsutil -m cp NA12878* gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/
gsutil -m cp SYNDIP* gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/
```



## script to add "AS_MAX_VQSLOD" to VCFs
```
python add_max_as_vqslod.py NA12878.bq.all.chr20.vcf.gz | bgzip > NA12878.bq.all.chr20.maxas.vcf.gz && tabix NA12878.bq.all.chr20.maxas.vcf.gz
python add_max_as_vqslod.py NA12878.warp.chr20.vcf.gz | bgzip > NA12878.warp.chr20.maxas.vcf.gz && tabix NA12878.warp.chr20.maxas.vcf.gz

python add_max_as_vqslod.py SYNDIP.bq.all.chr20.vcf.gz | bgzip > SYNDIP.bq.all.chr20.maxas.vcf.gz && tabix SYNDIP.bq.all.chr20.maxas.vcf.gz
python add_max_as_vqslod.py SYNDIP.warp.chr20.vcf.gz | bgzip > SYNDIP.warp.chr20.maxas.vcf.gz && tabix SYNDIP.warp.chr20.maxas.vcf.gz

```

## Evaluate
```
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.bq.all.chr20.maxas.vcf.gz -t human_REF_SDF -o NA12878_bq_all_roc
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.warp.chr20.maxas.vcf.gz -t human_REF_SDF -o NA12878_warp_roc

rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b chm/CHM-eval.kit/full.38.vcf.gz -e chm/CHM-eval.kit/full.38.bed.gz -c SYNDIP.bq.all.chr20.maxas.vcf.gz -t human_REF_SDF -o chm_bq_all_roc
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b chm/CHM-eval.kit/full.38.vcf.gz -e chm/CHM-eval.kit/full.38.bed.gz -c SYNDIP.warp.chr20.maxas.vcf.gz -t human_REF_SDF -o chm_warp_roc
```
## View
```
rtg roc NA12878_*_roc/snp_roc.tsv.gz 
rtg roc NA12878_*_roc/non_snp_roc.tsv.gz 

rtg roc chm_*_roc/snp_roc.tsv.gz 
rtg roc chm_*_roc/non_snp_roc.tsv.gz 
```

## RESEARCH: How to create list of excess het sites from WARP
```
rm excess_het_sites.bed
WORKFLOW_ID="36e7547e-3253-4888-9b1c-2a7437401aee" for f in `gsutil ls gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-HardFilterAndMakeSitesOnlyVcf/shard-*/*.sites_only.variant_filtered.vcf.gz `; do 
    echo "Processing $f" 
    gsutil cat $f | gunzip | awk '{ if ($7 == "ExcessHet") print $1"\t"($2-1)"\t"$2}' >> excess_het_sites.bed
done

cat /Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai | head -25 | cut -f1 > /tmp/genome.txt
bedtools sort -g /tmp/genome.txt -i excess_het_sites.bed > excess_het_sites.sorted.bed

bedtools subtract -header -a ../NA12878.bq.chr20.vcf.gz -b excess_het_sites.sorted.bed > filtered.vcf
```
