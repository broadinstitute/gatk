# NOTES
 - excess het sites are PRESENT but FILTERED (ExcessHet) in WARP pipeline
 - excess het sites are NOT MARKED OR FILTERED IN BQ PIPELINE
 - python strips out INFO field annotations with no `=` part (no effect)
 
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
```

First, create a full cohort extract using the `gvs_tieout_acmg_v3` (baseline), or otherwise desired, filtering model.  Then subset down to NA12878/CHM and chr20 as for WARP
dataset="spec-ops-aou.gvs_tieout_acmg_v1"

gatk --java-options "-Xms2g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
  -R $REFERENCE \
  -O full.bq.all.chr20.vcf.gz \
  --local-sort-max-records-in-ram 1000000 \
  --sample-table ${dataset}.metadata  \
  --variant-filter-table ${dataset}.filter_set_info \
  --filter-set-name gvs_tieout_acmg_v3 \
  --project-id spec-ops-aou \
  --cohort-extract-table ${dataset}.exported_cohort_all_samples \
	--interval-set-rule INTERSECTION \
  -L gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list \
  -L chr20
```

## script to add "AS_MAX_VQSLOD" to VCFs
```
python add_max_as_vqslod.py NA12878.bq.chr20.vcf.gz | bgzip > NA12878.bq.chr20.maxas.vcf.gz && tabix NA12878.bq.chr20.maxas.vcf.gz
python add_max_as_vqslod.py NA12878.bq.all.chr20.vcf.gz | bgzip > NA12878.bq.all.chr20.maxas.vcf.gz && tabix NA12878.bq.all.chr20.maxas.vcf.gz
python add_max_as_vqslod.py NA12878.warp.chr20.vcf.gz | bgzip > NA12878.warp.chr20.maxas.vcf.gz && tabix NA12878.warp.chr20.maxas.vcf.gz

python add_max_as_vqslod.py SYNDIP.bq.chr20.vcf.gz | bgzip > SYNDIP.bq.chr20.maxas.vcf.gz && tabix SYNDIP.bq.chr20.maxas.vcf.gz
python add_max_as_vqslod.py SYNDIP.bq.all.chr20.vcf.gz | bgzip > SYNDIP.bq.all.chr20.maxas.vcf.gz && tabix SYNDIP.bq.all.chr20.maxas.vcf.gz
python add_max_as_vqslod.py SYNDIP.warp.chr20.vcf.gz | bgzip > SYNDIP.warp.chr20.maxas.vcf.gz && tabix SYNDIP.warp.chr20.maxas.vcf.gz
```

## Evaluate
```
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.bq.chr20.maxas.vcf.gz -t human_REF_SDF -o NA12878_bq_roc
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.bq.all.chr20.maxas.vcf.gz -t human_REF_SDF -o NA12878_bq_all_roc
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.warp.chr20.maxas.vcf.gz -t human_REF_SDF -o NA12878_warp_roc

rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b chm/CHM-eval.kit/full.38.vcf.gz -e chm/CHM-eval.kit/full.38.bed.gz -c SYNDIP.bq.chr20.maxas.vcf.gz -t human_REF_SDF -o chm_bq_roc
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

## RESEARCH: Create list of excess het sites from WARP
```
rm excess_het_sites.bed
WORKFLOW_ID="c3733f56-218d-4242-af11-6557f101c5e2"
for f in `gsutil ls gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/$WORKFLOW_ID/call-HardFilterAndMakeSitesOnlyVcf/shard-*/cacheCopy/*.sites_only.variant_filtered.vcf.gz `; do
    echo "Processing $f"
    gsutil cat $f | gunzip | awk '{ if ($7 == "ExcessHet") print $1"\t"($2-1)"\t"$2}' >> excess_het_sites.bed
done
cat /Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai | head -25 | cut -f1 > /tmp/genome.txt
bedtools sort -g /tmp/genome.txt -i excess_het_sites.bed > excess_het_sites.sorted.bed

bedtools subtract -header -a ../NA12878.bq.chr20.vcf.gz -b excess_het_sites.sorted.bed > filtered.vcf
```

## DEBUGGING
```
NOTE: WHERE TO GET WARP SITES ONLY INPUT
gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/6d57ab66-8668-495f-8947-3b0ee7387389/call-SitesOnlyGatherVcf/warp_tieout_acmg_cohort_v1.sites_only.vcf.gz

reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

#OUT=NA12878.bq.filtered.vcf.gz
#INTERVALS=gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list

OUT=debug.vcf
INTERVALS=chr20:64150204

gatk --java-options "-Xms2g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
  -R $reference \
  -O ${OUT} \
  --local-sort-max-records-in-ram 1000000 \
  --sample-table ${dataset}.sample_NA12878  \
  --variant-filter-table ${dataset}.filter_set_info \
  --filter-set-name gvs_tieout_acmg_v5_warp_model \
  --project-id spec-ops-aou \
  --cohort-extract-table ${dataset}.exported_cohort_NA12878_test \
	-L ${INTERVALS}
```

