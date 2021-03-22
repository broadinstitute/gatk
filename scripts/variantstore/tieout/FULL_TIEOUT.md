### NOTES
 - excess het sites are PRESENT but FILTERED (ExcessHet) in WARP pipeline
 - excess het sites are NOT MARKED OR FILTERED IN BQ PIPELINE
 - python strips out INFO field annotations with no `=` part (no effect)
 
# One time tasks
REFERENCE="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

gsutil -m cp gs://dsde-palantir/NISTTruthData/* .
gatk IntervalListToBed -I NIST_highconfidenceregions.interval_list -O NIST_highconfidenceregions.bed
rtg format --output human_REF_SDF $REFERENCE

# WARP -- take full results, and strip down to chr20
gsutil cp gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/NA12878.warp.vcf.gz .
bcftools view -O z NA12878.warp.vcf.gz chr20 > NA12878.warp.chr20.vcf.gz
tabix NA12878.warp.chr20.vcf.gz

# BIGQUERY
dataset="spec-ops-aou.gvs_tieout_acmg_v1"

gatk --java-options "-Xms2g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
  -R $REFERENCE \
  -O NA12878.bq.chr20.vcf.gz \
  --local-sort-max-records-in-ram 1000000 \
  --sample-table ${dataset}.sample_NA12878  \
  --variant-filter-table ${dataset}.filter_set_info \
  --filter-set-name gvs_tieout_acmg_v1 \
  --project-id spec-ops-aou \
  --cohort-extract-table ${dataset}.exported_cohort_NA12878_test \
	--interval-set-rule INTERSECTION \
  -L gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list \
  -L chr20


# script to add "AS_MAX_VQSLOD" to VCFs
python add_max_as_vqslod.py NA12878.bq.chr20.vcf.gz | bgzip > NA12878.bq.chr20.maxas.vcf.gz && tabix NA12878.bq.chr20.maxas.vcf.gz
python add_max_as_vqslod.py NA12878.warp.chr20.vcf.gz | bgzip > NA12878.warp.chr20.maxas.vcf.gz && tabix NA12878.warp.chr20.maxas.vcf.gz

# Evaluate
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.bq.chr20.maxas.vcf.gz -t human_REF_SDF -o bq_roc
rtg vcfeval --region chr20 --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD -b NIST_truth.vcf.gz -e NIST_highconfidenceregions.bed -c NA12878.warp.chr20.maxas.vcf.gz -t human_REF_SDF -o warp_roc

# View
rtg roc bq_roc/*_roc.tsv.gz warp_roc/*_roc.tsv.gz

# DEBUGGING
reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

OUT=NA12878.bq.filtered.vcf.gz
INTERVALS=gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list

OUT=debug.vcf
INTERVALS=chr20:64150204

gatk --java-options "-Xms2g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  ExtractCohort --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
  -R $reference \
  -O ${OUT} \
  --local-sort-max-records-in-ram 1000000 \
  --sample-table ${dataset}.sample_NA12878  \
  --variant-filter-table ${dataset}.filter_set_info \
  --filter-set-name gvs_tieout_acmg_v1 \
  --project-id spec-ops-aou \
  --cohort-extract-table ${dataset}.exported_cohort_NA12878_test \
	-L ${INTERVALS}

