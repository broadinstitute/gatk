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
WORKFLOW_ID=49f5c0a5-dbfb-4b05-b293-11ddadeae8e5
gsutil -m cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-FinalGatherVcf/warp_tieout_acmg_cohort_v1.vcf.gz* .
bcftools view -O z warp_tieout_acmg_cohort_v1.vcf.gz chr20 > warp_tieout_acmg_cohort_v1.chr20.vcf.gz
tabix warp_tieout_acmg_cohort_v1.chr20.vcf.gz

INPUT_VCF=../warp_tieout_acmg_cohort_v1.chr20.vcf.gz
SOURCE=warp
gatk SelectVariants -V ${INPUT_VCF} --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name CHMI_CHMI3_WGS1 --select-type-to-exclude NO_VARIATION -O SYNDIP.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202194294 --select-type-to-exclude NO_VARIATION -O BI_HG002.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202243693 --select-type-to-exclude NO_VARIATION -O BI_HG003.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 573673 --select-type-to-exclude NO_VARIATION -O UW_HG002.${SOURCE}.chr20.vcf.gz

gsutil -m cp *chr20.vcf.gz* gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/warp/
```

## BIGQUERY 

First, create a full cohort extract (as described in README.md) using the `gvs_tieout_acmg_v3` (baseline), or otherwise desired, filtering model.  

```
~/gatk SelectVariants -V gvs.bq.all.vcf.gz --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.bq.all.vcf.gz
bcftools view -O z NA12878.bq.all.vcf.gz chr20 > NA12878.bq.all.chr20.vcf.gz
tabix NA12878.bq.all.chr20.vcf.gz

INPUT_VCF=../gvs_tieout_acmg_v2.bq.all.vcf.gz
SOURCE=bq
gatk SelectVariants -V ${INPUT_VCF} --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name CHMI_CHMI3_WGS1 --select-type-to-exclude NO_VARIATION -O SYNDIP.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202194294 --select-type-to-exclude NO_VARIATION -O BI_HG002.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202243693 --select-type-to-exclude NO_VARIATION -O BI_HG003.${SOURCE}.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 573673 --select-type-to-exclude NO_VARIATION -O UW_HG002.${SOURCE}.chr20.vcf.gz

#bcftools view -e 'GT="ref"' -a -O z -s SM-G947Y ${INPUT_VCF} > NA12878.${SOURCE}.chr20.vcf.gz && tabix NA12878.${SOURCE}.chr20.vcf.gz
#bcftools view -e 'GT="ref"' -a -O z -s CHMI_CHMI3_WGS1 ${INPUT_VCF} > SYNDIP.${SOURCE}.chr20.vcf.gz && tabix SYNDIP.${SOURCE}.chr20.vcf.gz
#bcftools view -e 'GT="ref"' -a -O z -s 1202194294 ${INPUT_VCF} > BI_HG002.${SOURCE}.chr20.vcf.gz && tabix BI_HG002.${SOURCE}.chr20.vcf.gz
#bcftools view -e 'GT="ref"' -a -O z -s 1202243693 ${INPUT_VCF} > BI_HG003.${SOURCE}.chr20.vcf.gz && tabix BI_HG003.${SOURCE}.chr20.vcf.gz
#bcftools view -e 'GT="ref"' -a -O z -s 573673  ${INPUT_VCF} > UW_HG002.${SOURCE}.chr20.vcf.gz && tabix UW_HG002.${SOURCE}.chr20.vcf.gz

```

## script to add "AS_MAX_VQSLOD" to VCFs
```
for source in bq warp
do
  for sample in NA12878 SYNDIP BI_HG002 BI_HG003 UW_HG002
  do
    echo "Processing ${sample} for ${source}"
    python ../add_max_as_vqslod.py ${sample}.${source}.chr20.vcf.gz | bgzip > ${sample}.${source}.chr20.maxas.vcf.gz && tabix ${sample}.${source}.chr20.maxas.vcf.gz
  done
done

```

## Evaluate
```

for source in bq warp
do
    rtg vcfeval --roc-subset snp,indel --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -e ../truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
    -c NA12878.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o NA12878_${source}_roc

    rtg vcfeval --roc-subset snp,indel --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c BI_HG002.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o BI_HG002_${source}_roc

    rtg vcfeval --roc-subset snp,indel --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c UW_HG002.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o UW_HG002_${source}_roc

    rtg vcfeval --roc-subset snp,indel --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c BI_HG003.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o BI_HG003_${source}_roc

    rtg vcfeval --roc-subset snp,indel --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/CHM.full.38.vcf.gz -e ../truth/CHM.full.38.bed.gz -c SYNDIP.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o syndip_${source}_roc
done


for source in bq warp
do
    rtg vcfeval --roc-subset snp,indel --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -e ../truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
    -c NA12878.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o NA12878_${source}_roc_filtered

    rtg vcfeval --roc-subset snp,indel --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c BI_HG002.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o BI_HG002_${source}_roc_filtered
    rtg vcfeval --roc-subset snp,indel --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c UW_HG002.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o UW_HG002_${source}_roc_filtered
    rtg vcfeval --roc-subset snp,indel --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c BI_HG003.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o BI_HG003_${source}_roc_filtered
    rtg vcfeval --roc-subset snp,indel --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/CHM.full.38.vcf.gz -e ../truth/CHM.full.38.bed.gz -c SYNDIP.${source}.chr20.maxas.vcf.gz -t ../human_REF_SDF -o syndip_${source}_roc_filtered
done

ONLY-EH
sample=BI_HG002
source=bq
python ../add_max_as_vqslod.py ${sample}.${source}.chr20.vcf.gz loose | bgzip > ${sample}.${source}.chr20.maxas.onlyeh.vcf.gz && tabix ${sample}.${source}.chr20.maxas.onlyeh.vcf.gz
rtg vcfeval --roc-subset snp,indel --all-records --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c BI_HG002.${source}.chr20.maxas.onlyeh.vcf.gz -t ../human_REF_SDF -o BI_HG002_${source}_roc_onlyeh
rtg vcfeval --roc-subset snp,indel               --vcf-score-field=INFO.MAX_AS_VQSLOD --region chr20 -b ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e ../truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c BI_HG002.${source}.chr20.maxas.onlyeh.vcf.gz -t ../human_REF_SDF -o BI_HG002_${source}_roc_onlyeh_filt


for f in $(ls -1 *_roc_onlyeh/indel*.tsv.gz); do
    d=$(cat $f | gunzip | tail -1 | cut -f3,5,6,7)
    s=$(echo $f | cut -d"/" -f1 | sed s/_roc_filtered//g)
    echo -e "$s\t$d"
done


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
