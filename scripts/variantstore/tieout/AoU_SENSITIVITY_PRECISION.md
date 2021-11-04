## Prerequisites

1. Use conda to create a fresh environment to add these new tools to it:
 ```
 conda create --name gvs python=3.8
 conda activate gvs
 conda install -c bioconda samtools=1.9 --force-reinstall
 conda install -c bioconda bcftools
 conda install -c bioconda rtg-tools
```

2. Download reference genome and obtain truth sample VCFs
```
# create SDF for reference (used by rtg eval)
# use local path to GRC38 reference fasta
REFERENCE="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"
rtg format --output human_REF_SDF $REFERENCE

# copy down truth data locally
mkdir -p truth && gsutil -m cp gs://broad-dsp-spec-ops/gvs/truth/* truth/
```

## Generate gVCFs for Chromosome 20

Once you have a full cohort extract you want to evaluate, we inspect the shards to find the set of shards from `GvsExtractCallset` that encompass chr20 entirely, e.g. looking at the following for different values of `shard-???`

```
gsutil cat alpha1_1000_???.vcf.gz | gunzip | grep -v "#" | cut -f1-8 | head
```

Once the range of shards is known we can gather those into a single chr20 gVCF for all samples with something like:

```
gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I aou_callset_1K_223.vcf.gz \
-I aou_callset_1K_224.vcf.gz \
-I aou_callset_1K_225.vcf.gz \
-I aou_callset_1K_226.vcf.gz \
-I aou_callset_1K_227.vcf.gz \
-I aou_callset_1K_228.vcf.gz \
-I aou_callset_1K_229.vcf.gz \
--output aou_callset_1K.chr20.vcf.gz
```

Now we can create single sample gVCFs for the control samples in AoU

```
INPUT_VCF=aou_callset_1K.chr20.vcf.gz
SOURCE=aou-bq
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name BI_HG-002 --select-type-to-exclude NO_VARIATION -O BI_HG002.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name BI_HG-003 --select-type-to-exclude NO_VARIATION -O BI_HG003.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name UW_HG-002 --select-type-to-exclude NO_VARIATION -O UW_HG002.${SOURCE}.chr20.vcf.gz
```

## Run Script to Add "AS_MAX_VQSLOD" to VCFs

The first line should contain only the control sample names in your callset.

```
for sample in NA12878 SYNDIP BI_HG002 BI_HG003 UW_HG002
do
  python ../add_max_as_vqslod.py ${sample}.gvs.chr20.vcf.gz | bgzip > ${sample}.gvs.chr20.maxas.vcf.gz && tabix ${sample}.gvs.chr20.maxas.vcf.gz
done
```

## Evaluate

First create the ROC curves using only the PASSing records
```
SOURCE="gvs"
BASE_CMD="rtg vcfeval --region chr20 --roc-subset snp,indel  --vcf-score-field=INFO.MAX_AS_VQSLOD -t human_REF_SDF"
SUFFIX="_roc_filtered"

${BASE_CMD} -b truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
-e truth/HG001.gvs.evaluation.bed \
-c NA12878.${SOURCE}.chr20.maxas.vcf.gz -o NA12878_${SOURCE}${SUFFIX}

${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c BI_HG002.${SOURCE}.chr20.maxas.vcf.gz -o BI_HG002_${SOURCE}${SUFFIX}
${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c UW_HG002.${SOURCE}.chr20.maxas.vcf.gz -o UW_HG002_${SOURCE}${SUFFIX}
${BASE_CMD} -b truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG003.gvs.evaluation.bed -c BI_HG003.${SOURCE}.chr20.maxas.vcf.gz -o BI_HG003_${SOURCE}${SUFFIX}
${BASE_CMD} -b truth/CHM.full.38.vcf.gz -e truth/CHM.gvs.evaluation.bed -c SYNDIP.${SOURCE}.chr20.maxas.vcf.gz -o syndip_${SOURCE}${SUFFIX}
```

Then do the same thing but use all records
```
SOURCE="gvs"
BASE_CMD="rtg vcfeval --region chr20 --all-records --roc-subset snp,indel  --vcf-score-field=INFO.MAX_AS_VQSLOD -t human_REF_SDF"
SUFFIX="_roc_all"

${BASE_CMD} -b truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
-e truth/HG001.gvs.evaluation.bed \
-c NA12878.${SOURCE}.chr20.maxas.vcf.gz -o NA12878_${SOURCE}${SUFFIX}

${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c BI_HG002.${SOURCE}.chr20.maxas.vcf.gz -o BI_HG002_${SOURCE}${SUFFIX}
${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c UW_HG002.${SOURCE}.chr20.maxas.vcf.gz -o UW_HG002_${SOURCE}${SUFFIX}
${BASE_CMD} -b truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG003.gvs.evaluation.bed -c BI_HG003.${SOURCE}.chr20.maxas.vcf.gz -o BI_HG003_${SOURCE}${SUFFIX}
${BASE_CMD} -b truth/CHM.full.38.vcf.gz -e truth/CHM.gvs.evaluation.bed -c SYNDIP.${SOURCE}.chr20.maxas.vcf.gz -o syndip_${SOURCE}${SUFFIX}
```

## Gather Precision and Sensitivity numbers

The columns here are: filename, type, FPs (false positives), FNs (false negatives), precision, sensitivity

```
for type in "snp" "indel"
do
    for suffix in "_all" "_filtered"
    do
      for f in $(ls -1 *${suffix}/${type}*.tsv.gz); do
        d=$(cat $f | gunzip | tail -1 | cut -f3,5,6,7)
        s=$(echo $f | cut -d"/" -f1 )
        echo -e "$s\t$type\t$d"
      done
    done
done

```
Paste the output from this call into a TSV file for delivery.

## Optional: View ROCs

You can use wildcards to pull in multiple datasets into the graphical viewer. For example

```
rtg roc NA12878_*_roc*/snp_roc.tsv.gz 
rtg roc NA12878_*_roc*/indel_roc.tsv.gz 
```
