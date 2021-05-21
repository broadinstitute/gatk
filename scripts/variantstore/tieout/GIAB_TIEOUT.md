# NOTES
 - Previously there was an attempt to tie out the feature extract input to VQSR between WARP and GVS.  However, that was abandoned in favor of the truth-sample based comparison here due to non-determinism and bugs in GenomicsDB that were hard to compensate for.  However, the docs and scripts for doing this are in the branch `gvs_vqsr_input_tieout` under `scripts/variantstore/tieout/feature_extract`
 
## Prerequisites

 - rtg
 - bcftools
 - tabix
 - python 3.7+
 - gatk
 
 
 We suggest using conda to create a fresh environment to add these new tools to
 ```
 conda create --name gvs python=3.8
 conda activate gvs
 conda install -c bioconda samtools=1.9 --force-reinstall
 conda install -c bioconda bcftools
 conda install -c bioconda rtg-tools
```

## One time tasks

```
# create SDF for reference (used by rtg eval)
# use local path to GRC38 reference fasta
REFERENCE="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"
rtg format --output human_REF_SDF $REFERENCE

# copy down truth data locally
mkdir -p truth && gsutil -m cp gs://broad-dsp-spec-ops/gvs/truth/* truth/
```

## Obtain Truth sample VCFs 

First, create a full cohort extract (as described in README.md) using the desired filter_set_name.  Assuming this is in a single gathered VCF of `gvs.vcf.gz`

```
# subset to chr20
bcftools view -O z gvs.vcf.gz chr20 > gvs.chr20.vcf.gz
tabix gvs.chr20.vcf.gz

# extract each of the samples
# NOTE: excluding NO_VARIATION just means to drop sites that are reference in the selected sample
INPUT_VCF="gvs.chr20.vcf.gz"

gatk SelectVariants -V ${INPUT_VCF} --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.gvs.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name CHMI_CHMI3_WGS1 --select-type-to-exclude NO_VARIATION -O SYNDIP.gvs.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202194294 --select-type-to-exclude NO_VARIATION -O BI_HG002.gvs.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202243693 --select-type-to-exclude NO_VARIATION -O BI_HG003.gvs.chr20.vcf.gz
gatk SelectVariants -V ${INPUT_VCF} --sample-name 573673 --select-type-to-exclude NO_VARIATION -O UW_HG002.gvs.chr20.vcf.gz

```

## script to add "AS_MAX_VQSLOD" to VCFs
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

The do the same thing but use all records
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

The columns here are: filename, type, FPs, FNs, precision, sensitivity

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
## View ROCs

You can use wildcards to pull in multiple datasets into the graphical viewer. For example

```
rtg roc NA12878_*_roc*/snp_roc.tsv.gz 
rtg roc NA12878_*_roc*/indel_roc.tsv.gz 
```

## Appendix A: Processing WARP results against gvs_sci_tieout_acmg_cohort workspace

Starting with the output of FinalGatherVcf run from the WDLs and inputs in the `wdl/legacy` directory.

If this was run on the specops cromwell you can get that final VCF with:

```
WORKFLOW_ID=49f5c0a5-dbfb-4b05-b293-11ddadeae8e5
gsutil -m cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-FinalGatherVcf/warp_tieout_acmg_cohort_v1.vcf.gz* .
```

Next is to subset to chr20, and the select out the 5 truth samples into their own VCFs.  Assuming a gathered initial VCF of "warp_tieout_acmg_cohort_v1.vcf.gz"

```
bcftools view -O z warp_tieout_acmg_cohort_v1.vcf.gz chr20 > warp_tieout_acmg_cohort_v1.chr20.vcf.gz
tabix warp_tieout_acmg_cohort_v1.chr20.vcf.gz

INPUT_VCF="warp_tieout_acmg_cohort_v1.chr20.vcf.gz"
SOURCE="warp"
~/gatk SelectVariants -V ${INPUT_VCF} --sample-name SM-G947Y --select-type-to-exclude NO_VARIATION -O NA12878.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} --sample-name CHMI_CHMI3_WGS1 --select-type-to-exclude NO_VARIATION -O SYNDIP.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202194294 --select-type-to-exclude NO_VARIATION -O BI_HG002.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} --sample-name 1202243693 --select-type-to-exclude NO_VARIATION -O BI_HG003.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} --sample-name 573673 --select-type-to-exclude NO_VARIATION -O UW_HG002.${SOURCE}.chr20.vcf.gz
```

From here it's a similar process to the evaluation of the GVS-based data

## Appendix B: AoU Tieout

The processing for extracting these samples from a large AoU joint callset is slightly different.  First we inspect the shards to find the set of shards that encompass chr20 entirely.

E.g. looking at the following for different values of `shard-???`

```
gsutil cat gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/8331ab64-7bb5-4b1c-9c6f-06765bfd6b15/GvsExtractCallset/e045f320-4f58-440c-9680-6c557c6538ee/call-ExtractTask/shard-223/alpha1_1000_*.vcf.gz | gunzip | grep -v "#" | cut -f1-8 | head
```

Once the range of shards is known we can gather those into a single chr20 gVCF for all samples with something like:

```
gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-223/alpha1_1000_223.vcf.gz \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-224/alpha1_1000_224.vcf.gz \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-225/attempt-2/alpha1_1000_225.vcf.gz \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-226/alpha1_1000_226.vcf.gz \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-227/alpha1_1000_227.vcf.gz \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-228/alpha1_1000_228.vcf.gz \
-I gs://fc-secure-daa074a0-884f-4b67-a79b-cb5cd1556c9f/e73809fa-9ef2-4f3b-965c-846eff4b905e/GvsExtractCallset/48c56780-2520-444f-b1b3-83091a74bb08/call-ExtractTask/shard-229/alpha1_1000_229.vcf.gz \
--output aou.alpha1v2.chr19-21.vcf.gz
```

Now we can create single sample gVCFs for the 3 GIAB samples in AoU

```
INPUT_VCF=aou.alpha1v2.chr19-21.vcf.gz
SOURCE=aou-bq
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name BI_HG-002 --select-type-to-exclude NO_VARIATION -O BI_HG002.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name BI_HG-003 --select-type-to-exclude NO_VARIATION -O BI_HG003.${SOURCE}.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name UW_HG-002 --select-type-to-exclude NO_VARIATION -O UW_HG002.${SOURCE}.chr20.vcf.gz
```

These can be analyzed similar to the initial GVS tieout
