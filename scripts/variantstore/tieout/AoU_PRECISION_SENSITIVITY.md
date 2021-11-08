# Generating Callset Precision and Sensitivity Values
## Prerequisites (hopefully only do once)
1. Use `conda` to create a fresh environment to add these new tools to it:
 ```
 conda create --name gvs python=3.8
 conda activate gvs
 conda install -c bioconda samtools=1.9 --force-reinstall
 conda install -c bioconda bcftools
 conda install -c bioconda bedtools
 conda install -c bioconda rtg-tools
```
2. Download reference genome from `gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta` and pass it in to `rtg`:
```
rtg format --output human_REF_SDF /path/to/Homo_sapiens_assembly38.fasta
```
3. Copy truth sample VCFs to your local environment. Make sure you have truth files for all your control samples (e.g. `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` for `UW_HG-002` and `BI_HG-002`); you will need these for the "Evaluate" step.  If you are missing any see Appendix A and follow the steps there to add and prepare them.
```
mkdir -p truth && gsutil -m cp gs://broad-dsp-spec-ops/gvs/truth/* truth/
```
4. Complete (at least) through the training step (`GvsCreateFilterSet`) of the [GVS pipeline](../AOU_DELIVERABLES.md) with all the samples (control and non-control).  Then run the callset extract step with a view of the sample-mapping table that's controls-only but with the training (i.e. "filter_set") from all the samples.
## Generate gVCFs for Chromosome 20
Once you have the "controls only" run complete, inspect the shards to find the set of shards from `GvsExtractCallset` that encompass chr20 entirely, e.g. looking at different values of `shard-???` where the output of:
```
gsutil cat alpha1_1000_???.vcf.gz | gunzip | grep -v "#" | cut -f1-8 | head
```
starts with "chr20". Once the range of shards is known, gather those into a single chr20 gVCF for all samples with something like:
```
gatk --java-options -Xms6g GatherVcfsCloud --ignore-safety-checks --gather-type BLOCK \
-I aou_callset_1K_223.vcf.gz \
-I aou_callset_1K_224.vcf.gz \
-I aou_callset_1K_225.vcf.gz \
-I aou_callset_1K_226.vcf.gz \
-I aou_callset_1K_227.vcf.gz \
-I aou_callset_1K_228.vcf.gz \
-I aou_callset_1K_229.vcf.gz \
--output chr20.vcf.gz
```
And create an index for that VCF:
```
tabix chr20.vcf.gz
```
Now create single sample gVCFs for the control samples; in this example the sample names for the controls are "BI_HG-002", "UW_HG-002" and "BI_HG-003":
```
INPUT_VCF=chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name BI_HG-002 --select-type-to-exclude NO_VARIATION -O BI_HG-002.chr20.vcf.gz
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name UW_HG-002 --select-type-to-exclude NO_VARIATION -O UW_HG-002.chr20.vcf.gz```
~/gatk SelectVariants -V ${INPUT_VCF} -L chr20 --sample-name BI_HG-003 --select-type-to-exclude NO_VARIATION -O BI_HG-003.chr20.vcf.gz
```
## Run Script to Add "AS_MAX_VQSLOD" to VCFs
The first line should contain only the control sample names in your callset (same list as previous step).
```
for sample in BI_HG-002 UW_HG-002 BI_HG-003
do
  python ../add_max_as_vqslod.py ${sample}.chr20.vcf.gz | bgzip > ${sample}.chr20.maxas.vcf.gz && tabix ${sample}.chr20.maxas.vcf.gz
done
```
## Evaluate
First create the ROC curves using only the PASSing records.  For each sample, you will need to find the corresponding files under `truth/` and pass those into the `BASE_CMD`.

```
BASE_CMD="rtg vcfeval --region chr20 --roc-subset snp,indel  --vcf-score-field=INFO.MAX_AS_VQSLOD -t human_REF_SDF"
SUFFIX="_roc_filtered"

${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c BI_HG-002.chr20.maxas.vcf.gz -o BI_HG-002_aou-bq${SUFFIX}
${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c UW_HG-002.chr20.maxas.vcf.gz -o UW_HG-002_aou-bq${SUFFIX}
${BASE_CMD} -b truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG003.gvs.evaluation.bed -c BI_HG-003.chr20.maxas.vcf.gz -o BI_HG-003_aou-bq${SUFFIX}
```
Then do the same thing but use all records.
```
SOURCE="gvs"
BASE_CMD="rtg vcfeval --region chr20 --all-records --roc-subset snp,indel  --vcf-score-field=INFO.MAX_AS_VQSLOD -t human_REF_SDF"
SUFFIX="_roc_all"

${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c BI_HG-002.chr20.maxas.vcf.gz -o BI_HG-002_aou-bq${SUFFIX}
${BASE_CMD} -b truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG002.gvs.evaluation.bed -c UW_HG-002.chr20.maxas.vcf.gz -o UW_HG-002_aou-bq${SUFFIX}
${BASE_CMD} -b truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e truth/HG003.gvs.evaluation.bed -c BI_HG-003.chr20.maxas.vcf.gz -o BI_HG-003_aou-bq${SUFFIX}
```
## Gather Precision and Sensitivity numbers
Run the following script to display the data for the output file:
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
The columns here are: sample name, variant type, FPs (false positives), FNs (false negatives), precision, and sensitivity. Copy and paste the output from this call into a TSV file.
## Optional: View ROCs
You can use wildcards to pull in multiple datasets into the graphical viewer. For example
```
rtg roc NA12878_*_roc*/snp_roc.tsv.gz 
rtg roc NA12878_*_roc*/indel_roc.tsv.gz 
```
## Appendix A: New Control Samples
"Truth" versions of control samples are from the [Genome in a Bottle Consortium (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) and [releases of their data can be found here](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/) organized by source.
1. Download the `.bed`, `.vcf` and index files associated with the control sample.
2. Convert the calling region GVS is using to BED format
```
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list .
gatk IntervalListToBed -I wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list -O wgs_calling_regions.hg38.noCentromeres.noTelomeres.bed
rm wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list
```
3. Then intersect with each of the truth datasets (this example is for HG-001):
```
bedtools intersect -a wgs_calling_regions.hg38.noCentromeres.noTelomeres.bed -b HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed > HG001.gvs.evaluation.bed
bedtools intersect -a wgs_calling_regions.hg38.noCentromeres.noTelomeres.bed -b HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > HG002.gvs.evaluation.bed
bedtools intersect -a wgs_calling_regions.hg38.noCentromeres.noTelomeres.bed -b HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > HG003.gvs.evaluation.bed
bedtools intersect -a wgs_calling_regions.hg38.noCentromeres.noTelomeres.bed -b CHM.full.38.bed.gz > CHM.gvs.evaluation.bed
```
4. Copy the GIAB data and the new `evaluation.bed` file to the `gs://broad-dsp-spec-ops/gvs/truth/` directory so the files will be available for future callsets.
