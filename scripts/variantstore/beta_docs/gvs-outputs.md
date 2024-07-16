# Working with the outputs of the Genomic Variant Store

## Content in the VCF
This is the full list of the filters, format, comment, and info fields in the GVS VCFs. Note that they are different from the ones in the GATK WARP Joint Calling VCFs. Additionally, some are new with the introduction of VETS. Read more about VETS [here](https://github.com/broadinstitute/gatk/blob/gvs_0.5.5/scripts/variantstore/docs/release_notes/VETS_Release.pdf).

```
##FILTER=<ID=ExcessHet,Description="Site has excess het value larger than the threshold">
##FILTER=<ID=NO_HQ_GENOTYPES,Description="Site has no high quality variant genotypes">
##FILTER=<ID=EXCESS_ALLELES,Description="Site has an excess of alternate alleles based on the input threshold">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype Filter Field">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AS_QUALapprox,Number=1,Type=String,Description="Allele-specific QUAL approximations">
##INFO=<ID=QUALapprox,Number=1,Type=Integer,Description="Sum of PL[0] values; used to approximate the QUAL score">
##COMMENT=<ID=high_CALIBRATION_SENSITIVITY_SNP,Description="Sample Genotype FT filter value indicating that the genotyped allele failed SNP model calibration sensitivity cutoff (0.997)">
##COMMENT=<ID=high_CALIBRATION_SENSITIVITY_INDEL,Description="Sample Genotype FT filter value indicating that the genotyped allele failed INDEL model calibration sensitivity cutoff (0.99)">
```

Exomes will additionally have:

```
##FILTER=<ID=OUTSIDE_OF_TARGETS,Description="Outside of sequencing target intervals">
```

## Notable content NOT in the VCF
### DP
DP is not in GVS VCFs. It can be approximated from SUM of AD and in most cases should be equal, unless there were some low AF alternates that were dropped in reblocking. AD is a 2 or 3 length vector (only 3 if ½ genotype) and is only included for variant sites.

### PL: phred-scaled genotype likelihoods
PLs (phred-scaled genotype likelihoods) are expected by some downstream tools, but they can (mostly) be derived from the GQ value. Removing PLs from VCFs has enabled us to reduce file sizes by over an order of magnitude because the typical PL representation leads to an explosion of VCF sites at highly multi-alleleic sites.

Note that normally the array for PL is length "G" where G is the number of genotypes. Instead, we output "." which is length 1, to represent the missing PL format field.

PL can be approximated as [0, GQ, ~2*GQ] for biallelic sites.

### Missing (`./.`) and `GQ0 0/0`
Reblocking, a step required before running GVS, combines `./.` and `0/0`. So GVS VCFs will output `./.` and `0/0 GQ0` calls as `./.`.

## Applying the GVS Joint Calling Filters
During joint calling using GVS, filters are evaluated at the site and allele level. All data is retained, but sites failing the site filters will be indicated by one of the site level filters from the table below appearing in the VCF FILTER field. Genotypes that have been filtered by the Genotype filters in the table below will be tagged in the FORMAT FT field. Variants that pass site level filters will have “.”.

| Filter Names                       | Type of Filter | Definition                                                                                                                    | Details                                                                                                                                                                                                                                                                                                                                                                                                                               |
|------------------------------------|----------------|-------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| NO_HQ_GENOTYPES                    | Site Level     | Site has *no* high quality variant genotypes                                                                                  | No high-quality genotype (GQ>=20, DP>=10, and AB>=0.2 for heterozygotes) called for the variant. If there is one genotype at the variant site, the filter will not be applied and the variant site will pass. Allele Balance (AB) is calculated for each heterozygous variant as the number of bases supporting the least-represented allele over the total number of base observations.  In other words, min(AD)/DP for diploid GTs. |
| ExcessHet                          | Site level     | Site has excess het value larger than the threshold                                                                           | Genotypes with this filter show a higher proportion of heterozygotes than expected under Hardy-Weinberg equilibrium, which can be indicative of mapping errors in the reads contributing to the variant. The filter is applied to anything more extreme than a z-score of -4.5 (p-value of 3.4e-06), which phred-scaled is 54.69.                                                                                                     |
| LowQual                            | Site level     | Low Quality                                                                                                                   | QUALapprox is too low (lower than 60 for SNPs; lower than 69 for Indels). QUAL tells you how confident we are that there is some kind of variation at a given site. The variation may be present in one or more samples.                                                                                                                                                                                                              |
| EXCESS_ALLELES                     | Site level     | Site has an excess of alternate alleles based on the input threshold                                                          | When a site has more than 100 alternate alleles, this filter will be present.                                                                                                                                                                                                                                                                                                                                                         |
| high_CALIBRATION_SENSITIVITY_SNP   | Genotype level | Sample Genotype FT filter value indicating that the genotyped allele failed SNP model calibration sensitivity cutoff (0.997)  | The VETS filtering score. This is a ##COMMENT in the VCF because this is not applied as a site-level filter. GVS is allele-specific. The FT for an allele will have this FT tag if it hits the CALIBRATION_SENSITIVITY threshold. A site can pass with failing genotypes, however we recommend filtering these genotypes. See example variants below.                                                                                 |
| high_CALIBRATION_SENSITIVITY_INDEL | Genotype level | Sample Genotype FT filter value indicating that the genotyped allele failed INDEL model calibration sensitivity cutoff (0.99) | The VETS filtering score. This is a ##COMMENT in the VCF because this is not applied as a site-level filter. GVS is allele-specific. The FT for an allele will have this FT tag if it hits the CALIBRATION_SENSITIVITY threshold. A site can pass with failing genotypes, however we recommend filtering these genotypes. See example variants below.                                                                                 |
| OUTSIDE_OF_TARGETS                 | Site level     | Outside of sequencing target intervals                                                                                        | Exome only. The site is not within the target intervals of the exome assay. We recommend filtering out these sites as we cannot stand behind the quality of sites called outside of the targets.                                                                                                                                                                                                                                      |

To post process to keep only data that was not filtered at the site level, you could use this command:

```
bcftools view -f 'PASS,.' -Oz -o output.vcf.gz
```

You can also use bcftools to filter based on the genotype level filters based on your analysis needs:

Filter out spanning deletions and variants with an AC of 0
```
bcftools view -e 'ALT[0]="*" || AC=0' input.vcf.gz -Oz -o output.vcf.gz
```

Ensure that we respect the FT tag to exclude alleles with high_CALIBRATION_SENSITIVITY_SNP and INDEL and recalculate AC, AN, AF
```
bcftools filter -i "FORMAT/FT='PASS,.'" --set-GTs . input.vcf.gz -Oz -o output.vcf.gz
```


### Example sites and how to interpret them
#### Passing site
```
chr20   66369   .       C       T       .       .       AC=2;AF=1.00;AN=2;AS_QUALapprox=0|982;CALIBRATION_SENSITIVITY=0.2492;QUALapprox=982;SCORE=-0.3836      GT:AD:GQ:RGQ    ./.         ./.     1/1:0,29:87:982
```

#### A site with failing calibration sensitivity SNPs for all sample genotypes
```
chr20   62898   .       T       A       .       .       AS_QUALapprox=0|933;CALIBRATION_SENSITIVITY=0.9999;QUALapprox=324;SCORE=-0.7456        GT:AD:FT:GQ:PGT:PID:RGQ 0/1:26,11:high_CALIBRATION_SENSITIVITY_SNP:99:0|1:62818_T_G:384     0/1:23,7:high_CALIBRATION_SENSITIVITY_SNP:99:0|1:62818_T_G:225  0/1:18,9:high_CALIBRATION_SENSITIVITY_SNP:99:0|1:62818_T_G:324
```

#### Multi-allelic site with one passing and two failed sample genotypes
```
chr20   245443  .       GTATATATA       GTATATA,G       .       .       AC=0,1;AF=0.00,0.500;AN=2;AS_QUALapprox=0|339|238;CALIBRATION_SENSITIVITY=0.9957,0.9402;QUALapprox=238;SCORE=-0.6448,-0.5305     GT:AD:FT:GQ:RGQ 0/2:11,0,9:PASS:99:238  0/1:10,5,0:high_CALIBRATION_SENSITIVITY_INDEL:99:264    0/1:5,3,0:high_CALIBRATION_SENSITIVITY_INDEL:75:75
```

#### A failing site due to site level filters
```
chr1	10177	.	A	C	.	LowQual;NO_HQ_GENOTYPES	AC=1;AF=0.500;AN=2;AS_QUALapprox=0|26;CALIBRATION_SENSITIVITY=.;QUALapprox=26;SCORE=.	GT:AD:GQ:PGT:PID:RGQ	0/1:1,2:26:0|1:10145_AAC_A:26	./.	./.
```

---

**Note**
The AC and AN numbers that we calculate for sites with passing genotypes are calculated based on following the allele filtering in the VCF.

Currently, there will be no AC or AN output when all alleles fail filtering. For sites with only site level filters, we do calculate AC and AN.

---

## Working with Sharded VCFs
The VCFs and indices that come out of GVS are sharded based on the size of the callset. They follow approximately this schema:

**Genomes**

| Number of Samples | Number of Shards |
|-------------------|----------------- |
| 0-4,999           | 1 shard per chromosome |
| 5,000-19,999      | 2,000 shards  |
| 20,000-25,000     | 10,000 shards |

**Exomes**

| Number of Samples | Number of Shards |
|-------------------|----------------- |
| 0-4,999           | 1 shard per chromosome |
| 5,000-19,999      | 1,000 shards  |
| 20,000-49,999     | 2,500 shards |
| > 50,000          | 7,500 shards |

No shard will ever contain data from more than one chromosome, so for callsets over 5k samples, the number of shards will be close to these but not exactly.

### Merging Shards
You can merge the VCF shards using [GatherVCFsCloud from GATK](https://app.terra.bio/#workspaces/help-gatk/GATK4-Germline-Preprocessing-VariantCalling-JointCalling/workflows/broad-firecloud-dsde/Optional-Gatk-GatherVCFsCloud).

Other tools, in order of our recommendation, that are available for working with VCFs include:
- [Hail](https://www.hail.is/) (see below).
- `merge` from [bcftools](https://samtools.github.io/bcftools/).
- `MergeVcfs` from [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037226612-MergeVcfs-Picard).

#### Merging GVS outputs in Hail
Here is a block of Hail code that will convert a GVS VCF to a Hail Matrix Table (MT).

```
mt = hl.import_vcf(vcf_bgz, force_bgz=True)
mt.write(output_gs_url)
```

## Interval lists
The interval lists are named consistently with the vcfs: 00000000.vcf.gz.interval-list will go with 00000000.vcf.gz and 00000000.vcf.gz.tbi
