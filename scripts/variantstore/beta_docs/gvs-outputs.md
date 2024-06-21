# Working with the outputs of the Genomic Variant Store

## Applying the GVS Joint Calling Filters
During joint calling using GVS, filters are tagged onto the data based on the filtering model. However, the VCFs we output from GVS are not hard filtered, which means we include the filters in the FT flag and have a recommended filtering approach. Here, we outline that approach.

```
[TODO] Insert filtering best practices here
```

---

**Warning:**
The AC and AN numbers we calculate in the VCFs are based on following the filtering recommendations. This follows GATK best practices.

---

## Working with Sharded VCFs
The VCFs and indices that come out of GVS are sharded based on the size of the callset. They follow approximately this schema:

**Genomes**
0-4999 samples: 1 shard per chr
5,000-19,999 samples: 2,000 shards
20,000-25,000 samples = 10k shards

**Exomes**
0-4999 samples: 1 shard per chr
5000-19,999 samples: 1,000 shards
20,000-49,999 samples: 2,500 shards
>50,000 samples: 7,500 shards

No shard will ever contain data from more than one chromosome, so for callsets over 5k samples, the number of shards will be close to these but not exactly.

You can merge the shards using MergeVcfs from [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037226612-MergeVcfs-Picard).

## Interval lists
The interval lists are named consistently with the vcfs: 00000000.vcf.gz.interval-list will go with 00000000.vcf.gz and 00000000.vcf.gz.tbi

## Content in the VCF [In progress]
These are the filters and format and info fields in the GVS VCFs. Note that they are different from the ones in the GATK WARP Joint Calling VCFs. Additionally, some are new with the introduction of VETS. Read more about VETS [here](https://github.com/broadinstitute/gatk/blob/gvs_0.5.5/scripts/variantstore/docs/release_notes/VETS_Release.pdf). 

```
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=ExcessHet,Description="Site has excess het value larger than the threshold">
##FILTER=<ID=NO_HQ_GENOTYPES,Description="Site has no high quality variant genotypes">
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
##COMMENT=<ID=high_CALIBRATION_SENSITIVITY_SNP,Description="Site failed SNP model calibration sensitivity cutoff (0.997)">
##COMMENT=<ID=high_CALIBRATION_SENSITIVITY_INDEL,Description="Site failed SNP model calibration sensitivity cutoff (0.999)">
```

Exomes will additionally have:

```
##FILTER=<ID=OUTSIDE_OF_TARGETS,Description="Doesn't overlap the user-defined interval list">
```

## Notable content NOT in the VCF
### DP
DP is not in GVS VCFs. It can be approximated from SUM of AD and in most cases should be equal, unless there were some low AF alternates that were dropped in reblocking. AD is a 2 or 3 length vector (only 3 if ½ genotype) and is only included for variant sites.

### PL: phred-scaled genotype likelihoods
PLs (phred-scaled genotype likelihoods) are critically important for downstream analysis, but they can (mostly) be derived from the GQ value. Removing PLs from VCFs has enabled us to reduce file sizes by over an order of magnitude because the typical PL representation leads to an explosion of VCF sites at highly multi-alleleic sites.

Note that normally the array for PL is length "G" where G is the number of genotypes. Instead, we output "." which is length 1, to represent the missing PL format field. 

PL can be approximated as [0, GQ, ~2*GQ]