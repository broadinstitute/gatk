### This is the truth data preparation part of the CRSP sensitivity validation

#Conceptual Overview
# To measure sensitivity, we sequence a pool 5, 10 or 20 normal samples. The pool contains a variety of allele fractions,
# like a real tumor sample. These normal samples were also sequenced individually as part of HapMap, so we have "truth" vcfs
# for them.  We calculate sensitivity by comparing the calls to the truth data.  For a variety of reasons we don't call
# against a matched normal.

#Workflow Steps
# 1.  Restrict a huge HapMap wgs VCF to given lists of samples and intervals.
# 2.  Annotate this VCF using a bam sequenced from a pool derived from the given samples
# 3.  Run Mutect in tumor-only mode on this pooled bam.
# 4.  Compare Mutect calls to the truth data and output a table of true positives and false negatives along with
#     annotations from the truth VCF prepared in steps 1 and 2.

# Here we implement steps 1 and 2.

# Uses select_variants to create a VCF that is subsetted by sample and by interval list
task CreateSubVcf {
    File gatk_jar
    File hapmap
    File hapmap_idx        # only here so the index file is copied along, or else SelectVariants will create it.
    File sample_file
    File intervals            # the primary intervals to restrict to
    File dbsnp                # to filter out false positives
    File dbsnp_idx

    command {
        # subsampling and restriction to biallelics and intervals
        java -jar ${gatk_jar} SelectVariants -V ${hapmap} -O sub.vcf \
            -restrictAllelesTo BIALLELIC --sample_file ${sample_file} -L ${intervals} -excludeNonVariants

        #SNPs, filtered for overlap with dbSNP
        java -jar ${gatk_jar} SelectVariants -V sub.vcf -O snps.vcf -L ${dbsnp} -selectType SNP

        #indels
        java -jar ${gatk_jar} SelectVariants -V sub.vcf -O indels.vcf -selectType INDEL

        #recombine SNPs and indel
        java -jar ${gatk_jar} MergeVcfs -I snps.vcf -I indels.vcf -O subsampled_hapmap.vcf
    }

    output {
        File output_vcf = "subsampled_hapmap.vcf"
        File output_vcf_idx = "subsampled_hapmap.vcf.idx"
    }
}

task RemoveNearbyIndels {
    File gatk_jar
    File input_vcf
    File input_vcf_idx
    Int min_indel_spacing

    command {
        java -jar ${gatk_jar} RemoveNearbyIndels -V ${input_vcf} -O close_indel_removed.vcf -minIndelSpacing ${min_indel_spacing}
    }

    output {
        File output_vcf = "close_indel_removed.vcf"
        File output_vcf_idx = "close_indel_removed.vcf.idx"
    }
}

task AnnotateVcfWithBamDepth {
    File gatk_jar
    File input_vcf
    File input_vcf_idx
    File bam
    File bam_idx

    command {
        java -jar ${gatk_jar} AnnotateVcfWithBamDepth -V ${input_vcf} -I ${bam} -O "bam_depth.vcf"
    }

    output {
        File output_vcf = "bam_depth.vcf"
        File output_vcf_idx = "bam_depth.vcf"
    }
}

task CalculateMixingFractions {
    File gatk_jar
    File input_vcf
    File input_vcf_idx
    File bam
    File bam_idx

    command {
        java -jar ${gatk_jar} CalculateMixingFractions -V ${input_vcf} -I ${bam} -O "mixing_fractions.table"
    }

    output {
        File output_table = "mixing_fractions.table"
    }
}

task AnnotateVcfWithExpectedAlleleFraction {
    File gatk_jar
    File input_vcf
    File input_vcf_idx
    File mixing_fractions_table

    command {
        java -jar ${gatk_jar} AnnotateVcfWithExpectedAlleleFraction -V ${input_vcf} -O truth.vcf \
            --mixingFractions  ${mixing_fractions_table}
    }

    output {
        File output_vcf = "truth.vcf"
        File output_vcf_idx = "truth.vcf.idx"
    }
}

workflow HapmapSensitivityTruth {
    File gatk_jar
    File hapmap
    File hapmap_idx
    File pooled_bam
    File pooled_bam_idx
    File sample_file
    File intervals
    File dbsnp
    File dbsnp_idx
    Int min_indel_spacing


    call CreateSubVcf {
        input:
            gatk_jar = gatk_jar,
            hapmap = hapmap,
            hapmap_idx = hapmap_idx,
            sample_file = sample_file,
            intervals = intervals,
            dbsnp = dbsnp,
            dbsnp_idx = dbsnp_idx,
    }

    call RemoveNearbyIndels {
        input:
            gatk_jar = gatk_jar,
            input_vcf = CreateSubVcf.output_vcf,
            input_vcf_idx = CreateSubVcf.output_vcf_idx,
            min_indel_spacing = min_indel_spacing
    }

    call AnnotateVcfWithBamDepth {
        input:
            gatk_jar = gatk_jar,
            input_vcf = RemoveNearbyIndels.output_vcf,
            input_vcf_idx = RemoveNearbyIndels.output_vcf_idx,
            bam = pooled_bam,
            bam_idx = pooled_bam_idx
    }

    call CalculateMixingFractions {
        input:
            gatk_jar = gatk_jar,
            input_vcf = CreateSubVcf.output_vcf,
            input_vcf_idx = CreateSubVcf.output_vcf_idx,
            bam = pooled_bam,
            bam_idx = pooled_bam_idx
    }

    call AnnotateVcfWithExpectedAlleleFraction {
        input:
            gatk_jar = gatk_jar,
            input_vcf = AnnotateVcfWithBamDepth.output_vcf,
            input_vcf_idx = AnnotateVcfWithBamDepth.output_vcf_idx,
            mixing_fractions_table = CalculateMixingFractions.output_table,
    }

    output {
        File mixing_fractions_table = CalculateMixingFractions.output_table
        File output_vcf = AnnotateVcfWithExpectedAlleleFraction.output_vcf
    }

}

