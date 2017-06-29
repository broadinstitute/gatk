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

# Here we implement steps 1 and 2 for several replicates of a single plex, eg several replicates of the 10-plex pool.

# Uses select_variants to create a VCF that is subsetted by sample and by interval list
task CreateSubVcf {
    File gatk_jar
    File hapmap
    File hapmap_idx           # only here so the index file is copied along, or else SelectVariants will create it.
    File sample_file
    File? intervals           # the primary intervals to restrict to
    File common_vcf         # hapmap variants confirmed by HaplotypeCaller to reduce false positives
    File common_vcf_idx

    command {
        # subsampling and restriction to biallelics and intervals
        java -jar ${gatk_jar} SelectVariants -V ${hapmap} -O sub.vcf \
            -restrictAllelesTo BIALLELIC --sample_name ${sample_file} \
            ${"-L " + intervals} -L ${common_vcf} --interval_set_rule INTERSECTION \
            -maxIndelSize 10 \
            --excludeNonVariants

         #remove NEGATIVE_TRAIN_SITE variants and re-index
         grep -v NEGATIVE_TRAIN_SITE sub.vcf > subsampled_hapmap.vcf
         java -jar ${gatk_jar} IndexFeatureFile -F subsampled_hapmap.vcf
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
    Int max_depth   #ignore sites with depth greater than this because they are alignment artifacts

    command {
        java -jar ${gatk_jar} AnnotateVcfWithBamDepth -V ${input_vcf} -I ${bam} -O "bam_depth.vcf"
        java -jar ${gatk_jar} SelectVariants -V bam_depth.vcf --select "BAM_DEPTH < ${max_depth}" -O truth.vcf
        java -jar ${gatk_jar} IndexFeatureFile -F truth.vcf
    }

    output {
        File output_vcf = "truth.vcf"
        File output_vcf_idx = "truth.vcf.idx"
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
        java -jar ${gatk_jar} AnnotateVcfWithExpectedAlleleFraction -V ${input_vcf} -O af_exp.vcf \
            --mixingFractions  ${mixing_fractions_table}
    }

    output {
        File output_vcf = "af_exp.vcf"
        File output_vcf_idx = "af_exp.vcf.idx"
    }
}

workflow HapmapSensitivityTruth {
    File gatk_jar
    File hapmap
    File hapmap_idx
    File pooled_bams_list
    Array[Array[String]] bams_and_index = read_tsv(pooled_bams_list)
    File sample_file
    File? intervals
    Int min_indel_spacing
    File common_vcf         # common variants eg gnomad variants with AF > 0.001
    File common_vcf_idx
    Int max_depth


    call CreateSubVcf {
        input:
            gatk_jar = gatk_jar,
            hapmap = hapmap,
            hapmap_idx = hapmap_idx,
            sample_file = sample_file,
            intervals = intervals,
            common_vcf = common_vcf,
            common_vcf_idx = common_vcf_idx
    }

    call RemoveNearbyIndels {
        input:
            gatk_jar = gatk_jar,
            input_vcf = CreateSubVcf.output_vcf,
            input_vcf_idx = CreateSubVcf.output_vcf_idx,
            min_indel_spacing = min_indel_spacing
    }

    scatter(row in bams_and_index) {
        File bam = row[0]
        File index = row[1]

        call CalculateMixingFractions {
            input:
                gatk_jar = gatk_jar,
                input_vcf = CreateSubVcf.output_vcf,
                input_vcf_idx = CreateSubVcf.output_vcf_idx,
                bam = bam,
                bam_idx = index
        }

        call AnnotateVcfWithExpectedAlleleFraction {
            input:
                gatk_jar = gatk_jar,
                input_vcf = RemoveNearbyIndels.output_vcf,
                input_vcf_idx = RemoveNearbyIndels.output_vcf_idx,
                mixing_fractions_table = CalculateMixingFractions.output_table,
        }

        call AnnotateVcfWithBamDepth {
                input:
                    gatk_jar = gatk_jar,
                    input_vcf = AnnotateVcfWithExpectedAlleleFraction.output_vcf,
                    input_vcf_idx = AnnotateVcfWithExpectedAlleleFraction.output_vcf_idx,
                    bam = bam,
                    bam_idx = index,
                    max_depth = max_depth
            }
    }

    output {
        Array[File] mixing_fractions_table = CalculateMixingFractions.output_table
        Array[File] output_vcf = AnnotateVcfWithBamDepth.output_vcf
    }
}

