### This is the truth data preprocessing part of the CRSP sensitivity validation

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

# Here we implement the subsampling part of step 1, along with some other preprocessing such as removing pairs of
# neighboring indels.

workflow PreprocessHapmap {
    # inputs
    File hapmap
    File hapmap_idx
    File five_plex_samples
    File ten_plex_samples
    File twenty_plex_samples
    Int  min_indel_spacing
    File gnomad         # common variants eg gnomad variants with AF > 0.001
    File gnomad_idx

    File? gatk_override

    # runtime
    String gatk_docker

    call Subsample as SubsampleFive {
        input: gatk_override = gatk_override, gatk_docker = gatk_docker, hapmap = hapmap, hapmap_idx = hapmap_idx, samples = five_plex_samples, gnomad = gnomad, gnomad_idx = gnomad_idx
    }

    call Subsample as SubsampleTen {
        input: gatk_override = gatk_override, gatk_docker = gatk_docker,  hapmap = hapmap, hapmap_idx = hapmap_idx, samples = ten_plex_samples, gnomad = gnomad, gnomad_idx = gnomad_idx
    }

    call Subsample as SubsampleTwenty {
         input: gatk_override = gatk_override, gatk_docker = gatk_docker,  hapmap = hapmap, hapmap_idx = hapmap_idx, samples = twenty_plex_samples, gnomad = gnomad, gnomad_idx = gnomad_idx
    }

    call RemoveNearbyIndels as RemoveFive {
        input: gatk_override = gatk_override, gatk_docker = gatk_docker,  input_vcf = SubsampleFive.output_vcf, input_vcf_idx = SubsampleFive.output_vcf_idx,
            min_indel_spacing = min_indel_spacing, name = "five_plex"
    }

    call RemoveNearbyIndels as RemoveTen {
        input: gatk_override = gatk_override, gatk_docker = gatk_docker,  input_vcf = SubsampleTen.output_vcf, input_vcf_idx = SubsampleTen.output_vcf_idx,
            min_indel_spacing = min_indel_spacing, name = "ten_plex"
    }

    call RemoveNearbyIndels as RemoveTwenty {
        input: gatk_override = gatk_override, gatk_docker = gatk_docker,  input_vcf = SubsampleTwenty.output_vcf, input_vcf_idx = SubsampleTwenty.output_vcf_idx,
            min_indel_spacing = min_indel_spacing, name = "twenty_plex"
    }

    output {
        File five_plex_preprocessed = RemoveFive.output_vcf
        File five_plex_preprocessed_idx = RemoveFive.output_vcf
        File ten_plex_preprocessed = RemoveTen.output_vcf
        File ten_plex_preprocessed_idx = RemoveTen.output_vcf
        File twenty_plex_preprocessed = RemoveTwenty.output_vcf
        File twenty_plex_preprocessed_idx = RemoveTwenty.output_vcf
    }
}

task Subsample {
    # inputs
    File hapmap
    File hapmap_idx
    File samples
    File gnomad           # common variants, to reduce false positives eg mapping artifacts in the original Hapmap vcf
    File gnomad_idx

    File? gatk_override

    # runtime
    String gatk_docker

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        # subsampling and restriction to biallelics
        gatk --java-options "-Xmx4g" SelectVariants -V ${hapmap} -O sub.vcf \
            -restrict-alleles-to BIALLELIC \
            --sample-name ${samples} \
            -L ${gnomad} \
            -max-indel-size 10 \
            --exclude-non-variants

         #remove NEGATIVE_TRAIN_SITE variants and re-index
         grep -v NEGATIVE_TRAIN_SITE sub.vcf > subsampled.vcf
         gatk --java-options "-Xmx4g" IndexFeatureFile -F subsampled.vcf
    }

    runtime {
        docker: "${gatk_docker}"
        preemptible: 2
    }

    output {
        File output_vcf = "subsampled.vcf"
        File output_vcf_idx = "subsampled.vcf.idx"
    }
}

task RemoveNearbyIndels {
    # inputs
    File input_vcf
    File input_vcf_idx
    Int min_indel_spacing
    String name

    File? gatk_override

    # runtime
    String gatk_docker

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx4g" RemoveNearbyIndels -V ${input_vcf} -O ${name}.vcf -min-indel-spacing ${min_indel_spacing}
    }

    runtime {
        docker: "${gatk_docker}"
        preemptible: 2
    }

    output {
        File output_vcf = "${name}.vcf"
        File output_vcf_idx = "${name}.vcf.idx"
    }
}

