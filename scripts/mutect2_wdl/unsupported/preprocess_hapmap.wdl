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
    File gatk
    File hapmap
    File hapmap_idx
    File five_plex_samples
    File ten_plex_samples
    File twenty_plex_samples
    Int  min_indel_spacing
    File gnomad         # common variants eg gnomad variants with AF > 0.001
    File gnomad_idx

    call Subsample as SubsampleFive {
        input: gatk = gatk, hapmap = hapmap, hapmap_idx = hapmap_idx, samples = five_plex_samples, gnomad = gnomad, gnomad_idx = gnomad_idx
    }

    call Subsample as SubsampleTen {
        input: gatk = gatk, hapmap = hapmap, hapmap_idx = hapmap_idx, samples = ten_plex_samples, gnomad = gnomad, gnomad_idx = gnomad_idx
    }

    call Subsample as SubsampleTwenty {
         input: gatk = gatk, hapmap = hapmap, hapmap_idx = hapmap_idx, samples = twenty_plex_samples, gnomad = gnomad, gnomad_idx = gnomad_idx
    }

    call RemoveNearbyIndels as RemoveFive {
        input: gatk = gatk, input_vcf = SubsampleFive.output_vcf, input_vcf_idx = SubsampleFive.output_vcf_idx,
            min_indel_spacing = min_indel_spacing, name = "five_plex"
    }

    call RemoveNearbyIndels as RemoveTen {
        input: gatk = gatk, input_vcf = SubsampleTen.output_vcf, input_vcf_idx = SubsampleTen.output_vcf_idx,
            min_indel_spacing = min_indel_spacing, name = "ten_plex"
    }

    call RemoveNearbyIndels as RemoveTwenty {
        input: gatk = gatk, input_vcf = SubsampleTwenty.output_vcf, input_vcf_idx = SubsampleTwenty.output_vcf_idx,
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
    File gatk
    File hapmap
    File hapmap_idx
    File samples
    File gnomad           # common variants, to reduce false positives eg mapping artifacts in the original Hapmap vcf
    File gnomad_idx

    command {
        # subsampling and restriction to biallelics
        java -jar ${gatk} SelectVariants -V ${hapmap} -O sub.vcf \
            -restrictAllelesTo BIALLELIC \
            --sample_name ${samples} \
            -L ${gnomad} \
            -maxIndelSize 10 \
            --excludeNonVariants

         #remove NEGATIVE_TRAIN_SITE variants and re-index
         grep -v NEGATIVE_TRAIN_SITE sub.vcf > subsampled.vcf
         java -jar ${gatk} IndexFeatureFile -F subsampled.vcf
    }

    output {
        File output_vcf = "subsampled.vcf"
        File output_vcf_idx = "subsampled.vcf.idx"
    }
}

task RemoveNearbyIndels {
    File gatk
    File input_vcf
    File input_vcf_idx
    Int min_indel_spacing
    String name

    command {
        java -jar ${gatk} RemoveNearbyIndels -V ${input_vcf} -O ${name}.vcf -minIndelSpacing ${min_indel_spacing}
    }

    output {
        File output_vcf = "${name}.vcf"
        File output_vcf_idx = "${name}.vcf.idx"
    }
}

