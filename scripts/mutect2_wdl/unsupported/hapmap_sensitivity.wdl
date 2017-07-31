### This is the concordance part of the CRSP sensitivity validation

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

# Here we implement these steps, except for the subsampling in step 1, for several replicate bams of a *single* plex,
# e.g. 4 different 10-plex bams.  Note that each replicate has its own truth vcf because although the truth variants are identical,
# the bam-derived annotations differ slightly.

import "mutect2.wdl" as MutectSingleSample

workflow HapmapSensitivity {
    File gatk
  	Int scatter_count
  	File bam_list
  	Array[Array[String]] replicates = read_tsv(bam_list)
  	File ref_fasta
  	File ref_fasta_index
  	File ref_dict
  	File? intervals
  	File? pon
  	File? pon_index
  	Boolean is_run_orientation_bias_filter
    Array[String] artifact_modes
    File picard_jar
    String? m2_extra_args
    String? m2_extra_filtering_args
    String prefix   #a prefix string like "5plex"
    File python_script
    Int max_depth
    File preprocessed_hapmap
    File preprocessed_hapmap_idx

    call RestrictIntervals {
        input: gatk = gatk, vcf = preprocessed_hapmap, vcf_idx = preprocessed_hapmap_idx, intervals = intervals
    }

    scatter (row in replicates) {
        File bam = row[0]
        File index = row[1]
        String sample = row[2]

        call MixingFractions {
            input: gatk = gatk, vcf = RestrictIntervals.output_vcf, vcf_idx = RestrictIntervals.output_vcf_idx, bam = bam, bam_idx = index
        }

        call ExpectedAlleleFraction {
            input: gatk = gatk, vcf = RestrictIntervals.output_vcf, vcf_idx = RestrictIntervals.output_vcf_idx, mixing_fractions = MixingFractions.mixing
        }

        call BamDepth {
            input: gatk = gatk, vcf = ExpectedAlleleFraction.output_vcf, vcf_idx = ExpectedAlleleFraction.output_vcf_idx,
                bam = bam, bam_idx = index, max_depth = max_depth
        }

        call MutectSingleSample.Mutect2 {
            input:
                gatk4_jar = "OVERRIDDEN",
                scatter_count = scatter_count,
                tumor_bam = bam,
                tumor_bam_index = index,
                tumor_sample_name = sample,
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                pon = pon,
                pon_index = pon_index,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator = false,
		m2_docker = "ubuntu:16.04",
                oncotator_docker = "ubuntu:16.04",
		gatk4_jar_override = gatk,
                preemptible_attempts = 2,
                artifact_modes = artifact_modes,
                picard_jar = picard_jar,
                m2_extra_args = m2_extra_args,
                m2_extra_filtering_args = m2_extra_filtering_args
        }

        call Concordance {
            input: gatk = gatk, intervals = intervals,
                  truth = BamDepth.output_vcf,
                  truth_idx = BamDepth.output_vcf_idx,
                  eval = Mutect2.filtered_vcf,
                  eval_idx = Mutect2.filtered_vcf_index
        }

        call ConvertToTable {
            input: gatk = gatk, input_vcf = Concordance.tpfn, input_vcf_idx = Concordance.tpfn_idx
        }
    } #done with scatter over replicates

    call CombineTables as SensitivityTables {
        input: input_tables = ConvertToTable.table, prefix = prefix
    }

    call AnalyzeSensitivity {
        input: input_table = SensitivityTables.table, python_script = python_script, prefix = prefix
    }

    call CombineTables as SummaryTables {
        input: input_tables = Concordance.summary, prefix = prefix
    }

    call Jaccard as JaccardSNP {
        input: gatk = gatk, calls = Mutect2.filtered_vcf, calls_idx = Mutect2.filtered_vcf_index, prefix = prefix, type = "SNP"
    }

    call Jaccard as JaccardINDEL {
        input: gatk = gatk, calls = Mutect2.filtered_vcf, calls_idx = Mutect2.filtered_vcf_index, prefix = prefix, type = "INDEL"
    }

    output {
        File snp_table = AnalyzeSensitivity.snp_table
        File snp_plot = AnalyzeSensitivity.snp_plot
        File indel_table = AnalyzeSensitivity.indel_table
        File indel_plot = AnalyzeSensitivity.indel_plot
        File summary = SummaryTables.table
        File raw_table = SensitivityTables.table
        File snp_jaccard_table = JaccardSNP.table
        File indel_jaccard_table = JaccardINDEL.table
        Array[File] tpfn = Concordance.tpfn
        Array[File] tpfn_idx = Concordance.tpfn_idx
        Array[File] ftnfn = Concordance.ftnfn
        Array[File] ftnfn_idx = Concordance.ftnfn_idx
    }
} #end of workflow

#### Tasks for making truth
task RestrictIntervals {
    File gatk
    File vcf
    File vcf_idx
    File? intervals

    command {
        # subsampling and restriction to biallelics and intervals
        java -jar ${gatk} SelectVariants -V ${vcf} -O restricted.vcf \
            ${"-L " + intervals} \
    }

    output {
        File output_vcf = "restricted.vcf"
        File output_vcf_idx = "restricted.vcf.idx"
    }
}

task BamDepth {
    File gatk
    File vcf
    File vcf_idx
    File bam
    File bam_idx
    Int  max_depth   #ignore sites with depth greater than this because they are alignment artifacts

    command {
        java -jar ${gatk} AnnotateVcfWithBamDepth -V ${vcf} -I ${bam} -O "bam_depth.vcf"
        java -jar ${gatk} SelectVariants -V bam_depth.vcf --select "BAM_DEPTH < ${max_depth}" -O truth.vcf
    }

    output {
        File output_vcf = "truth.vcf"
        File output_vcf_idx = "truth.vcf.idx"
    }
}

task MixingFractions {
    File gatk
    File vcf
    File vcf_idx
    File bam
    File bam_idx

    command {
        java -jar ${gatk} CalculateMixingFractions -V ${vcf} -I ${bam} -O "mixing.table"
    }

    output { File mixing = "mixing.table" }
}

task ExpectedAlleleFraction {
    File gatk
    File vcf
    File vcf_idx
    File mixing_fractions

    command {
        java -jar ${gatk} AnnotateVcfWithExpectedAlleleFraction -V ${vcf} -O af_exp.vcf --mixingFractions  ${mixing_fractions}
    }

    output {
        File output_vcf = "af_exp.vcf"
        File output_vcf_idx = "af_exp.vcf.idx"
    }
}

### Tasks for analysing sensitivity

task ConvertToTable {
  File gatk
  File input_vcf
  File input_vcf_idx

  command {
      java -jar ${gatk} VariantsToTable -V ${input_vcf} -F STATUS -F BAM_DEPTH -F AF_EXP -F TYPE -O "result.table"
  }

  runtime { memory: "5 GB" }

  output { File table = "result.table" }
}

task CombineTables {
    Array[File] input_tables
    String prefix

    command {
        for file in ${sep=' ' input_tables}; do
           head -n 1 $file > header
           tail -n +2 $file >> body
        done

        cat header body > "${prefix}_combined.txt"
    }

    runtime { memory: "5 GB" }

    output { File table = "${prefix}_combined.txt" }
}

task AnalyzeSensitivity {
    File input_table
    File python_script
    String prefix

    command {
        . /broad/software/scripts/useuse
        use Python-2.7
	python ${python_script} ${input_table}
        mv "SNP_sensitivity.tsv" "${prefix}_SNP_sensitivity.tsv"
        mv "SNP_sensitivity.png" "${prefix}_SNP_sensitivity.png"
        mv "Indel_sensitivity.tsv" "${prefix}_Indel_sensitivity.tsv"
        mv "Indel_sensitivity.png" "${prefix}_Indel_sensitivity.png"
    }

    runtime {
            continueOnReturnCode: [0,1]
            memory: "5 GB"
    }

    output {
        File snp_table = "${prefix}_SNP_sensitivity.tsv"
        File snp_plot = "${prefix}_SNP_sensitivity.png"
        File indel_table = "${prefix}_Indel_sensitivity.tsv"
        File indel_plot = "${prefix}_Indel_sensitivity.png"
    }
}

#Make Jaccard index table for SNVs or indels from an array of called vcfs
task Jaccard {
    File gatk
    Array[File] calls
    Array[File] calls_idx
    String prefix
    String type #SNP or INDEL

    command {
        result="${prefix}_${type}_jaccard.txt"
        touch $result

        count=0
        for vcf in ${sep = ' ' calls}; do
            ((count++))
            java -jar ${gatk} SelectVariants -V $vcf -selectType ${type} -O ${type}_only_$count.vcf
        done

        for file1 in ${type}_only*.vcf; do
            column=0
            for file2 in ${type}_only*.vcf; do
            ((column++))
            if [ $column != 1 ]; then
                printf "\t" >> $result
            fi

            if [ $file1 == $file2 ]; then
                printf 1.0000 >> $result
            else
                java -jar ${gatk} SelectVariants -V $file1 --concordance $file2 -O overlap.vcf
                overlap=`grep -v '#' overlap.vcf | wc -l`

                num1=`grep -v '#' $file1 | wc -l`
                num2=`grep -v '#' $file2 | wc -l`
                just1=$((num1 - overlap))
                just2=$((num2 - overlap))

                total=$((overlap + just1 + just2))
                jaccard=`echo "$overlap / $total" | bc -l`

                printf "%0.6f" $jaccard >> $result
            fi
            done
            printf "\n" >> $result
        done

    }

    runtime { memory: "5 GB" }

    output { File table = "${prefix}_${type}_jaccard.txt" }
}

task Concordance {
      File gatk
      File? intervals
      File truth
      File truth_idx
      File eval
      File eval_idx

      command {
          java -jar ${gatk} Concordance ${"-L " + intervals} \
            -truth ${truth} -eval ${eval} \
            -tpfn "tpfn.vcf" \
            -ftnfn "ftnfn.vcf" \
            -summary summary.tsv
      }

      runtime { memory: "5 GB" }

      output {
            File tpfn = "tpfn.vcf"
            File tpfn_idx = "tpfn.vcf.idx"
            File ftnfn = "ftnfn.vcf"
            File ftnfn_idx = "ftnfn.vcf.idx"
            File summary = "summary.tsv"
      }
}
