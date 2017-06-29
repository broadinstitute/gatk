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

# Here we implement steps 3 and 4 for several replicate bams of a *single* plex, e.g. 4 different 10-plex bams.  Note that each
# replicate has its own truth vcf because although the truth variants are identical, the bam-derived annotations may
# differ slightly.

#Inputs are identical to mutect2_multi_sample_concordance.wdl except for the addition of a python script for generating plots

import "mutect2_multi_sample_concordance.wdl" as concordance


task ConvertToTable {
  String gatk4_jar
  File? gatk4_jar_override
  File input_vcf
  File input_vcf_idx

  command {
  # Use GATK Jar override if specified
            GATK_JAR=${gatk4_jar}
            if [[ "${gatk4_jar_override}" == *.jar ]]; then
                GATK_JAR=${gatk4_jar_override}
            fi

        java -jar $GATK_JAR VariantsToTable -V ${input_vcf} -F STATUS -F BAM_DEPTH -F AF_EXP -F TYPE -O "result.table"
  }

  runtime {
          memory: "5 GB"
      }

  output {
    File table = "result.table"
  }
}

task AnalyzeSensitivity {
    File input_table
    File python_sensitivity_script
    String prefix

    command {
        python ${python_sensitivity_script} ${input_table}
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

    runtime {
            memory: "5 GB"
        }

    output {
        File table = "${prefix}_combined.txt"
    }
}

task Jaccard {
    String gatk4_jar
    File? gatk4_jar_override
    Array[File] tpfp_vcfs
    Array[File] tpfp_vcfs_idx
    String prefix

    #type is SNP or INDEL
    String type

    command {

        # Use GATK Jar override if specified
        GATK_JAR=${gatk4_jar}
        if [[ "${gatk4_jar_override}" == *.jar ]]; then
          GATK_JAR=${gatk4_jar_override}
        fi

        result="${prefix}_${type}_jaccard.txt"
        touch $result

        count=0
        for vcf in ${sep = ' ' tpfp_vcfs}; do
            ((count++))
            java -jar $GATK_JAR SelectVariants -V $vcf -selectType ${type} -O ${type}_only_$count.vcf
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
                java -jar $GATK_JAR SelectVariants -V $file1 --concordance $file2 -O overlap.vcf
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

    runtime {
        memory: "5 GB"
    }

    output {
        File table = "${prefix}_${type}_jaccard.txt"
    }
}

workflow HapmapSensitivity {
   # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
  	String gatk4_jar
  	Int scatter_count
  	File bam_list
  	File ref_fasta
  	File ref_fasta_index
  	File ref_dict
  	File? intervals
  	File? pon
  	File? pon_index
  	Boolean is_run_orientation_bias_filter
    String m2_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    Array[String] artifact_modes
    File picard_jar
    File truth_list
    String? m2_extra_args
    String? m2_extra_filtering_args
    #a prefix string like "5plex"
    String prefix

    File python_sensitivity_script

    Array[Array[String]] truth_files = read_tsv(truth_list)

  call concordance.Mutect2_Multi_Concordance {
    input:
     # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
    	gatk4_jar = gatk4_jar,
    	scatter_count = scatter_count,
    	pair_list = bam_list,
    	intervals = intervals,
    	ref_fasta = ref_fasta,
    	ref_fasta_index = ref_fasta_index,
    	ref_dict = ref_dict,
    	pon = pon,
    	pon_index = pon_index,
    	is_run_orientation_bias_filter = is_run_orientation_bias_filter,
        m2_docker = m2_docker,
        gatk4_jar_override = gatk4_jar_override,
        preemptible_attempts = preemptible_attempts,
        artifact_modes = artifact_modes,
        picard_jar = picard_jar,
        truth_list = truth_list,
        m2_extra_args = m2_extra_args,
        m2_extra_filtering_args = m2_extra_filtering_args
  }

  scatter(n in range(length(read_tsv(bam_list)))) {
    call ConvertToTable {
        input:
              gatk4_jar = gatk4_jar,
              gatk4_jar_override = gatk4_jar_override,
              input_vcf = Mutect2_Multi_Concordance.tpfn[n],
              input_vcf_idx = Mutect2_Multi_Concordance.tpfn_idx[n]
      }
  }

  call CombineTables as SensitivityTables { input: input_tables = ConvertToTable.table, prefix = prefix }

  call AnalyzeSensitivity {
             input:
              input_table = SensitivityTables.table,
              python_sensitivity_script = python_sensitivity_script,
              prefix = prefix
  }

  call CombineTables as SummaryTables { input: input_tables = Mutect2_Multi_Concordance.summary, prefix = prefix }

  call Jaccard as JaccardSNP {
    input:
        gatk4_jar = gatk4_jar,
        gatk4_jar_override = gatk4_jar_override,
        tpfp_vcfs = Mutect2_Multi_Concordance.tpfp,
        tpfp_vcfs_idx = Mutect2_Multi_Concordance.tpfp_idx,
        prefix = prefix,
        type = "SNP"
  }

  call Jaccard as JaccardINDEL {
      input:
          gatk4_jar = gatk4_jar,
          gatk4_jar_override = gatk4_jar_override,
          tpfp_vcfs = Mutect2_Multi_Concordance.tpfp,
          tpfp_vcfs_idx = Mutect2_Multi_Concordance.tpfp_idx,
          prefix = prefix,
          type = "INDEL"
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
    Array[File] tpfn = Mutect2_Multi_Concordance.tpfn
    Array[File] tpfn_idx = Mutect2_Multi_Concordance.tpfn_idx
    Array[File] ftnfn = Mutect2_Multi_Concordance.ftnfn
    Array[File] ftnfn_idx = Mutect2_Multi_Concordance.ftnfn_idx
  }
}