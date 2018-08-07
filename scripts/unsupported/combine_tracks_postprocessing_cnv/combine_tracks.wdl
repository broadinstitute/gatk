
# A postprocessing workflow for evaluating CNV
workflow CombineTracksWorkflow {
	File tumor_called_seg
	File matched_normal_called_seg
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File centromere_tracks_seg
    File gistic_blacklist_tracks_seg
    File? gatk4_jar_override
    Array[String] columns_of_interest
    Int? germline_tagging_padding
    String group_id
    String gatk_docker
    call CombineTracks {
        input:
            tumor_called_seg = tumor_called_seg,
            matched_normal_called_seg = matched_normal_called_seg,
            centromere_tracks_seg = centromere_tracks_seg,
            columns_of_interest = columns_of_interest,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            germline_tagging_padding = germline_tagging_padding,
            gistic_blacklist_tracks_seg = gistic_blacklist_tracks_seg,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker
    }

    call IGVConvert as IGVConvertNormal {
        input:
            COMMENTCHAR="@",
            INPUT=matched_normal_called_seg,
            VALUE=basename(matched_normal_called_seg),
            FIELD="SAMPLE",
            OUTPUT = basename(matched_normal_called_seg) + ".igv.seg",
            PRE_POST="PRE"
    }
    call IGVConvert as IGVConvertTumor {
        input:
            COMMENTCHAR="@",
            INPUT=tumor_called_seg,
            VALUE=basename(tumor_called_seg),
            FIELD="SAMPLE",
            OUTPUT = basename(tumor_called_seg) + ".igv.seg",
            PRE_POST="PRE"
    }
    call IGVConvert as IGVConvertTumorOutput {
        input:
            COMMENTCHAR="@",
            INPUT=CombineTracks.germline_tagged_with_tracks_seg,
            VALUE=basename(tumor_called_seg),
            FIELD="SAMPLE",
            OUTPUT = basename(tumor_called_seg) + ".tagged.igv.seg",
            PRE_POST="PRE"
    }
    call PruneGermlineTagged {
        input:
            germline_tagged_seg = IGVConvertTumorOutput.outFile
    }

    call MergeSegmentByAnnotation {
        input:
            seg_file = PruneGermlineTagged.tumor_with_germline_pruned_seg,
            annotations = ["MEAN_LOG2_COPY_RATIO"],
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk_docker = gatk_docker,
            gatk4_jar_override = gatk4_jar_override
    }

    output {
        File cnv_postprocessing_tumor_igv_compat = IGVConvertTumor.outFile
        File cnv_postprocessing_normal_igv_compat = IGVConvertNormal.outFile
        File cnv_postprocessing_tumor_with_tracks_pruned_seg = PruneGermlineTagged.tumor_with_germline_pruned_seg
        File cnv_postprocessing_tumor_with_tracks_pruned_merged_seg = MergeSegmentByAnnotation.cnv_merged_seg
        File cnv_postprocessing_tumor_with_tracks_tagged_seg = CombineTracks.germline_tagged_with_tracks_seg
    }

}

task CombineTracks {
	File tumor_called_seg
	File matched_normal_called_seg
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File centromere_tracks_seg
    File gistic_blacklist_tracks_seg
    Array[String] columns_of_interest

	File? gatk4_jar_override
	
	String output_name = basename(tumor_called_seg)

    Int? germline_tagging_padding

    # Runtime parameters
    Int? mem_gb
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu 
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
	Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

	command <<<
	    
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

		echo "======= Germline Tagging "
		gatk --java-options "-Xmx${command_mem}m" \
			TagGermlineEvents \
            --segments ${tumor_called_seg} --called-matched-normal-seg-file ${matched_normal_called_seg} \
            -O ${output_name}.germline_tagged.seg -R ${ref_fasta} \
            ${"--endpoint-padding " + germline_tagging_padding}

        echo "======= Centromeres "
    	gatk --java-options "-Xmx${command_mem}m" \
    	     CombineSegmentBreakpoints \
            --segments ${output_name}.germline_tagged.seg --segments ${centromere_tracks_seg}  \
            --columns-of-interest ${sep=" --columns-of-interest " columns_of_interest} \
            -O ${output_name}.centro.seg -R ${ref_fasta}

        echo "======= GISTIC blacklist "
    	gatk --java-options "-Xmx${command_mem}m" \
    	     CombineSegmentBreakpoints \
            --segments ${output_name}.centro.seg --segments ${gistic_blacklist_tracks_seg}  \
            --columns-of-interest ${sep=" --columns-of-interest " columns_of_interest} \
            --columns-of-interest ID \
            -O ${output_name}.final.seg -R ${ref_fasta}
	>>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        # Note that the space before SSD and HDD should be included.
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: 20
    }

    output {
        File germline_tagged_with_tracks_seg = "${output_name}.final.seg"
        File germline_tagged_with_centro_track_seg = "${output_name}.centro.seg"
        File germline_tagged_seg = "${output_name}.germline_tagged.seg"
    }
}

task PruneGermlineTagged {
    File germline_tagged_seg

    String output_filename = basename(germline_tagged_seg) + ".pruned.seg"
    # Ugh ... this command that I am writing here is definitely not ready for prime time
    command <<<
    set -e
        python <<EOF
import pandas
import os.path
tumor_tagged = "${germline_tagged_seg}"

tumor_tagged_df = pandas.read_csv(tumor_tagged, delimiter="\t")
tumor_tagged_pruned_df = tumor_tagged_df[(tumor_tagged_df["POSSIBLE_GERMLINE"] == "0") & (tumor_tagged_df["type"] != "centromere") & (tumor_tagged_df["ID"].isna())]
output_filename = "${output_filename}"
print(output_filename)
tumor_tagged_pruned_df.to_csv(output_filename, sep="\t", index=False)

EOF

    >>>

    runtime {
        docker: "amancevice/pandas"
        memory: "2000 MB"
        disks: "local-disk 100 HDD"
        preemptible: 3
        cpu: 1
    }

    output {
        File tumor_with_germline_pruned_seg = "${output_filename}"
    }
}

task IGVConvert {

    #Inputs and constants defined here

    File  INPUT
    String FIELD
    String VALUE
    String OUTPUT
    String? PRE_POST
    String? COMMENTCHAR

    String PRE=select_first([PRE_POST, "POST"])
    String COMMENT=select_first([COMMENTCHAR, "#"])

    # Runtime parameters
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_gb
    Boolean use_ssd = false

    command <<<
# Modified from original task written by Chip Stewart
grep -v "^"${COMMENT} ${INPUT} > tmp.tsv
head -1 tmp.tsv > tmp.header.txt
cat tmp.header.txt | while IFS=$'\n\r' read -r line
do
		echo ${FIELD} > tmp.x.tsv
done
sed 1,1d tmp.tsv > tmp.rest.txt



cat tmp.rest.txt | while IFS=$'\n\r' read -r line
do
   	echo ${VALUE} >> tmp.x.tsv
done

if [ ${PRE} = "PRE" ]; then
    echo expression evaluated as true
    paste tmp.x.tsv tmp.tsv  > ${OUTPUT}.tmp
else
    echo expression evaluated as false
    paste tmp.tsv tmp.x.tsv > ${OUTPUT}.tmp
fi
head -1 ${OUTPUT}.tmp > tmp_header2.txt

tr "\t" "\n" < tmp_header2.txt | grep -n MEAN_LOG2_COPY_RATIO | cut -f1 -d: > tmp_col_num
COL_NUM=`cat tmp_col_num`
echo $COL_NUM

cut -f$COL_NUM  ${OUTPUT}.tmp > col_data
cut --complement -f $COL_NUM ${OUTPUT}.tmp > ${OUTPUT}.tmp_header2
paste ${OUTPUT}.tmp_header2 col_data >${OUTPUT}
    >>>

    output {
        File outFile="${OUTPUT}"
    }

    runtime {
        docker : "ubuntu:16.04"
        memory: select_first([mem_gb, 2000]) + " GB"
        # Note that the space before SSD and HDD should be included.
        disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: select_first([boot_disk_gb, 20])
    }
}

task MergeSegmentByAnnotation {
	File seg_file
	Array[String] annotations
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai

	File? gatk4_jar_override

	String output_name = basename(seg_file)

    Int? max_merge_distance

    # Runtime parameters
    Int? mem_gb
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
	Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

	command <<<

        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

		echo "======= Merging "
		gatk --java-options "-Xmx${command_mem}m" \
			MergeAnnotatedRegionsByAnnotation \
            --segments ${seg_file} \
            ${"--max-merge-distance " + max_merge_distance} \
            --annotations-to-match ${sep=" --annotations-to-match " annotations} \
            -O ${output_name}.merged.seg -R ${ref_fasta}

	>>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        # Note that the space before SSD and HDD should be included.
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: 20
    }

    output {
        File cnv_merged_seg = "${output_name}.merged.seg"
    }
}