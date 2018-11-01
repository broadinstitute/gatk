# A postprocessing workflow for GATK CNV ModelSegments.
#
# THIS CANNOT BE RUN IN TUMOR-ONLY MODE.  A MATCHED NORMAL IS REQUIRED.
#
# Internally, we have shown that this workflow improves sensitivity and precision.  Additionally, this workflow is quite
#  cheap to run.  Note that auxiliary files are required and without these exact files, performance may suffer.
#
# UNSUPPORTED and we may modify this WDL before we fully support it.
#
# This will also generate ABSOLUTE (https://software.broadinstitute.org/cancer/cga/absolute) compatible files.
#  Currently, the balanced-segment calling is based on hard thresholds.
#
# This workflow:
#  - Does conversion of output seg files to a format compatible with IGV.  This will be removed soon, since the GATK tools
#     have been updated to produce files that can be easily read by IGV.
#  - Tags possible germline events that appear in the tumor results.  This is done by comparing breakpoints and
#     reciprocal overlaps.  This algorithm takes into account any hypersegmentation.
#  - Tags tumor events that appear in the given blacklists (centromere and GISTIC blacklist)
#  - Removes any events that have been tagged in the previous two steps.  Any resulting gaps are imputed, so long as
#     the segment mean would remain constant and the distance is less than max_merge_distance (default 1Mb as of
#     this writing)
#
#  Outputs:
#  - seg file that is has germline events removed (and gaps merged).  Same columns as a called seg file.
#      This should be considered the final output.  IGV compatible.
#  - seg file that can be used as input to ABSOLUTE.  This uses an identical format as AllelicCapSeg.  This is NOT compatible with IGV.
#      Please note that the balanced-segment calling in this file has not been evaluated heavily.
#  - seg file that can be used as input to GISTIC2.  Note that since this workflow is meant for a single pair, users will need to aggregate the GISTIC2 seg files from this workflow for all pairs.
#  - skew file for the skew parameter in ABSOLUTE
#  - skew as a float for the skew parameter in ABSOLUTE
#  - Intermediate files for browsing various steps of the workflow.
#
# Important notes:
#  - The probe and het counts become inaccurate in the results.  These need to be corrected, but that is not trivial.
#  - This workflow has not been tested with hg38
#  - Evaluation of germline tagging and blacklists on hg38 is still pending.
#  - Evaluation of the conversion to ACS format for ABSOLUTE is still pending.
#  - Performance increases (both sensitivity and precision) over this workflow when using a blacklist during PoN creation and running of case samples.
#
#  A blacklist for that can be found at:
#   - hg19: gs://gatk-best-practices/somatic-b37/CNV_and_centromere_blacklist.hg19.list
#   - hg38: gs://gatk-best-practices/somatic-hg38/CNV_and_centromere_blacklist.hg38liftover.seg
#
# Do not attempt to use the above blacklists for this workflow.  See below:
#
# Auxiliary files for this workflow (hg19):
#  - centromere_tracks_seg: gs://gatk-best-practices/somatic-b37/final_centromere_hg19.seg
#  - gistic_blacklist_tracks_seg:  gs://gatk-best-practices/somatic-b37/CNV.hg19.bypos.v1.CR1_event_added.mod.seg
#
# Auxiliary files for this workflow (hg38) -- these are untested and the gistic list is a liftover:
#  - centromere_tracks_seg:  gs://gatk-best-practices/somatic-hg38/final_centromere_hg38.seg
#  - gistic_blacklist_tracks_seg:  gs://gatk-best-practices/somatic-hg38/CNV.hg38liftover.bypos.v1.CR1_event_added.mod.seg
#
#  Dev note:  FYI...This workflow is easily tested on a laptop.
#
workflow CombineTracksWorkflow {
    File tumor_called_seg
    File tumor_modeled_seg
    File af_param
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
    Int? max_merge_distance
    Array[String]? annotations_on_which_to_merge
    Int? min_hets_acs_results
    Float? maf90_threshold

    Array[String]? annotations_on_which_to_merge_final = select_first([annotations_on_which_to_merge,
    ["MEAN_LOG2_COPY_RATIO", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_90",
        "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_50", "MINOR_ALLELE_FRACTION_POSTERIOR_90"]])

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
            PRE_POST="PRE",
            SEGMENT_MEAN_COL = "MEAN_LOG2_COPY_RATIO"
    }

    call IGVConvert as IGVConvertTumor {
        input:
            COMMENTCHAR="@",
            INPUT=tumor_called_seg,
            VALUE=basename(tumor_called_seg),
            FIELD="SAMPLE",
            OUTPUT = basename(tumor_called_seg) + ".igv.seg",
            PRE_POST="PRE",
            SEGMENT_MEAN_COL = "MEAN_LOG2_COPY_RATIO"
    }

    call PrepareForACSConversion {
        input:
            called_seg = CombineTracks.germline_tagged_with_tracks_seg,
            modeled_seg = tumor_modeled_seg,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker
    }

    call IGVConvert as IGVConvertTumorOutput {
        input:
            COMMENTCHAR="@",
            INPUT=PrepareForACSConversion.model_and_calls_merged_gatk_seg,
            VALUE=basename(PrepareForACSConversion.model_and_calls_merged_gatk_seg),
            FIELD="SAMPLE",
            OUTPUT = basename(PrepareForACSConversion.model_and_calls_merged_gatk_seg) + ".tagged.igv.seg",
            PRE_POST="PRE",
            SEGMENT_MEAN_COL = "MEAN_LOG2_COPY_RATIO"
    }

    call FilterGermlineTagged {
        input:
            germline_tagged_seg = PrepareForACSConversion.model_and_calls_merged_gatk_seg,
            docker = gatk_docker
    }

    call MergeSegmentByAnnotation {
        input:
            seg_file = FilterGermlineTagged.tumor_with_germline_filtered_seg,
            annotations = annotations_on_which_to_merge_final,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk_docker = gatk_docker,
            gatk4_jar_override = gatk4_jar_override,
            max_merge_distance = max_merge_distance
    }

    call PrototypeACSConversion {
        input:
            model_seg = MergeSegmentByAnnotation.cnv_merged_seg,
            af_param = af_param,
            docker = gatk_docker,
            min_hets_acs_results = min_hets_acs_results,
            maf90_threshold = maf90_threshold
    }

    call IGVConvert as IGVConvertMergedTumorOutput {
        input:
            COMMENTCHAR="@",
            INPUT=MergeSegmentByAnnotation.cnv_merged_seg,
            VALUE=basename(tumor_called_seg),
            FIELD="SAMPLE",
            OUTPUT = basename(tumor_called_seg) + ".pruned_merged.igv.seg",
            PRE_POST="PRE",
            SEGMENT_MEAN_COL = "MEAN_LOG2_COPY_RATIO"
    }

    call Gistic2Convert {
        input:
            input_file = IGVConvertMergedTumorOutput.outFile,
            docker = gatk_docker
    }


    output {
        File cnv_postprocessing_tumor_igv_compat = IGVConvertTumor.outFile
        File cnv_postprocessing_normal_igv_compat = IGVConvertNormal.outFile
        File cnv_postprocessing_tumor_with_tracks_filtered_seg = FilterGermlineTagged.tumor_with_germline_filtered_seg
        File cnv_postprocessing_tumor_with_tracks_filtered_merged_seg = IGVConvertMergedTumorOutput.outFile
        File cnv_postprocessing_tumor_with_tracks_filtered_merged_seg_ms_format = MergeSegmentByAnnotation.cnv_merged_seg
        File cnv_postprocessing_tumor_with_tracks_tagged_seg = CombineTracks.germline_tagged_with_tracks_seg
        File cnv_postprocessing_tumor_acs_seg = PrototypeACSConversion.cnv_acs_conversion_seg
        File cnv_postprocessing_tumor_acs_skew = PrototypeACSConversion.cnv_acs_conversion_skew
        Float cnv_postprocessing_tumor_acs_skew_float = PrototypeACSConversion.cnv_acs_conversion_skew_float
        String cnv_postprocessing_tumor_acs_skew_string = PrototypeACSConversion.cnv_acs_conversion_skew_string
        File cnv_postprocessing_tumor_with_tracks_filtered_merged_seg_gistic2 = Gistic2Convert.output_file_gistic2
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
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
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

task FilterGermlineTagged {
    File germline_tagged_seg
    String docker

    String output_filename = basename(germline_tagged_seg) + ".pruned.seg"
    command <<<
    set -e
        python <<EOF
import pandas
import os.path
tumor_tagged = "${germline_tagged_seg}"

tumor_tagged_df = pandas.read_csv(tumor_tagged, delimiter='\t', comment="@")
tumor_tagged_df["POSSIBLE_GERMLINE"] = tumor_tagged_df["POSSIBLE_GERMLINE"].astype('str')
tumor_tagged_pruned_df = tumor_tagged_df[((tumor_tagged_df["POSSIBLE_GERMLINE"] == "0.0") | (tumor_tagged_df["POSSIBLE_GERMLINE"] == "0") ) & (tumor_tagged_df["type"] != "centromere") & (tumor_tagged_df["ID"].isna())]

output_filename = "${output_filename}"
print(output_filename)
tumor_tagged_pruned_df.to_csv(output_filename, sep="\t", index=False)

EOF

    >>>

    runtime {
        docker: docker
        memory: "2000 MB"
        disks: "local-disk 100 HDD"
        preemptible: 3
        cpu: 1
    }

    output {
        File tumor_with_germline_filtered_seg = "${output_filename}"
    }
}

task IGVConvert {

    #Inputs and constants defined here

    File  INPUT
    String FIELD
    String VALUE
    String OUTPUT
    String SEGMENT_MEAN_COL
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

tr "\t" "\n" < tmp_header2.txt | grep -n ${SEGMENT_MEAN_COL} | cut -f1 -d: > tmp_col_num
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
        memory: select_first([mem_gb, 2]) + " GB"
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
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
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

# TODO: No non-trivial heredocs in WDL.  Add this to the script directory and call via anaconda (future release)
# TODO: This is a hard threholding algorithm.  Better approaches exist if this does not meet needs.
# TODO: The min hets hard thresholding should eventually be removed.  This is to mitigate some hypersegmentation when rerunning GATK CNV is too expensive.
task PrototypeACSConversion {
    File model_seg
    File af_param
    String docker
    Float? maf90_threshold
    String output_filename = basename(model_seg) + ".acs.seg"
    String output_skew_filename = output_filename + ".skew"
    Int? min_hets_acs_results
    Int min_hets_acs_results_final = select_first([min_hets_acs_results, 0])

    command <<<
        set -e
        python <<EOF
import sys
import re
import pandas as pd
import numpy as np
from collections import defaultdict
import scipy
from scipy import special as sp
import os.path

model_segments_seg_input_file = "${model_seg}"
model_segments_af_param_input_file = "${af_param}"
alleliccapseg_seg_output_file = "${output_filename}"
alleliccapseg_skew_output_file = "${output_skew_filename}"

HAM_FIST_THRESHOLD=${default="0.47" maf90_threshold}

# regular expression for matching sample name from header comment line
sample_name_header_regexp = "^@RG.*SM:(.*)[\t]*.*$"

#define AllelicCapSeg columns
alleliccapseg_seg_columns = [
    'Chromosome',
    'Start.bp',
    'End.bp',
    'n_probes',
    'length',
    'n_hets',
    'f',
    'tau',
    'sigma.tau',
    'mu.minor',
    'sigma.minor',
    'mu.major',
    'sigma.major',
    'SegLabelCNLOH']

def read_sample_name(input_file, max_scan_lines=10000):
    with open(input_file, 'r') as f:
        for _ in range(max_scan_lines):
            line = f.readline()
            match = re.search(sample_name_header_regexp, line, re.M)
            if match is None:
                continue
            groups = match.groups()
            return groups[0]
    raise Exception("Sample name could not be found in \"{0}\"".format(input_file))

#read GATK ModelSegments files and perform some basic checks
model_segments_seg_pd = pd.read_csv(model_segments_seg_input_file,
                                    sep='\t', comment='@', na_values='NA')
model_segments_af_param_pd = pd.read_csv(model_segments_af_param_input_file, sep='\t', comment='@')

def simple_determine_allelic_fraction(model_segments_seg_pd):
    result = model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_50']
    result[model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'] > HAM_FIST_THRESHOLD] = 0.5
    return result

def convert_model_segments_to_alleliccapseg(model_segments_seg_pd,
                                            model_segments_af_param_pd):
    alleliccapseg_seg_pd = pd.DataFrame(columns=alleliccapseg_seg_columns)

    #The following conversions are trivial.
    alleliccapseg_seg_pd['Chromosome'] = model_segments_seg_pd['CONTIG']
    alleliccapseg_seg_pd['Start.bp'] = model_segments_seg_pd['START']
    alleliccapseg_seg_pd['End.bp'] = model_segments_seg_pd['END']
    alleliccapseg_seg_pd['n_probes'] = model_segments_seg_pd['NUM_POINTS_COPY_RATIO_1']
    alleliccapseg_seg_pd['length'] = alleliccapseg_seg_pd['End.bp'] - alleliccapseg_seg_pd['Start.bp']
    alleliccapseg_seg_pd['n_hets'] = model_segments_seg_pd['NUM_POINTS_ALLELE_FRACTION']

    #ModelSegments estimates posterior credible intervals, while AllelicCapSeg performs maximum a posteriori (MAP) estimation.
    #The copy-ratio and allele-fraction models fit by both also differ.
    # We will attempt a rough translation of the model fits here.

    alleliccapseg_seg_pd['f'] = simple_determine_allelic_fraction(model_segments_seg_pd)

    alleliccapseg_seg_pd['tau'] = 2. * 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_50']
    alleliccapseg_seg_pd['sigma.tau'] = 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_90'] - 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_10']
    sigma_f = (model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'].values - model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_10'].values) / 2.
    sigma_mu = np.sqrt(sigma_f**2 + alleliccapseg_seg_pd['sigma.tau']**2) #we propagate errors in the products f * tau and (1 - f) * tau in the usual way
    alleliccapseg_seg_pd['mu.minor'] = alleliccapseg_seg_pd['f'] * alleliccapseg_seg_pd['tau']
    alleliccapseg_seg_pd['sigma.minor'] = sigma_mu
    alleliccapseg_seg_pd['mu.major'] = (1. - alleliccapseg_seg_pd['f']) * alleliccapseg_seg_pd['tau']
    alleliccapseg_seg_pd['sigma.major'] = sigma_mu

    #For whatever reason, AllelicCapSeg attempts to call CNLOH.  Documentation is spotty, but it seems like it attempts
    # to distinguish between three states ("0 is flanked on both sides, 1 is one side, 2 is no cn.loh").
    # Let's just set everything to 2 for now.
    # Hopefully, ABSOLUTE is robust to this ...
    alleliccapseg_seg_pd['SegLabelCNLOH'] = 2

    #One important caveat: for segments with less than 10 hets, AllelicCapSeg also tries to call whether a segment is "split" or not.
    #  This script will attempt to call "split" on all segments.
    # ACS performs a simple hypothesis test on the alternate-allele fractions to see if
    # a unimodal distribution peaked at 0.5 is supported over a bimodal distribution peaked at f and 1 - f.
    # If the former is supported, then AllelicCapSeg ignores the MAP estimate of f and simply sets it to be 0.5.
    # ABSOLUTE may actually be rather sensitive to this.  Again, let's ignore for now, and we can later port this
    # statistical test if necessary.

    #Finally, I believe that ABSOLUTE requires the value of the "skew" parameter from the AllelicCapSeg
    #allele-fraction model.  This parameter is supposed to allow the model to account for reference bias,
    #  but the model likelihood that AllelicCapSeg uses is not valid over the entire range of the skew parameter.
    # We corrected this during the development of AllelicCNV and retain the same corrected model in ModelSegments.
    # We will try to transform the relevant parameter in the corrected model back to a "skew",
    # but this operation is ill defined.  Luckily, for WGS, the reference bias is typically negligible.
    model_segments_reference_bias = model_segments_af_param_pd[
        model_segments_af_param_pd['PARAMETER_NAME'] == 'MEAN_BIAS']['POSTERIOR_50']
    alleliccapseg_skew = 2. / (1. + model_segments_reference_bias)

    # If a row has less than X (set by user) hets, then assume zero
    filter_rows = alleliccapseg_seg_pd['n_hets'] < ${min_hets_acs_results_final}
    # mu.minor  sigma.minor  mu.major  sigma.major
    alleliccapseg_seg_pd.ix[filter_rows, 'n_hets'] = 0
    alleliccapseg_seg_pd.ix[filter_rows, 'f'] = np.NaN
    alleliccapseg_seg_pd.ix[filter_rows, 'mu.minor'] = np.NaN
    alleliccapseg_seg_pd.ix[filter_rows, 'sigma.minor'] = np.NaN
    alleliccapseg_seg_pd.ix[filter_rows, 'mu.major'] = np.NaN
    alleliccapseg_seg_pd.ix[filter_rows, 'sigma.major'] = np.NaN


    return alleliccapseg_seg_pd, alleliccapseg_skew


#do the conversion
alleliccapseg_seg_pd, alleliccapseg_skew = convert_model_segments_to_alleliccapseg(model_segments_seg_pd,
                                                                                   model_segments_af_param_pd)

#write the results
alleliccapseg_seg_pd.to_csv(alleliccapseg_seg_output_file, sep='\t', index=False, na_rep='NaN')
np.savetxt(alleliccapseg_skew_output_file, alleliccapseg_skew)

EOF
    >>>

    runtime {
        docker: docker
        memory: "2000 MB"
        disks: "local-disk 100 HDD"
        preemptible: 3
        cpu: 1
    }

    output {
        File cnv_acs_conversion_seg = "${output_filename}"
        File cnv_acs_conversion_skew = "${output_skew_filename}"
        Float cnv_acs_conversion_skew_float = read_float(output_skew_filename)
        String cnv_acs_conversion_skew_string = read_string(output_skew_filename)
    }
}

# Merge the model segments and the called file.
task PrepareForACSConversion {
	File called_seg
	File modeled_seg
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai

    Array[String] columns_of_interest = [
     "NUM_POINTS_COPY_RATIO", "NUM_POINTS_ALLELE_FRACTION",
     "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_90",
     "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_50", "MINOR_ALLELE_FRACTION_POSTERIOR_90",
     "CALL", "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO", "POSSIBLE_GERMLINE", "type", "ID"
    ]

	File? gatk4_jar_override

	String output_name = basename(modeled_seg)

    # Runtime parameters
    Int? mem_gb
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
	Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

	command <<<

        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

		echo "======= Merging GATK Model Seg and GATK Segment caller file "
    	gatk --java-options "-Xmx${command_mem}m" \
    	     CombineSegmentBreakpoints \
            --segments ${called_seg} --segments ${modeled_seg}  \
            --columns-of-interest ${sep=" --columns-of-interest " columns_of_interest} \
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
        File model_and_calls_merged_gatk_seg = "${output_name}.final.seg"
    }
}

# Must be python3 in the docker image
task Gistic2Convert {
    File input_file
    String docker
    String output_file = basename(input_file) + ".gistic2.seg"

    command <<<
        set -e
        python <<EOF
import csv
input_file = "${input_file}"
output_file = "${output_file}"

"""
  The column headers are:

(1)  Sample           (sample name)
(2)  Chromosome  (chromosome number)
(3)  Start Position  (segment start position, in bases)
(4)  End Position   (segment end position, in bases)
(5)  Num markers      (number of markers in segment)
(6)  Seg.CN       (log2() -1 of copy number)
"""

if __name__ == "__main__":
    with open(input_file, 'r') as tsvinfp, open(output_file, 'w') as tsvoutfp:
        tsvin = csv.DictReader(tsvinfp, delimiter='\t')
        tsvout = csv.writer(tsvoutfp, delimiter="\t")
        for r in tsvin:
            int_ify_num_points = r["NUM_POINTS_COPY_RATIO_1"].replace(".0", "")
            outrow = [r["SAMPLE"], r["CONTIG"], r["START"], r["END"], int_ify_num_points, r["MEAN_LOG2_COPY_RATIO"]]
            print(outrow)
            tsvout.writerow(outrow)

EOF
    >>>

    runtime {
        docker: docker
        memory: "2000 MB"
        disks: "local-disk 100 HDD"
        preemptible: 3
        cpu: 1
    }
    output {
        File output_file_gistic2 = "${output_file}"
    }
}