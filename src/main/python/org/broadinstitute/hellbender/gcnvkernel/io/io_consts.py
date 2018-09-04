# interval list .tsv file column names
contig_column_name = "CONTIG"
start_column_name = "START"
end_column_name = "END"
count_column_name = "COUNT"

# prefix for saving posteriors for multiple samples
sample_folder_prefix = "SAMPLE_"

# log copy number posterior matrix column name prefix for each integer copy number state
copy_number_column_prefix = "COPY_NUMBER_"

# ploidy prior table header column names
ploidy_prior_contig_name_column = "CONTIG_NAME"
ploidy_prior_prefix = "PLOIDY_PRIOR_"

# column names for ploidy and depth .tsv outputs
sample_name_column_name = "SAMPLE_NAME"
global_read_depth_column_name = "GLOBAL_READ_DEPTH"
average_ploidy_column_name = "AVERAGE_PLOIDY"
ploidy_column_name = "PLOIDY"
ploidy_gq_column_name = "PLOIDY_GQ"

# column names for copy-number segments file
call_copy_number_column_name = "CALL_COPY_NUMBER"
num_points_column_name = "NUM_POINTS"
quality_some_called_column_name = "QUALITY_SOME_CALLED"
quality_all_called_column_name = "QUALITY_ALL_CALLED"
quality_start_column_name = "QUALITY_START"
quality_end_column_name = "QUALITY_END"

# column name for baseline copy-number file
baseline_copy_number_column_name = "BASELINE_COPY_NUMBER"

# regular expression for matching sample name from header comment line
sample_name_header_regexp = "^@RG.*SM:(.*)[\t]*.*$"

# prefix for adding sample name as a header comment line
sample_name_sam_header_prefix = "RG\tID:GATKCopyNumber\tSM:"

default_comment_char = "@"
default_delimiter_char = "\t"

# default file names for loading and saving models, posteriors, and configurations
default_sample_read_depth_tsv_filename = "global_read_depth.tsv"
default_sample_name_txt_filename = "sample_name.txt"
default_sample_contig_ploidy_tsv_filename = "contig_ploidy.tsv"
default_copy_number_log_posterior_tsv_filename = "log_q_c_tc.tsv"
default_copy_number_log_emission_tsv_filename = "log_c_emission_tc.tsv"
default_class_log_posterior_tsv_filename = "log_q_tau_tk.tsv"
default_baseline_copy_number_tsv_filename = "baseline_copy_number_t.tsv"
default_copy_number_segments_tsv_filename = "copy_number_segments.tsv"

default_denoising_config_json_filename = "denoising_config.json"
default_calling_config_json_filename = "calling_config.json"
default_ploidy_config_json_filename = "ploidy_config.json"
default_gcnvkernel_version_json_filename = "gcnvkernel_version.json"

default_interval_list_filename = "interval_list.tsv"
default_contig_ploidy_prior_tsv_filename = "contig_ploidy_prior.tsv"

default_adamax_m_filename = "adamax_m.npy"
default_adamax_u_filename = "adamax_u.npy"
default_adamax_res_filename = "adamax_res.npy"
