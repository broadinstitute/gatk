version 1.0

# Run PrintReads (GATK Version 4.1.5.0-27-g4c9f5fa-SNAPSHOT)
#
# Print reads in the SAM/BAM/CRAM file
#
#  General Workflow (non-tool) Arguments
#    gatk                                              Location of gatk used to run this workflow
#
#  Required Tool Arguments
#    I                                                  BAM/SAM/CRAM file containing reads                          
#    O                                                  Write output to this file                                   
#
#  Optional Tool Arguments
#    arguments_file                                     read one or more arguments files and add them to the command line
#    cloud_index_prefetch_buffer                        Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudP
#    cloud_prefetch_buffer                              Size of the cloud-only prefetch buffer (in MB; 0 to disable).
#    disable_bam_index_caching                          If true, don't cache bam indexes, this will reduce memory requirements but may h
#    disable_sequence_dictionary_validation             If specified, do not check the sequence dictionaries from our inputs for compati
#    gcs_max_retries                                    If the GCS bucket channel errors out, how many times it will attempt to re-initi
#    gcs_project_for_requester_pays                     Project to bill when accessing "requester pays" buckets. If unset, these buckets
#    help                                               display the help message                                    
#    interval_merging_rule                              Interval merging rule for abutting intervals                
#    intervals                                          One or more genomic intervals over which to operate         
#    reference                                          Reference sequence                                          
#    sites_only_vcf_output                              If true, don't emit genotype fields when writing vcf file output.
#    version                                            display the version number for this tool                    
#
#  Optional Common Arguments
#    add_output_sam_program_record                      If true, adds a PG tag to created SAM/BAM/CRAM files.       
#    add_output_vcf_command_line                        If true, adds a command line header line to created VCF files.
#    create_output_bam_index                            If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.
#    create_output_bam_md5                              If true, create a MD5 digest for any BAM/SAM/CRAM file created
#    create_output_variant_index                        If true, create a VCF index when writing a coordinate-sorted VCF file.
#    create_output_variant_md5                          If true, create a a MD5 digest any VCF file created.        
#    disable_read_filter                                Read filters to be disabled before analysis                 
#    disable_tool_default_read_filters                  Disable all tool default read filters (WARNING: many tools will not function cor
#    exclude_intervals                                  One or more genomic intervals to exclude from processing    
#    gatk_config_file                                   A configuration file to use with the GATK.                  
#    interval_exclusion_padding                         Amount of padding (in bp) to add to each interval you are excluding.
#    interval_padding                                   Amount of padding (in bp) to add to each interval you are including.
#    interval_set_rule                                  Set merging approach to use for combining interval inputs   
#    lenient                                            Lenient processing of VCF files                             
#    QUIET                                              Whether to suppress job-summary info on System.err.         
#    read_filter                                        Read filters to be applied before analysis                  
#    read_index                                         Indices to use for the read inputs. If specified, an index must be provided for 
#    read_validation_stringency                         Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The 
#    seconds_between_progress_updates                   Output traversal statistics every time this many seconds elapse
#    sequence_dictionary                                Use the given sequence dictionary as the master/canonical sequence dictionary.  
#    tmp_dir                                            Temp directory to use.                                      
#    use_jdk_deflater                                   Whether to use the JdkDeflater (as opposed to IntelDeflater)
#    use_jdk_inflater                                   Whether to use the JdkInflater (as opposed to IntelInflater)
#    verbosity                                          Control verbosity of logging.                               

workflow PrintReads {

  input {
    #GATK location
    String gatk

    # Required Arguments
    Array[String] I
    File O

    # Optional Tool Arguments
    Array[File]? arguments_file
    Int? cloud_index_prefetch_buffer
    Int? cloud_prefetch_buffer
    Boolean? disable_bam_index_caching
    Boolean? disable_sequence_dictionary_validation
    Int? gcs_max_retries
    String? gcs_project_for_requester_pays
    Boolean? help
    String? interval_merging_rule
    Array[String]? intervals
    File? reference
    File? referenceDictionary
    File? referenceIndex
    Boolean? sites_only_vcf_output
    Boolean? version

    # Optional Common Arguments
    Boolean? add_output_sam_program_record
    Boolean? add_output_vcf_command_line
    Boolean? create_output_bam_index
    Boolean? create_output_bam_md5
    Boolean? create_output_variant_index
    Boolean? create_output_variant_md5
    Array[String]? disable_read_filter
    Boolean? disable_tool_default_read_filters
    Array[String]? exclude_intervals
    String? gatk_config_file
    Int? interval_exclusion_padding
    Int? interval_padding
    String? interval_set_rule
    Boolean? lenient
    Boolean? QUIET
    Array[String]? read_filter
    Array[String]? read_index
    String? read_validation_stringency
    Float? seconds_between_progress_updates
    String? sequence_dictionary
    File? tmp_dir
    Boolean? use_jdk_deflater
    Boolean? use_jdk_inflater
    String? verbosity

  }

  call PrintReadsTask {

    input:

    #GATK location
    gatk                                               = gatk,

    # Required Arguments
    I                                                  = I,
    O                                                  = O,

    # Optional Tool Arguments
    arguments_file                                     = arguments_file,
    cloud_index_prefetch_buffer                        = cloud_index_prefetch_buffer,
    cloud_prefetch_buffer                              = cloud_prefetch_buffer,
    disable_bam_index_caching                          = disable_bam_index_caching,
    disable_sequence_dictionary_validation             = disable_sequence_dictionary_validation,
    gcs_max_retries                                    = gcs_max_retries,
    gcs_project_for_requester_pays                     = gcs_project_for_requester_pays,
    help                                               = help,
    interval_merging_rule                              = interval_merging_rule,
    intervals                                          = intervals,
    reference                                          = reference,
    referenceDictionary                                = referenceDictionary,
    referenceIndex                                     = referenceIndex,
    sites_only_vcf_output                              = sites_only_vcf_output,
    version                                            = version,

    # Optional Common Arguments
    add_output_sam_program_record                      = add_output_sam_program_record,
    add_output_vcf_command_line                        = add_output_vcf_command_line,
    create_output_bam_index                            = create_output_bam_index,
    create_output_bam_md5                              = create_output_bam_md5,
    create_output_variant_index                        = create_output_variant_index,
    create_output_variant_md5                          = create_output_variant_md5,
    disable_read_filter                                = disable_read_filter,
    disable_tool_default_read_filters                  = disable_tool_default_read_filters,
    exclude_intervals                                  = exclude_intervals,
    gatk_config_file                                   = gatk_config_file,
    interval_exclusion_padding                         = interval_exclusion_padding,
    interval_padding                                   = interval_padding,
    interval_set_rule                                  = interval_set_rule,
    lenient                                            = lenient,
    QUIET                                              = QUIET,
    read_filter                                        = read_filter,
    read_index                                         = read_index,
    read_validation_stringency                         = read_validation_stringency,
    seconds_between_progress_updates                   = seconds_between_progress_updates,
    sequence_dictionary                                = sequence_dictionary,
    tmp_dir                                            = tmp_dir,
    use_jdk_deflater                                   = use_jdk_deflater,
    use_jdk_inflater                                   = use_jdk_inflater,
    verbosity                                          = verbosity,

  }

  output {
    # Workflow Outputs                                  
    File PrintReadsO = PrintReadsTask.PrintReadsTask_O
  }
}

task PrintReadsTask {

  input {
    String gatk
    Array[String] I
    File O
    Array[File]? arguments_file
    Int? cloud_index_prefetch_buffer
    Int? cloud_prefetch_buffer
    Boolean? disable_bam_index_caching
    Boolean? disable_sequence_dictionary_validation
    Int? gcs_max_retries
    String? gcs_project_for_requester_pays
    Boolean? help
    String? interval_merging_rule
    Array[String]? intervals
    File? reference
    File? referenceDictionary
    File? referenceIndex
    Boolean? sites_only_vcf_output
    Boolean? version
    Boolean? add_output_sam_program_record
    Boolean? add_output_vcf_command_line
    Boolean? create_output_bam_index
    Boolean? create_output_bam_md5
    Boolean? create_output_variant_index
    Boolean? create_output_variant_md5
    Array[String]? disable_read_filter
    Boolean? disable_tool_default_read_filters
    Array[String]? exclude_intervals
    String? gatk_config_file
    Int? interval_exclusion_padding
    Int? interval_padding
    String? interval_set_rule
    Boolean? lenient
    Boolean? QUIET
    Array[String]? read_filter
    Array[String]? read_index
    String? read_validation_stringency
    Float? seconds_between_progress_updates
    String? sequence_dictionary
    File? tmp_dir
    Boolean? use_jdk_deflater
    Boolean? use_jdk_inflater
    String? verbosity

  }

  command <<<
    ~{gatk} PrintReads \
    --I ~{sep=' --I ' I} \
    --O ~{sep=' --O ' O} \
    ~{true='--arguments_file ' false='' defined(arguments_file)}~{sep=' --arguments_file ' arguments_file} \
    ~{true='--cloud_index_prefetch_buffer ' false='' defined(cloud_index_prefetch_buffer)}~{sep=' --cloud_index_prefetch_buffer ' cloud_index_prefetch_buffer} \
    ~{true='--cloud_prefetch_buffer ' false='' defined(cloud_prefetch_buffer)}~{sep=' --cloud_prefetch_buffer ' cloud_prefetch_buffer} \
    ~{true='--disable_bam_index_caching ' false='' defined(disable_bam_index_caching)}~{sep=' --disable_bam_index_caching ' disable_bam_index_caching} \
    ~{true='--disable_sequence_dictionary_validation ' false='' defined(disable_sequence_dictionary_validation)}~{sep=' --disable_sequence_dictionary_validation ' disable_sequence_dictionary_validation} \
    ~{true='--gcs_max_retries ' false='' defined(gcs_max_retries)}~{sep=' --gcs_max_retries ' gcs_max_retries} \
    ~{true='--gcs_project_for_requester_pays ' false='' defined(gcs_project_for_requester_pays)}~{sep=' --gcs_project_for_requester_pays ' gcs_project_for_requester_pays} \
    ~{true='--help ' false='' defined(help)}~{sep=' --help ' help} \
    ~{true='--interval_merging_rule ' false='' defined(interval_merging_rule)}~{sep=' --interval_merging_rule ' interval_merging_rule} \
    ~{true='--intervals ' false='' defined(intervals)}~{sep=' --intervals ' intervals} \
    ~{true='--reference ' false='' defined(reference)}~{sep=' --reference ' reference} \
    ~{true='--sites_only_vcf_output ' false='' defined(sites_only_vcf_output)}~{sep=' --sites_only_vcf_output ' sites_only_vcf_output} \
    ~{true='--version ' false='' defined(version)}~{sep=' --version ' version} \
    ~{true='--add_output_sam_program_record ' false='' defined(add_output_sam_program_record)}~{sep=' --add_output_sam_program_record ' add_output_sam_program_record} \
    ~{true='--add_output_vcf_command_line ' false='' defined(add_output_vcf_command_line)}~{sep=' --add_output_vcf_command_line ' add_output_vcf_command_line} \
    ~{true='--create_output_bam_index ' false='' defined(create_output_bam_index)}~{sep=' --create_output_bam_index ' create_output_bam_index} \
    ~{true='--create_output_bam_md5 ' false='' defined(create_output_bam_md5)}~{sep=' --create_output_bam_md5 ' create_output_bam_md5} \
    ~{true='--create_output_variant_index ' false='' defined(create_output_variant_index)}~{sep=' --create_output_variant_index ' create_output_variant_index} \
    ~{true='--create_output_variant_md5 ' false='' defined(create_output_variant_md5)}~{sep=' --create_output_variant_md5 ' create_output_variant_md5} \
    ~{true='--disable_read_filter ' false='' defined(disable_read_filter)}~{sep=' --disable_read_filter ' disable_read_filter} \
    ~{true='--disable_tool_default_read_filters ' false='' defined(disable_tool_default_read_filters)}~{sep=' --disable_tool_default_read_filters ' disable_tool_default_read_filters} \
    ~{true='--exclude_intervals ' false='' defined(exclude_intervals)}~{sep=' --exclude_intervals ' exclude_intervals} \
    ~{true='--gatk_config_file ' false='' defined(gatk_config_file)}~{sep=' --gatk_config_file ' gatk_config_file} \
    ~{true='--interval_exclusion_padding ' false='' defined(interval_exclusion_padding)}~{sep=' --interval_exclusion_padding ' interval_exclusion_padding} \
    ~{true='--interval_padding ' false='' defined(interval_padding)}~{sep=' --interval_padding ' interval_padding} \
    ~{true='--interval_set_rule ' false='' defined(interval_set_rule)}~{sep=' --interval_set_rule ' interval_set_rule} \
    ~{true='--lenient ' false='' defined(lenient)}~{sep=' --lenient ' lenient} \
    ~{true='--QUIET ' false='' defined(QUIET)}~{sep=' --QUIET ' QUIET} \
    ~{true='--read_filter ' false='' defined(read_filter)}~{sep=' --read_filter ' read_filter} \
    ~{true='--read_index ' false='' defined(read_index)}~{sep=' --read_index ' read_index} \
    ~{true='--read_validation_stringency ' false='' defined(read_validation_stringency)}~{sep=' --read_validation_stringency ' read_validation_stringency} \
    ~{true='--seconds_between_progress_updates ' false='' defined(seconds_between_progress_updates)}~{sep=' --seconds_between_progress_updates ' seconds_between_progress_updates} \
    ~{true='--sequence_dictionary ' false='' defined(sequence_dictionary)}~{sep=' --sequence_dictionary ' sequence_dictionary} \
    ~{true='--tmp_dir ' false='' defined(tmp_dir)}~{sep=' --tmp_dir ' tmp_dir} \
    ~{true='--use_jdk_deflater ' false='' defined(use_jdk_deflater)}~{sep=' --use_jdk_deflater ' use_jdk_deflater} \
    ~{true='--use_jdk_inflater ' false='' defined(use_jdk_inflater)}~{sep=' --use_jdk_inflater ' use_jdk_inflater} \
    ~{true='--verbosity ' false='' defined(verbosity)}~{sep=' --verbosity ' verbosity} \
  >>>

  runtime {
    memory: "1GB"
  }

  output {
    # Task Outputs                                      
    File PrintReadsTask_O = "${O}"
  }
 }

