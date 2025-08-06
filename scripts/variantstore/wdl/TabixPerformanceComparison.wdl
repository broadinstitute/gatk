version 1.0

workflow TabixPerformanceComparison {
    input {
        File vcf_file
        String? tabix_path_1
        String? tabix_path_2
        
        # Runtime parameters
        String docker = "us.gcr.io/broad-gatk/gatk:latest"
        Int memory_gb = 4
        Int disk_size_gb = 50
        Int cpu_count = 1
        Boolean use_preemptibles = true
    }
    
    parameter_meta {
        vcf_file: "VCF file to index with tabix (supports .vcf, .vcf.gz, .vcf.bz2 formats)"
        tabix_path_1: "Google Cloud Storage path to first tabix binary (optional, uses system default if not provided)"
        tabix_path_2: "Google Cloud Storage path to second tabix binary (optional, uses system default if not provided)"
        docker: "Docker image to use"
        memory_gb: "Memory in GB"
        disk_size_gb: "Disk size in GB"
        cpu_count: "Number of CPUs"
        use_preemptibles: "Use preemptible instances"
    }
    
    meta {
        description: "Compare performance of two tabix versions by timing indexing operations in parallel"
        author: "GATK Team"
    }
    
    # Always run both tabix tasks in parallel
    call TabixIndexTask as TabixTask1 {
        input:
            vcf_file = vcf_file,
            tabix_path = tabix_path_1,
            task_name = "tabix_1",
            docker = docker,
            memory_gb = memory_gb,
            disk_size_gb = disk_size_gb,
            cpu_count = cpu_count,
            use_preemptibles = use_preemptibles
    }
    
    call TabixIndexTask as TabixTask2 {
        input:
            vcf_file = vcf_file,
            tabix_path = tabix_path_2,
            task_name = "tabix_2",
            docker = docker,
            memory_gb = memory_gb,
            disk_size_gb = disk_size_gb,
            cpu_count = cpu_count,
            use_preemptibles = use_preemptibles
    }
    
    output {
        File tabix1_timing_log = TabixTask1.timing_log
        File tabix1_index_file = TabixTask1.index_file
        Float tabix1_execution_time = TabixTask1.execution_time
        
        File tabix2_timing_log = TabixTask2.timing_log
        File tabix2_index_file = TabixTask2.index_file
        Float tabix2_execution_time = TabixTask2.execution_time
    }
}

task TabixIndexTask {
    input {
        File vcf_file
        String? tabix_path
        String task_name
        
        # Runtime parameters
        String docker
        Int memory_gb
        Int disk_size_gb
        Int cpu_count
        Boolean use_preemptibles
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # Set up tabix command - copy custom binary if provided
        if [[ "~{defined(tabix_path)}" == "true" ]]; then
            echo "Custom tabix path provided: ~{select_first([tabix_path, ""])}"
            echo "Copying custom tabix binary..."
            gsutil cp "~{select_first([tabix_path, ""])}" ./custom_tabix
            chmod +x ./custom_tabix
            TABIX_CMD="./custom_tabix"
            echo "Custom tabix binary copied and made executable"
        else
            echo "Using system default tabix"
            TABIX_CMD="tabix"
        fi

        # Use input VCF directly without copying
        # Determine the file to index based on the input file
        INPUT_FILE="~{vcf_file}"
        VCF_TO_INDEX="$INPUT_FILE"
        
        # Check if the file is compressed and needs decompression/recompression for tabix
        if [[ "$INPUT_FILE" == *.vcf.gz ]]; then
            echo "Input is already bgzip compressed VCF, using directly"
            VCF_TO_INDEX="$INPUT_FILE"
        elif [[ "$INPUT_FILE" == *.vcf.bz2 ]]; then
            echo "Input is bzip2 compressed, decompressing and recompressing with bgzip..."
            bzcat "$INPUT_FILE" | bgzip -c > input_recompressed.vcf.gz
            VCF_TO_INDEX="input_recompressed.vcf.gz"
        elif [[ "$INPUT_FILE" == *.vcf ]]; then
            echo "Input is uncompressed VCF, compressing with bgzip..."
            bgzip -c "$INPUT_FILE" > input_compressed.vcf.gz
            VCF_TO_INDEX="input_compressed.vcf.gz"
        else
            echo "Unsupported file format. Expected .vcf, .vcf.gz, or .vcf.bz2"
            exit 1
        fi
        
        # Display tabix version information
        echo "=== Tabix Version Information ===" > ~{task_name}_timing.log
        echo "Tabix command: $TABIX_CMD" >> ~{task_name}_timing.log
        if command -v "$TABIX_CMD" &> /dev/null; then
            echo "Tabix version:" >> ~{task_name}_timing.log
            "$TABIX_CMD" --version >> ~{task_name}_timing.log 2>&1 || echo "Version command failed" >> ~{task_name}_timing.log
        else
            echo "ERROR: Tabix command '$TABIX_CMD' not found" >> ~{task_name}_timing.log
            exit 1
        fi
        echo "" >> ~{task_name}_timing.log
        
        # Get file size information
        echo "=== File Information ===" >> ~{task_name}_timing.log
        echo "VCF file size: $(stat -c%s $VCF_TO_INDEX) bytes" >> ~{task_name}_timing.log
        echo "VCF file: $(ls -lh $VCF_TO_INDEX)" >> ~{task_name}_timing.log
        echo "" >> ~{task_name}_timing.log
        
        # Time the tabix indexing operation
        echo "=== Timing Information ===" >> ~{task_name}_timing.log
        echo "Starting tabix indexing at: $(date)" >> ~{task_name}_timing.log
        
        START_TIME=$(date +%s.%N)
        
        # Run tabix indexing with timing
        /usr/bin/time -v "$TABIX_CMD" -p vcf $VCF_TO_INDEX 2>&1 | tee -a ~{task_name}_timing.log
        
        END_TIME=$(date +%s.%N)
        EXECUTION_TIME=$(echo "$END_TIME - $START_TIME" | bc -l)
        
        echo "Finished tabix indexing at: $(date)" >> ~{task_name}_timing.log
        echo "Total execution time: ${EXECUTION_TIME} seconds" >> ~{task_name}_timing.log
        echo "" >> ~{task_name}_timing.log
        
        # Verify index was created
        if [[ -f "${VCF_TO_INDEX}.tbi" ]]; then
            echo "=== Index File Information ===" >> ~{task_name}_timing.log
            echo "Index file created successfully: ${VCF_TO_INDEX}.tbi" >> ~{task_name}_timing.log
            echo "Index file size: $(stat -c%s ${VCF_TO_INDEX}.tbi) bytes" >> ~{task_name}_timing.log
            echo "Index file: $(ls -lh ${VCF_TO_INDEX}.tbi)" >> ~{task_name}_timing.log
            
            # Copy index file to output
            cp "${VCF_TO_INDEX}.tbi" ~{task_name}_output.tbi
        else
            echo "ERROR: Index file was not created" >> ~{task_name}_timing.log
            exit 1
        fi
        
        # Save execution time for WDL output
        echo $EXECUTION_TIME > execution_time.txt
        
        echo "=== Task completed successfully ===" >> ~{task_name}_timing.log
    >>>
    
    output {
        File timing_log = "~{task_name}_timing.log"
        File index_file = "~{task_name}_output.tbi"
        Float execution_time = read_float("execution_time.txt")
        File monitoring_log = "monitoring.log"
    }
    
    runtime {
        docker: docker
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        cpu: cpu_count
        preemptible: if use_preemptibles then 3 else 0
        bootDiskSizeGb: 15
    }
}