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
        vcf_file: "VCF file to index with tabix"
        tabix_path_1: "Path to first tabix binary (optional, uses system default if not provided)"
        tabix_path_2: "Path to second tabix binary (optional, uses system default if not provided)"
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
    
    String vcf_basename = basename(vcf_file, ".vcf")
    String tabix_cmd = if defined(tabix_path) then select_first([tabix_path]) else "tabix"
    
    command <<<
        set -euo pipefail
        
        # Copy input VCF to working directory
        cp ~{vcf_file} ./input.vcf
        
        # Compress the VCF if it's not already compressed
        if [[ "~{vcf_file}" != *.gz ]]; then
            echo "Compressing VCF file..."
            bgzip -c input.vcf > input.vcf.gz
            VCF_TO_INDEX="input.vcf.gz"
        else
            mv input.vcf input.vcf.gz
            VCF_TO_INDEX="input.vcf.gz"
        fi
        
        # Display tabix version information
        echo "=== Tabix Version Information ===" > ~{task_name}_timing.log
        echo "Tabix command: ~{tabix_cmd}" >> ~{task_name}_timing.log
        if command -v ~{tabix_cmd} &> /dev/null; then
            echo "Tabix version:" >> ~{task_name}_timing.log
            ~{tabix_cmd} --version >> ~{task_name}_timing.log 2>&1 || echo "Version command failed" >> ~{task_name}_timing.log
        else
            echo "ERROR: Tabix command '~{tabix_cmd}' not found" >> ~{task_name}_timing.log
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
        /usr/bin/time -v ~{tabix_cmd} -p vcf $VCF_TO_INDEX 2>&1 | tee -a ~{task_name}_timing.log
        
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