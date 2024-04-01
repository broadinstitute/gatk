version 1.0

workflow RetryWithMoreMemoryWorkflow {
    call RetryWithMoreMemory
}

task RetryWithMoreMemory {
    command <<<
        # dump everything
        set

        # Fake being out of memory
        echo "OutOfMemory" 1>&2

        exit 137
    >>>
    runtime {
        docker: "ubuntu:latest"
        memory: "7 GiB"
        maxRetries: 3
    }
}
