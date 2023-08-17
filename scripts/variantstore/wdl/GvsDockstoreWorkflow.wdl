version 1.0

workflow GvsDockstoreWorkflow {
    call Hello
}

task Hello {
    command <<<
        echo "hello from Dockstore workflow"
    >>>
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        File out = read_string(stdout())
    }
}
