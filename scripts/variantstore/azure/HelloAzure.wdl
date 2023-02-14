version 1.0

workflow HelloAzure {
    input {
        String salutation
    }
    call Hello {
        input:
            salutation = salutation
    }

    output {
        String salutation = Hello.salutation
    }
}

task Hello {
    input {
        String salutation
    }
    command <<<
        echo '~{salutation}'
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-09"
    }
    output {
        String salutation = read_string(stdout())
    }
}
