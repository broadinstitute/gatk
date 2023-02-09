version 1.0

workflow CromwellOnAzureSqlDatabase {
    call HelloAzure
}

task HelloAzure {
    command <<<
        echo 'Hello Azure!'
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-09"
    }
}
