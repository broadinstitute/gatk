version 1.0

workflow myWorkflow {
    input {
        String sample_id
        File gvcf_gz
        File pyscript
    }

    call Calculate { input: sample_id = sample_id, gvcf_gz = gvcf_gz, pyscript = pyscript }

    output {
        File output_stats = Calculate.output_stats
    }
}

task Calculate {
    input {
        String sample_id
        File gvcf_gz
        File pyscript
    }
    parameter_meta {gvcf_gz: {localization_optional: true}}

    command {
        cat ${gvcf_gz} | gunzip | python ${pyscript} ${sample_id} > ${sample_id}.gqstats.txt
    }
    runtime {
        docker: "python:3.7"
        memory: "4 GB"
        cpu: 2
    }
    output {
        File output_stats = "${sample_id}.gqstats.txt"
    }
}
