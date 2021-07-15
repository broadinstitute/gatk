version 1.0

import "GvsWarpTasks.wdl" as Tasks

workflow GvsValidateVatTable {
    input {
        String query_project_id
        String default_dataset
        String vat_table_name
        String? service_account_json_path
    }

    String fq_vat_table = "~{query_project_id}.~{default_dataset}.~{vat_table_name}"

    call EnsureVatTableHasVariants {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            service_account_json_path = service_account_json_path
    }
}

task EnsureVatTableHasVariants {
    input {
        String query_project_id
        String fq_vat_table
        String? service_account_json_path
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project_id}
        fi
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT COUNT (DISTINCT vid) AS count FROM ~{fq_vat_table}' > bq_variant_count.csv

        NUMVARS=$(python3 -c "csvObj=open('bq_variant_count.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

        if [[ $NUMVARS =~ ^[0-9]+$ ]]; then
            echo "The VAT table ~{fq_vat_table} has $NUMVARS variants in it." > validation_results.txt
        else
            echo "Something went wrong. The attempt to count the variants returned: " + $(cat bq_variant_count.csv) > validation_results.txt
        fi
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        String result = read_string('validation_results.txt')
    }
}
