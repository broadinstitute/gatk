version 1.0

workflow ImportGenomes {

  input {
    String project_id
    String table_id
    String dataset_name
    String storage_location
    String datatype
    String superpartitioned
    File pet_schema
    File? service_account_json
  }

  String docker="us.gcr.io/broad-gatk/gatk:4.1.7.0"

  call LoadTable {
    input:
      project_id = project_id,
      table_id = table_id,
      dataset_name = dataset_name,
      storage_location = storage_location,
      datatype = datatype,
      superpartitioned = superpartitioned,
      schema = pet_schema,
      service_account_json = service_account_json,
      docker = docker
  }
}

task LoadTable {
  meta {
    volatile: true
  }

  input {
    String project_id
    String table_id
    String dataset_name
    String storage_location
    String datatype
    String superpartitioned
    File schema
    File? service_account_json

    String docker
  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gcloud auth activate-service-account --key-file='~{service_account_json}'
      gcloud config set project ~{project_id}
    fi

    DIR="~{storage_location}/~{datatype}_tsvs/"
    # check for existence of the correct lockfile
    #LOCKFILE="~{storage_location}/LOCKFILE"
    # EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE}) || { echo "Error retrieving lockfile from ${LOCKFILE}" 1>&2 ; exit 1; }
    # CURRENT_RUN_ID="{run_uuid}"

    #if [ "${EXISTING_LOCK_ID}" != "${CURRENT_RUN_ID}" ]; then
    #echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
    #exit 1
    #fi

    #DIR="~{storage_location}/~{datatype}_tsvs/"

    printf -v PADDED_TABLE_ID "%03d" ~{table_id}

    # even for non-superpartitioned tables (e.g. metadata), the TSVs do have the suffix
    FILES="~{datatype}_${PADDED_TABLE_ID}_*"

    NUM_FILES=$(gsutil ls "${DIR}${FILES}" | wc -l)

    if [ ~{superpartitioned} = "true" ]; then
      TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
    else
      TABLE="~{dataset_name}.${PREFIX}~{datatype}"
    fi

    if [ $NUM_FILES -gt 0 ]; then
        # get list of of pet files and their byte sizes
        echo "Getting the file sizes(bytes) and path to each file."
        gsutil du ${DIR}pet_001_* | tr " " "\t" | tr -s "\t" > ~{datatype}_du.txt

        # get total memory in bytes
        echo "Calculating total files' size(bytes)."
        TOTAL_FILE_SIZE=$(awk '{print $1}' OFS="\t" ~{datatype}_du.txt| paste -sd+ - | bc)
        echo "$TOTAL_FILE_SIZE"

        # get number of iterations to loop through file - round up to get full set of files
        num_sets=$(((TOTAL_FILE_SIZE+16492674416639)/16492674416640))
        echo "$num_sets"

        echo "Starting chunking of files."
        for set in $(seq 1 $num_sets)
        do
        # write set of data totaling 16000000000000 bytes to file labeled by set #
        awk '{s+=$1}{print $1"\t"$2"\t"s}' ~{datatype}_du.txt | awk '$3 < 16000000000000 {print $1"\t"$2}' > "${set}"_files_to_load.txt
        # subtract the files put into set file from full list
        awk 'NR==FNR{a[$2]=$0;next}{$3=a[$2]}1' OFS="\t" "${set}"_files_to_load.txt ~{datatype}_du.txt | awk 'NF<3' | cut -f1-2 \
        > tmp_~{datatype}_du.txt && mv tmp_~{datatype}_du.txt ~{datatype}_du.txt

        echo "Copying chunk $set data into separate directory."
        cut -f2 "${set}"_files_to_load.txt | gsutil -m cp -I "${DIR}set_${set}/" 2> gsutil_copy.log

        echo "Running BigQuery load for set $set."
        #bq_job_id=$(bq load --nosync --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" \
        #            "$TABLE" "${DIR}set_${set}/${FILES}" ~{schema} | tr ":" "\t" | cut -f2)
        bq load --nosync --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" \
          "$TABLE" "${DIR}set_${set}/${FILES}" ~{schema} > status_bq_submission

        bq_job_id=$(sed 's/.*://' status_bq_submission)
        # add job ID as key and gs path to the data set uploaded as value
        echo -e "${bq_job_id}\t${set}\t${DIR}set_${set}/" >> bq_load_details.txt

        bq wait "$bq_job_id"
        done

        # for each bq job submitted, run the bq wait command and capture output to log file
        # while IFS="\t" read -r job_id
        # do
        # bq wait $(echo "$job_id" | cut -f1) > bq_wait_status
        # wait_status=$(sed '5q;d' bq_wait_status | tr " " "\t" | tr -s "\t" | cut -f3) 
        # echo "$wait_status"
        # done < bq_load_details.txt >> bq_wait_details.txt

        # paste bq_load_details.txt bq_wait_details.txt > bq_final_job_statuses.txt
        
        # gsutil -m mv "${DIR}set_${set}/${FILES}" "${DIR}set_${set}/done/"
    
    else
        echo "no ${FILES} files to process in $DIR"
    fi
  
>>>

runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 0
    cpu: 1
  }

  output {
    File load_details = "bq_load_details.txt"
    # File wait_details = "bq_wait_details.txt"
    # File final_job_statuses = "bq_final_job_statuses.txt"
  }
}

