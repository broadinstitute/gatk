version 1.0

workflow GvsExtractAvroFilesForHail {
    input {
        String project_id
        String dataset
        String filter_set_name
    }

    call OutputPath { input: go = true }

    call ExtractFromNonSuperpartitionedTables {
        input:
            project_id = project_id,
            dataset = dataset,
            filter_set_name = filter_set_name,
            avro_sibling = OutputPath.out
    }

    call CountSamples {
        input:
            project_id = project_id,
            dataset = dataset
    }

    # Superpartitions have max size 4000. The inner '- 1' is so the 4000th (and multiples of 4000) sample lands in the
    # appropriate partition, the outer '+ 1' is to iterate over the correct number of partitions.
    scatter (i in range(((CountSamples.num_samples - 1) / 4000) + 1)) {
        call ExtractFromSuperpartitionedTables {
            input:
                project_id = project_id,
                dataset = dataset,
                filter_set_name = filter_set_name,
                avro_sibling = OutputPath.out,
                table_index = i + 1 #  'i' is 0-based so add 1.
        }
    }

    call GenerateHailScript {
        avro_prefix = ExtractFromNonSuperpartitionedTables.output_prefix
    }
    output {
        String output_prefix = ExtractFromNonSuperpartitionedTables.output_prefix
    }
}


task OutputPath {
    meta {
        description: "Does nothing but produce the cloud path to its stdout."
    }
    input {
        Boolean go = true
    }
    command <<<
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


task CountSamples {
    meta {
        description: "Counts the number of samples in the sample_info table efficiently."
    }
    input {
        String project_id
        String dataset
    }
    command <<<
        python3 <<FIN

        from google.cloud import bigquery

        client = bigquery.Client(project="~{project_id}")
        sample_info_table_id = f'~{project_id}.~{dataset}.sample_info'
        sample_info_table = client.get_table(sample_info_table_id)

        print(str(sample_info_table.num_rows))

        FIN
    >>>

    output {
        Int num_samples = read_int(stdout())
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}


task ExtractFromNonSuperpartitionedTables {
    meta {
        description: "Extracts from the non-superpartitioned tables: sample_info, filter_set_info, filter_set_sites"
    }
    input {
        String project_id
        String dataset
        String filter_set_name
        String avro_sibling
    }
    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        echo $avro_prefix > "avro_prefix.out"

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/sample_mapping/sample_mapping_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT sample_id, sample_name, '40',
            'gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list' as intervals_file
            FROM \`~{project_id}.~{dataset}.sample_info\`
            ORDER BY sample_id
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vqsr_filtering_data/vqsr_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, type as model, ref, alt, vqslod, yng_status
            FROM \`~{project_id}.~{dataset}.filter_set_info\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/site_filtering_data/site_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, filters
            FROM \`~{project_id}.~{dataset}.filter_set_sites\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vqsr_tranche_data/vqsr_tranche_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT model, truth_sensitivity, min_vqslod, filter_name
            FROM \`~{project_id}.~{dataset}.filter_set_tranches\`
            WHERE filter_set_name = '~{filter_set_name}'
        "
    >>>

    output {
        Boolean done = true
        String output_prefix = read_string("avro_prefix.out")
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
    }
}


task ExtractFromSuperpartitionedTables {
    meta {
        description: "Extracts from the superpartitioned tables: vet_<table index>, ref_ranges_<table index>"
    }
    input {
        String project_id
        String dataset
        String filter_set_name
        String avro_sibling
        Int table_index
    }
    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
        table_index: "1-based index for the superpartitioned ref_ranges and vet tables to be extracted by this call."
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        str_table_index=$(printf "%03d" ~{table_index})

        # These bq exports error out if there are any objects at the sibling level to where output files would be written
        # so an extra layer of `vet_${str_table_index}` is inserted here.
        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vets/vet_${str_table_index}/vet_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, sample_id, ref, REPLACE(alt,',<NON_REF>','') alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,',')[OFFSET(0)] as int64) as RGQ
            FROM \`~{project_id}.~{dataset}.vet_${str_table_index}\`
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/refs/ref_ranges_${str_table_index}/ref_ranges_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, sample_id, length, state
            FROM \`~{project_id}.~{dataset}.ref_ranges_${str_table_index}\`
            ORDER BY location
        "
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
    }
}

task GenerateHailScript {
    input {
        String avro_prefix
        String slashed_avro_prefix = avro_prefix + "/"
    }
    command <<<
        gsutil ls -r "${avro_prefix}" | grep '\.avro$' | sed 's!~{slashed_avro_prefix}!!' > avro_files.out

        python <<FIN
            print """
            import hail as hl

            rg38 = hl.get_reference('GRCh38')
            rg38.add_sequence(
                'gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz',
                'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai'
            )

            hl.import_gvs(
            """

            def superpartition(name, files):
                for sp in

            f = open('avro_files.out', 'r')
            files=defaultdict(list)
            for avro in f.readlines():
                parts=avro.split('/')
                if len(parts) == 1:
                    parts[part[0]].append(avro)
                else:
                    m=re.search(r"^(?P<kind>[^/]+)/(?P=kind)_(?P<super>[^/]+)/.*", avro)
                    if not m:
                        raise ValueError(f"'{avro}' does not look like a path to a superpartitioned Avro file!")
                    kind=m.group('kind')
                    superpartition=m.group('super')
                    zero_based_index=int(superpartition) - 1
                    if len(files[kind]) == zero_based_index:
                        files[kind].append([])
                    files[kind][zero_based_index].append(avro)


        FIN


            refs=[
                $ref_ranges
            ],
            vets=[
                $vets
            ],
            sample_mapping=[$sample_mapping],
            site_filtering=[$site_filtering],
            vqsr_filtering=[$vqsr_filtering],
            vqsr_tranches=[$vqsr_tranches],
            final_path=$final_vds_path,
            tmp_dir=$temp_vds_path,
            reference_genome=rg38,
        )

        vds = hl.vds.read_vds($final_vds_path)

        mt = hl.vds.to_dense_mt(vds)
        fail_case = 'FAIL'
        mt = mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))
        hl.export_vcf(mt, $final_vcf_path)
        FIN
    >>>
    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
    }
}
