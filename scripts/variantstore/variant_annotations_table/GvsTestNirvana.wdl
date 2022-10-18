version 1.0

workflow GvsTestNirvana {
    input {
        Int start_index
        Int end_index
        String prefix_path
        String vcf_template
        File nirvana_references = "gs://gvs-internal/bigquery-jointcalling/VAT/Nirvana/Nirvana-references-2022-10-07.tgz"
        File ref = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        Int scatter_width = 100
    }

    if (false) {
        call DownloadNirvanaReferences {
        }
    }

    scatter (i in range(scatter_width)) {
        call CreateMiniSitesOnlyVcfs {
            input:
                shard_index = i,
                scatter_width = scatter_width,
                start_index = start_index,
                end_index = end_index,
                prefix_path = prefix_path,
                vcf_template = vcf_template,
        }
    }

    Array[File] mini_sites_only_vcfs = flatten(CreateMiniSitesOnlyVcfs.sites_only_vcfs)
    Array[File] mini_sites_only_vcf_indexes = flatten(CreateMiniSitesOnlyVcfs.sites_only_vcf_indexes)

    call RunNirvana {
        input:
            mini_sites_only_vcfs = mini_sites_only_vcfs,
            mini_sites_only_vcf_indexes = mini_sites_only_vcf_indexes,
            references = nirvana_references,
    }

    output {
        Array[File] sites_only_vcfs = mini_sites_only_vcfs
        Array[File] sites_only_vcf_indexes = mini_sites_only_vcf_indexes
        Array[File] nirvana_jsons = RunNirvana.jsons
        Array[File] nirvana_json_indexes = RunNirvana.json_indexes
    }
}

task CreateMiniSitesOnlyVcfs {
    input {
        Boolean go = true
        String prefix_path
        String vcf_template
        Int shard_index
        Int scatter_width
        Int start_index
        Int end_index
        Int mini_sites_cpus = 4
        Int preemptible = 2
    }
    meta {
        description: "Makes mini (same sharding as GVS output shards) sites only vcfs. Drop 'genotypes' (samples) ASAP so the remainder of the pipeline flows more smoothly."
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o xtrace -o pipefail

        # Strip any trailing slash on the prefix for consistency
        prefix_path="~{prefix_path}"
        prefix_path=${prefix_path%/}

        # Localize and process only the vcfs and indexes for this shard.
        touch localize-manifest.txt
        for i in $(seq ~{start_index + shard_index} ~{scatter_width} ~{end_index})
        do
            vcf=$(printf "${prefix_path}/~{vcf_template}" $i)
            echo $vcf >> localize-manifest.txt
            echo ${vcf}.tbi >> localize-manifest.txt
        done

        # Localize for this shard with `gcloud storage`
        # https://cloud.google.com/sdk/gcloud/reference/storage/cp
        cat localize-manifest.txt | gcloud storage cp -I .

        # [W::hts_idx_load3] The index file is older than the data file
        touch *.tbi

        # Extract the base names (remove bucket and "directory" structure as well as .vcf.gz extension).
        sed -n -E 's;.*/(.*)\.vcf\.gz$;\1;p' localize-manifest.txt > vcf_base_names.txt

        cat vcf_base_names.txt | xargs -I % -n 1 -P ~{mini_sites_cpus} bash -c '
            bcftools view --no-update --drop-genotypes --output-type b %.vcf.gz |
            bcftools annotate --remove FORMAT/AD --output-type b |
            bcftools norm --multiallelics - --output-type z --output sites-%.vcf.gz;
            bcftools index --tbi sites-%.vcf.gz
        '
    >>>
    output {
        Array[File] sites_only_vcfs = glob("sites-*.vcf.gz")
        Array[File] sites_only_vcf_indexes = glob("sites-*.tbi")
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2022-10-12-alpine"
        cpu: "~{mini_sites_cpus}"
        disks: "local-disk 500 HDD"
        preemptible: "~{preemptible}"
    }
}


task RunNirvana {
    input {
        # These vcfs and indexes are coming from many different shards with 2 globs each, so Cromwell would localize
        # these inputs in 2 separate batches per shard. It's a LOT more efficient to just do the localization ourselves
        # in one batch, plus we can use `gcloud storage` to do that. Also localize the references ourselves.
        Array[File] mini_sites_only_vcfs
        Array[File] mini_sites_only_vcf_indexes
        File references
        File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
        Int nirvana_memory_gib = 30
    }
    parameter_meta {
        mini_sites_only_vcfs:
        {
            localization_optional: true
        }
        mini_sites_only_vcf_indexes:
        {
            localization_optional: true
        }
        references:
        {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o xtrace -o pipefail

        bash ~{monitoring_script} &

        # manifest of files to be localized explicitly
        touch manifest.txt
        for file in ~{references} ~{sep=' ' mini_sites_only_vcfs} ~{sep=' ' mini_sites_only_vcf_indexes}
        do
            echo $file >> manifest.txt
        done

        # https://cloud.google.com/sdk/gcloud/reference/storage/cp
        cat manifest.txt | gcloud storage cp -I .

        tar xfz $(basename ~{references})

        for vcf in *.vcf.gz
        do
            # strip the first dot and everything after
            base=${vcf%%.*}
            /root/.dotnet/dotnet /Nirvana/Nirvana.dll \
                -c references/Cache/GRCh38/Both \
                -r references/References/Homo_sapiens.GRCh38.Nirvana.dat \
                --sd references/SupplementaryAnnotation/GRCh38 \
                -i $vcf -o $base
        done
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19"
        disks: "local-disk 2000 HDD"
        memory: "~{nirvana_memory_gib} GiB"
    }
    output {
        Array[File] jsons = glob("*.json.gz")
        Array[File] json_indexes = glob("*.json.gz.jsi")
    }
}

task DownloadNirvanaReferences {
    input {
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o xtrace -o pipefail

        mkdir references

        /usr/bin/dotnet /Nirvana/Downloader.dll --ga GRCh38 -o references
        tar cfz Nirvana-references-2022-10-07.tgz references
    >>>
    output {
        File references = "Nirvana-references-2022-10-07.tgz"
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:vs_500_nirvana_3_18_1"
        disks: "local-disk 1000 HDD"
        memory: "3 GiB"
    }
}
