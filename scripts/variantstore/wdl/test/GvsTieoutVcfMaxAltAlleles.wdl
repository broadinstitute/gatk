version 1.0

import "../GvsUtils.wdl" as Utils

workflow GvsTieoutVcfMaxAltAlleles {
    input {
        String unfiltered_output_gcs_path
        String filtered_output_gcs_path
        Int max_alt_alleles
        String? variants_docker
    }

    if (!defined(variants_docker)) {
        call Utils.GetToolVersions
    }

    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])

    call Tieout {
        input:
            unfiltered_output_gcs_path = unfiltered_output_gcs_path,
            filtered_output_gcs_path = filtered_output_gcs_path,
            max_alt_alleles = max_alt_alleles,
            variants_docker = effective_variants_docker,
    }
}

task Tieout {
    input {
        String unfiltered_output_gcs_path
        String filtered_output_gcs_path
        Int max_alt_alleles
        String variants_docker
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace
        . /localvenv/bin/activate

        mkdir unfiltered filtered compare

        cd unfiltered
        gcloud storage cp '~{unfiltered_output_gcs_path}/*' .

        cd ..
        cd filtered
        gcloud storage cp '~{filtered_output_gcs_path}/*' .

        cd ..

        for unfiltered in unfiltered/*.vcf.gz
        do
          base="$(basename $unfiltered)"
          filtered="filtered/${base}"
          out=${base%-*}
          bcftools isec $filtered $unfiltered -p compare/${out}
        done

        # Identify the 'isec' output files that represent entries which are "private" to either VCF file being compared.
        find . -name README.txt | xargs grep -h private | awk '{print $1}' > private_files.txt

        # Turn off errexit, the grep below is expected to return non-zero
        set +o errexit

        # Find all the private files, print out all the ALTs, look for any lines that do *not* have at least two commas.
        # All ALT entries are expected to have two or more commas because the maximum number of ALT alleles was specified
        # as two. So the unfiltered file would have entries like ALT1,ALT2,...,ALTN while the filtered file would not,
        # and these should be the only differences between the two files.
        cat private_files.txt | xargs -n 1 bcftools query -f '%ALT\n' | grep -E -v ',.*,' > unexpected.txt

        rc=$?
        if [[ $rc -eq 0 ]]
        then
            echo "Unexpected ALTs found, see 'unexpected' output for details." 1>&2
            exit 1
        fi

        set -o errexit

        # Similar to above, except making sure that we did observe some filtered sites. At the time of this
        # writing there are 4001 sites filtered for > 2 entries; conservatively setting the threshold at 1000.
        filtered_count=$(cat private_files.txt | xargs -n 1 bcftools query -f '%ALT\n' | grep -E ',.*,' | wc -l)
        if [[ $filtered_count -lt 1000 ]]
        then
            echo "Unexpectedly found fewer than 1000 filtered ALTs, see 'tarball' output for details." 1>&2
            exit 1
        fi

        tar cfvz tieout.tgz compare
    >>>
    output {
        File tarball = "tieout.tgz"
        File unexpected = "unexpected.txt"
    }
    runtime {
        docker: variants_docker
    }
}
