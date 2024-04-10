version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsTieoutPgenToVcf {
    input {
        String extract_task_call_root
        String merge_pgen_workflow_root
    }

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
        File? none = ""
    }

    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = none,
    }

    call Tieout {
        input:
            gatk_docker = GetToolVersions.gatk_docker,
            extract_task_call_root = extract_task_call_root,
            merge_pgen_workflow_root = merge_pgen_workflow_root,
    }
}

task Tieout {
    input {
        String gatk_docker
        String extract_task_call_root
        String merge_pgen_workflow_root
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        mkdir tieout
        cd tieout

        # Localize all the sharded VCFs and their indexes.
        mkdir vcf
        cd vcf
        gcloud storage cp -R '~{extract_task_call_root}/**/*-quickit.vcf.gz' .
        gcloud storage cp -R '~{extract_task_call_root}/**/*-quickit.vcf.gz.tbi' .

        # Get the names of the samples, sort them, and stash this in a file in the "tieout" directory.
        # PGEN appears to order samples lexicographically, while VCF extract orders them randomly.
        # We will need to explicitly order the PGEN and VCF data similarly for our comparisons.
        bcftools query --list-samples 0000000000-quickit.vcf.gz | sort > ../samples.txt

        # Process the VCFs:
        # - Write an appropriate header
        # - Group by chromosome
        # - Include only the data we want to use in our PGEN comparison
        # - Assemble everything into a VCF
        # - bgzip

        # Write an appropriate header
        # Take all the ## header lines that aren't contigs
        zgrep -E '^##' 0000000000-quickit.vcf.gz | grep -E -v '^##contig=' > vcf_header.txt

        # Now just the contigs we care about
        zgrep -E '^##' 0000000000-quickit.vcf.gz | grep -E '^##contig=<ID=chr[12]?[0-9],|^##contig=<ID=chr[XY],' >> vcf_header.txt

        # The TSV header minus the sample names
        echo -n '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT' >> vcf_header.txt

        # The sample names in the correct order
        for sample in $(cat ../samples.txt
        do
            echo -n "\t${sample}" >> vcf_header.txt
        done

        # Finally a newline
        echo >> vcf_header.txt

        # Group by chromosome
        for vcf in 0000000*quickit.vcf.gz
        do
            echo -n "$vcf,"
            zgrep -E -v '^#' $vcf | head -1 | awk '{print $1}'
        done > vcf_to_chr.csv

        while IFS=, read -r vcf chr
        do
            echo "vcf: $vcf / chr: $chr"
            # Include only the data we want to use in our PGEN comparison
            bcftools query -S ../samples.txt -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT[\t%GT]\n' $vcf >> vcf_body_${chr}.vcf
        done < vcf_to_chr.csv

        # Assemble everything into a VCF
        # bgzip
        for chr in $(seq 1 22) X Y
        do
            out=vcf_compare_chr${chr}.vcf
            cat vcf_header.txt > $out
            cat vcf_body_chr${chr}.vcf >> $out
            bgzip $out
        done

        cd ../..
        tar cfz tieout.tgz tieout
    >>>
    output {
        File out = "tieout.tgz"
    }
    runtime {
        docker: gatk_docker
    }
}
