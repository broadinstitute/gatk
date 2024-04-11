version 1.0


workflow GvsTieoutPgenToVcf {
    input {
        String extract_task_call_root
        String merge_pgen_workflow_root
    }

    call Tieout {
        input:
            extract_task_call_root = extract_task_call_root,
            merge_pgen_workflow_root = merge_pgen_workflow_root,
    }
}

task Tieout {
    input {
        String extract_task_call_root
        String merge_pgen_workflow_root
        String plink_download_url = "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240318.zip"
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        mkdir tieout
        cd tieout

        # Localize all the sharded VCFs and their indexes.
        # Note everything here is specific to 'quickit' but probably shouldn't be.
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
        for sample in $(cat ../samples.txt)
        do
            # literal tab below as \t is not working in WDL though it does work locally
            echo -n "   ${sample}" >> vcf_header.txt
        done

        # Finally a newline
        echo >> vcf_header.txt

        # Group by chromosome
        # Expecting some SIGPIPEs with the `zgrep` trying to push data to a `head -1` that stopped listening, so turn off pipefail for now.
        set +o pipefail
        for vcf in 0000000*quickit.vcf.gz
        do
            echo -n "$vcf,"
            zgrep -E -v '^#' $vcf | head -1 | awk '{print $1}'
        done > vcf_to_chr.csv
        set -o pipefail

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

        mkdir pgen
        cd pgen
        for kind in pgen psam pvar.zst
        do
            # Note the 'kind' var is intentionally outside the single quotes because we do not want the globs expanded
            # but we do want the 'kind' var expanded.
            gcloud storage cp -R '~{merge_pgen_workflow_root}/**/call-FinalMerge/quickit.chr*.'${kind} .
        done

        # Download and "install" plink2. Building this on an Apple Silicon Mac was completely straightforward, but
        # the GATK Docker image was another story...
        curl -O --location '~{plink_download_url}'
        unzip $(basename '~{plink_download_url}')
        mv plink2 /usr/bin

        for chr in $(seq 1 22) X Y
        do
            plink2 export --pfile quickit.chr${chr} --export vcf-4.2 --out pgen_compare_chr${chr}
        done

        cd ..

        # compare
        mkdir compare
        cd compare
        for chr in $(seq 1 22) X Y
        do
            bcftools isec vcf/vcf_compare_chr${chr}.vcf pgen/pgen_compare_chr${chr}.vcf -p chr${chr}
        done

        cd ../..
        tar cfz tieout.tgz tieout
    >>>
    output {
        File out = "tieout.tgz"
    }
    runtime {
        docker: "broadinstitute/gatk:4.5.0.0"
    }
}