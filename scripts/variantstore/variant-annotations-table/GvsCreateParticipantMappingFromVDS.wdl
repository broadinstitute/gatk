version 1.0

# Participant mapping built from the VDS (VS-2013). hail_create_participant_mapping_inputs.py writes the carrier pairs
# and the VDS keys that are not left-aligned; this workflow gives those keys their VAT VIDs.
import "../wdl/GvsUtils.wdl" as Utils

workflow GvsCreateParticipantMappingFromVDS {
    input {
        File keys_to_normalize_vcf
        String reference_name = "hg38"

        String? git_branch_or_tag
        String? basic_docker
        String? variants_docker
    }

    parameter_meta {
        keys_to_normalize_vcf: {
            help: "keys_to_normalize.vcf.bgz from hail_create_participant_mapping_inputs.py: sites-only, the VDS key in ID."
        }
    }

    if (!defined(basic_docker) || !defined(variants_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = effective_basic_docker,
    }

    call NormalizeKeys {
        input:
            keys_to_normalize_vcf = keys_to_normalize_vcf,
            ref = GetReference.reference.reference_fasta,
            variants_docker = effective_variants_docker,
    }

    output {
        File key_to_vid = NormalizeKeys.key_to_vid
    }
}

task NormalizeKeys {
    input {
        File keys_to_normalize_vcf
        File ref
        String variants_docker
    }

    Int disk_size = ceil(size(ref, "GB") * 2 + size(keys_to_normalize_vcf, "GB") * 5) + 20

    # Writes key_to_vid.tsv, one headerless `vds_key<TAB>vid` row per input key.
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # The same normalization RemoveDuplicatesFromSitesOnlyVCF in GvsCreateVATfromVDS.wdl applies, against the same
        # reference, so each key gets the VID the VAT gave its allele. ID passes through untouched.
        bcftools norm -m- --check-ref w -f ~{ref} ~{keys_to_normalize_vcf} -O v -o normalized.vcf

        # VID format: contig without the chr prefix.
        bcftools query -f '%ID\t%CHROM-%POS-%REF-%ALT\n' normalized.vcf \
            | awk -F'\t' -v OFS='\t' '{ sub(/^chr/, "", $2); print }' > key_to_vid.tsv

        # Every key is biallelic, so norm -m- neither splits nor drops; anything else means the input is not what the
        # Hail script writes.
        input_keys=$(bcftools view -H ~{keys_to_normalize_vcf} | wc -l)
        output_rows=$(wc -l < key_to_vid.tsv)
        distinct_keys=$(cut -f1 key_to_vid.tsv | sort -u | wc -l)
        if [[ ${input_keys} -ne ${output_rows} || ${input_keys} -ne ${distinct_keys} ]]
        then
            echo "Expected one row per key: ${input_keys} keys in, ${output_rows} rows and ${distinct_keys} distinct keys out." >&2
            exit 1
        fi

        # The Hail script flags a key only when its indel can move left, so every one must come out at a new VID.
        unmoved=$(awk -F'\t' '$1 == $2' key_to_vid.tsv | wc -l)
        if [[ ${unmoved} -ne 0 ]]
        then
            echo "${unmoved} flagged keys normalized to themselves, e.g.:" >&2
            awk -F'\t' '$1 == $2' key_to_vid.tsv | head -5 >&2
            exit 1
        fi
    >>>

    runtime {
        docker: variants_docker
        maxRetries: 3
        memory: "4 GB"
        preemptible: 3
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File key_to_vid = "key_to_vid.tsv"
    }
}
