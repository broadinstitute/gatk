version 1.0

import "GvsCreateVATfromVDS.wdl" as GvsCreateVATfromVDS

workflow RemoveDuplicatesFromSitesOnlyVCF {
    input {
        File sites_only_vcf
        File ref
        String variants_docker
    }

    call GvsCreateVATfromVDS.RemoveDuplicatesFromSitesOnlyVCF {
        input:
            sites_only_vcf = sites_only_vcf,
            ref = ref,
            variants_docker = variants_docker
    }

    output {
        File track_dropped = RemoveDuplicatesFromSitesOnlyVCF.track_dropped
        File output_vcf = RemoveDuplicatesFromSitesOnlyVCF.output_vcf
        File monitoring_log = RemoveDuplicatesFromSitesOnlyVCF.monitoring_log
        File filtered_synonyms = RemoveDuplicatesFromSitesOnlyVCF.filtered_synonyms
    }
}
