version 1.0

import "GvsCreateVATFromVDS.wdl"

workflow repro {
    File input_vcf = "gs://fc-34a42a2b-29aa-43f8-a6dd-2a50f7f3408b/submissions/bef81325-9321-4b4f-a741-b5d151122ceb/GvsCreateVATfromVDS/35a55d99-84e3-4cd4-ad57-1929d99f350c/call-StripCustomAnnotationsFromSitesOnlyVCF/shard-44/0000000044-scattered.sites-only-vcf-c05be276-20b4.unannotated.sites_only.vcf"
    String output_annotated_file_name = "0000000044-scattered.sites-only-vcf-c05be276-20b4_annotated"
    File custom_annotations_file = "gs://fc-34a42a2b-29aa-43f8-a6dd-2a50f7f3408b/submissions/bef81325-9321-4b4f-a741-b5d151122ceb/GvsCreateVATfromVDS/35a55d99-84e3-4cd4-ad57-1929d99f350c/call-StripCustomAnnotationsFromSitesOnlyVCF/shard-44/0000000044-scattered.sites-only-vcf-c05be276-20b4.custom_annotations.tsv"
    String cromwell_root = "/mnt/disks/cromwell_root"

    call GvsCreateVATFromVDS.AnnotateVCF {
        input:
            input_vcf = input_vcf,
            output_annotated_file_name = output_annotated_file_name,
            custom_annotations_file = custom_annotations_file,
            cromwell_root = cromwell_root
    }

    output {
        Boolean done = true
    }
}



