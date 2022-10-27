version 1.0

import "GvsUnifiedHail.wdl" as Unified
import "GvsUtils.wdl" as Utils

workflow GvsQuickstartHailIntegration {

    input {
        String branch_name
        String expected_output_prefix = "gs://gvs-internal-quickstart/integration/2022-07-05/"

        Array[String] external_sample_names = [
                                              "ERS4367795",
                                              "ERS4367796",
                                              "ERS4367797",
                                              "ERS4367798",
                                              "ERS4367799",
                                              "ERS4367800",
                                              "ERS4367801",
                                              "ERS4367803",
                                              "ERS4367804",
                                              "ERS4367805"
                                              ]

        Array[File] input_vcfs = [
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00420.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00423.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00427.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00429.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00444.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00447.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00450.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz"
                                 ]

        Array[File] input_vcf_indexes = [
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00420.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00423.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00427.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00429.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00444.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00447.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00450.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi"
                                        ]

        String drop_state = "FORTY"

        Array[File] tieout_vcfs = [
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-0/0000000000-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-1/0000000001-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-2/0000000002-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-3/0000000003-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-4/0000000004-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-5/0000000005-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-6/0000000006-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-7/0000000007-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-8/0000000008-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-9/0000000009-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-10/0000000010-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-11/0000000011-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-12/0000000012-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-13/0000000013-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-14/0000000014-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-15/0000000015-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-16/0000000016-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-17/0000000017-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-18/0000000018-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-19/0000000019-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-20/0000000020-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-21/0000000021-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-22/0000000022-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-23/0000000023-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-24/0000000024-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-25/0000000025-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-26/0000000026-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-27/0000000027-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-28/0000000028-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-29/attempt-2/0000000029-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-30/0000000030-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-31/0000000031-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-32/0000000032-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-33/0000000033-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-34/0000000034-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-35/0000000035-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-36/0000000036-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-37/0000000037-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-38/0000000038-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-39/0000000039-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-40/0000000040-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-41/0000000041-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-42/0000000042-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-43/0000000043-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-44/0000000044-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-45/0000000045-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-46/0000000046-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-47/0000000047-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-48/0000000048-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-49/0000000049-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-50/0000000050-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-51/0000000051-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-52/0000000052-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-53/0000000053-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-54/0000000054-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-55/0000000055-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-56/0000000056-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-57/0000000057-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-58/0000000058-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-59/0000000059-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-60/0000000060-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-61/0000000061-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-62/0000000062-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-63/0000000063-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-64/0000000064-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-65/0000000065-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-66/0000000066-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-67/0000000067-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-68/0000000068-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-69/0000000069-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-70/0000000070-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-71/0000000071-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-72/0000000072-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-73/0000000073-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-74/0000000074-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-75/0000000075-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-76/0000000076-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-77/0000000077-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-78/0000000078-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-79/0000000079-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-80/0000000080-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-81/0000000081-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-82/0000000082-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-83/0000000083-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-84/0000000084-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-85/0000000085-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-86/0000000086-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-87/0000000087-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-88/0000000088-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-89/0000000089-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-90/0000000090-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-91/0000000091-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-92/0000000092-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-93/0000000093-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-94/0000000094-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-95/0000000095-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-96/0000000096-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-97/0000000097-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-98/0000000098-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-99/0000000099-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-100/0000000100-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-101/0000000101-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-102/0000000102-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-103/0000000103-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-104/0000000104-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-105/0000000105-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-106/0000000106-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-107/0000000107-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-108/0000000108-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-109/0000000109-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-110/0000000110-quickit.vcf.gz",
                                  "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-111/0000000111-quickit.vcf.gz"
                                  ]

        Array[File] tieout_vcf_indexes = [
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-0/0000000000-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-1/0000000001-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-2/0000000002-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-3/0000000003-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-4/0000000004-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-5/0000000005-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-6/0000000006-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-7/0000000007-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-8/0000000008-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-9/0000000009-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-10/0000000010-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-11/0000000011-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-12/0000000012-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-13/0000000013-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-14/0000000014-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-15/0000000015-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-16/0000000016-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-17/0000000017-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-18/0000000018-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-19/0000000019-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-20/0000000020-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-21/0000000021-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-22/0000000022-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-23/0000000023-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-24/0000000024-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-25/0000000025-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-26/0000000026-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-27/0000000027-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-28/0000000028-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-29/attempt-2/0000000029-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-30/0000000030-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-31/0000000031-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-32/0000000032-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-33/0000000033-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-34/0000000034-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-35/0000000035-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-36/0000000036-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-37/0000000037-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-38/0000000038-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-39/0000000039-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-40/0000000040-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-41/0000000041-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-42/0000000042-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-43/0000000043-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-44/0000000044-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-45/0000000045-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-46/0000000046-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-47/0000000047-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-48/0000000048-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-49/0000000049-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-50/0000000050-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-51/0000000051-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-52/0000000052-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-53/0000000053-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-54/0000000054-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-55/0000000055-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-56/0000000056-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-57/0000000057-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-58/0000000058-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-59/0000000059-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-60/0000000060-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-61/0000000061-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-62/0000000062-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-63/0000000063-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-64/0000000064-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-65/0000000065-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-66/0000000066-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-67/0000000067-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-68/0000000068-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-69/0000000069-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-70/0000000070-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-71/0000000071-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-72/0000000072-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-73/0000000073-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-74/0000000074-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-75/0000000075-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-76/0000000076-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-77/0000000077-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-78/0000000078-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-79/0000000079-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-80/0000000080-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-81/0000000081-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-82/0000000082-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-83/0000000083-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-84/0000000084-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-85/0000000085-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-86/0000000086-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-87/0000000087-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-88/0000000088-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-89/0000000089-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-90/0000000090-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-91/0000000091-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-92/0000000092-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-93/0000000093-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-94/0000000094-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-95/0000000095-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-96/0000000096-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-97/0000000097-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-98/0000000098-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-99/0000000099-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-100/0000000100-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-101/0000000101-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-102/0000000102-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-103/0000000103-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-104/0000000104-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-105/0000000105-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-106/0000000106-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-107/0000000107-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-108/0000000108-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-109/0000000109-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-110/0000000110-quickit.vcf.gz.tbi",
                                         "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/1cfae8ea-2e3b-4ba2-b6db-802349ef69bc/GvsQuickstartIntegration/b270ddee-1bc8-4488-bc60-79d0cb82797e/call-GvsUnified/GvsUnified/1a15df10-9020-4c50-aa61-55f241ebb874/call-GvsExtractCallset/GvsExtractCallset/ce4cced0-9863-4887-8a35-a85c8bf7b89a/call-ExtractTask/shard-111/0000000111-quickit.vcf.gz.tbi"
                                         ]
    }
    String project_id = "gvs-internal"

#    call Utils.BuildGATKJarAndCreateDataset {
#        input:
#            branch_name = branch_name,
#            dataset_prefix = "quickit"
#    }

    call Unified.GvsUnifiedHail {
        input:
            call_set_identifier = branch_name,
            dataset_name = "quickit_vs_639_hail_testing_spike_80af46a",
            project_id = project_id,
            external_sample_names = external_sample_names,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            filter_set_name = "quickit",
            drop_state = "NONE",
            tieout_vcfs = tieout_vcfs,
            tieout_vcf_indexes = tieout_vcf_indexes,
    }

#    call AssertIdenticalOutputs {
#        input:
#            expected_output_prefix = expected_output_prefix,
#            actual_vcfs = GvsUnified.output_vcfs
#    }
#
#    call AssertCostIsTrackedAndExpected {
#        input:
#            go = GvsUnified.done,
#            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
#            project_id = project_id,
#            expected_output_csv = expected_output_prefix + "cost_observability_expected.csv"
#    }
#
#    call AssertTableSizesAreExpected {
#        input:
#            go = GvsUnified.done,
#            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
#            project_id = project_id,
#            expected_output_csv = expected_output_prefix + "table_sizes_expected.csv"
#    }
#
    output {
        # String vds_output_path = GvsUnifiedHail.vds_output_path
        Boolean done = true
    }
}


task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        Array[File] actual_vcfs
    }

    command <<<
        set -o errexit
        set -o nounset
        set -o pipefail

        failures=()

        # Where the current set of expected results lives in the cloud
        expected_prefix="~{expected_output_prefix}"
        # Remove a trailing slash if there is one
        expected_prefix=${expected_prefix%/}

        # Download all the expected data
        mkdir expected
        cd expected
        # Make a FOFN for more efficient downloading
        for file in ~{sep= ' ' actual_vcfs}; do
            echo "$expected_prefix/$(basename $file)" >> expected_fofn.txt
        done
        # Download and unzip all the expected data
        cat expected_fofn.txt | gsutil -m cp -I .
        gzip -d *.gz
        cd ..

        # Also unzip actual result data
        gzip -d ~{sep= ' ' actual_vcfs}

        # Headers first, these can yield useful diagnostics when there are mismatches.
        for file in ~{sep=' ' actual_vcfs}; do
          unzipped=${file%.gz}
          expected="expected/$(basename $unzipped)"
          set +o errexit
          cmp <(grep '^#' $unzipped) <(grep '^#' $expected)
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
            # If there is a mismatch add it to a list of failures but keep on looking for mismatches.
            failures+=( $unzipped )
          fi
        done

        if [[ ${#failures[@]} -ne 0 ]]; then
          echo "Error: headers for the following files do not match:"
          for failure in ${failures[@]}; do
            echo $(basename $failure)
            expected="expected/$(basename $failure)"
            diff $failure $expected
          done
          exit 1
        fi

        # If the headers all matched look for any mismatches in overall file content.
        fail=0
        for file in ~{sep=' ' actual_vcfs}; do
          unzipped=${file%.gz}
          expected="expected/$(basename $unzipped)"
          set +o errexit
          cmp $unzipped $expected
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
            echo "Error: file contents of expected and actual do not match: $(basename $unzipped)"
            fail=1
          fi
        done

        if [[ $fail -ne 0 ]]; then
          exit 1
        fi

        echo "All vcfs compared and matched!"
    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:402.0.0-alpine"
        disks: "local-disk 500 HDD"
    }

    output {
        File fofn = "expected/expected_fofn.txt"
        Boolean done = true
    }
}

task AssertCostIsTrackedAndExpected {
    meta {
        # we want to check the database each time this runs
        volatile: true
    }

    input {
        Boolean go = true
        String dataset_name
        String project_id
        File expected_output_csv
    }

    command <<<
        set -o errexit
        set -o nounset
        set -o pipefail

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
            "SELECT call, step, event_key, sum(event_bytes) \
              FROM \`~{dataset_name}.cost_observability\` \
              GROUP BY call, step, event_key \
              ORDER BY call, step, event_key" > cost_observability_output.csv

        # Put the exit code in a file because we are using a subshell (while) later and changes to the variable *in* the subshell are lost
        echo "0" > ret_val.txt

        paste cost_observability_output.csv ~{expected_output_csv} | while  IFS=$'\t' read observed expected; do
        IFS=, read -ra OBS <<< "$observed"
        IFS=, read -ra EXP <<< "$expected"
        if [[ "${#OBS[@]}" -ne 4  || "${#EXP[@]}" -ne 4 ]]; then
          echo "Unexpected number of rows found in the input files"
          exit 1
        fi

        OBS_KEY=${OBS[0]}.${OBS[1]}.${OBS[2]}
        EXP_KEY=${EXP[0]}.${EXP[1]}.${EXP[2]}
        if [[ "$OBS_KEY" != "$EXP_KEY" ]]; then
          echo "Mismatched keys in results files - were these sorted properly?"
          exit 1
        fi

        if [[ "$OBS_KEY" == "call.step.event_key" ]]; then
          # Skip the header
          continue;
        fi

        OBS_BYTES=${OBS[3]}
        EXP_BYTES=${EXP[3]}

        TOLERANCE=0

        # For these two costs, there is non-determinism in the pipeline - we allow a % difference
        if [[ $OBS_KEY == "ExtractFilterTask.GvsCreateFilterSet.BigQuery Query Scanned" ]]; then
          TOLERANCE=0.015   # 1.5% tolerance  (Note - have seen as high as: 0.0109429)
        elif [[ $OBS_KEY == "ExtractTask.GvsCreateCallset.Storage API Scanned" ]]; then
          TOLERANCE=0.01   # 1% tolerance  (Note - have seen as high as: 0.00608656)
        elif [[ $OBS_KEY == "ExtractFilterTask.GvsCreateFilterSet.Storage API Scanned" ]]; then
          TOLERANCE=0.05  # 5% tolerance (Note - have seen as high as: 0.0281223)
        fi

        if [[ $OBS_BYTES -ne $EXP_BYTES ]]; then
          echo "The bytes observed ($OBS_BYTES) for '$OBS_KEY' differ from those expected ($EXP_BYTES)"

          if [[ $OBS_BYTES -ge $EXP_BYTES ]]; then
            DIFF_FOUND=$(echo $OBS_BYTES $EXP_BYTES | awk '{print ($1-$2)/$1}')
          else
            DIFF_FOUND=$(echo $EXP_BYTES $OBS_BYTES | awk '{print ($1-$2)/$1}')
          fi

          if ! awk "BEGIN{ exit ($DIFF_FOUND > $TOLERANCE) }"
          then
            echo "FAIL!!! The relative difference between these is $DIFF_FOUND, which is greater than the allowed tolerance ($TOLERANCE)"
            echo "1" > ret_val.txt
          else
            echo "However, the relative difference between these is $DIFF_FOUND, which is below the allowed tolerance ($TOLERANCE)"
          fi
        fi
        done

        RET_VAL=`cat ret_val.txt`
        exit $RET_VAL

    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:402.0.0-alpine"
        disks: "local-disk 10 HDD"
    }

    output {
        File cost_observability_output_csv = "cost_observability_output.csv"
    }
}

task AssertTableSizesAreExpected {
    meta {
        # we want to check the database each time this runs
        volatile: true
    }

    input {
        Boolean go = true
        String dataset_name
        String project_id
        File expected_output_csv
    }

    command <<<
        set -o errexit
        set -o nounset
        set -o pipefail

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
            "SELECT 'vet_total' AS total_name, sum(total_billable_bytes) AS total_bytes FROM \
            \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` WHERE table_name LIKE 'vet_%' \
            UNION ALL \
            SELECT 'ref_ranges_total' AS total_name, sum(total_billable_bytes) AS total_bytes \
            FROM \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` \
            WHERE table_name LIKE 'ref_ranges_%' ORDER BY total_name" > table_size_output.csv

        set +o errexit
        diff -w table_size_output.csv ~{expected_output_csv} > differences.txt
        set -o errexit

        if [[ -s differences.txt ]]; then
            echo "Differences found:"
            cat differences.txt
            exit 1
        fi
    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:402.0.0-alpine"
        disks: "local-disk 10 HDD"
    }

    output {
        File table_size_output_csv = "table_size_output.csv"
        File differences = "differences.txt"
    }
}
