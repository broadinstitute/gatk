package org.broadinstitute.hellbender.tools.spark.sv.integration;


import org.broadinstitute.hellbender.CommandLineProgramTest;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SVIntegrationTestDataProvider extends CommandLineProgramTest {

    public static final String THIS_TEST_FOLDER = getTestDataDir() + "/spark/sv/integration";

    public static final File reference = new File(b37_reference_20_21);
    public static final File reference_2bit = new File(b37_2bit_reference_20_21);
    public static final File reference_fai = new File(b37_reference_20_21+".fai");
    public static final File reference_dict = new File(b37_reference_20_21.replace(".fasta", ".dict"));
    public static final String ALIGNER_INDEX_IMG = largeFileTestDir + "human_g1k_v37.20.21.fasta.img";

    // inputs to tests
    public static final String TEST_BAM = "src/test/resources/large/SVIntegrationTest.bam";
    public static final String TEST_CONTIG_SAM = THIS_TEST_FOLDER + "/inputs/hg19_DEL_contigAssemblies.sam";
    public static final String EXTERNAL_CNV_CALLS = THIS_TEST_FOLDER + "/inputs/hg19_DEL_cnv_calls.vcf";
    public static final String KMER_KILL_LIST = THIS_TEST_FOLDER + "/inputs/dummy.kill.kmers";


    // expected outputs
    public static final String EXPECTED_SIMPLE_DEL_VCF = THIS_TEST_FOLDER + "/outputs/hg19_DEL.vcf";
    public static final String EXPECTED_SIMPLE_INV_VCF = THIS_TEST_FOLDER + "/outputs/hg19_INV.vcf";
    public static final String EXPECTED_ALIGNED_CONTIGS = THIS_TEST_FOLDER + "/outputs/hg19_DEL_actualContigAssemblies.sam";
}
