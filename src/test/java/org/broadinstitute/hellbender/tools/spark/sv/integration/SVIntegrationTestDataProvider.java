package org.broadinstitute.hellbender.tools.spark.sv.integration;


import org.broadinstitute.hellbender.CommandLineProgramTest;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

class SVIntegrationTestDataProvider extends CommandLineProgramTest {

    static final File reference = new File(b37_reference_20_21);
    static final File reference_2bit = new File(b37_2bit_reference_20_21);
    static final File reference_fai = new File(b37_reference_20_21+".fai");
    static final File reference_dict = new File(b37_reference_20_21.replace(".fasta", ".dict"));

    private static final String THIS_TEST_FOLDER = getTestDataDir() + "/spark/sv/integration";
    static final String TEST_BAM_LEFT = THIS_TEST_FOLDER + "/hg19_DEL_leftReads.bam";
    static final String TEST_BAM_RIGHT = THIS_TEST_FOLDER + "/hg19_DEL_rightReads.bam";
    static final String KMER_KILL_LIST = THIS_TEST_FOLDER + "/dummy.kill.kmers";
    static final String ALIGNER_INDEX_IMG = largeFileTestDir + "human_g1k_v37.20.21.fasta.img";
    static final String TEST_CONTIG_SAM = THIS_TEST_FOLDER + "/hg19_DEL_contigAssemblies.sam";

    static final String EXPECTED_SIMPLE_DEL_VCF = THIS_TEST_FOLDER + "/hg19_DEL.vcf";
    static final String EXPECTED_SIMPLE_INV_VCF = THIS_TEST_FOLDER + "/hg19_INV.vcf";

    @SuppressWarnings("unchecked")
    static final List<String> dummyExpectedFileNames = new ArrayList<>(Collections.EMPTY_LIST);
}
