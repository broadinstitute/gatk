package org.broadinstitute.hellbender.tools.funcotator;

import java.io.File;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;

/**
 * A class to hold the constants for the Funcotator Tests.
 * Created by jonn on 11/1/17.
 */
public class FuncotatorTestConstants {
    public static final String FUNCOTATOR_TEST_DIR = GATKBaseTest.toolsTestDir + "funcotator" + File.separator;

    public static final String HG19_CHR19_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "GRCh37.p13.chr19.fasta";
    public static final String HG19_CHR3_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "GRCh37.p13.chr3.fasta";

    public static final String MUC16_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_TEST_DIR + "gencode.v19.MUC16.gtf";
    public static final String MUC16_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_TEST_DIR + "gencode.v19.MUC16_transcript.fasta";
    public static final String MUC_16_TRANSCRIPT = "ENST00000397910.4";

    public static final String PIK3CA_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_TEST_DIR + "gencode.v19.PIK3CA.gtf";
    public static final String PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_TEST_DIR + "gencode.v19.PIK3CA_transcript.fasta";

    public static final String GTF_CHR3_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "gencode.v19.chr_patch_hapl_scaff.chr3.gtf";

    public static final String GENCODE_TRANSCRIPT_FASTA_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "gencode.v19.pc_transcripts.fasta";

    public static final String VARIANT_FILE_HG19_CHR3 = FUNCOTATOR_TEST_DIR + "singleSnpTest_chr3_hg19.vcf";
}
