package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;

/**
 * A class to hold the constants for the Funcotator Tests.
 * Created by jonn on 11/1/17.
 */
public class FuncotatorTestConstants {

    public static final String DUMMY_DATA_SOURCES_TAR_GZ             = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "dummyDataSources.tar.gz";
    public static final String DUMMY_DATA_SOURCES_TAR_GZ_SHA256_FILE = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "dummyDataSources.sha256";
    public static final String DUMMY_DATA_SOURCES_FOLDER             = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "dummyDataSources";

    public static final double FUNCOTATOR_DOUBLE_COMPARISON_EPSILON = 0.0001;

    public static final String FUNCOTATOR_TEST_DIR = GATKBaseTest.toolsTestDir + "funcotator" + File.separator;

    public static final String HG19_CHR19_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "GRCh37.p13.chr19.tar.gz";
    public static final String HG19_CHR3_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "GRCh37.p13.chr3.tar.gz";
    public static final String HG19_3_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "b37.3.tar.gz";
    public static final String HG38_3_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "hg38.3.tar.gz";
    public static final String HG19_2_REFERENCE_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "b37.2.tar.gz";

    public static final String FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "funcotator_dataSources" + File.separator;

    public static final String REFERENCE_VERSION_HG19 = "hg19";
    public static final String REFERENCE_VERSION_HG38 = "hg38";

    public static final SimpleInterval MUC16_POSITION = new SimpleInterval("chr19", 8959520, 9092018);
    public static final String MUC16_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode_muc16" + File.separator + "hg19" + File.separator + "gencode.v19.MUC16.gtf";
    public static final String MUC16_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode_muc16" + File.separator + "hg19" + File.separator + "gencode.v19.MUC16_transcript.fasta";
    public static final String MUC16_TRANSCRIPT = "ENST00000397910.4";
    public static final String MUC16_PATHOLOGICAL_TRANSCRIPT = "ENST00000599436.1";
    public static final String MUC16_GENCODE_NON_BASIC_ANNOTATIONS_FILE_NAME = FUNCOTATOR_TEST_DIR + "gencode.v19.MUC16.non-basic.gtf";

    public static final SimpleInterval PIK3CA_POSITION = new SimpleInterval("chr3", 178866311, 178957881);
    public static final String PIK3CA_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode_pik3ca" + File.separator + "hg19" + File.separator + "gencode.v19.PIK3CA.gtf";
    public static final String PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode_pik3ca" + File.separator + "hg19" + File.separator + "gencode.v19.PIK3CA_transcript.fasta";
    public static final String PIK3CA_TRANSCRIPT = "ENST00000263967.3";

    public static final String GTF_CHR3_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "gencode.v19.chr_patch_hapl_scaff.chr3.gtf";

    public static final String GENCODE_TRANSCRIPT_FASTA_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "gencode.v19.pc_transcripts.fasta";

    public static final String VARIANT_FILE_HG19_CHR3 = FUNCOTATOR_TEST_DIR + "snpTest_chr3_hg19.vcf";
    public static final String VARIANT_FILE_HG19_CHR19 = FUNCOTATOR_TEST_DIR + "snpTest_chr19_hg19.vcf";

    // ----------------------------------------------------------------------
    // Data source variables:

    public static final String XSV_CSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_testFile.csv";

    public static final String XSV_TSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_TSV_testFile.csv";
    public static final String XSV_PIPESV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_PIPESV_testFile.xsv";
    public static final String XSV_DEADBEEFSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_DEADBEEFSV_testFile.csv";
    public static final String XSV_CSV_PIK3CA_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_PIK3CA.csv";
    public static final String XSV_CSV_MUC16_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_MUC16.csv";

    public static final String XSV_LOCATABLE_TEST_FILE1_PATH = FUNCOTATOR_TEST_DIR + "xsv_locatable_test.csv";
    public static final String XSV_LOCATABLE_TEST_FILE2_PATH = FUNCOTATOR_TEST_DIR + "xsv_locatable_test2.csv";
    public static final String XSV_LOCATABLE_TEST_FILE3_PATH = FUNCOTATOR_TEST_DIR + "xsv_locatable_test3.tsv";

    public static final String COSMIC_TEST_DB = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "cosmic" + File.separator + "hg19" + File.separator + "CosmicTest.db";

    public static final String DBSNP_HG19_SNIPPET_FILE_PATH = FUNCOTATOR_TEST_DIR + "dbSNP_hg19_snippet.vcf";

    // ----------------------------------------------------------------------
    // Integration Test Variables:
    public static final String MUC16_MNP_FILE_BASE_NAME    = FUNCOTATOR_TEST_DIR + "MUC16_MNP";
    public static final String PIK3CA_SNP_FILE_BASE_NAME   = FUNCOTATOR_TEST_DIR + "PIK3CA_SNPS";
    public static final String PIK3CA_INDEL_FILE_BASE_NAME = FUNCOTATOR_TEST_DIR + "PIK3CA_INDELS";
}
