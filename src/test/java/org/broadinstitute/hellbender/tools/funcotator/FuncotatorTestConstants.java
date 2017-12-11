package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.GATKBaseTest;

import java.io.File;

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
    public static final String MUC16_TRANSCRIPT = "ENST00000397910.4";

    public static final String PIK3CA_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_TEST_DIR + "gencode.v19.PIK3CA.gtf";
    public static final String PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_TEST_DIR + "gencode.v19.PIK3CA_transcript.fasta";
    public static final String PIK3CA_TRANSCRIPT = "ENST00000263967.3";

    public static final String GTF_CHR3_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "gencode.v19.chr_patch_hapl_scaff.chr3.gtf";

    public static final String GENCODE_TRANSCRIPT_FASTA_FILE_NAME = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "gencode.v19.pc_transcripts.fasta";

    public static final String VARIANT_FILE_HG19_CHR3 = FUNCOTATOR_TEST_DIR + "snpTest_chr3_hg19.vcf";
    public static final String VARIANT_FILE_HG19_CHR19 = FUNCOTATOR_TEST_DIR + "snpTest_chr19_hg19.vcf";

    // ----------------------------------------------------------------------
    // Data source variables:

    public static final String XSV_CSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_testFile.csv";
    public static final String XSV_TSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_TSV_testFile.csv";
    public static final String XSV_DEADBEEFSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_DEADBEEFSV_testFile.csv";
    public static final String XSV_CSV_PIK3CA_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_PIK3CA.csv";
    public static final String XSV_CSV_MUC16_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_MUC16.csv";

    public static final String DATA_SOURCES_FOLDER_PATH = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator + "dataSources" + File.separator;
    public static final String HGNC_HG19_TSV_PATH = DATA_SOURCES_FOLDER_PATH + "hgnc" + File.separator + "hg19" + File.separator + "hgnc_download_Nov302017.tsv";
    public static final String SIMPLE_UNIPROT_HG19_TSV_PATH = DATA_SOURCES_FOLDER_PATH + "simple_uniprot" + File.separator + "hg19" + File.separator + "simple_uniprot_Dec012014.tsv";
}
