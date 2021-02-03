package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;

/**
 * A class to hold the constants for the Funcotator Tests.
 * Created by jonn on 11/1/17.
 */
public class FuncotatorTestConstants {

    public static final String FUNCOTATOR_LARGE_FILES_DIR = GATKBaseTest.largeFileTestDir + "funcotator" + File.separator;
    public static final String FUNCOTATOR_TEST_DIR        = GATKBaseTest.toolsTestDir + "funcotator" + File.separator;

    public static final String REFERENCE_VERSION_HG19 = "hg19";
    public static final String REFERENCE_VERSION_HG38 = "hg38";

    public static final String REFERENCE_VERSION_ECOLI = "ASM584v2";

    // ----------------------------------------------------------------------
    // Data source variables:

    /**
     * This folder is the top-level folder containing data sources for all funcotator tests.
     * This can be treated as a normal data sources folder for the purposes of funcotator, with the special properties
     * that for Gencode and other large interval-based data sources, the data sources themselves have been trimmed to
     * the extents of any and all variants involved in the testing of funcotator.
     * This is to say that when creating new test data, you may very well have to regenerate the large interval-based
     * data sources so that they have data which will overlap the new variants you add.
     *
     * For gencode you can use the following script as a starting point (and please read it before you run it), though
     * there will be a bit of manual work as well (I did not have the will or the time to automate everything - Jonn Smith):
     *     GATK_DEVELOPMENT_TOP_DIRECTORY/scripts/funcotator/testing/getGencodeGenesForVcfVariants.sh
     */
    public static final String FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER        = FUNCOTATOR_LARGE_FILES_DIR + "funcotator_dataSources" + File.separator;
    /** Local folder containing data sources that point to the cloud. */
    public static final String FUNCOTATOR_DATA_SOURCES_LOCAL_CLOUD_FOLDER = FUNCOTATOR_LARGE_FILES_DIR + "funcotator_dataSources_cloud" + File.separator;
    /** Local folder containing local data sources and one that points to gnomAD on the cloud. */
    public static final String FUNCOTATOR_DATA_SOURCES_LOCAL_CLOUD_GNOMAD_FOLDER = FUNCOTATOR_LARGE_FILES_DIR + "funcotator_dataSources_cloud_gnomad" + File.separator;
    /** Cloud-based folder containing data sources that point to the cloud. */
    public static final String FUNCOTATOR_DATA_SOURCES_REMOTE_CLOUD_FOLDER = "gs://hellbender/test/resources/large/funcotatorDataSourceCollection/funcotator_dataSources_cloud" + File.separator;

    public static final String DUMMY_DATA_SOURCES_TAR_GZ             = FUNCOTATOR_LARGE_FILES_DIR + "dummyDataSources.tar.gz";
    public static final String DUMMY_DATA_SOURCES_TAR_GZ_SHA256_FILE = FUNCOTATOR_LARGE_FILES_DIR + "dummyDataSources.sha256";
    public static final String DUMMY_DATA_SOURCES_FOLDER             = FUNCOTATOR_LARGE_FILES_DIR + "dummyDataSources";

    // MT Info:
    public static final String MT_TRANSCRIPT = "ENST00000361567.2";

    // Gencode main GTF file:
    public static final String GENCODE_DATA_SOURCE_GTF_PATH_HG19   = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode" + File.separator + REFERENCE_VERSION_HG19 + File.separator + "gencode.v19.testVariantSubset.gtf";
    public static final String GENCODE_DATA_SOURCE_FASTA_PATH_HG19 = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode" + File.separator + REFERENCE_VERSION_HG19 + File.separator + "gencode.v19.testVariantSubset.pc_transcripts.fa";
    public static final String GENCODE_DATA_SOURCE_GTF_PATH_HG38   = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode" + File.separator + REFERENCE_VERSION_HG38 + File.separator + "gencode.v28.regressionTestVariantSet.gtf";
    public static final String GENCODE_DATA_SOURCE_FASTA_PATH_HG38 = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode" + File.separator + REFERENCE_VERSION_HG38 + File.separator + "gencode.v28.regressionTestVariantSet.pc_transcripts.fa";

    public static final String XSV_CSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_CSV_testFile.csv";

    public static final String XSV_TSV_FILE_PATH        = FUNCOTATOR_TEST_DIR + "xsv_TSV_testFile.csv";
    public static final String XSV_PIPESV_FILE_PATH     = FUNCOTATOR_TEST_DIR + "xsv_PIPESV_testFile.xsv";
    public static final String XSV_DEADBEEFSV_FILE_PATH = FUNCOTATOR_TEST_DIR + "xsv_DEADBEEFSV_testFile.csv";
    public static final String XSV_CSV_PIK3CA_PATH      = FUNCOTATOR_TEST_DIR + "xsv_CSV_PIK3CA.csv";
    public static final String XSV_CSV_MUC16_PATH       = FUNCOTATOR_TEST_DIR + "xsv_CSV_MUC16.csv";

    public static final String XSV_LOCATABLE_TEST_FILE1_DATA_PATH   = FUNCOTATOR_TEST_DIR + "xsv_locatable_test.csv";
    public static final String XSV_LOCATABLE_TEST_FILE1_CONFIG_PATH = FUNCOTATOR_TEST_DIR + "xsv_locatable_test.config";
    public static final String XSV_LOCATABLE_TEST_FILE2_DATA_PATH   = FUNCOTATOR_TEST_DIR + "xsv_locatable_test2.csv";
    public static final String XSV_LOCATABLE_TEST_FILE2_CONFIG_PATH = FUNCOTATOR_TEST_DIR + "xsv_locatable_test2.config";
    public static final String XSV_LOCATABLE_TEST_FILE3_DATA_PATH   = FUNCOTATOR_TEST_DIR + "xsv_locatable_test3.tsv";
    public static final String XSV_LOCATABLE_TEST_FILE3_CONFIG_PATH = FUNCOTATOR_TEST_DIR + "xsv_locatable_test3.config";

    public static final String COSMIC_TEST_DB = FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "cosmic" + File.separator + "hg19" + File.separator + "CosmicTest.db";

    public static final String DBSNP_HG19_SNIPPET_FILE_PATH = FUNCOTATOR_TEST_DIR + "dbSNP_hg19_snippet.vcf";
    public static final String DBSNP_HG19_SNIPPET_WITH_FILTERS_FILE_PATH = FUNCOTATOR_TEST_DIR + "dbSNP_hg19_snippet_with_filters.vcf";

    public static final String GENCODE_TRANSCRIPT_FASTA_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "gencode.v19.pc_transcripts.fasta";

    public static final double FUNCOTATOR_DOUBLE_COMPARISON_EPSILON = 0.0001;

    public static final String HG19_CHR19_REFERENCE_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "GRCh37.p13.chr19.tar.gz";
    public static final String HG19_CHR17_REFERENCE_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "GRCh37.p13.chr17.tar.gz";
    public static final String HG19_CHR3_REFERENCE_FILE_NAME  = FUNCOTATOR_LARGE_FILES_DIR + "GRCh37.p13.chr3.tar.gz";
    public static final String HG19_3_REFERENCE_FILE_NAME     = FUNCOTATOR_LARGE_FILES_DIR + "b37.3.tar.gz";
    public static final String HG38_3_REFERENCE_FILE_NAME     = FUNCOTATOR_LARGE_FILES_DIR + "hg38.3.tar.gz";
    public static final String HG19_2_REFERENCE_FILE_NAME     = FUNCOTATOR_LARGE_FILES_DIR + "b37.2.tar.gz";

    public static final String ECOLI_REFERENCE_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "e.coli_K12_MG1655.NC_000913.3.fasta";
    public static final String ECOLI_VCF_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "e_coli_K12_MG1655.vcf";

    // A MUC16-only datasource that contains all transcripts:
    public static final String MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "pik3ca_muc16_all_transcripts_ds" + File.separator + "gencode_muc16" + File.separator + "hg19" + File.separator + "gencode.v19.MUC16.gtf";
    public static final String MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_LARGE_FILES_DIR + "pik3ca_muc16_all_transcripts_ds" + File.separator + "gencode_muc16" + File.separator + "hg19" + File.separator + "gencode.v19.MUC16_transcript.fasta";

    // A PIK3CA-only datasource that contains all transcripts:
    public static final String PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "pik3ca_muc16_all_transcripts_ds" + File.separator + "gencode_pik3ca" + File.separator + "hg19" + File.separator + "gencode.v19.PIK3CA.gtf";
    public static final String PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_LARGE_FILES_DIR + "pik3ca_muc16_all_transcripts_ds" + File.separator + "gencode_pik3ca" + File.separator + "hg19" + File.separator + "gencode.v19.PIK3CA_transcript.fasta";


    // A TP53-only datasource for hg19 that contains all transcripts with subunits in numerical order, not the order
    // in which they are transcribed.  This is for testing the fix for
    // https://github.com/broadinstitute/gatk/issues/7051 and ensuring it does not revert:
    public static final String TP53_REVERSE_ORDER_GENCODE_ANNOTATIONS_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "small_ds_order_bugfix_7051_tests" + File.separator + "gencode_tp53" + File.separator + "hg19" + File.separator + "gencode.v34lift37.annotation.REORDERED.TP53.gtf";
    public static final String TP53_REVERSE_ORDER_GENCODE_TRANSCRIPT_FASTA_FILE = FUNCOTATOR_LARGE_FILES_DIR + "small_ds_order_bugfix_7051_tests" + File.separator + "gencode_tp53" + File.separator + "hg19" + File.separator + "gencode.v34lift37.pc_transcripts.TP53.fa";

    /*
     * MUC16 info:
     */

    // The full span of the MUC16 gene on HG19 (including non-basic transcripts):
    public static final SimpleInterval MUC16_HG19_GENE_POSITION = new SimpleInterval("chr19", 8959520, 9092018);

    // The union of the spans of just the *basic* transcripts for the MUC16 gene on HG19:
    // (MUC16 basic transcripts: chr19:8959520-9092018, chr19:8959522-9003586, chr19:8973992-8977690)
    public static final SimpleInterval MUC16_HG19_BASIC_TRANSCRIPTS_SPAN = new SimpleInterval("chr19", 8959520, 9092018);

    public static final String         MUC16_TRANSCRIPT                              = "ENST00000397910.4";
    public static final String         MUC16_PATHOLOGICAL_TRANSCRIPT                 = "ENST00000599436.1";
    public static final String         MUC16_GENCODE_NON_BASIC_ANNOTATIONS_FILE_NAME = FUNCOTATOR_TEST_DIR + "gencode.v19.MUC16.non-basic.gtf";

    /*
     * PIK3CA info:
     */

    // The full span of the PIK3CA gene on HG19 (including non-basic transcripts):
    public static final SimpleInterval PIK3CA_HG19_GENE_POSITION = new SimpleInterval("chr3", 178865902, 178957881);

    // The union of the spans of just the *basic* transcripts for the PIK3CA gene on HG19:
    // (PIK3CA basic transcripts: chr3: 178866311-178957881)
    public static final SimpleInterval PIK3CA_HG19_BASIC_TRANSCRIPTS_SPAN = new SimpleInterval("chr3", 178866311, 178957881);

    public static final String         PIK3CA_TRANSCRIPT = "ENST00000263967.3";
    // Position on HG19 of the PIK3CA_TRANSCRIPT
    public static final SimpleInterval PIK3CA_TRANSCRIPT_POSITION = new SimpleInterval("chr3", 178866311, 178957881);


    public static final String GTF_CHR3_FILE_NAME = FUNCOTATOR_LARGE_FILES_DIR + "gencode.v19.chr_patch_hapl_scaff.chr3.gtf";

    public static final String VARIANT_FILE_HG19_CHR3  = FUNCOTATOR_TEST_DIR + "snpTest_chr3_hg19.vcf";
    public static final String VARIANT_FILE_HG19_CHR19 = FUNCOTATOR_TEST_DIR + "snpTest_chr19_hg19.vcf";

    public static final String BAD_DATA_ONLY_ONE_AD_FIELD_VALUE = FUNCOTATOR_TEST_DIR + "badDataOneAlleleDepthValue_hg38.vcf";

    // ----------------------------------------------------------------------
    // Integration Test Variables:
    public static final String MUC16_MNP_FILE_BASE_NAME    = FUNCOTATOR_TEST_DIR + "MUC16_MNP";
    public static final String PIK3CA_SNP_FILE_BASE_NAME   = FUNCOTATOR_TEST_DIR + "PIK3CA_SNPS";
    public static final String PIK3CA_INDEL_FILE_BASE_NAME = FUNCOTATOR_TEST_DIR + "PIK3CA_INDELS";

    // ----------------------------------------------------------------------
    // Data Validation Test Variables:

    // The following data sets were produced by annotating two WGS VCF files:
    //
    //  B37/HG19: 0816201804HC0_R01C01.vcf
    //      HG38: hg38_trio.vcf
    //
    // and then taking the first 20 (or first 100 in the case of the `LARGE DATA` set) variants from each variant
    // classification and creating a new VCF file (without any annotations) from them.
    //
    // HG19 set 2 was created by doing the same process after lifting over the `hg38_trio.vcf` data to hg19.
    //
    // A script was used to do so, and is checked in in the following path:
    //     GATK_DEVELOPMENT_TOP_DIR/scripts/funcotator/testing/getRepresentativeListOfFuncotations.sh
    public static final String VALIDATION_TEST_DATA_DIR                                             = FUNCOTATOR_TEST_DIR + "validationTestData" + File.separator;
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1                     = VALIDATION_TEST_DATA_DIR + "regressionTestVariantSet1.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_1_EXPECTED_OUTPUT     = VALIDATION_TEST_DATA_DIR + "regressionTestVariantSet1_expected.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_2                     = VALIDATION_TEST_DATA_DIR + "regressionTestVariantSet2.vcf";
    public static final String SINGLE_LINE                                                          = VALIDATION_TEST_DATA_DIR + "hashSetOrderingIssue.vcf";
    public static final String SINGLE_LINE_EXPECTED                                                 = VALIDATION_TEST_DATA_DIR + "hashSetOrderingIssue_expected.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_DATA_SET_2_EXPECTED_OUTPUT     = VALIDATION_TEST_DATA_DIR + "regressionTestVariantSet2_expected.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG38                                = VALIDATION_TEST_DATA_DIR + "regressionTestVariantSetHG38.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_EXPECTED_OUTPUT                     = VALIDATION_TEST_DATA_DIR + "regressionTestVariantSetHG38_expected.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_LARGE_DATA_SET                 = VALIDATION_TEST_DATA_DIR + "regressionTestHg19Large.vcf";
    public static final String NON_TRIVIAL_DATA_VALIDATION_TEST_HG19_LARGE_DATA_SET_EXPECTED_OUTPUT = VALIDATION_TEST_DATA_DIR + "regressionTestHg19Large_expected.vcf";
}
