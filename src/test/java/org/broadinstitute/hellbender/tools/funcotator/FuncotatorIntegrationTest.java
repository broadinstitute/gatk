package org.broadinstitute.hellbender.tools.funcotator;

import avro.shaded.com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.CustomMafFuncotationCreator;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRendererConstants;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription;

/**
 * An integration test for the {@link Funcotator} tool.
 * Created by jonn on 8/29/17.
 */
public class FuncotatorIntegrationTest extends CommandLineProgramTest {

    //TODO: SUPER IMPORTANT! Add checks on the output file to make sure it's equal to an expected output file!

    // Temp directory in which to place output files.
    private static final File tmpOutDir;

    // Whether to do debug output (i.e. leave output around).
    // This should always be false when checked in.
    private static final boolean doDebugTests = false;
    private static final String LARGE_DATASOURCES_FOLDER = "funcotator_dataSources_latest";
    private static final String GERMLINE_DATASOURCES_FOLDER = "funcotator_dataSources_germline_latest";

    private static final String XSV_CLINVAR_MULTIHIT_TEST_VCF = toolsTestDir + "funcotator" + File.separator + "clinvar_hg19_multihit_test.vcf";
    private static final String DS_XSV_CLINVAR_TESTS          = largeFileTestDir + "funcotator" + File.separator + "small_ds_clinvar_hg19" + File.separator;
    private static final String NOT_M2_TEST_HG19 = toolsTestDir + "funcotator/NotM2_test_custom_maf_fields.vcf";
    private static final String M2_TEST_HG19 = toolsTestDir + "funcotator/M2_test_custom_maf_fields.vcf";
    private static final String NOT_M2_TEST_HG19_TUMOR_ONLY = toolsTestDir + "funcotator/NotM2_test_custom_maf_fields_tumor_only.vcf";
    private static final String M2_TEST_HG19_TUMOR_ONLY = toolsTestDir + "funcotator/M2_test_custom_maf_fields_tumor_only.vcf";
    private static final String THREE_SAMPLE_SOMATIC = toolsTestDir + "funcotator/Three_sample_somatic.vcf";
    private static final String PIK3CA_VCF_HG19          = toolsTestDir + "funcotator" + File.separator + "0816201804HC0_R01C01.pik3ca.vcf";
    private static final String PIK3CA_VCF_HG38          = toolsTestDir + "funcotator" + File.separator + "hg38_trio.pik3ca.vcf";
    private static final String PIK3CA_VCF_HG19_SNPS     = toolsTestDir + "funcotator" + File.separator + "PIK3CA_SNPS_3.vcf";
    private static final String PIK3CA_VCF_HG19_INDELS   = toolsTestDir + "funcotator" + File.separator + "PIK3CA_INDELS_3.vcf";
    private static final String MUC16_VCF_HG19           = toolsTestDir + "funcotator" + File.separator + "MUC16_MNP.vcf";
    private static final String PIK3CA_VCF_HG19_ALTS     = toolsTestDir + "funcotator" + File.separator + "PIK3CA_3_miss_clinvar_alt_only.vcf";
    private static final String SPANNING_DEL_VCF         = toolsTestDir + "funcotator" + File.separator + "spanning_del.vcf";
    private static final String DS_PIK3CA_DIR            = largeFileTestDir + "funcotator" + File.separator + "small_ds_pik3ca" + File.separator;
    private static final String DS_MUC16_DIR             = largeFileTestDir + "funcotator" + File.separator + "small_ds_muc16" + File.separator;
    private static final String MAF_TEST_CONFIG          = toolsTestDir + "funcotator" + File.separator + "maf.config";
    private static final String XSV_CLINVAR_COL_TEST_VCF = toolsTestDir + "funcotator" + File.separator + "clinvar_hg19_column_test.vcf";
    private static final String DS_XSV_CLINVAR_COL_TEST  = largeFileTestDir + "funcotator" + File.separator + "small_ds_clinvar_hg19" + File.separator;
    private static final String EMPTY_VCF  = publicTestDir + File.separator + "empty.vcf";
    private static final String PIK3CA_DBSNP_DS          = toolsTestDir + "funcotator" + File.separator + "small_pik3ca_dbsnp_ds";
    private static final String MAF_DBSNP_TEST           = toolsTestDir + "funcotator" + File.separator + "maf_dbsnp_test_input.vcf";

    private static final List<String> VCF_FIELDS_GENCODE_19_DS = Arrays.asList("Gencode_19_hugoSymbol","Gencode_19_ncbiBuild","Gencode_19_chromosome","Gencode_19_start","Gencode_19_end","Gencode_19_variantClassification","Gencode_19_variantType","Gencode_19_refAllele","Gencode_19_tumorSeqAllele1","Gencode_19_tumorSeqAllele2","Gencode_19_genomeChange","Gencode_19_annotationTranscript","Gencode_19_transcriptStrand","Gencode_19_transcriptExon","Gencode_19_transcriptPos","Gencode_19_cDnaChange","Gencode_19_codonChange","Gencode_19_proteinChange","Gencode_19_gcContent","Gencode_19_referenceContext","Gencode_19_otherTranscripts");//,"Achilles_Top_Genes","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_HGMD_ID","ClinVar_SYM","ClinVar_TYPE","ClinVar_ASSEMBLY","ClinVar_rs","Cosmic_overlapping_mutations","CosmicFusion_fusion_genes","CosmicFusion_fusion_id","CosmicTissue_total_alterations_in_gene","CosmicTissue_tissue_types_affected","DNARepairGenes_Activity_linked_to_OMIM","DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView","DNARepairGenes_Accession_number_linked_to_NCBI_Entrez","Familial_Cancer_Genes_Syndrome","Familial_Cancer_Genes_Synonym","Familial_Cancer_Genes_Reference","Gencode_XHGNC_hgnc_id","Gencode_XRefSeq_mRNA_id","Gencode_XRefSeq_prot_acc","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","Oreganno_Build","Oreganno_ID","Oreganno_Values","Simple_Uniprot_uniprot_entry_name","Simple_Uniprot_DrugBank","Simple_Uniprot_alt_uniprot_accessions","Simple_Uniprot_uniprot_accession","Simple_Uniprot_GO_Biological_Process","Simple_Uniprot_GO_Cellular_Component","Simple_Uniprot_GO_Molecular_Function","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID");
    private static final List<String> VCF_FIELDS_GENCODE_28_DS = Arrays.asList("Gencode_28_hugoSymbol","Gencode_28_ncbiBuild","Gencode_28_chromosome","Gencode_28_start","Gencode_28_end","Gencode_28_variantClassification","Gencode_28_variantType","Gencode_28_refAllele","Gencode_28_tumorSeqAllele1","Gencode_28_tumorSeqAllele2","Gencode_28_genomeChange","Gencode_28_annotationTranscript","Gencode_28_transcriptStrand","Gencode_28_transcriptExon","Gencode_28_transcriptPos","Gencode_28_cDnaChange","Gencode_28_codonChange","Gencode_28_proteinChange","Gencode_28_gcContent","Gencode_28_referenceContext","Gencode_28_otherTranscripts");//,"Achilles_Top_Genes","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_HGMD_ID","ClinVar_SYM","ClinVar_TYPE","ClinVar_ASSEMBLY","ClinVar_rs","Cosmic_overlapping_mutations","CosmicFusion_fusion_genes","CosmicFusion_fusion_id","CosmicTissue_total_alterations_in_gene","CosmicTissue_tissue_types_affected","DNARepairGenes_Activity_linked_to_OMIM","DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView","DNARepairGenes_Accession_number_linked_to_NCBI_Entrez","Familial_Cancer_Genes_Syndrome","Familial_Cancer_Genes_Synonym","Familial_Cancer_Genes_Reference","Gencode_XHGNC_hgnc_id","Gencode_XRefSeq_mRNA_id","Gencode_XRefSeq_prot_acc","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","Oreganno_Build","Oreganno_ID","Oreganno_Values","Simple_Uniprot_uniprot_entry_name","Simple_Uniprot_DrugBank","Simple_Uniprot_alt_uniprot_accessions","Simple_Uniprot_uniprot_accession","Simple_Uniprot_GO_Biological_Process","Simple_Uniprot_GO_Cellular_Component","Simple_Uniprot_GO_Molecular_Function","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID");
    private static final List<String> MAF_FIELDS_GENCODE_DS = Arrays.asList(MafOutputRendererConstants.FieldName_Hugo_Symbol, MafOutputRendererConstants.FieldName_NCBI_Build, MafOutputRendererConstants.FieldName_Chromosome,
            MafOutputRendererConstants.FieldName_Start_Position, MafOutputRendererConstants.FieldName_End_Position, MafOutputRendererConstants.FieldName_Variant_Classification,  MafOutputRendererConstants.FieldName_Variant_Type,
            MafOutputRendererConstants.FieldName_Reference_Allele, MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, MafOutputRendererConstants.FieldName_Genome_Change, MafOutputRendererConstants.FieldName_Annotation_Transcript, MafOutputRendererConstants.FieldName_Transcript_Strand, MafOutputRendererConstants.FieldName_Transcript_Exon, MafOutputRendererConstants.FieldName_Transcript_Position, MafOutputRendererConstants.FieldName_cDNA_Change, MafOutputRendererConstants.FieldName_Codon_Change, MafOutputRendererConstants.FieldName_Protein_Change, MafOutputRendererConstants.FieldName_gc_content, MafOutputRendererConstants.FieldName_ref_context, MafOutputRendererConstants.FieldName_Other_Transcripts);

    private static String hg38Chr3Ref;
    private static String b37Chr3Ref;
    private static String b37Chr2Ref;
    private static String hg19Chr3Ref;
    private static String hg19Chr19Ref;

    static {
        if (!doDebugTests) {
            tmpOutDir = createTempDir("funcotatorTmpFolder");
        } else {
            tmpOutDir = new File("funcotatorTmpFolder" + File.separator);
            if (!tmpOutDir.mkdirs() && !tmpOutDir.exists()) {
                throw new GATKException("Error making output folder for test: " + FuncotatorIntegrationTest.class.getName());
            }
        }

        hg38Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg38Chr3Ref();
        b37Chr3Ref = FuncotatorReferenceTestUtils.retrieveB37Chr3Ref();
        b37Chr2Ref = FuncotatorReferenceTestUtils.retrieveB37Chr2Ref();
        hg19Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref();
        hg19Chr19Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref();
    }

    //==================================================================================================================
    // Helper methods to create output files and maybe leave them around to debug the test output.

    private static File getOutputFile(final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType) {
        return getOutputFile("funcotator_tmp_out", outputFormatType.toString().toLowerCase());
    }

    private static File getOutputFile(final String outfileBaseName,
                                      final String outFileExtension) {
        final File outputFile;
        if (!doDebugTests) {
            outputFile = createTempFile(tmpOutDir + File.separator + outfileBaseName, "." + outFileExtension);
        } else {
            outputFile = new File(tmpOutDir, outfileBaseName + "." + outFileExtension);
        }
        return outputFile;
    }

    private static void addManualAnnotationsToArguments(final ArgumentsBuilder arguments) {

        // ================================================================================
        // Annotation Defaults:
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "dbSNP_RS:0");
        arguments.addArgumentWithValueThatIncludesWhitespace(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "dbSNP_Val_Status:No Value");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Center:broad.mit.edu");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "source:WES");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "normal_barcode:normal_sample");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "tumor_barcode:tumor_sample");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "NCBI_Build:37");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Strand:+");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "status:Somatic");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "phase:Phase_I");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "sequencer:Illumina");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Tumor_Validation_Allele1:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Tumor_Validation_Allele2:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Validation_Allele1:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Validation_Allele2:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Verification_Status:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Validation_Status:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Validation_Method:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Score:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "BAM_file:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Seq_Allele1:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Seq_Allele2:");

        // ================================================================================
        // Annotation Overrides:
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME, "Oreganno_Build:BUILDED_GOOD_REAL_BIG");
    }

    private ArgumentsBuilder createBaselineArgumentsForFuncotator(final String variantFileName,
                                                                  final File outputFile,
                                                                  final String referenceFileName,
                                                                  final String dataSourcesPath,
                                                                  final String refVer,
                                                                  final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                                                  final boolean shouldValidateSeqDicts) {

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(variantFileName));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(referenceFileName));
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, dataSourcesPath);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, refVer);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.addArgument("verbosity", "INFO");

        if ( !shouldValidateSeqDicts ) {
            // Disable the sequence dictionary check for the tests:
            arguments.addBooleanArgument(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);
        }

        return arguments;
    }

    //==================================================================================================================

    @DataProvider
    Iterator<Object[]> provideForIntegrationTest() {

        final ArrayList<Object[]> testCases = new ArrayList<>();

        // TODO: Fix this set of tests!  THEY DON'T ALL PASS!  (Issue: https://github.com/broadinstitute/gatk/issues/4344)
        // VCF SNP / MNP / INDEL tests
//        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values() ) {
//            for ( final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values() ) {
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//            }
//        }

        // Basic Test Cases:
        for (final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            for (final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values()) {
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                                hg19Chr3Ref,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3,
                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                                hg19Chr19Ref,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19,
                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
            }
        }

        return testCases.iterator();
    }

    @DataProvider
    public Object[][] provideForLargeDataValidationTest() {
        return new Object[][]{
                {
                        "M2_01115161-TA1-filtered.vcf",
                        "Homo_sapiens_assembly19.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER
                },
                {
                        "C828.TCGA-D3-A2JP-06A-11D-A19A-08.3-filtered.PASS.vcf",
                        "Homo_sapiens_assembly19.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER
                },
                {
                        "hg38_test_variants.vcf",
                        "Homo_sapiens_assembly38.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER
                },
                {
                        "sample21.trimmed.vcf",
                        "Homo_sapiens_assembly38.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER
                },
                {
                        "0816201804HC0_R01C01.vcf",
                        "Homo_sapiens_assembly19.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER
                },
                {
                        "hg38_trio.vcf",
                        "Homo_sapiens_assembly38.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER
                },
                {
                        "0816201804HC0_R01C01.vcf",
                        "Homo_sapiens_assembly19.fasta",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        getFuncotatorLargeDataValidationTestInputPath() + GERMLINE_DATASOURCES_FOLDER
                },
        };
    }

    //==================================================================================================================

    // This test is to make sure we don't create a bunch of temp files anywhere.
    // It will force anyone who changes the outputToTmpDir flag to make it true when they check in this test file.
    @Test(groups = {"funcotatorValidation"})
    public void metaTestEnsureTempDirs() {
        Assert.assertEquals(doDebugTests, false);
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void testFuncotatorWithoutValidatingResults(final String dataSourcesPath,
                                                       final String refVer,
                                                       final String referenceFileName,
                                                       final String variantFileName,
                                                       final String transcriptName,
                                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                                       final String expectedOutputFilePath) throws IOException {

        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(variantFileName, outputFile, referenceFileName, dataSourcesPath, refVer, outputFormatType, false);

        runCommandLine(arguments);
    }

    private void validateFuncotationsOnVcf(final Iterable<VariantContext> vcfIterable, final String[] funcotationFieldNames) {
        for (final VariantContext vc : vcfIterable ) {
            final String funcotation = vc.getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, "");

            Assert.assertNotEquals(funcotation, "");

            final String rawFuncotations = funcotation.substring(1,funcotation.length()-1);

            Assert.assertEquals(StringUtils.countMatches(rawFuncotations, VcfOutputRenderer.FIELD_DELIMITER), funcotationFieldNames.length - 1);

            // This is here to make sure we can create the FuncotationMap object without exploding.
            // It serves as a secondary check.
            final FuncotationMap funkyMap = FuncotationMap.createAsAllTableFuncotationsFromVcf(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, funcotationFieldNames,
                    funcotation, vc.getAlternateAllele(0), "VCF");
        }
    }

    @Test(enabled = doDebugTests,
          groups = {"funcotatorValidation"},
          dataProvider = "provideForLargeDataValidationTest")
    public void largeDataValidationTest(final String inputVcfName,
                                        final String referencePath,
                                        final String referenceVersion,
                                        final String dataSourcesPath) throws IOException {

        // Get our main test folder path from our environment:
        final String testFolderInputPath = getFuncotatorLargeDataValidationTestInputPath();

        long startTime = 0, endTime = 0;
        final long overallStartTime = System.nanoTime();

        final String outFileBaseName = inputVcfName + ".funcotator";

        for (final FuncotatorArgumentDefinitions.OutputFormatType outFormat : FuncotatorArgumentDefinitions.OutputFormatType.values()) {

            startTime = System.nanoTime();

            final File outputFile;
            if (outFormat == FuncotatorArgumentDefinitions.OutputFormatType.VCF) {
                outputFile = getOutputFile(outFileBaseName, outFormat.toString().toLowerCase());
            } else {
                outputFile = getOutputFile(outFileBaseName + ".maf", "tsv");
            }

            final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                    testFolderInputPath + inputVcfName,
                    outputFile,
                    testFolderInputPath + referencePath,
                    getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER,
                    referenceVersion,
                    outFormat,
                    true);

            // Add our manual annotations to the arguments:
            addManualAnnotationsToArguments(arguments);

            // Run the tool with our args:
            runCommandLine(arguments);

            endTime = System.nanoTime();

            System.out.println("  Elapsed Time (" + outFormat.toString() + "): " + (endTime - startTime) / 1e9 + "s");

            // Perform a simple validation if the file was a VCF:
            if ( outFormat == FuncotatorArgumentDefinitions.OutputFormatType.VCF) {

                try ( final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(outputFile.getAbsolutePath()) ) {
                    Assert.assertTrue(vcfReader.getHeader() instanceof VCFHeader, "Header is not a VCFHeader!");
                    final VCFHeader vcfHeader = (VCFHeader)vcfReader.getHeader();

                    final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
                    final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

                    validateFuncotationsOnVcf(vcfReader, funcotationFieldNames);

                }
            }
        }

        System.out.println("Total Elapsed Time: " + (endTime - overallStartTime) / 1e9 + "s");
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void exhaustiveArgumentTest(final String dataSourcesPath,
                                       final String refVer,
                                       final String referenceFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                       final String expectedOutputFilePath) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + transcriptName + "." + outputFormatType.toString().toLowerCase();

        final File outputFile = getOutputFile(outFileName.substring(0, outFileName.length() - 4), outFileName.substring(outFileName.length() - 3));

        // Set up our transcript of interest in a file:
        final File transcriptIdFile = getSafeNonExistentFile("TranscriptIdFile.txt");
        try (final BufferedWriter bw = new BufferedWriter(new FileWriter(transcriptIdFile))) {
            bw.write(transcriptName);
        }

        // Required Args:
        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                variantFileName,
                outputFile,
                referenceFileName,
                dataSourcesPath,
                refVer,
                outputFormatType,
                false);

        // Transcript selection:
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.BEST_EFFECT.toString());
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME, transcriptIdFile.toString());

        // Annotation Defaults and Overrides:
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "GARBAGEDAY:SPUMONI");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME, "Oreganno_Build:BUILDED_GOOD_REAL_BIG");

        // Run the beast:
        runCommandLine(arguments);

        // Only test for content-correctness if the output file was specified:
        if (expectedOutputFilePath != null) {
            // Get the expected output file:
            final File expectedOutputFile = new File(expectedOutputFilePath);

            // Make sure that the actual and expected output files are the same:
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedOutputFile, "#");
        }
    }

    /**
     * Test that we can annotate a b37 (here it is clinvar) datasource when GENCODE is using "chr*" and the datasource is not.
     */
    @Test
    public void testCanAnnotateMixedContigHg19Clinvar() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        final int NUM_VARIANTS = 21;
        final int NUM_CLINVAR_HITS = 4;
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS, "Found unexpected number of variants!");

        // Look for "MedGen" to know that we have a clinvar hit.
        Assert.assertEquals(variantContexts.stream()
                .filter(vc -> StringUtils.contains(vc.getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, ""), "MedGen"))
                .count(), NUM_CLINVAR_HITS, "Found unexpected number of ClinVar hits!");
    }

    @Test
    public void testXsvLocatableAnnotationsHaveOnlyOneEntryForMultiHitLocations() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                XSV_CLINVAR_MULTIHIT_TEST_VCF,
                outputFile,
                b37Chr2Ref,
                DS_XSV_CLINVAR_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final VCFInfoHeaderLine funcotationHeaderLine = vcfInfo.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final int EXPECTED_NUM_VARIANTS = 1;
        Assert.assertEquals(vcfInfo.getRight().size(), EXPECTED_NUM_VARIANTS, "Found more than " + EXPECTED_NUM_VARIANTS + " variants!");

        validateFuncotationsOnVcf(vcfInfo.getRight(), funcotationFieldNames);
    }

    @Test
    public void testXsvLocatableAnnotationsHaveCorrectColsForOnlyOnePositionSpecified() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                XSV_CLINVAR_COL_TEST_VCF,
                outputFile,
                b37Chr2Ref,
                DS_XSV_CLINVAR_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final VCFInfoHeaderLine funcotationHeaderLine = vcfInfo.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        final String[] funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final int EXPECTED_NUM_VARIANTS = 10;
        Assert.assertEquals(vcfInfo.getRight().size(), EXPECTED_NUM_VARIANTS);

        validateFuncotationsOnVcf(vcfInfo.getRight(), funcotationFieldNames);
    }

    @Test
    public void testCanAnnotateHg38ClinvarAndGencodeV28() {
        // Clinvar datasource did  go through one round of preprocessing to make contig names "1" --> "chr1" (for example).  This is an issue with ClinVar, not GATK.
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG38,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);

        runCommandLine(arguments);

        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        final int NUM_VARIANTS = 100;
        final int NUM_CLINVAR_HITS = 1; // This was verified with SelectVariants
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS);

        // Look for "MedGen" to know that we have a clinvar hit.
        Assert.assertEquals(variantContexts.stream()
                .filter(vc -> StringUtils.contains(vc.getAttributeAsString("FUNCOTATION", ""), "MedGen"))
                .count(), NUM_CLINVAR_HITS);
    }

    @DataProvider(name = "provideForMafVcfConcordance")
    final Object[][] provideForMafVcfConcordance() {
        return new Object[][]{
                {PIK3CA_VCF_HG19_SNPS, b37Chr3Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList("Gencode_19_proteinChange"), Collections.singletonList(MafOutputRendererConstants.FieldName_Protein_Change), DS_PIK3CA_DIR, true, 15},
                {PIK3CA_VCF_HG19_INDELS, b37Chr3Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList("Gencode_19_proteinChange"), Collections.singletonList(MafOutputRendererConstants.FieldName_Protein_Change), DS_PIK3CA_DIR, true, 57},
                {MUC16_VCF_HG19, hg19Chr19Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList("Gencode_19_proteinChange"), Collections.singletonList(MafOutputRendererConstants.FieldName_Protein_Change), DS_MUC16_DIR, false, 2057},
                {
                        PIK3CA_VCF_HG38,
                        hg38Chr3Ref,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        VCF_FIELDS_GENCODE_28_DS,
                        MAF_FIELDS_GENCODE_DS,
                        DS_PIK3CA_DIR,
                        false,
                        104,
                },
                {PIK3CA_VCF_HG19_INDELS, b37Chr3Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, VCF_FIELDS_GENCODE_19_DS, MAF_FIELDS_GENCODE_DS, DS_PIK3CA_DIR, true, 57},
        };
    }

    private void createConfigFileForMAF(final File mafConfigFile) {
        try ( final PrintWriter printWriter = new PrintWriter(mafConfigFile) ) {
            printWriter.println("contig_column = " + MafOutputRendererConstants.FieldName_Chromosome);
            printWriter.println("start_column = " + MafOutputRendererConstants.FieldName_Start_Position);
            printWriter.println("end_column = " + MafOutputRendererConstants.FieldName_End_Position);
            printWriter.println("xsv_delimiter = \\t");
            printWriter.println("name = ");
        }
        catch (final FileNotFoundException ex) {
            throw new GATKException("Could not create the tmp config file to test maf/vcf concorance: " + mafConfigFile.toURI().toString(), ex);
        }
    }

    /**
     * Make sure that VCFs and MAFs have exactly the same annotation strings.  This test does not look for
     *  multiallelics.
     */
    @Test(dataProvider = "provideForMafVcfConcordance")
    public void testVcfMafConcordance(final String inputVcf,
                                      final String inputRef,
                                      final String funcotatorRef,
                                      final List<String> annotationsToCheckVcf,
                                      final List<String> annotationsToCheckMaf,
                                      final String datasourceDir,
                                      final boolean forceB37Hg19Conversion,
                                      final int gtNumVariants) {

        // ===========================================================================================
        // Run in VCF Mode:
        // ===========================================
        final FuncotatorArgumentDefinitions.OutputFormatType vcfOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File vcfOutputFile = getOutputFile(vcfOutputFormatType);

        final ArgumentsBuilder argumentsVcf = createBaselineArgumentsForFuncotator(
                inputVcf,
                vcfOutputFile,
                inputRef,
                datasourceDir,
                funcotatorRef,
                vcfOutputFormatType,
                false);

        argumentsVcf.addBooleanArgument(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        if ( forceB37Hg19Conversion ) {
            // We need this argument since we are testing on a subset of b37
            argumentsVcf.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        }

        runCommandLine(argumentsVcf);

        // ===========================================================================================
        // Run in MAF Mode:
        // ===========================================

        final FuncotatorArgumentDefinitions.OutputFormatType mafOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File mafOutputFile = getOutputFile(mafOutputFormatType);

        final ArgumentsBuilder argumentsMaf = createBaselineArgumentsForFuncotator(
                inputVcf,
                mafOutputFile,
                inputRef,
                datasourceDir,
                funcotatorRef,
                mafOutputFormatType,
                false);

        argumentsMaf.addBooleanArgument(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);
        argumentsMaf.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        if ( forceB37Hg19Conversion ) {
            // We need this argument since we are testing on a subset of b37
            argumentsMaf.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        }

        runCommandLine(argumentsMaf);

        // TODO: Create another config file for this MAF that has good column headers in it.
        final File mafConfigFile = createTempFile(mafOutputFile.getName(), ".config");
        createConfigFileForMAF(mafConfigFile);

        // ===========================================================================================
        // Read in and validate VCF:
        // ===========================================

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(vcfOutputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        Assert.assertEquals(variantContexts.stream().map(vc -> vc.getAlleles().size() - 1).mapToInt(Integer::intValue).sum(), gtNumVariants);
        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.hasAttribute(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME)));

        // ===========================================================================================
        // Read in and validate MAF:
        // ===========================================

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(mafOutputFile.toPath(), mafConfigFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(), gtNumVariants);

        // Some errors manifest as all of the variant classifications being IGR.  Check to make sure that is not the case.
        Assert.assertTrue(maf.getRecords().stream()
                .anyMatch(v -> !v.getAnnotationValue(MafOutputRendererConstants.FieldName_Variant_Classification).equals("IGR")), "Output produced only IGR annotations!");

        // ===========================================================================================
        // Compare VCF and MAF Annotations:
        // ===========================================

        // Get the annotation fields:
        final String[] funcotationKeys = extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        for (int i = 0; i < annotationsToCheckMaf.size(); i++) {
            final String annotationToCheckVcf = annotationsToCheckVcf.get(i);
            final String annotationToCheckMaf = annotationsToCheckMaf.get(i);

            // Have to get the contig / start / end from the interval itself for MAF:
            List<String> mafFieldValues;
            if ( annotationsToCheckMaf.get(i).equals(MafOutputRendererConstants.FieldName_Chromosome) ) {
                mafFieldValues = maf.getRecords().stream().map(AnnotatedInterval::getContig).collect(Collectors.toList());
            }
            else if ( annotationsToCheckMaf.get(i).equals(MafOutputRendererConstants.FieldName_Start_Position) ) {
                mafFieldValues = maf.getRecords().stream().map(AnnotatedInterval::getStart).map(x -> new Integer(x)).map(Object::toString).collect(Collectors.toList());
            }
            else if ( annotationsToCheckMaf.get(i).equals(MafOutputRendererConstants.FieldName_End_Position) ) {
                mafFieldValues = maf.getRecords().stream().map(AnnotatedInterval::getEnd).map(x -> new Integer(x)).map(Object::toString).collect(Collectors.toList());
            }
            else {
                mafFieldValues = maf.getRecords().stream().map(v -> v.getAnnotationValue(annotationToCheckMaf)).collect(Collectors.toList());
            }
            mafFieldValues = mafFieldValues.stream().map(val -> MafOutputRenderer.mafTransformInvert(annotationToCheckMaf, val, funcotatorRef)).collect(Collectors.toList());

            // Note that we assume that each variant context has one allele and one transcript.  This is true due to the
            // datasources and input VCF.
            // Don't try to refactor this for-loop to a stream here.
            final List<String> vcfFieldValues = new ArrayList<>();
            for (final VariantContext v: variantContexts) {
                final Map<Allele, FuncotationMap> alleleFuncotationMapMap = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                        funcotationKeys, v, "Gencode_19_annotationTranscript", "TEST");
                for (final Allele alternateAllele : v.getAlternateAlleles() ) {
                    final FuncotationMap funcotationMap = alleleFuncotationMapMap.get(alternateAllele);
                    vcfFieldValues.add(funcotationMap.getFieldValue(funcotationMap.getTranscriptList().get(0), annotationToCheckVcf, alternateAllele));
                }
            }
            Assert.assertEquals(mafFieldValues, vcfFieldValues, "Failed matching (VCF: " + annotationToCheckVcf + " , MAF: " + annotationToCheckMaf + ")");
        }
    }

    @Test
    public void testVcfDatasourceAccountsForAltAlleles() {
        final FuncotatorArgumentDefinitions.OutputFormatType vcfOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File vcfOutputFile = getOutputFile(vcfOutputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19_ALTS,
                vcfOutputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                vcfOutputFormatType,
                false);

        // We need this argument since we are testing on a subset of b37
        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(vcfOutputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        Assert.assertTrue(variantContexts.size() > 0);
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationKeys = extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        // The first variant context should have clinvar annotations, since it hit on the alt allele.  None of the rest.
        // This test assumes that each test variant context has only one alt allele.
        // The rest should not have any clinvar hits.
        for (int i = 0; i < variantContexts.size(); i++) {
            final String gtString = (i == 0) ? "MedGen:C0027672,SNOMED_CT:699346009" : "";
            final Map<Allele, FuncotationMap> alleleToFuncotationMap = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(funcotationKeys,
                    variantContexts.get(i), "Gencode_19_annotationTranscript", "TEST");
            Assert.assertEquals(alleleToFuncotationMap.entrySet().size(), 1, "Found more than 1 alternate allele!");

            final FuncotationMap funcotationMap = alleleToFuncotationMap.values().iterator().next();
            Assert.assertEquals(funcotationMap.getTranscriptList().size(), 1, "Found more than 1 funcotation!");
            Assert.assertTrue(funcotationMap.getTranscriptList().stream().noneMatch(k -> k.equals(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY)), "Found a funcotation with an unknown transcript name: " + FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY);
            Assert.assertTrue(funcotationMap.getTranscriptList().stream().noneMatch(StringUtils::isEmpty), "Found a funcotation with an empty transcript!");
            final List<Funcotation> funcotations = funcotationMap.get(funcotationMap.getTranscriptList().get(0));
            Assert.assertEquals(funcotations.size(), 1, "Found more than one funcotation in the funcotation map!");
            final Funcotation funcotation = funcotations.get(0);

            Assert.assertEquals(funcotation.getField("dummy_ClinVar_VCF_CLNDISDB"), FuncotatorUtils.sanitizeFuncotationForVcf(gtString), "Field (dummy_ClinVar_VCF_CLNDISDB) was unsanititzed: " + funcotation.getField("dummy_ClinVar_VCF_CLNDISDB"));
        }
    }

    @Test
    public void testCanAnnotateSpanningDeletions() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                SPANNING_DEL_VCF,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        Assert.assertTrue(variantContexts.size() > 0);
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationKeys = extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        final int NUM_VARIANTS = 10;
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS);
        Assert.assertEquals(variantContexts.stream().filter(v -> v.hasAllele(Allele.SPAN_DEL)).count(), NUM_VARIANTS);

        for (final VariantContext vc : variantContexts) {
            final Map<Allele, FuncotationMap> alleleToFuncotationMap = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(funcotationKeys,
                    vc, "Gencode_28_annotationTranscript", "TEST");

            // Make sure that all spanning deletions are funcotated with could not determine and that the others are something else.
            for (final Allele allele : alleleToFuncotationMap.keySet()) {
                final List<String> txIds = alleleToFuncotationMap.get(allele).getTranscriptList();
                for (final String txId : txIds) {
                    if (allele.equals(Allele.SPAN_DEL)) {
                        Assert.assertEquals(alleleToFuncotationMap.get(allele).get(txId).get(0).getField("Gencode_28_variantClassification"), GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE.toString());
                        Assert.assertEquals(alleleToFuncotationMap.get(allele).get(txId).get(0).getField("Gencode_28_variantType"), GencodeFuncotation.VariantType.NA.toString());
                    } else {
                        Assert.assertNotEquals(alleleToFuncotationMap.get(allele).get(txId).get(0).getField("Gencode_28_variantClassification"), GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE.toString());
                    }
                }
            }
        }
    }

    @Test
    public void testNoSpanningDeletionWriteWithMAF() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                SPANNING_DEL_VCF,
                outputFile,
                hg38Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                outputFormatType,
                false);

        runCommandLine(arguments);

        // There should only be one variant in the MAF, not two.
        //  TODO:  This input VCF has an "END" field which (currently) throws off the reading of an AnnotatedIntervalCollection.  The MAF_TEST_CONFIG is a workaround/hack.  See issue https://github.com/broadinstitute/gatk/issues/4897
        final AnnotatedIntervalCollection annotatedIntervalCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), Paths.get(MAF_TEST_CONFIG), null);
        Assert.assertEquals(annotatedIntervalCollection.getRecords().size(), 10);
    }

    @Test
    public void testVCFToMAFPreservesFields() {

        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        // Disable the sequence dictionary check for the tests:
        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertTrue(maf.getRecords().size() > 0);
        Assert.assertTrue(maf.getRecords().stream().allMatch(r -> r.hasAnnotation("ILLUMINA_BUILD")));
        Assert.assertTrue(maf.getRecords().stream().allMatch(r -> r.getAnnotationValue("ILLUMINA_BUILD").startsWith("37")));

        // Needs to get aliases from the MAF, since AF (and maybe more) has its name changed.  So create a dummy
        //  MafOutputRenderer that mimics the one that is used in the command line invocation above and get the aliases.
        final File dummyOutputFile = getOutputFile(outputFormatType);
        final MafOutputRenderer dummyMafOutputRenderer = new MafOutputRenderer(dummyOutputFile.toPath(), Collections.emptyList(), new VCFHeader(), new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), "b37");
        final Map<String, Set<String>> mafAliasMap = dummyMafOutputRenderer.getReverseOutputFieldNameMap();

        // Get all of the alias lists
        final Pair<VCFHeader, List<VariantContext>> vcf  = VariantContextTestUtils.readEntireVCFIntoMemory(PIK3CA_VCF_HG19);
        final Set<String> vcfHeaderInfoSet = vcf.getLeft().getInfoHeaderLines().stream()
                .map(VCFCompoundHeaderLine::getID)
                .map(s -> mafAliasMap.getOrDefault(s, Collections.singleton(s)))
                .flatMap(Set::stream)
                .collect(Collectors.toSet());
        Assert.assertTrue(vcfHeaderInfoSet.size() > 0);
        Assert.assertTrue(maf.getAnnotations().containsAll(vcfHeaderInfoSet));
    }

    @Test
    public void testVCFToVCFPreservesFields() {

        final Set<String> PIK3CA_VCF_HG19_INPUT_FIELDS = new HashSet<>(Arrays.asList("FUNCOTATION", "AC", "AF", "ALLELE_A", "ALLELE_B", "AN", "BEADSET_ID", "DP", "GC_SCORE", "ILLUMINA_BUILD", "ILLUMINA_CHR", "ILLUMINA_POS", "ILLUMINA_STRAND", "N_AA", "N_AB", "N_BB", "PROBE_A", "PROBE_B", "SOURCE", "devR_AA", "devR_AB", "devR_BB", "devTHETA_AA", "devTHETA_AB", "devTHETA_BB", "devX_AA", "devX_AB", "devX_BB", "devY_AA", "devY_AB", "devY_BB", "meanR_AA", "meanR_AB", "meanR_BB", "meanTHETA_AA", "meanTHETA_AB", "meanTHETA_BB", "meanX_AA", "meanX_AB", "meanX_BB", "meanY_AA", "meanY_AB", "meanY_BB", "refSNP", "zthresh_X", "zthresh_Y"));

        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                PIK3CA_VCF_HG19,
                outputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        // We need this argument since we are testing on a subset of b37
        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);

        final Pair<VCFHeader, List<VariantContext>> vcf  = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());

        final Set<String> vcfFields = vcf.getLeft().getInfoHeaderLines().stream()
                .map(VCFCompoundHeaderLine::getID)
                .collect(Collectors.toSet());

        final Set<String> missingFields = Sets.difference(PIK3CA_VCF_HG19_INPUT_FIELDS, vcfFields);
        final Set<String> additionalFields = Sets.difference(vcfFields, PIK3CA_VCF_HG19_INPUT_FIELDS);

        Assert.assertTrue(missingFields.size() == 0, "Fields were missing in the output: " + missingFields.stream().collect(Collectors.joining(", ")));
        Assert.assertTrue(additionalFields.size() == 0, "Fields were added in the output: " + additionalFields.stream().collect(Collectors.joining(", ")));

        final List<VariantContext> variantContexts = vcf.getRight();

        Assert.assertTrue(variantContexts.size() > 0);
        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.hasAttribute("ILLUMINA_BUILD")));
        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.getAttributeAsString("ILLUMINA_BUILD", "").startsWith("37")));
        final VCFInfoHeaderLine funcotationHeaderLine = vcf.getLeft().getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        // Make sure that none of the input fields appear in the funcotation header.
        final Set<String> funcotationKeys = new HashSet<>(Arrays.asList(extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription())));
        Assert.assertEquals(Sets.intersection(funcotationKeys, PIK3CA_VCF_HG19_INPUT_FIELDS).size(), 0);
    }

    @DataProvider
    public Object[][] provideTNVcfs() {
        // These two VCFs are exactly the same, except how the sample names are handled.
        return new Object[][] {
                {NOT_M2_TEST_HG19},
                {M2_TEST_HG19}
        };
    }
    @Test(dataProvider = "provideTNVcfs")
    public void testMafCustomCountFields(final String tnVcf) {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                tnVcf,
                outputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Painstakingly lifted from the input VCF file
        final int[] tRefCounts = new int[]{42, 44, 117, 16, 24, 23};
        final int[] tAltCounts = new int[]{3, 2, 4, 6, 8, 7};
        final int[] nRefCounts = new int[]{43, 45, 145, 22, 45, 13};
        final int[] nAltCounts = new int[]{1, 1, 4, 0, 0, 0};
        final double[] tumorF = new double[]{0.075, 0.049, 0.034, 0.278, 0.259, 0.19};

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_alt_count)),
                   tAltCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_ref_count)),
                    tRefCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_ref_count)),
                    nRefCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_alt_count)),
                    nAltCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Double.parseDouble(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_tumor_f)),
                    tumorF[i], "Record did not match GT at entry " + i));
    }
    @DataProvider
    public Object[][] provideTOnlyVcfs() {
        // These two VCFs are exactly the same, except how the sample names are handled.
        return new Object[][] {
                {NOT_M2_TEST_HG19_TUMOR_ONLY},
                {M2_TEST_HG19_TUMOR_ONLY}
        };
    }
    @Test(dataProvider = "provideTOnlyVcfs")
    public void testMafCustomCountFieldsTumorOnly(final String tnVcf) {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                tnVcf,
                outputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Painstakingly lifted from the input VCF file
        final int[] tRefCounts = new int[]{42, 44, 117, 16, 24, 23};
        final int[] tAltCounts = new int[]{3, 2, 4, 6, 8, 7};
        final double[] tumorF = new double[]{0.075, 0.049, 0.034, 0.278, 0.259, 0.19};

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_alt_count)),
                   tAltCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Integer.parseInt(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_t_ref_count)),
                    tRefCounts[i], "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_ref_count),
                    "", "Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_n_alt_count),
                        "","Record did not match GT at entry " + i));

        IntStream.range(0, maf.getRecords().size()).boxed()
                .forEach(i -> Assert.assertEquals(Double.parseDouble(maf.getRecords().get(i).getAnnotationValue(MafOutputRendererConstants.FieldName_tumor_f)),
                    tumorF[i], "Record did not match GT at entry " + i));
    }

    /**
     * Test that an input file with more than one possible tumor-normal pairing fails.
     */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testMoreThanOneTNPair() {
        runSomaticVcf(THREE_SAMPLE_SOMATIC);
    }

    private void runSomaticVcf(final String vcfFile) {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                vcfFile,
                outputFile,
                b37Chr3Ref,
                DS_PIK3CA_DIR,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                false);

        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());

        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testSequenceDictionaryCheck() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        // The input VCF and Reference file are incompatible because
        // the reference file dictionary has only chromosome 2 and the
        // input VCF has a dictionary that contains all contigs for HG19.
        // Therefore the reference dictionary is NOT a superset of the input VCF dictionary.

        final ArgumentsBuilder arguments = createBaselineArgumentsForFuncotator(
                XSV_CLINVAR_COL_TEST_VCF,
                outputFile,
                b37Chr2Ref,
                DS_XSV_CLINVAR_TESTS,
                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                outputFormatType,
                true);

        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION, true);

        runCommandLine(arguments);
    }

    @Test
    public void testNoVariantsProduceMaf() {
        // Make sure that a MAF is actually produced (not a blank file).  Testing https://github.com/broadinstitute/gatk/issues/4937
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(EMPTY_VCF));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(b37Chr3Ref));
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG19);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());
        arguments.addBooleanArgument(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);
        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(),  0);
        assertCustomFieldsArePresent(maf);


    }

    private void assertCustomFieldsArePresent(final AnnotatedIntervalCollection maf) {
        // Double-check that the custom MAF fields are present.
        Assert.assertTrue(CustomMafFuncotationCreator.COUNT_FIELD_NAMES.stream()
                .allMatch(f -> maf.getAnnotations().contains(f)));

        Assert.assertTrue(maf.getAnnotations().contains(MafOutputRendererConstants.FieldName_dbSNP_Val_Status));
    }

    @Test
    public void testEnsureDbSnpInMaf() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(MAF_DBSNP_TEST));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(b37Chr3Ref));
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, PIK3CA_DBSNP_DS);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG19);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.CANONICAL.toString());
        arguments.addBooleanArgument(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);
        runCommandLine(arguments);

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(),  5);
        assertCustomFieldsArePresent(maf);

        final List<String> gtDbSnpValStatus = Arrays.asList("",
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG,
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG,
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG,
                CustomMafFuncotationCreator.BY_FREQ + CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER + CustomMafFuncotationCreator.BY_1KG);
        final List<String> guessDbSnpValStatus = maf.getRecords().stream().map(r -> r.getAnnotationValue(MafOutputRendererConstants.FieldName_dbSNP_Val_Status))
                .collect(Collectors.toList());
        Assert.assertEquals(guessDbSnpValStatus, gtDbSnpValStatus);
    }
}

