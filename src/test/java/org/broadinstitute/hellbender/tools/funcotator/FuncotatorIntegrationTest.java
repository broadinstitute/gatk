package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * An integration test for the {@link Funcotator} tool.
 * Created by jonn on 8/29/17.
 */
public class FuncotatorIntegrationTest extends CommandLineProgramTest {

    //TODO: Add checks on the output file.

    //==================================================================================================================

    @DataProvider
    Object[][] provideForIntegrationTest() {
        return new Object[][] {
//                {FuncotatorTestConstants.GTF_CHR3_FILE_NAME, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, 0},
//                {FuncotatorTestConstants.GTF_CHR3_FILE_NAME, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, 1},
//                {FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, 0},
//                {FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, 1},
//                {FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_PATHOLOGICAL_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, 1},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID}
        };
    }

    //==================================================================================================================

    @Test(dataProvider = "provideForIntegrationTest")
    public void basicMarbleRoll(final String dataSourcesPath,
                                final FuncotatorArgumentDefinitions.ReferenceVersionType refVer,
                                final String referenceFileName,
                                final String variantFileName,
                                final String transcriptName,
                                final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType) throws IOException {

        final File outputFile = createTempFile("funcotator_tmp_out", ".vcf");
        final List<String> arguments = new ArrayList<>();

        arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
        arguments.add(dataSourcesPath);
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add(refVer.toString());
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());

        runCommandLine(arguments);
    }

//    @Test()
//    public void spotCheck() throws IOException {
//
//        final File outputFile = createTempFile("funcotator_tmp_out", ".vcf");
//        final List<String> arguments = new ArrayList<>();
//
//        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
//        arguments.add("/Users/jonn/Development/oncotator_testing/BENCHMARK_INPUT.funcotator.vcf");
//
//        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
//        arguments.add("/Users/jonn/Development/references/GRCh37.p13.genome.fasta");
//
//        arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
//        arguments.add(FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER);
//        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
//        arguments.add(FuncotatorArgumentDefinitions.ReferenceVersionType.hg19.toString());
//        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
//        arguments.add(outputFile.getAbsolutePath());
//
//        runCommandLine(arguments);
//    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void exhaustiveArgumentTest(final String dataSourcesPath,
                                       final FuncotatorArgumentDefinitions.ReferenceVersionType refVer,
                                       final String referenceFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + transcriptName + ".vcf";

        final File outputFile = createTempFile(outFileName.substring(0,outFileName.length()-4), outFileName.substring(outFileName.length()-4));

        final List<String> arguments = new ArrayList<>();

        // Required Args:
        arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
        arguments.add(dataSourcesPath);
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add(refVer.toString());
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());

        // Transcript selection:
        arguments.add("--" + FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME);
        arguments.add(FuncotatorArgumentDefinitions.TranscriptSelectionMode.BEST_EFFECT.toString());

        arguments.add("--" + FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME);
        arguments.add(transcriptName);

        // Annotation Defaults and Overrides:
        arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME);
        arguments.add("GARBAGEDAY:SPUMONI");

        arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME);
        arguments.add("Gencode_hugoSymbol:Freddie Mercury");

        // Run the beast:
        runCommandLine(arguments);
    }
}
