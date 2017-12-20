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
    Object[][] provideDataForBasicMarbleRoll() {
        return new Object[][] {
                {FuncotatorTestConstants.GTF_CHR3_FILE_NAME, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, 0},
                {FuncotatorTestConstants.GTF_CHR3_FILE_NAME, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, 1},
                {FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, 0},
                {FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, 1},
                {FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_PATHOLOGICAL_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, 1},
        };
    }

    //==================================================================================================================

    // TODO: REFACTOR ALL THIS - IT'S KIND OF GROSS:

    @Test(dataProvider = "provideDataForBasicMarbleRoll")
    public void basicMarbleRoll(final String gtfFileName,
                                final String referenceFileName,
                                final String fastaFileName,
                                final String variantFileName,
                                final String transcriptName,
                                final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                final int xsvMatchColumn) throws IOException {
        final File outputFile = createTempFile("funcotator_tmp_out", ".vcf");
        final List<String> arguments = new ArrayList<>();

        arguments.add("-" + FuncotatorArgumentDefinitions.GTF_FILE_ARG_SHORT_NAME);
        arguments.add(gtfFileName);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + FuncotatorArgumentDefinitions.GENCODE_FASTA_ARG_NAME);
        arguments.add(fastaFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());

        runCommandLine(arguments);
    }

    @Test(dataProvider = "provideDataForBasicMarbleRoll")
    public void exhaustiveArgumentTest(final String gtfFileName,
                                       final String referenceFileName,
                                       final String fastaFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                       final int xsvMatchColumn) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + xsvMatchColumn + "_" + transcriptName + ".vcf";

        final File outputFile = createTempFile(outFileName.substring(0,outFileName.length()-4), outFileName.substring(outFileName.length()-4));

        final List<String> arguments = new ArrayList<>();

        // Required Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.GTF_FILE_ARG_SHORT_NAME);
        arguments.add(gtfFileName);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + FuncotatorArgumentDefinitions.GENCODE_FASTA_ARG_NAME);
        arguments.add(fastaFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());

        // Optional Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_SHORT_NAME);
        arguments.add(FuncotatorArgumentDefinitions.TranscriptSelectionMode.BEST_EFFECT.toString());

        arguments.add("-" + FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_SHORT_NAME);
        arguments.add(transcriptName);

        // Arbitrary XSV Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_INPUT_ARG_SHORT_NAME);
        arguments.add(FuncotatorTestConstants.XSV_CSV_PIK3CA_PATH);
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME);
        arguments.add(",");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_SHORT_NAME);
        arguments.add(Integer.toString(xsvMatchColumn));
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_SHORT_NAME);
        arguments.add(xsvMatchType.toString());
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_NAME_ARG_SHORT_NAME);
        arguments.add("PIK3CA_XSV_INPUT");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_SHORT_NAME);
        arguments.add("true");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_VERSION_ARG_SHORT_NAME);
        arguments.add("0.001.TEST");

        // Arbitrary XSV Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_INPUT_ARG_SHORT_NAME);
        arguments.add(FuncotatorTestConstants.XSV_CSV_MUC16_PATH);
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME);
        arguments.add(",");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_SHORT_NAME);
        arguments.add(Integer.toString(xsvMatchColumn));
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_SHORT_NAME);
        arguments.add(xsvMatchType.toString());
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_NAME_ARG_SHORT_NAME);
        arguments.add("MUC16_XSV_INPUT");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_SHORT_NAME);
        arguments.add("true");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_VERSION_ARG_SHORT_NAME);
        arguments.add("0.002.TEST");

        // Locatable XSV Data Source:
        arguments.add("-" + FuncotatorArgumentDefinitions.LOCATABLE_XSV_IN_ARG_SHORT_NAME);
        arguments.add(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE1_PATH);

        // --------------------------------------------------------------------
        // Packaged Test Data Source Arguments:

        // hgnc:
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_INPUT_ARG_SHORT_NAME);
        arguments.add(FuncotatorTestConstants.HGNC_HG19_TSV_PATH);
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME);
        arguments.add("\t");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_SHORT_NAME);
        arguments.add(Integer.toString(1));
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_SHORT_NAME);
        arguments.add(SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME.toString());
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_NAME_ARG_SHORT_NAME);
        arguments.add("hgnc");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_SHORT_NAME);
        arguments.add("true");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_VERSION_ARG_SHORT_NAME);
        arguments.add("0.003.HGNC");

        // simple_uniprot:
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_INPUT_ARG_SHORT_NAME);
        arguments.add(FuncotatorTestConstants.SIMPLE_UNIPROT_HG19_TSV_PATH);
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME);
        arguments.add("\t");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_SHORT_NAME);
        arguments.add(Integer.toString(0));
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_SHORT_NAME);
        arguments.add(SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME.toString());
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_NAME_ARG_SHORT_NAME);
        arguments.add("simple_uniprot");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_SHORT_NAME);
        arguments.add("true");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_VERSION_ARG_SHORT_NAME);
        arguments.add("0.004.SIMPLE_UNIPROT");

        // Override Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_SHORT_NAME);
        arguments.add("GARBAGEDAY:SPUMONI");

        arguments.add("-" + FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_SHORT_NAME);
        arguments.add("Gencode_hugoSymbol:Freddie Mercury");

        // Run the beast:
        runCommandLine(arguments);
    }
}
