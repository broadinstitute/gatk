package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.XSV.SimpleKeyXsvFuncotationFactory;
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
                {FuncotatorTestConstants.GTF_CHR3_FILE_NAME, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3},
        };
    }

    //==================================================================================================================

    // TODO: REFACTOR ALL THIS:

    @Test(dataProvider = "provideDataForBasicMarbleRoll")
    public void basicMarbleRoll(final String gtfFileName, final String referenceFileName, final String fastaFileName, final String variantFileName) throws IOException {
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
    public void exhaustiveArgumentTest(final String gtfFileName, final String referenceFileName, final String fastaFileName, final String variantFileName) throws IOException {
//        final File outputFile = createTempFile("funcotator_tmp_out", ".vcf");
        final File outputFile = new File("funcotator_tmp_out.vcf");
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
        arguments.add("ENST00000263967.3");

        // XSV Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_INPUT_ARG_SHORT_NAME);
        arguments.add(FuncotatorTestConstants.XSV_CSV_PIK3CA_PATH);
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME);
        arguments.add(",");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_SHORT_NAME);
        arguments.add("0");
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_SHORT_NAME);
        arguments.add(SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME.toString());
        arguments.add("-" + FuncotatorArgumentDefinitions.XSV_NAME_ARG_SHORT_NAME);
        arguments.add("TEST_XSV_INPUT");

        // Override Args:
        arguments.add("-" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_SHORT_NAME);
        arguments.add("GARBAGEDAY:SPUMONI");

        arguments.add("-" + FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_SHORT_NAME);
        arguments.add("Gencode_hugoSymbol:Freddie Mercury");

        runCommandLine(arguments);
    }
}
