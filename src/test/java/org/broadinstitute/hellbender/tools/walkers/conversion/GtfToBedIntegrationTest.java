package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;

public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    final File outputFile = createTempFile("outputBed", ".bed");

    private static final File mapk1 = new File(ConversionTestUtils.getMapk1Gtf());
    private static final File decoys = new File(ConversionTestUtils.getDecoySampleGtf());
    private static final File manyTranscripts = new File(ConversionTestUtils.getManyTranscriptsGtf());
    private static final File mouse = new File(ConversionTestUtils.getMouseGtf());
    private static final File mouseDictionary = new File(ConversionTestUtils.getMouseDict());
    private static final File dictionary = new File(ConversionTestUtils.getReferenceDict());
    private static final File expectedDecoysGeneBed = new File(ConversionTestUtils.getDecoySamplesGeneBed());
    private static final File expectedDecoysTranscriptBed = new File(ConversionTestUtils.getDecoySamplesTranscriptBed());
    private static final File expectedMapk1GeneBed = new File(ConversionTestUtils.getMapk1GeneBed());
    private static final File expectedMapk1TranscriptBed = new File(ConversionTestUtils.getMapk1TranscriptBed());
    private static final File expectedManyTranscriptsBed = new File (ConversionTestUtils.getManyTranscriptsBed());
    private static final File expectedNotBasicBed = new File(ConversionTestUtils.getNotBasicBed());
    private static final File expectedMouseBed = new File(ConversionTestUtils.getMouseBed());

    // tests specifically mapk1 gene
    @Test
    public void testMapk1Gene() throws IOException {
        runMapk1(false);
    }

    // tests transcripts in mapk1 gene
    @Test
    public void testMapk1Transcript() throws IOException {
        runMapk1(true);
    }

    // tests a sample of decoy genes (gene) - decoy = any gene that doesn't start with "chr"
    @Test
    public void testDecoyGenes() throws IOException {
        runDecoySample(false);
    }

    // tests a sample of decoy genes (transcript)
    @Test
    public void testDecoyTranscripts() throws IOException {
        runDecoySample(true);
    }

    // tests a gene with many transcripts and a gene in another chr
    @Test
    public void testManyTranscripts() throws IOException {
        runManyTranscripts();
    }

    @Test
    public void testNotBasic() throws IOException {
        runNotBasic();
    }

    @Test
    public void testMouse() throws IOException {
        runMouse();
    }

    private void runMapk1(boolean transcript) throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add(GtfToBed.INPUT_LONG_NAME, mapk1)
                .add(GtfToBed.SEQUENCE_DICTIONARY_LONG_NAME, dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .add(GtfToBed.SORT_BY_BASIC_LONG_NAME, true)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);

        if (transcript) {
            Assert.assertEquals(countLines(outputFile), 3);
            IntegrationTestSpec.assertEqualTextFiles(expectedMapk1TranscriptBed, outputFile);
        } else {
            Assert.assertEquals(countLines(outputFile), 1);
            IntegrationTestSpec.assertEqualTextFiles(expectedMapk1GeneBed, outputFile);
        }
    }

    private void runDecoySample(boolean transcript) throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add(GtfToBed.INPUT_LONG_NAME, decoys)
                .add(GtfToBed.SEQUENCE_DICTIONARY_LONG_NAME, dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .add(GtfToBed.SORT_BY_BASIC_LONG_NAME, true)
                .addOutput(outputFile);
        countLines(outputFile);
        runCommandLine(argsBuilder);

        if (transcript) {
            Assert.assertEquals(countLines(outputFile), 20);
            IntegrationTestSpec.assertEqualTextFiles(expectedDecoysTranscriptBed, outputFile);

        } else {
            Assert.assertEquals(countLines(outputFile), 19);
            IntegrationTestSpec.assertEqualTextFiles(expectedDecoysGeneBed, outputFile);
        }
    }

    private void runManyTranscripts() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add(GtfToBed.INPUT_LONG_NAME, manyTranscripts)
                .add(GtfToBed.SEQUENCE_DICTIONARY_LONG_NAME, dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, true)
                .add(GtfToBed.SORT_BY_BASIC_LONG_NAME, true)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        IntegrationTestSpec.assertEqualTextFiles(expectedManyTranscriptsBed, outputFile);

    }

    private void runNotBasic() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add(GtfToBed.INPUT_LONG_NAME, manyTranscripts)
                .add(GtfToBed.SEQUENCE_DICTIONARY_LONG_NAME, dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, true)
                .add(GtfToBed.SORT_BY_BASIC_LONG_NAME, false)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        IntegrationTestSpec.assertEqualTextFiles(expectedNotBasicBed, outputFile);
    }

    private void runMouse() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add(GtfToBed.INPUT_LONG_NAME, mouse)
                .add(GtfToBed.SEQUENCE_DICTIONARY_LONG_NAME, mouseDictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, true)
                .add(GtfToBed.SORT_BY_BASIC_LONG_NAME, false)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        IntegrationTestSpec.assertEqualTextFiles(expectedMouseBed, outputFile);
    }

    private int countLines(File file) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        int lines = 0;
        while (reader.readLine() != null) {
            lines++;
        }
        return lines;
    }
}
