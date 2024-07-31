package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;


public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir);
    private static final File input = new File(TEST_SUB_DIR, "/funcotator/small_pik3ca_dbsnp_ds/gencode_pik3ca/hg19/gencode.v19.PIK3CA.gtf");

    private static final File mapk1 = new File(ConversionTestUtils.getMapk1Gtf());
    private static final File decoys = new File(ConversionTestUtils.getDecoySampleGtf());
    private static final File dictionary = new File(ConversionTestUtils.getReferenceDict());
    private static final File decoysGeneBed = new File(ConversionTestUtils.getDecoySamplesGeneBed());
    private static final File decoysTranscriptBed = new File(ConversionTestUtils.getDecoySamplesTranscriptBed());
    private static final File mapk1GeneBed = new File(ConversionTestUtils.getMapk1GeneBed());
    private static final File mapk1TranscriptBed = new File(ConversionTestUtils.getMapk1TranscriptBed());


    @Test
    public void testGene() throws IOException {
        runGtfToBed(false);
    }

    @Test
    public void testTranscript() throws IOException {
        runGtfToBed(true);
    }

    @Test
    public void testMapk1Gene() throws IOException {
        runMapk1(false);
    }

    @Test
    public void testMapk1Transcript() throws IOException {
        runMapk1(true);
    }

    @Test
    public void testDecoyGenes() throws IOException {
        runDecoySample(false);
    }

    @Test
    public void testDecoyTranscripts() throws IOException {
        runDecoySample(true);
    }

    private void runGtfToBed(boolean transcript) throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", input)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        if (transcript) {
            Assert.assertEquals(countLines(outputFile), 1);
        } else {
            Assert.assertEquals(countLines(outputFile), 1);
        }
    }

    private void runMapk1(boolean transcript) throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", mapk1)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        if (transcript) {
            Assert.assertEquals(countLines(outputFile), 3);
            Assert.assertTrue(compareFiles(mapk1TranscriptBed, outputFile));
        } else {
            Assert.assertEquals(countLines(outputFile), 1);
            Assert.assertTrue(compareFiles(mapk1GeneBed, outputFile));
        }
    }

    private void runDecoySample(boolean transcript) throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", decoys)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        countLines(outputFile);
        runCommandLine(argsBuilder);
        if (transcript) {
            Assert.assertEquals(countLines(outputFile), 20);
            Assert.assertTrue(compareFiles(decoysTranscriptBed, outputFile));

        } else {
            Assert.assertEquals(countLines(outputFile), 19);
            Assert.assertTrue(compareFiles(decoysGeneBed, outputFile));

        }
    }

    private int countLines(File file) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        int lines = 0;
        while (reader.readLine() != null) {
            lines++;
        }
        return lines;
    }

    private boolean compareFiles(File file1, File file2) throws IOException {
        BufferedReader reader1 = new BufferedReader(new FileReader(file1));
        BufferedReader reader2 = new BufferedReader(new FileReader(file2));

        String line1 = reader1.readLine();
        String line2 = reader2.readLine();
        int lineNum = 1;

        while (line1 != null || line2 != null) {
            if (line1 == null || line2 == null) {
                System.out.println("File are different lenghts");
                return false;
            } else if (!line1.equals(line2)) {
                System.out.println("Difference at line " + lineNum);
                System.out.println("File1: " + line1);
                System.out.println("File2: " + line2);
                return false;
            }
            line1 = reader1.readLine();
            line2 = reader2.readLine();
            lineNum++;
        }
        System.out.println("Files are identical");
        return true;
    }
}
