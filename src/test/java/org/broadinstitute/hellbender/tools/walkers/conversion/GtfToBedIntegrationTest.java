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
    private static final File manyTranscripts = new File(ConversionTestUtils.getManyTranscriptsGtf());
    private static final File dictionary = new File(ConversionTestUtils.getReferenceDict());
    private static final File decoysGeneBed = new File(ConversionTestUtils.getDecoySamplesGeneBed());
    private static final File decoysTranscriptBed = new File(ConversionTestUtils.getDecoySamplesTranscriptBed());
    private static final File mapk1GeneBed = new File(ConversionTestUtils.getMapk1GeneBed());
    private static final File mapk1TranscriptBed = new File(ConversionTestUtils.getMapk1TranscriptBed());
    private static final File manyTranscriptsBed = new File (ConversionTestUtils.getManyTranscriptsBed());

    // tests any Gtf file (gene)
    @Test
    public void testGene() throws IOException {
        runGtfToBed(false);
    }

    // tests any Gtf file (transcript)
    @Test
    public void testTranscript() throws IOException {
        runGtfToBed(true);
    }

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


    private void runGtfToBed(boolean transcript) throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", input)
                .add("SD", dictionary)
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
                .add("SD", dictionary)
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
                .add("SD", dictionary)
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

    private void runManyTranscripts() throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", manyTranscripts)
                .add("SD", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, true)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Assert.assertTrue(compareFiles( manyTranscriptsBed, outputFile));

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
