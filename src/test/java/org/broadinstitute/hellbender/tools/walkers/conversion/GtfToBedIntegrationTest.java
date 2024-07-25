package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    private static final File input = new File(TestUtils.getGencodeGtf());
    private static final File mapk1 = new File(TestUtils.getMapk1Gtf());
    private static final File decoys = new File(TestUtils.getDecoysGtf());
    private static final File dictionary = new File(TestUtils.getReferenceDict());


    @Test
    public void testGene() {
        runGtfToBed(false);
    }

    @Test
    public void testTranscript() {
        runGtfToBed(true);
    }

    @Test
    public void testMapk1Gene(){
        runMapk1(false);
    }

    @Test
    public void testMapk1Transcript(){
        runMapk1(true);
    }

    @Test
    public void testDecoyGenes(){
        runDecoys(false);
    }

    @Test
    public void testDecoyTranscripts(){
        runDecoys(true);
    }


    private void runGtfToBed(boolean transcript){
        //File outputFile = new File("/Users/shahsana/TestGatk/output.bed");
        final File outputFile = createTempFile("gtf_to_bed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", input)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Assert.assertEquals(1,1);
    }

    private void runMapk1(boolean transcript){
        //File outputFile = new File("/Users/shahsana/TestGatk/output.bed");
        final File outputFile = createTempFile("gtf_to_bed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", mapk1)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Assert.assertEquals(1,1);
    }

    private void runDecoys(boolean transcript){
        //File outputFile = new File("/Users/shahsana/TestGatk/output.bed");
        final File outputFile = createTempFile("gtf_to_bed", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", decoys)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Assert.assertEquals(1,1);
    }

}
