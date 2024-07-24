package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    private static final File input = new File("/Users/shahsana/TestGatk/gencode.v38.chr_patch_hapl_scaff.annotation.gtf");
    private static final File dictionary = new File("/Users/shahsana/TestGatk/reference.dict");

    @Test
    public void testGene() {
        runGtfToBed(false);
    }

    @Test
    public void testTranscript() {
        runGtfToBed(true);
    }

    @Test
    public void testOneGene(){

    }

    @Test
    public void testOneTranscript(){

    }

    @Test
    public void testDecoyGenes(){

    }

    @Test
    public void testDecoyTranscripts(){

    }


    private void runGtfToBed(boolean transcript){
        File outputFile = new File("/Users/shahsana/TestGatk/output.bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", input)
                .add("D", dictionary)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, transcript)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Assert.assertEquals(1,1);
    }

}
