package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class RevertBaseQualityScoresIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testRevertQualityScores() throws IOException {
        final File inputBam = new File(getTestDataDir(), "has.original.quals.bam");
        Assert.assertFalse(hasOriginalQualScores(inputBam));

        final File outputBam = createTempFile("testRevertQualityScores", ".bam");
        final String[] args = new String[] {
                "--input" , inputBam.getAbsolutePath(),
                "--output", outputBam.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertTrue(hasOriginalQualScores(outputBam));
    }

    private static boolean hasOriginalQualScores( final File bam ) throws IOException {
        try ( final SamReader reader = SamReaderFactory.makeDefault().open(bam) ) {
            for ( SAMRecord read : reader ) {
                final byte[] originalBaseQualities = read.getOriginalBaseQualities();
                Assert.assertNotNull(originalBaseQualities);
                final byte[] baseQualities = read.getBaseQualities();
                Assert.assertEquals(originalBaseQualities.length, baseQualities.length);
                for (int i = 0; i < originalBaseQualities.length; ++i) {
                    if (originalBaseQualities[i] != baseQualities[i]) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    @Test(expectedExceptions = UserException.class)
    public void testNoOriginalQuals() throws IOException {
        // This tool should fail when original quality scores are missing.
        final File inputBam = new File(getTestDataDir(), "without.original.quals.bam");

        final File outputBam = createTempFile("testRevertQualityScores", ".bam");
        final String[] args = new String[] {
                "--input" , inputBam.getAbsolutePath(),
                "--output", outputBam.getAbsolutePath()
        };
        runCommandLine(args);
    }

    //Regression test for https://github.com/broadinstitute/gatk/issues/1473
    @Test
    public void testReadThatsEntirelyInsertion() throws IOException {
        final File inputBam = new File(getTestDataDir(), "badRead.sam");
        Assert.assertFalse(hasOriginalQualScores(inputBam));

        final File outputBam = createTempFile("testRevertQualityScores", ".bam");
        final String[] args = new String[] {
                "--input" , inputBam.getAbsolutePath(),
                "--output", outputBam.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertTrue(hasOriginalQualScores(outputBam));
    }
}
