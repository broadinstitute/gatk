package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class RevertQualityScoresIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testRevertQualityScores() throws IOException {
        final File inputBam = new File(getTestDataDir(), "BQSR/expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam");
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
        boolean hasOriginalQuals = true;
        try ( final SamReader reader = SamReaderFactory.makeDefault().open(bam) ) {
            for ( SAMRecord read : reader ) {
                final byte[] originalBaseQualities = read.getOriginalBaseQualities();
                final byte[] baseQualities = read.getBaseQualities();
                Assert.assertEquals(originalBaseQualities.length, baseQualities.length);
                for (int i = 0; i < originalBaseQualities.length; ++i) {
                    if (originalBaseQualities[i] != baseQualities[i]) {
                        hasOriginalQuals = false;
                    }
                }
            }
        }
        return hasOriginalQuals;
    }

    @Test(expectedExceptions = GATKException.class)
    public void testNoOriginalQuals() throws IOException {
        // This tool should fail when original quality scores are missing.
        final File inputBam = new File(getTestDataDir(), "count_reads.bam");

        final File outputBam = createTempFile("testRevertQualityScores", ".bam");
        final String[] args = new String[] {
                "--input" , inputBam.getAbsolutePath(),
                "--output", outputBam.getAbsolutePath()
        };
        runCommandLine(args);
    }
}