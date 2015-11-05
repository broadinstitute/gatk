package org.broadinstitute.hellbender.tools.walkers;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class UnmarkDuplicatesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testUnmarkDuplicates() throws IOException {
        final File inputBam = new File(getTestDataDir(), "walkers/UnmarkDuplicates/allDuplicates.bam");
        Assert.assertEquals(getDuplicateCountForBam(inputBam), 11, "Wrong number of duplicates in original input file");

        final File outputBam = createTempFile("testUnmarkDuplicates", ".bam");
        final String[] args = new String[] {
                "--input" , inputBam.getAbsolutePath(),
                "--output", outputBam.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertEquals(getDuplicateCountForBam(outputBam), 0, "Duplicates still present after UnmarkDuplicates");
    }

    private static int getDuplicateCountForBam( final File bam ) throws IOException {
        int duplicateCount = 0;
        try ( final SamReader reader = SamReaderFactory.makeDefault().open(bam) ) {
            for ( SAMRecord read : reader ) {
                if ( read.getDuplicateReadFlag() ) {
                    ++duplicateCount;
                }
            }
        }

        return duplicateCount;
    }
}
