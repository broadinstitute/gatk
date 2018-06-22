package org.broadinstitute.hellbender.tools.walkers;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class UnmarkDuplicatesIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "testUnmarkDuplicatesProvider")
    private Object[][] testUmiSetsDataProvider() {
        return new Object[][]{{
                "walkers/UnmarkDuplicates/allDuplicates.bam",
                11,
                null
        }, {
                "walkers/UnmarkDuplicates/allDuplicates.cram",
                11,
                new File("src/test/resources/hg19mini.fasta")
        } };
    }

    @Test(dataProvider = "testUnmarkDuplicatesProvider")
    public void testUnmarkDuplicates(String samFile, int duplicateCount, File referenceFile) throws IOException {
        final File inputBam = new File(getTestDataDir(), samFile);
        Assert.assertEquals(getDuplicateCountForBam(inputBam, referenceFile), duplicateCount, "Wrong number of duplicates in original input file");

        final File outputBam = createTempFile("testUnmarkDuplicates", ".bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputBam).addOutput(outputBam);

        if(referenceFile != null) {
            args.addReference(referenceFile);
        }

        runCommandLine(args);

        Assert.assertEquals(getDuplicateCountForBam(outputBam, referenceFile), 0, "Duplicates still present after UnmarkDuplicates");
    }

    private static int getDuplicateCountForBam(final File bam, final File referenceFile) throws IOException {
        int duplicateCount = 0;
        final SamReaderFactory factory = SamReaderFactory.makeDefault();
        if(referenceFile != null) {
            factory.referenceSequence(referenceFile);
        }
        try ( final SamReader reader = factory.open(bam) ) {
            for ( SAMRecord read : reader ) {
                if ( read.getDuplicateReadFlag() ) {
                    ++duplicateCount;
                }
            }
        }

        return duplicateCount;
    }
}
