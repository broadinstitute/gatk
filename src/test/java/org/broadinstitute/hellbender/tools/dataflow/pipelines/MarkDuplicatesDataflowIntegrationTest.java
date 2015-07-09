package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.picard.sam.markduplicates.MarkDuplicatesIntegrationTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;


public class MarkDuplicatesDataflowIntegrationTest extends CommandLineProgramTest{

    @DataProvider(name = "md")
    public Object[][] md(){
        return new Object[][]{
          // Duplicates here are number of duplicate reads. To get number of
          // duplicate pairs, divide by 2.
          //{new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"), 20, 0},
          //{new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.bam"), 90, 6},
          //{new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.markedDups.bam"), 90, 6},  //90 total reads, 6 dups
          {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR, "optical_dupes.bam"), 4, 2},
          //{new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR, "optical_dupes_casava.bam"), 4, 2},
        };
    }

    @Test(groups = "dataflow", dataProvider = "md")
    public void testMarkDuplicatesDataflowIntegrationTestLocal(final File input, final long totalExpected, final long dupsExpected) throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(input.getPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        File outputFile = createTempFile("markdups", ".bam");
        args.add(outputFile.getAbsolutePath());

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        int totalReads = 0;
        int duplicateReads = 0;
        try ( final ReadsDataSource outputReads = new ReadsDataSource(outputFile) ) {
            for ( GATKRead read : outputReads ) {
                ++totalReads;

                if ( read.isDuplicate() ) {
                    ++duplicateReads;
                }
            }
        }

        Assert.assertEquals(totalReads, totalExpected, "Wrong number of reads in output BAM");
        Assert.assertEquals(duplicateReads, dupsExpected, "Wrong number of duplicate reads in output BAM");
        // Add assertion about optical duplicates
    }

}
