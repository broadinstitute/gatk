package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.picard.sam.markduplicates.MarkDuplicatesIntegrationTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Objects;


public class MarkDuplicatesDataflowIntegrationTest extends CommandLineProgramTest{

    @DataProvider(name = "md")
    public Object[][] md(){
        return new Object[][]{
                {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"), 20, 0},
                {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.bam"), 90, 6},
                {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.markedDups.bam"), 90, 6},  //90 total reads, 6 dups
        };
    }

    @Test(groups = "dataflow", dataProvider = "md")
    public void testMarkDuplicatesDataflowIntegrationTestLocal(final File input, final long totalExpected, final long dupsExpected) throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(input.getPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        File placeHolder = createTempFile("markdups", ".txt");
        args.add(placeHolder.getPath());

        runCommandLine(args.getArgsArray());
        File outputFile = findDataflowOutput(placeHolder);

        Assert.assertTrue(outputFile.exists());
        Assert.assertEquals(new XReadLines(outputFile).readLines().size(), totalExpected);

        Assert.assertEquals(new XReadLines(outputFile).readLines().stream().filter(line -> line.contains("duplicateFragment\":true")).count(), dupsExpected);
    }
}