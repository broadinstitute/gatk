package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.IntegerReadCountFileComparator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Integration test for {@link CollectFragmentCounts}.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class CollectFragmentCountsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/collectfragmentcounts";

    private static final File NA12878_BAM = new File(TEST_SUB_DIR, "exome-read-counts-NA12878.bam");
    private static final File NA12878_FRAGMENT_COUNT_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "fragment-counts-NA12878-output.tsv");
    private static final File SITES_FILE = new File(TEST_SUB_DIR, "read-counts-test.interval_list");

    @DataProvider(name = "testData")
    public Object[][] testData() {
        return new Object[][] {
                {NA12878_BAM, NA12878_FRAGMENT_COUNT_EXPECTED_OUTPUT}
        };
    }

    @Test(dataProvider = "testData")
    public void test(final File inputBAMFile, final File expectedOutputFile) {
        final File actualOutputFile = createTempFile("fragment-counts-test-file", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBAMFile.getAbsolutePath(),
                "-L", SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, actualOutputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        IntegerReadCountFileComparator.AssertEquals(expectedOutputFile, actualOutputFile, null);
    }

    @DataProvider(name = "artificialReadsData")
    public Object[][] artificialReadsTestData() {
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final List<GATKRead> readPairFirstInPairIsLeft = ArtificialReadUtils.createPair(samHeader, "firstReadPair", 101, 1000, 1300, true, false);
        final GATKRead firstInPairPositive = readPairFirstInPairIsLeft.get(0);

        final List<GATKRead> readPairFirstInPairIsRight = ArtificialReadUtils.createPair(samHeader, "secondReadPair", 101, 800, 1000, false, true);
        final GATKRead firstInPairNegative = readPairFirstInPairIsRight.get(1);

        return new Object[][] {
                {firstInPairPositive, 1200},
                {firstInPairNegative, 950}
        };
    }

    @Test(dataProvider = "artificialReadsData")
    public void testFragmentCenterComputation(final GATKRead read, final int expectedFragmentCenter) {

        final int actualFragmentCenter = CollectFragmentCounts.ReadOrientation.getCenterOfFragment(read);
        Assert.assertEquals(actualFragmentCenter, expectedFragmentCenter);
    }

}