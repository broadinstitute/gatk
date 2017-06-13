package org.broadinstitute.hellbender.tools.walkers.indels;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

public class RealignerTargetCreatorIntegrationTest extends CommandLineProgramTest {

    private final static File TEMP_DIR = createTempDir("RealignerTargetCreatorIntegrationTest");

    private final static File b37_reference_20_21_file = new File(b37_reference_20_21);
    private final static File NA12878_20_21_WGS_bam_file = new File(NA12878_20_21_WGS_bam);

    private GenomeLocParser b37GenomeLocParser;

    @BeforeClass
    public void initGenomeLocParserb37() throws IOException {
        try(CachingIndexedFastaSequenceFile reader = new CachingIndexedFastaSequenceFile(b37_reference_20_21_file)) {
            b37GenomeLocParser = new GenomeLocParser(reader);
        }
    }

    @DataProvider(name = "RealignerTargetCreator_Arguments")
    public Object[][] argumentsForTesting() {
        final String interval_20_9m_11m = "20:9,000,000-10,500,000";
        return new Object[][]{
                {"test_mismatch_0.15", new ArgumentsBuilder()
                        .addArgument("mismatchFraction", "0.15")
                        .addArgument("intervals", interval_20_9m_11m)
                }
        };
    }

    @Test(dataProvider = "RealignerTargetCreator_Arguments")
    public void testIntervalListOutput(String testName, ArgumentsBuilder argsToTest) throws Exception {
        final File expectedOutput = getTestFile(testName + ".interval_list");
        Assert.assertTrue(expectedOutput.exists(), "Test file not found: " + expectedOutput.getAbsolutePath());

        final File actualOutput = new File(TEMP_DIR, expectedOutput.getName());
        final ArgumentsBuilder arguments = new ArgumentsBuilder(argsToTest.getArgsArray())
                .addArgument("verbosity", "DEBUG")
                .addReference(b37_reference_20_21_file).addInput(NA12878_20_21_WGS_bam_file)
                .addOutput(actualOutput);

        runCommandLine(arguments);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test(dataProvider = "RealignerTargetCreator_Arguments")
    public void testTargetListOutputAgainstIntervalList(String testName, ArgumentsBuilder argsToTest) throws Exception {
        final File expectedOutput = getTestFile(testName + ".interval_list");
        Assert.assertTrue(expectedOutput.exists(), "Test file not found: " + expectedOutput.getAbsolutePath());

        final File actualOutput = new File(TEMP_DIR, testName + ".targets");
        final ArgumentsBuilder arguments = new ArgumentsBuilder(argsToTest.getArgsArray())
                .addReference(b37_reference_20_21_file).addInput(NA12878_20_21_WGS_bam_file)
                .addOutput(actualOutput);

        runCommandLine(arguments);

        final List<Interval> targetListResult = IntervalUtils.intervalFileToList(b37GenomeLocParser, actualOutput.getAbsolutePath())
                .stream().map(t -> new Interval(t.getContig(), t.getStart(), t.getStop()))
                .collect(Collectors.toList());

        final List<Interval> expectedListResult = IntervalList.fromFile(expectedOutput).getIntervals();

        Assert.assertFalse(targetListResult.isEmpty());
        Assert.assertFalse(expectedListResult.isEmpty());
        Assert.assertEquals(targetListResult, expectedListResult);
    }
}
