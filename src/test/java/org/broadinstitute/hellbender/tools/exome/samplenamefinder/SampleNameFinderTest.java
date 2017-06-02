package org.broadinstitute.hellbender.tools.exome.samplenamefinder;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class SampleNameFinderTest extends BaseTest {
    static final File INPUT_DIR = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/");
    static final File INPUT_FILE = new File(INPUT_DIR, "HCC1143_reduced_log.tsv");
    static final File TWO_SAMPLE_INPUT_FILE = new File(INPUT_DIR, "HCC1143_short_2samples.tsv");
    static final File DUPE_SAMPLE_INPUT_FILE = new File(INPUT_DIR, "HCC1143_short_dupe_sample.tsv");
    static final File NO_SAMPLE_INPUT_FILE = new File(INPUT_DIR, "HCC1143_short_no_samples.tsv");


    @Test
    public void testSampleNameFinderEasy() {
        final List<String> guess = SampleNameFinder.determineSampleNamesFromReadCountsFile(INPUT_FILE);
        Assert.assertEquals(guess.size(), 1);
        Assert.assertEquals(guess.get(0), "HCC1143");
    }

    @Test
    public void testSampleNameFinder2Sample() {
        final List<String> guess = SampleNameFinder.determineSampleNamesFromReadCountsFile(TWO_SAMPLE_INPUT_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testSampleNameFinderDupe() {
        final List<String> guess = SampleNameFinder.determineSampleNamesFromReadCountsFile(DUPE_SAMPLE_INPUT_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testSampleNameFinderNoSample() {
        final List<String> guess = SampleNameFinder.determineSampleNamesFromReadCountsFile(NO_SAMPLE_INPUT_FILE);
    }
}
