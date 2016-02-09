package org.broadinstitute.hellbender.tools.exome.samplenamefinder;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class SampleNameFinderTest extends BaseTest {
    static final File INPUT_DIR = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/");
    static final File INPUT_FILE = new File(INPUT_DIR, "HCC1143_reduced_log.tsv");

    @Test
    public void testSampleNameFinder() {

        final List<String> guess = SampleNameFinder.determineSampleNamesFromTargetCoverageFile(INPUT_FILE);
        Assert.assertEquals(guess.size(), 1);
        Assert.assertEquals(guess.get(0), "HCC1143");

    }
}
