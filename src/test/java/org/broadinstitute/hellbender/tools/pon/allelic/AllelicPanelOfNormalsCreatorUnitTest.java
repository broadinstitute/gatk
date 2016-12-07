package org.broadinstitute.hellbender.tools.pon.allelic;

import htsjdk.samtools.util.Log;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Tests for {@link AllelicPanelOfNormalsCreator}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormalsCreatorUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final FileFilter PULLDOWN_FILTER = new WildcardFileFilter("allelic-pon-test-pulldown-*tsv");
    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_SUB_DIR).listFiles(PULLDOWN_FILTER));

    private static final File PON_EXPECTED_SITE_FREQUENCY_50 = new File(TEST_SUB_DIR, "allelic-pon-test-pon-freq-50.tsv");
    private static final File PON_EXPECTED_SITE_FREQUENCY_75 = new File(TEST_SUB_DIR, "allelic-pon-test-pon-freq-75.tsv");

    @DataProvider(name = "dataCreate")
    public Object[][] dataCreate() {
        return new Object[][]{
                {0.5, AllelicPanelOfNormals.read(PON_EXPECTED_SITE_FREQUENCY_50)},
                {0.75, AllelicPanelOfNormals.read(PON_EXPECTED_SITE_FREQUENCY_75)}
        };
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyPulldownList() {
        new AllelicPanelOfNormalsCreator(new ArrayList<>());
    }

    @Test(dataProvider = "dataCreate")
    public void testCreate(final double siteFrequencyThreshold, final AllelicPanelOfNormals expected) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicPanelOfNormalsCreator allelicPoNCreator = new AllelicPanelOfNormalsCreator(PULLDOWN_FILES);
        final AllelicPanelOfNormals result = allelicPoNCreator.create(siteFrequencyThreshold);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(result, expected);
    }
}