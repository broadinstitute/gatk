package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.copynumber.annotation.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.annotation.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.annotation.AnnotationSet;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Integration tests for {@link AnnotateIntervals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AnnotateIntervalsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/";
    private static final File INTERVALS_FILE = new File(TEST_SUB_DIR, "annotate-intervals-test.interval_list");
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);

    /**
     * Test that intervals are sorted according to {@link LocatableCollection#LEXICOGRAPHICAL_ORDER_COMPARATOR}
     * and adjacent intervals are not merged.  GC content truth was taken from AnnotateTargets (a previous version of the tool).
     */
    @Test
    public void test() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument("L",  INTERVALS_FILE.getAbsolutePath())
                .addArgument("imr", IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(outputFile);

        final AnnotatedIntervalCollection expected = new AnnotatedIntervalCollection(Arrays.asList(
                new AnnotatedInterval(new SimpleInterval("20", 1000001,	1001000), new AnnotationSet(0.49)),
                new AnnotatedInterval(new SimpleInterval("20", 1001001,	1002000), new AnnotationSet(0.483)),
                new AnnotatedInterval(new SimpleInterval("20", 1002001,	1003000), new AnnotationSet(0.401)),
                new AnnotatedInterval(new SimpleInterval("20", 1003001,	1004000), new AnnotationSet(0.448)),
                new AnnotatedInterval(new SimpleInterval("21", 1,	100), new AnnotationSet(Double.NaN)),
                new AnnotatedInterval(new SimpleInterval("21", 101,	200), new AnnotationSet(Double.NaN))));
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalSetRule() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument("L", INTERVALS_FILE.getAbsolutePath())
                .addArgument("isr", IntervalSetRule.INTERSECTION.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalExclusionPadding() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument("L", INTERVALS_FILE.getAbsolutePath())
                .addArgument("ixp", "1")
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalPadding() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument("L", INTERVALS_FILE.getAbsolutePath())
                .addArgument("ip", "1")
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalMergingRule() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument("L", INTERVALS_FILE.getAbsolutePath())
                .addArgument("imr", IntervalMergingRule.ALL.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }
}