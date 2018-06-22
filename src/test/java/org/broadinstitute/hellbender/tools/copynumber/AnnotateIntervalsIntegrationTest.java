package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotationSet;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
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
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File INTERVALS_FILE = new File(TEST_SUB_DIR, "annotate-intervals-test.interval_list");
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);

    private static final SAMSequenceDictionary SEQUENCE_DICTIONARY = ReferenceDataSource.of(REFERENCE_FILE.toPath()).getSequenceDictionary();
    private static final LocatableMetadata LOCATABLE_METADATA = new SimpleLocatableMetadata(SEQUENCE_DICTIONARY);

    /**
     * Test that intervals are sorted according to {@link #SEQUENCE_DICTIONARY}
     * and adjacent intervals are not merged.  GC content truth was taken from AnnotateTargets (a previous version of the tool).
     */
    @Test
    public void test() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(outputFile);

        final AnnotatedIntervalCollection expected = new AnnotatedIntervalCollection(
                LOCATABLE_METADATA,
                Arrays.asList(
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
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME, IntervalSetRule.INTERSECTION.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalExclusionPadding() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_EXCLUSION_PADDING_LONG_NAME,"1")
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalPadding() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME,"1")
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalMergingRule() {
        final File resultOutputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.ALL.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }
}