package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollectionUnitTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

/**
 * Integration tests for {@link AnnotateIntervals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AnnotateIntervalsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File INTERVALS_FILE = new File(TEST_SUB_DIR, "annotate-intervals-test.interval_list");
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);
    private static final File MAPPABILITY_TRACK_FILE = new File(TEST_SUB_DIR,
            "annotate-intervals-hg19-umap-k100-single-read-mappability-merged-20-21.bed.gz");
    private static final File SEGMENTAL_DUPLICATION_TRACK_FILE = new File(TEST_SUB_DIR,
            "annotate-intervals-hg19-segmental-duplication-20-21.bed.gz");
    private static final File SEGMENTAL_DUPLICATION_WITH_OVERLAPS_TRACK_FILE = new File(TEST_SUB_DIR,
            "annotate-intervals-hg19-segmental-duplication-20-21-with-overlaps.bed.gz");

    private static final SAMSequenceDictionary SEQUENCE_DICTIONARY = ReferenceDataSource.of(REFERENCE_FILE.toPath()).getSequenceDictionary();
    private static final LocatableMetadata LOCATABLE_METADATA = new SimpleLocatableMetadata(SEQUENCE_DICTIONARY);

    /**
     * Test case checks that intervals are sorted according to {@link #SEQUENCE_DICTIONARY} and
     * adjacent intervals are not merged.  This test case is also used in {@link AnnotatedIntervalCollectionUnitTest}.
     */
    public static final AnnotatedIntervalCollection EXPECTED_ALL_ANNOTATIONS = new AnnotatedIntervalCollection(
            LOCATABLE_METADATA,
            Arrays.asList(
                    new AnnotatedInterval(new SimpleInterval("20", 50001,	70000),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.3962),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.5),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.)))),
                    new AnnotatedInterval(new SimpleInterval("20", 1000001,	1001000),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.49),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 1.),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.)))),
                    new AnnotatedInterval(new SimpleInterval("20", 1001001,	1002000),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.483),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 1.),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.)))),
                    new AnnotatedInterval(new SimpleInterval("20", 1002001,	1003000),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.401),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 1.),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.)))),
                    new AnnotatedInterval(new SimpleInterval("20", 1003001,	1004000),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.448),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 1.),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.)))),
                    new AnnotatedInterval(new SimpleInterval("20", 1238001,	1239000),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.3),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 1.),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.218)))),
                    new AnnotatedInterval(new SimpleInterval("21", 1,	100),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, Double.NaN),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.0),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.)))),
                    new AnnotatedInterval(new SimpleInterval("21", 101,	200),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, Double.NaN),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.0),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.))))));

    @Test
    public void testGCContentOnly() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(outputFile);
        final AnnotatedIntervalCollection expected = AnnotatedIntervalCollectionUnitTest.subsetAnnotations(
                EXPECTED_ALL_ANNOTATIONS,
                Collections.singletonList(CopyNumberAnnotations.GC_CONTENT));
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @Test
    public void testMappability() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addFileArgument(AnnotateIntervals.MAPPABILITY_TRACK_PATH_LONG_NAME, MAPPABILITY_TRACK_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(outputFile);
        final AnnotatedIntervalCollection expected = AnnotatedIntervalCollectionUnitTest.subsetAnnotations(
                EXPECTED_ALL_ANNOTATIONS,
                Arrays.asList(
                        CopyNumberAnnotations.GC_CONTENT,
                        CopyNumberAnnotations.MAPPABILITY));
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @Test
    public void testSegmentalDuplicationContent() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addFileArgument(AnnotateIntervals.SEGMENTAL_DUPLICATION_TRACK_PATH_LONG_NAME, SEGMENTAL_DUPLICATION_TRACK_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(outputFile);
        final AnnotatedIntervalCollection expected = AnnotatedIntervalCollectionUnitTest.subsetAnnotations(
                EXPECTED_ALL_ANNOTATIONS,
                Arrays.asList(
                        CopyNumberAnnotations.GC_CONTENT,
                        CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT));
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @Test
    public void testAllAnnotations() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addFileArgument(AnnotateIntervals.MAPPABILITY_TRACK_PATH_LONG_NAME, MAPPABILITY_TRACK_FILE)
                .addFileArgument(AnnotateIntervals.SEGMENTAL_DUPLICATION_TRACK_PATH_LONG_NAME, SEGMENTAL_DUPLICATION_TRACK_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(outputFile);
        final AnnotatedIntervalCollection expected = EXPECTED_ALL_ANNOTATIONS;
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testSegmentalDuplicationContentWithOverlaps() {
        final File outputFile = createTempFile("annotate-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addFileArgument(AnnotateIntervals.SEGMENTAL_DUPLICATION_TRACK_PATH_LONG_NAME, SEGMENTAL_DUPLICATION_WITH_OVERLAPS_TRACK_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
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