package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.AnnotateIntervalsIntegrationTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public final class AnnotatedIntervalCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir + "copynumber/formats/collections");
    private static final File ANNOTATED_INTERVALS_ALL_ANNOTATIONS_FILE = new File(TEST_SUB_DIR,
            "annotated-intervals-all-annotations.tsv");
    private static final File ANNOTATED_INTERVALS_EXTRA_ANNOTATION_FILE = new File(TEST_SUB_DIR,
            "annotated-intervals-extra-annotation.tsv");
    private static final File ANNOTATED_INTERVALS_REPEATED_ANNOTATION_FILE = new File(TEST_SUB_DIR,
            "annotated-intervals-repeated-annotation.tsv");
    private static final File ANNOTATED_INTERVALS_GC_CONTENT_ONLY_FILE = new File(TEST_SUB_DIR,
            "annotated-intervals-gc-content-only.tsv");

    private static final AnnotatedIntervalCollection EXPECTED_ALL_ANNOTATIONS = AnnotateIntervalsIntegrationTest.EXPECTED_ALL_ANNOTATIONS;

    @Test
    public void testRead() {
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(ANNOTATED_INTERVALS_ALL_ANNOTATIONS_FILE);
        Assert.assertEquals(result, EXPECTED_ALL_ANNOTATIONS);
        Assert.assertNotSame(result, EXPECTED_ALL_ANNOTATIONS);
    }

    /**
     * Extra annotations not listed in {@link CopyNumberAnnotations} should be ignored.
     */
    @Test
    public void testReadExtraAnnotation() {
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(ANNOTATED_INTERVALS_EXTRA_ANNOTATION_FILE);
        Assert.assertEquals(result, EXPECTED_ALL_ANNOTATIONS);
        Assert.assertNotSame(result, EXPECTED_ALL_ANNOTATIONS);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadRepeatedAnnotation() {
        new AnnotatedIntervalCollection(ANNOTATED_INTERVALS_REPEATED_ANNOTATION_FILE);
    }

    @Test
    public void testWrite() throws IOException {
        final File outputFile = createTempFile("annotated-interval-collection-test-output", ".tsv");
        EXPECTED_ALL_ANNOTATIONS.write(outputFile);
        Assert.assertTrue(FileUtils.contentEquals(outputFile, ANNOTATED_INTERVALS_ALL_ANNOTATIONS_FILE));
    }

    @Test
    public void testReadGCContentOnly() {
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(ANNOTATED_INTERVALS_GC_CONTENT_ONLY_FILE);
        final AnnotatedIntervalCollection expected = AnnotatedIntervalCollectionUnitTest.subsetAnnotations(
                EXPECTED_ALL_ANNOTATIONS,
                Collections.singletonList(CopyNumberAnnotations.GC_CONTENT));
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    private static AnnotatedInterval subsetAnnotations(final AnnotatedInterval annotatedInterval,
                                                      final List<AnnotationKey<?>> annotationKeys) {
        final List<Pair<AnnotationKey<?>, Object>> subsetAnnotationEntries = new ArrayList<>();
        for (final AnnotationKey<?> annotationKey : annotationKeys) {
            subsetAnnotationEntries.add(Pair.of(
                    annotationKey,
                    annotatedInterval.getAnnotationMap().getValue(annotationKey)));
        }
        final AnnotationMap subsetAnnotationMap = new AnnotationMap(subsetAnnotationEntries);
        return new AnnotatedInterval(annotatedInterval.getInterval(), subsetAnnotationMap);
    }

    public static AnnotatedIntervalCollection subsetAnnotations(final AnnotatedIntervalCollection annotatedIntervals,
                                                                final List<AnnotationKey<?>> annotationKeys) {
        return new AnnotatedIntervalCollection(
                annotatedIntervals.getMetadata(),
                annotatedIntervals.getRecords().stream()
                        .map(i -> subsetAnnotations(i, annotationKeys))
                        .collect(Collectors.toList()));
    }
}