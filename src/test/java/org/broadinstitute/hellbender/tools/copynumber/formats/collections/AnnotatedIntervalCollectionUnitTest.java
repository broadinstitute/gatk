package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);

    private static final SAMSequenceDictionary SEQUENCE_DICTIONARY = ReferenceDataSource.of(REFERENCE_FILE.toPath()).getSequenceDictionary();
    private static final LocatableMetadata LOCATABLE_METADATA = new SimpleLocatableMetadata(SEQUENCE_DICTIONARY);

    private static final AnnotatedIntervalCollection EXPECTED_ALL_ANNOTATIONS = new AnnotatedIntervalCollection(
            LOCATABLE_METADATA,
            Arrays.asList(
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
        final AnnotatedIntervalCollection result = new AnnotatedIntervalCollection(ANNOTATED_INTERVALS_REPEATED_ANNOTATION_FILE);
        Assert.assertEquals(result, EXPECTED_ALL_ANNOTATIONS);
        Assert.assertNotSame(result, EXPECTED_ALL_ANNOTATIONS);
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