package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Unit tests for {@link TargetTableWriter}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetTableWriterUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullAnnotationSet() throws IOException {
        final File testFile = BaseTest.createTempFile("ttw-test", ".tsv");
        new TargetTableWriter(testFile, null).close();
    }

    @Test(dataProvider="correctAnnotationSets")
    public void testSimpleWrite(final Set<TargetAnnotation> annotationSet) throws IOException {
        final File testFile = BaseTest.createTempFile("ttw-test", ".tsv");
        final int TARGET_COUNT = 10;

        final TargetTableWriter subject = new TargetTableWriter(testFile, annotationSet);
        for (int i = 0; i < TARGET_COUNT; i++) {
            final Map<TargetAnnotation, String> annotationMap = new HashMap<>(annotationSet.size());
            for (final TargetAnnotation annotation : annotationSet) {
                annotationMap.put(annotation, annotation.name() + "_" + i);
            }
            final TargetAnnotationCollection annotationCollection = TargetAnnotationCollection.fromMap(annotationMap);
            subject.writeRecord(new Target("target_" + i, new SimpleInterval("1", (i + 1) * 100, (i + 1) * 100 + 50), annotationCollection));
        }
        subject.close();
        // Now check the contents are correct using a reader.
        final TargetTableReader reader = new TargetTableReader(testFile);
        for (int i = 0; i < TARGET_COUNT; i++) {
            final Target target = reader.readRecord();
            Assert.assertEquals(target.getName(), "target_" + i);
            Assert.assertEquals(target.getInterval(), new SimpleInterval("1", (i + 1) * 100, (i + 1) * 100 + 50));
            final TargetAnnotationCollection annotationCollection = target.getAnnotations();
            Assert.assertEquals(annotationCollection.size(), annotationSet.size());
            Assert.assertEquals(annotationCollection.annotationSet(), annotationSet);
            for (final TargetAnnotation annotation : annotationSet) {
                Assert.assertEquals(annotationCollection.get(annotation), annotation.name() + "_" + i);
            }
        }
        reader.close();
    }

    @Test(dataProvider="wrongTargetAnnotationSetPairs", expectedExceptions = NoSuchElementException.class)
    public void testWrongTargetAnnotations(final Set<TargetAnnotation> expectedSet, final Set<TargetAnnotation> actualAnnotationSet) throws IOException {
        final File testFile = BaseTest.createTempFile("ttw-test", ".tsv");
        final int TARGET_COUNT = 5;

        final TargetTableWriter subject = new TargetTableWriter(testFile, expectedSet);
        for (int i = 0; i < TARGET_COUNT; i++) {
            final Map<TargetAnnotation, String> annotationMap = new HashMap<>(expectedSet.size());
            for (final TargetAnnotation annotation : i == TARGET_COUNT - 1 ? actualAnnotationSet : expectedSet) {
                annotationMap.put(annotation, annotation.name() + "_" + i);
            }
            final TargetAnnotationCollection annotationCollection = TargetAnnotationCollection.fromMap(annotationMap);
            subject.writeRecord(new Target("target_" + i, new SimpleInterval("1", (i + 1) * 100, (i + 1) * 100 + 50), annotationCollection));
        }
        subject.close();
    }


    @Test(dataProvider="extraTargetAnnotationSetPairs")
    public void testExtraTargetAnnotations(final Set<TargetAnnotation> expectedSet, final Set<TargetAnnotation> actualAnnotationSet) throws IOException {
        final File testFile = BaseTest.createTempFile("ttw-test", ".tsv");
        final int TARGET_COUNT = 5;

        final TargetTableWriter subject = new TargetTableWriter(testFile, expectedSet);
        for (int i = 0; i < TARGET_COUNT; i++) {
            final Map<TargetAnnotation, String> annotationMap = new HashMap<>(expectedSet.size());
            for (final TargetAnnotation annotation : i == TARGET_COUNT - 1 ? actualAnnotationSet : expectedSet) {
                annotationMap.put(annotation, annotation.name() + "_" + i);
            }
            final TargetAnnotationCollection annotationCollection = TargetAnnotationCollection.fromMap(annotationMap);
            subject.writeRecord(new Target("target_" + i, new SimpleInterval("1", (i + 1) * 100, (i + 1) * 100 + 50), annotationCollection));
        }
        subject.close();
    }

    @DataProvider(name="correctAnnotationSets")
    public Object[][] correctAnnotationSets() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] { Stream.of(TargetAnnotation.values()).collect(Collectors.toSet()) });
        result.add(new Object[] { Collections.emptySet()});
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            result.add(new Object[] { Collections.singleton(annotation) });
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "wrongTargetAnnotationSetPairs")
    public Object[][] wrongTargetAnnotationSetPairs() {
        final List<Object[]> result = new ArrayList<>();
        final Set<TargetAnnotation> all = Stream.of(TargetAnnotation.values()).collect(Collectors.toSet());
        result.add(new Object[] { all, Collections.emptySet() });
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            result.add(new Object[] { all, Collections.singleton(annotation) });
        }
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            for (final TargetAnnotation annotation2 : TargetAnnotation.values()) {
                if (annotation2 == annotation) {
                    continue;
                }
                result.add(new Object[] { Collections.singleton(annotation), Collections.singleton(annotation2)});
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "extraTargetAnnotationSetPairs")
    public Object[][] extraTargetAnnotationSetPairs() {
        final List<Object[]> result = new ArrayList<>();
        final Set<TargetAnnotation> all = Stream.of(TargetAnnotation.values()).collect(Collectors.toSet());
        result.add(new Object[] { Collections.emptySet(), all });
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            result.add(new Object[] { Collections.singleton(annotation), all });
        }
        return result.toArray(new Object[result.size()][]);
    }
}
