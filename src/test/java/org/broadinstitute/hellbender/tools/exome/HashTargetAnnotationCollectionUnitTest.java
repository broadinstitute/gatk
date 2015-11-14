package org.broadinstitute.hellbender.tools.exome;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Unit tests for {@link HashTargetAnnotationCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HashTargetAnnotationCollectionUnitTest {

    @Test
    public void testCreateNoValueInstance() {
         final TargetAnnotationCollection subject = new HashTargetAnnotationCollection(Collections.emptyMap());
         Assert.assertEquals(subject.size(), 0);
         Assert.assertEquals(subject.annotationSet(), Collections.emptySet());
    }

    @Test(dataProvider = "annotations")
    public void testSingleAnnotationCreation(final TargetAnnotation annotation) {
         final TargetAnnotationCollection subject = new HashTargetAnnotationCollection(Collections.singletonMap(annotation, "text"));
         Assert.assertEquals(subject.size(), 1);
         Assert.assertTrue(subject.hasAnnotation(annotation));
         for (final TargetAnnotation annotation2 : TargetAnnotation.values()) {
             if (annotation != annotation2) {
                 Assert.assertFalse(subject.hasAnnotation(annotation2));
             }
         }
         Assert.assertEquals(subject.get(annotation), "text");
         Assert.assertEquals(subject.annotationSet(), Collections.singleton(annotation));
    }

    @Test
    public void testMultipleAnnotationCreation() {
        final Map<TargetAnnotation, String> values = Stream.of(TargetAnnotation.values()).collect(Collectors.toMap(a -> a, a -> "text"));
        final TargetAnnotationCollection subject = new HashTargetAnnotationCollection(values);
        // change values in input map to check that that does not affect the annotation collection.
        for (final TargetAnnotation annotation : values.keySet()) {
            values.put(annotation, "text2");
        }
        Assert.assertEquals(subject.size(), TargetAnnotation.values().length);
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            Assert.assertTrue(subject.hasAnnotation(annotation));
            Assert.assertEquals(subject.get(annotation), "text");
            try {
                subject.getDouble(annotation);
                Assert.fail("expected an exception");
            } catch (final IllegalStateException ex) {
                // good.
            } catch (final RuntimeException ex) {
                Assert.fail("Wrong exception type " + ex.getClass().getName() + ", expected is " + IllegalStateException.class);
            }
        }
        Assert.assertEquals(subject.annotationSet(), Stream.of(TargetAnnotation.values()).collect(Collectors.toSet()));
    }

    @Test
    public void testMultipleDoubleAnnotationCreation() {
        final Map<TargetAnnotation, String> values = Stream.of(TargetAnnotation.values()).collect(Collectors.toMap(a -> a, a -> "-1.1"));
        final TargetAnnotationCollection subject = new HashTargetAnnotationCollection(values);
        // change values in input map to check that that does not affect the annotation collection.
        for (final TargetAnnotation annotation : values.keySet()) {
            values.put(annotation, "text2");
        }
        Assert.assertEquals(subject.size(), TargetAnnotation.values().length);
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            Assert.assertTrue(subject.hasAnnotation(annotation));
            Assert.assertEquals(subject.get(annotation), "-1.1");
            Assert.assertEquals(subject.getDouble(annotation), -1.1);
        }
    }

    @Test(dataProvider = "annotations", expectedExceptions = NoSuchElementException.class)
    public void testAccessNoValueInstanceAnnotation(final TargetAnnotation annotation) {
        final TargetAnnotationCollection subject = new HashTargetAnnotationCollection(Collections.emptyMap());
        subject.get(annotation);
    }

    @DataProvider(name = "annotations")
    public Object[][] annotations() {
        return Stream.of(TargetAnnotation.values()).map(a -> new Object[] { a } ).toArray(Object[][]::new);
    }
}
