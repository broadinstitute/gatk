package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import com.google.common.collect.ImmutableSortedMap;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.collections.Sets;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

public final class AnnotatedIntervalUnitTest extends GATKBaseTest {
    private static final File TEST_FILE = new File(toolsTestDir,
            "copynumber/utils/annotatedinterval/annotated-interval-many-columns.seg");

    @Test
    public void basicTest() throws IOException {
        final Set<String> headersOfInterest = Sets.newHashSet(Arrays.asList("name", "learning_SAMPLE_0"));
        final List<AnnotatedInterval> annotatedIntervals =
                AnnotatedIntervalCollection.create(TEST_FILE.toPath(), headersOfInterest).getRecords();

        Assert.assertEquals(annotatedIntervals.size(), 15);
        Assert.assertTrue(annotatedIntervals.stream()
                .mapToInt(s -> s.getAnnotations().entrySet().size())
                .allMatch(i -> i == headersOfInterest.size()));
        Assert.assertTrue(annotatedIntervals.stream().allMatch(s -> s.getAnnotations().keySet().containsAll(headersOfInterest)));

        // Grab the first 15 and test values
        List<AnnotatedInterval> gtRegions = Arrays.asList(
                new AnnotatedInterval(new SimpleInterval("1", 30365, 30503), ImmutableSortedMap.of("name", "target_1_None", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 69088, 70010), ImmutableSortedMap.of("name", "target_2_None", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 367656, 368599), ImmutableSortedMap.of("name", "target_3_None", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 621093, 622036), ImmutableSortedMap.of("name", "target_4_None", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 861319, 861395), ImmutableSortedMap.of("name", "target_5_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 865532, 865718), ImmutableSortedMap.of("name", "target_6_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 866416, 866471), ImmutableSortedMap.of("name", "target_7_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 871149, 871278), ImmutableSortedMap.of("name", "target_8_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 874417, 874511), ImmutableSortedMap.of("name", "target_9_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 874652, 874842), ImmutableSortedMap.of("name", "target_10_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 876521, 876688), ImmutableSortedMap.of("name", "target_11_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 877513, 877633), ImmutableSortedMap.of("name", "target_12_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 877787, 877870), ImmutableSortedMap.of("name", "target_13_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 877936, 878440), ImmutableSortedMap.of("name", "target_14_SAMD11", "learning_SAMPLE_0", "2")),
                new AnnotatedInterval(new SimpleInterval("1", 878630, 878759), ImmutableSortedMap.of("name", "target_15_SAMD11", "learning_SAMPLE_0", "2"))
        );

        Assert.assertEquals(annotatedIntervals.subList(0, gtRegions.size()), gtRegions);
    }

    @Test
    public void basicTestWithAllColumnsFile() throws IOException {

        // If no columns of interest are given in a read call, the method will try to load all columns as "interesting".
        final List<AnnotatedInterval> annotatedIntervals =
                AnnotatedIntervalCollection.create(TEST_FILE.toPath(), null).getRecords();

        Assert.assertEquals(annotatedIntervals.size(), 15);
        Assert.assertTrue(annotatedIntervals.stream()
                .mapToInt(s -> s.getAnnotations().entrySet().size())
                .allMatch(i -> i == 101)); // The number of columns in the TEST_FILE (name, learning_SAMPLE_0...99
    }

    @Test(dataProvider = "equalsTests")
    public void testEquals(AnnotatedInterval test1, AnnotatedInterval test2, boolean gt) {
        Assert.assertEquals(test1.equals(test2), gt);
    }

    @DataProvider(name = "equalsTests")
    public Object [][] createEqualsTests() {

        return new Object[][] {
            {
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                true
            }, {
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                false
            }, {
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foobar1", "bar1")),
                false
            }, {
                new AnnotatedInterval( new SimpleInterval("1", 1000, 2000),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                false
            }, {
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                        ImmutableSortedMap.of("Foo", "bar")),
                false
            }
        };
    }
}
