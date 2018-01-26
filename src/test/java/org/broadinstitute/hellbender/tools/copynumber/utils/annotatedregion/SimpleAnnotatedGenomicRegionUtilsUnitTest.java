package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.ImmutableSortedMap;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.spark_project.guava.collect.Lists;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class SimpleAnnotatedGenomicRegionUtilsUnitTest extends GATKBaseTest {
    @Test(dataProvider = "mergeTests")
    public void testMerge(List<SimpleAnnotatedGenomicRegion> test1, List<SimpleAnnotatedGenomicRegion> gt) {
        List<SimpleAnnotatedGenomicRegion> mergeTestResults = SimpleAnnotatedGenomicRegionUtils.mergeRegions(test1,
                ReferenceUtils.loadFastaDictionary(new File(ReferenceUtils.getFastaDictionaryFileName(hg19MiniReference))),
                "__", l -> {});
        Assert.assertEquals(mergeTestResults, gt);
    }

    @DataProvider(name = "mergeTests")
    public Object [][] createMergeTests() {

        return new Object[][] {
                {
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                                new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1"))),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                        ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")))
                }, {
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 201, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1"))),
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 201, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")))
        }, {
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("2", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foobar1", "bar1"))),
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("2", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foobar1", "bar1")))
        }, {
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 190, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1"))),
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1__bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1")))
        }, {
                Lists.newArrayList(
                        // Same as previous test, but trying to see if input order hoses things.  Note that the output will
                        //  be sorted.
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 190, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1"))),
                Lists.newArrayList(
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1__bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1")))
        }
        };
    }

}
