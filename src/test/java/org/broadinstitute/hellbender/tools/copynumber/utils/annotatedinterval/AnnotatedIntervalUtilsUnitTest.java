package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import com.google.common.collect.ImmutableSortedMap;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AnnotatedIntervalUtilsUnitTest extends GATKBaseTest {
    @Test(dataProvider = "mergeTests")
    public void testMerge(List<AnnotatedInterval> test1, List<AnnotatedInterval> gt) {
        List<AnnotatedInterval> mergeTestResults = AnnotatedIntervalUtils.mergeRegions(test1,
                ReferenceUtils.loadFastaDictionary(new File(ReferenceUtils.getFastaDictionaryFileName(hg19MiniReference))),
                "__", l -> {});
        Assert.assertEquals(mergeTestResults, gt);
    }

    @DataProvider(name = "mergeTests")
    public Object [][] createMergeTests() {

        return new Object[][] {
            {
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1"))),
                Collections.singletonList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")))
            }, {
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new AnnotatedInterval( new SimpleInterval("1", 201, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1"))),
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new AnnotatedInterval( new SimpleInterval("1", 201, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")))
            }, {
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new AnnotatedInterval( new SimpleInterval("2", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foobar1", "bar1"))),
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new AnnotatedInterval( new SimpleInterval("2", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foobar1", "bar1")))
            }, {
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new AnnotatedInterval( new SimpleInterval("1", 190, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new AnnotatedInterval( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1"))),
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1__bar2")),
                        new AnnotatedInterval( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1")))
            }, {
                Arrays.asList(
                        // Same as previous test, but trying to see if input order hoses things.  Note that the output will
                        //  be sorted.
                        new AnnotatedInterval( new SimpleInterval("1", 190, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1")),
                        new AnnotatedInterval( new SimpleInterval("1", 100, 200),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar2")),
                        new AnnotatedInterval( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1"))),
                Arrays.asList(
                        new AnnotatedInterval( new SimpleInterval("1", 100, 300),
                                ImmutableSortedMap.of("Foo", "bar", "Foo1", "bar1__bar2")),
                        new AnnotatedInterval( new SimpleInterval("1", 301, 500),
                                ImmutableSortedMap.of("Foo2", "bar", "Foo3", "bar1")))
            }
        };
    }

    @DataProvider
    public Object[][] provideMergeByAnnotation() {
        return new Object[][] {
                {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 200001, 300000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 100000,
                        Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    )
                }, {
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 100000,
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 200001, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))

                    ), Arrays.asList("Foo", "Baz"), 100000,
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 200001, 300000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 1,
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 200002, 300000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 1,
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 200002, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 200001, 300000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcde", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 100000,
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd__abcde", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcde", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 100000,
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd__abcde", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                ImmutableSortedMap.of("Foo", "1234", "Bar", "abcde", "Baz", "wxyz"))
                    ), Arrays.asList("Foo", "Baz"), 10000,
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcde", "Baz", "wxyz"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcde", "Baz", "wxyzabc"))
                    ), Arrays.asList("Foo", "Baz"), 100000,
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcde", "Baz", "wxyzabc"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyzabc"))
                    ), Collections.singletonList("Foo"), 100000,
                    Collections.singletonList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz__wxyzabc"))
                    )
                }, {
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 200000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz")),
                            new AnnotatedInterval( new SimpleInterval("1", 250001, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyzabc")),
                            new AnnotatedInterval( new SimpleInterval("1", 350001, 300000),
                                    ImmutableSortedMap.of("Foo", "5678", "Bar", "abcd", "Baz", "wxyzabc"))
                    ), Collections.singletonList("Foo"), 100000,
                    Arrays.asList(
                            new AnnotatedInterval( new SimpleInterval("1", 100000, 300000),
                                    ImmutableSortedMap.of("Foo", "1234", "Bar", "abcd", "Baz", "wxyz__wxyzabc")),
                            new AnnotatedInterval( new SimpleInterval("1", 350001, 400000),
                                    ImmutableSortedMap.of("Foo", "5678", "Bar", "abcd", "Baz", "wxyzabc"))
                    )
                }
        };
    }

    @Test(dataProvider = "provideMergeByAnnotation")
    public void testMergeByAnnotation(List<AnnotatedInterval> initialRegions, List<String> annotationNames,
                                      int maxDistanceInBp, List<AnnotatedInterval> gtMergedRegions) {
        final List<AnnotatedInterval> mergedRegions = AnnotatedIntervalUtils.mergeRegionsByAnnotation(initialRegions,
                ReferenceUtils.loadFastaDictionary(new File(ReferenceUtils.getFastaDictionaryFileName(hg19MiniReference))),
                annotationNames, t -> {}, "__", maxDistanceInBp);
        Assert.assertEquals(mergedRegions, gtMergedRegions);
    }
}
