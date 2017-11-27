package org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Maps;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.spark_project.guava.collect.Lists;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class SimpleGermlineTaggerUnitTest extends GATKBaseTest {
    private static final String REF = hg19MiniReference;
    public static final String TEST_GERMLINE_TAGGING_ANNOTATION = "germline_tagging";

    @Test(dataProvider = "simpleTests")
    public void testSimpleTagging(List<SimpleAnnotatedGenomicRegion> tumorSegments, List<SimpleAnnotatedGenomicRegion> normalSegments,
                                  List<SimpleAnnotatedGenomicRegion> gt) {

        final List<SimpleAnnotatedGenomicRegion> testResult = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(tumorSegments,
                normalSegments, "call",
                ReferenceUtils.loadFastaDictionary(new File(ReferenceUtils.getFastaDictionaryFileName(REF))), TEST_GERMLINE_TAGGING_ANNOTATION, 10);

        Assert.assertEquals(testResult, gt);
    }

    @DataProvider(name = "simpleTests")
    public Object [][] createSimpleTests() {
        return new Object[][] {
                {
                        // Trivial case
                        Lists.newArrayList(
                                // Tumor segments are assumed to be mutable.
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0")))
                        ),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+")),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 500), ImmutableSortedMap.of("call", "0"))
                        ),
                        Lists.newArrayList(
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                        )
                },
                {
                        // Almost-trivial case
                        Lists.newArrayList(
                                // Tumor segments are assumed to be mutable.
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                        ),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+")),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 500), ImmutableSortedMap.of("call", "0"))
                        ),
                        Lists.newArrayList(
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                        )
                },
                {
                        // More than one tumor segment maps to one germline seg.  This example would indicate that a tumor sample had deleted a germline amp.
                        Lists.newArrayList(
                                // Tumor segments are assumed to be mutable.
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                        ),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                        ),
                        Lists.newArrayList(
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                        )
                },
               {
                        // Same as above, but using the padding.
                        Lists.newArrayList(
                                // Tumor segments are assumed to be mutable.
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                        ),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 99, 295), ImmutableSortedMap.of("call", "+")),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 296, 500), ImmutableSortedMap.of("call", "0"))
                        ),
                        Lists.newArrayList(
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                        )
                },
               {
                        // Same as above, but with deletions.
                        Lists.newArrayList(
                                // Tumor segments are assumed to be mutable.
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                        ),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 99, 295), ImmutableSortedMap.of("call", "-")),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 296, 500), ImmutableSortedMap.of("call", "0"))
                        ),
                        Lists.newArrayList(
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                        )
               }, {
                        // No tumor segment matches at the endpoint
                        Lists.newArrayList(
                                // Tumor segments are assumed to be mutable.
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 251, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                        ),
                        Lists.newArrayList(
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                        ),
                        Lists.newArrayList(
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 251, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                        )
                }
        };
    }
}
