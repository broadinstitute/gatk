package org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Maps;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class SimpleGermlineTaggerUnitTest extends GATKBaseTest {
    private static final String REF = hg19MiniReference;
    public static final String TEST_GERMLINE_TAGGING_ANNOTATION = "germline_tagging";

    @Test(dataProvider = "simpleTests")
    public void testSimpleTagging(List<AnnotatedInterval> tumorSegments, List<AnnotatedInterval> normalSegments,
                                  List<AnnotatedInterval> gt) {

        final List<AnnotatedInterval> testResult = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(tumorSegments,
                normalSegments, "call",
                ReferenceUtils.loadFastaDictionary(new File(ReferenceUtils.getFastaDictionaryFileName(REF))), TEST_GERMLINE_TAGGING_ANNOTATION, 10);

        Assert.assertEquals(testResult, gt);
    }

    @DataProvider(name = "simpleTests")
    public Object[][] createSimpleTests() {
        return new Object[][]{
                {
                    // Trivial case
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // Almost-trivial case
                    Arrays.asList(
                            // Tumor segments are assumed to be mutable.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // More than one tumor segment maps to one germline seg.  This example would indicate that a tumor sample had deleted a germline amp.
                    Arrays.asList(
                            // Tumor segments are assumed to be mutable.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // Same as above, but using the padding.
                    Arrays.asList(
                            // Tumor segments are assumed to be mutable.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 99, 295), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 296, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // Same as above, but with deletions.
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 99, 295), ImmutableSortedMap.of("call", "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 296, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // No tumor segment matches at the endpoint
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 400), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 400), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // Three tumor segment matches at the endpoint
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // Three tumor segment matches at the start, but just barely not at the end
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 315), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 315), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    // Three tumor segment matches at the endpoint
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 310), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 310), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 90, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 90, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 89, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 89, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    Arrays.asList(
                            // Small segment is within padding, so it gets included.
                            new AnnotatedInterval(new SimpleInterval("1", 95, 100), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 95, 100), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    Arrays.asList(
                            // The first segment should not get tagged.
                            new AnnotatedInterval(new SimpleInterval("1", 89, 100), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 89, 100), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }, {
                    Arrays.asList(
                            // The third segment should not get tagged.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("call", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 311), Maps.newTreeMap(ImmutableSortedMap.of("call", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("call", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("call", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("call", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 311), ImmutableSortedMap.of("call", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("call", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    )
            }
        };
    }
}
