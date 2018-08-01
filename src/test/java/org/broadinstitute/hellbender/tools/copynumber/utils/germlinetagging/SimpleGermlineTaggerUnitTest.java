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
                                  List<AnnotatedInterval> gt, int padding) {

        final List<AnnotatedInterval> testResult = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(tumorSegments,
                normalSegments, "CALL",
                ReferenceUtils.loadFastaDictionary(new File(ReferenceUtils.getFastaDictionaryFileName(REF))), TEST_GERMLINE_TAGGING_ANNOTATION, padding);

        Assert.assertEquals(testResult, gt);
    }

    @DataProvider(name = "simpleTests")
    public Object[][] createSimpleTests() {
        return new Object[][]{
                {
                    // Trivial case
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // Almost-trivial case
                    Arrays.asList(
                            // Tumor segments are assumed to be mutable.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // More than one tumor segment maps to one germline seg.  This example would indicate that a tumor sample had deleted a germline amp.
                    Arrays.asList(
                            // Tumor segments are assumed to be mutable.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // Same as above, but using the padding.
                    Arrays.asList(
                            // Tumor segments are assumed to be mutable.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 99, 295), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 296, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // Same as above, but with deletions.
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 99, 295), ImmutableSortedMap.of("CALL", "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 296, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 400), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // No tumor segment matches at the endpoint
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 400), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 400), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // Three tumor segment matches at the endpoint
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // Three tumor segment matches at the start, but just barely not at the end
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 315), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 315), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    // Three tumor segment matches at the endpoint
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 310), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 310), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 90, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 90, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 89, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 89, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    Arrays.asList(
                            // Small segment is within padding, so it gets included.
                            new AnnotatedInterval(new SimpleInterval("1", 95, 100), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            // The first segment should not be tagged as a germline event, since it mostly does not overlap a germline segment.
                            new AnnotatedInterval(new SimpleInterval("1", 95, 100), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    Arrays.asList(
                            // The first segment should not get tagged.
                            new AnnotatedInterval(new SimpleInterval("1", 89, 100), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 89, 100), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 101, 250), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 251, 305), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    Arrays.asList(
                            // The third segment should not get tagged.
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+"))),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 311), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "0"))),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), Maps.newTreeMap(ImmutableSortedMap.of("CALL", "+")))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 300), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 200), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 201, 300), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 301, 311), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 401, 500), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10
            }, {
                    Arrays.asList(
                            // The big segment should not get tagged
                            new AnnotatedInterval(new SimpleInterval("1", 100, 23000), ImmutableSortedMap.of("CALL", "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 23001, 25000), ImmutableSortedMap.of("CALL", "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 26001, 300500), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 17001, 31000), ImmutableSortedMap.of("CALL", "-"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 100, 23000), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0")),
                            new AnnotatedInterval(new SimpleInterval("1", 23001, 25000), ImmutableSortedMap.of("CALL", "-", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 26001, 300500), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10000
            }, {
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 10459001, 10472000), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 10472001, 11427000), ImmutableSortedMap.of("CALL", "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 11427001, 30709000), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 10459001, 10469000), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 10469001, 11489000), ImmutableSortedMap.of("CALL", "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 11489001, 77286000), ImmutableSortedMap.of("CALL", "0"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 10459001, 10472000), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 10472001, 11427000), ImmutableSortedMap.of("CALL", "-", TEST_GERMLINE_TAGGING_ANNOTATION, "-")),
                            new AnnotatedInterval(new SimpleInterval("1", 11427001, 30709000), ImmutableSortedMap.of("CALL", "0", TEST_GERMLINE_TAGGING_ANNOTATION, "0"))
                    ), 10000
            }, {
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 2602001, 2792000), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 2792001, 3101000), ImmutableSortedMap.of("CALL", "+"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 2539001, 2792000), ImmutableSortedMap.of("CALL", "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 2792001, 3101000), ImmutableSortedMap.of("CALL", "+"))
                    ),
                    Arrays.asList(
                            new AnnotatedInterval(new SimpleInterval("1", 2602001, 2792000), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+")),
                            new AnnotatedInterval(new SimpleInterval("1", 2792001, 3101000), ImmutableSortedMap.of("CALL", "+", TEST_GERMLINE_TAGGING_ANNOTATION, "+"))
                    ),10000
            }
        };
    }
}
