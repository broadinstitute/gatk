package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

public class SimpleAnnotatedGenomicRegionCollectionUnitTest  extends GATKBaseTest {

    private static final File TEST_FILE = new File(toolsTestDir,
            "copynumber/utils/combine-segment-breakpoints-with-legacy-header-learning-combined-copy-number.tsv");


    @Test
    public void basicTest() throws IOException {
        final Set<String> headersOfInterest = Sets.newHashSet(Arrays.asList("name", "learning_SAMPLE_0"));
        final SimpleAnnotatedGenomicRegionCollection simpleAnnotatedGenomicRegions =
                SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(TEST_FILE, headersOfInterest);

        Assert.assertEquals(simpleAnnotatedGenomicRegions.size(), 15);
        Assert.assertTrue(simpleAnnotatedGenomicRegions.getRecords().stream()
                .mapToInt(s -> s.getAnnotations().entrySet().size())
                .allMatch(i -> i == headersOfInterest.size()));
        Assert.assertTrue(simpleAnnotatedGenomicRegions.getRecords().stream().allMatch(s -> s.getAnnotations().keySet().containsAll(headersOfInterest)));

        // Grab the first 15 and test values
        List<SimpleAnnotatedGenomicRegion> gtRegions = Arrays.asList(
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 30365, 30503), ImmutableSortedMap.of("name", "target_1_None", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 69088, 70010), ImmutableSortedMap.of("name", "target_2_None", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 367656, 368599), ImmutableSortedMap.of("name", "target_3_None", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 621093, 622036), ImmutableSortedMap.of("name", "target_4_None", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 861319, 861395), ImmutableSortedMap.of("name", "target_5_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 865532, 865718), ImmutableSortedMap.of("name", "target_6_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 866416, 866471), ImmutableSortedMap.of("name", "target_7_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 871149, 871278), ImmutableSortedMap.of("name", "target_8_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 874417, 874511), ImmutableSortedMap.of("name", "target_9_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 874652, 874842), ImmutableSortedMap.of("name", "target_10_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 876521, 876688), ImmutableSortedMap.of("name", "target_11_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 877513, 877633), ImmutableSortedMap.of("name", "target_12_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 877787, 877870), ImmutableSortedMap.of("name", "target_13_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 877936, 878440), ImmutableSortedMap.of("name", "target_14_SAMD11", "learning_SAMPLE_0", "2")),
                new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 878630, 878759), ImmutableSortedMap.of("name", "target_15_SAMD11", "learning_SAMPLE_0", "2"))
        );

        Assert.assertEquals(simpleAnnotatedGenomicRegions.getRecords().subList(0, gtRegions.size()), gtRegions);
    }

    @Test
    public void testCreationFromList() {
        final Set<String> headersOfInterest = Sets.newHashSet(Arrays.asList("name", "learning_SAMPLE_0"));
        final List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions =
                SimpleAnnotatedGenomicRegion.readAnnotatedRegions(TEST_FILE, headersOfInterest);
        final SimpleAnnotatedGenomicRegionCollection collection = SimpleAnnotatedGenomicRegionCollection.create(simpleAnnotatedGenomicRegions,
                ReferenceUtils.loadFastaDictionary(new File(hg19_chr1_1M_dict)), Lists.newArrayList("name", "learning_SAMPLE_0"));

        Assert.assertEquals(collection.getRecords(), simpleAnnotatedGenomicRegions);
    }

    @Test
    public void basicTestWithAllColumnsFile() throws IOException {

        // If no columns of interest are given in a read call, the method will try to load all columns as "interesting".
        final SimpleAnnotatedGenomicRegionCollection simpleAnnotatedGenomicRegions =
                SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(TEST_FILE);

        Assert.assertEquals(simpleAnnotatedGenomicRegions.size(), 15);
        Assert.assertTrue(simpleAnnotatedGenomicRegions.getRecords().stream()
                .mapToInt(s -> s.getAnnotations().entrySet().size())
                .allMatch(i -> i == 101)); // The number of columns in the TEST_FILE (name, learning_SAMPLE_0...99
    }
}
