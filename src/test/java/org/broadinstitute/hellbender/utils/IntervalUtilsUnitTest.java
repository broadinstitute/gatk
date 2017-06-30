package org.broadinstitute.hellbender.utils;

import com.google.common.collect.Lists;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.SimpleFeature;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;


/**
 * test out the interval utility methods
 */
public final class IntervalUtilsUnitTest extends BaseTest {
    public static final String INTERVAL_TEST_DATA = TestResources.publicTestDir + "org/broadinstitute/hellbender/utils/interval/";
    public static final String emptyIntervals = TestResources.publicTestDir + "empty_intervals.list";
    private List<GenomeLoc> hg19ReferenceLocs;
    private List<GenomeLoc> hg19exomeIntervals;


    @BeforeClass
    public void init() {
        hg19ReferenceLocs = Collections.unmodifiableList(GenomeLocSortedSet.createSetFromSequenceDictionary(TestResources.getHg19Header().getSequenceDictionary()).toList()) ;
        hg19exomeIntervals = Collections.unmodifiableList(IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(TestResources.hg19MiniIntervalFile)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullArgSpanningInterval() throws Exception {
        IntervalUtils.getSpanningInterval(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSpanningInterval_differentConfigs() throws Exception {
        final SimpleInterval chr1_200_300 = new SimpleInterval("1:200-300");
        final SimpleInterval chr2_1_100 = new SimpleInterval("2:1-100");

        IntervalUtils.getSpanningInterval(Arrays.asList(chr1_200_300, chr2_1_100));
    }

    @Test
    public void testCompareLocatables() throws Exception {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("1", 1000));
        dict.addSequence(new SAMSequenceRecord("2", 1000));
        final SimpleInterval chr1_1_100 = new SimpleInterval("1:1-100");
        final SimpleInterval chr1_5_100 = new SimpleInterval("1:5-100");
        final SimpleInterval chr2_1_100 = new SimpleInterval("2:1-100");

        // equal intervals comparison return 0
        Assert.assertEquals(IntervalUtils.compareLocatables(chr1_1_100, chr1_1_100, dict), 0);
        Assert.assertEquals(IntervalUtils.compareLocatables(chr1_5_100, chr1_5_100, dict), 0);
        Assert.assertEquals(IntervalUtils.compareLocatables(chr2_1_100, chr2_1_100, dict), 0);
        // first < second return negative
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_1_100, chr1_5_100, dict) < 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_1_100, chr2_1_100, dict) < 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_5_100, chr2_1_100, dict) < 0);
        // first > second return positive
        Assert.assertTrue(IntervalUtils.compareLocatables(chr2_1_100, chr1_1_100, dict) > 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr2_1_100, chr1_5_100, dict) > 0);
        Assert.assertTrue(IntervalUtils.compareLocatables(chr1_5_100, chr1_1_100, dict) > 0);
    }

    @Test
    public void testSpanningInterval_nullIfEmptyInput() throws Exception {
        Assert.assertNull(IntervalUtils.getSpanningInterval(Collections.emptyList()));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSpanningInterval_nullElement() throws Exception {
        final SimpleInterval chr1_200_300 = new SimpleInterval("1:200-300");
        Assert.assertNull(IntervalUtils.getSpanningInterval(Arrays.asList(chr1_200_300, null)));
    }

    @DataProvider(name = "SpanningInterval")
    public Object[][] SpanningInterval() {
        final SimpleInterval chr1_1_100   = new SimpleInterval("1:1-100");
        final SimpleInterval chr1_1_200 = new SimpleInterval("1:1-200");
        final SimpleInterval chr1_1_300 = new SimpleInterval("1:1-300");
        final SimpleInterval chr1_100_200 = new SimpleInterval("1:100-200");
        final SimpleInterval chr1_200_300 = new SimpleInterval("1:200-300");
        return new Object[][]{
           {Arrays.asList(chr1_1_100),                             chr1_1_100, },
           {Arrays.asList(chr1_1_100, chr1_100_200),               chr1_1_200, },
           {Arrays.asList(chr1_1_100, chr1_100_200, chr1_200_300), chr1_1_300, },
           {Arrays.asList(chr1_1_100, chr1_200_300),               chr1_1_300, },
        };
    }

    @Test(dataProvider = "SpanningInterval")
    public void testSpanningInterval(final List<? extends Locatable> locs, final SimpleInterval expectedResult) throws Exception {
        Assert.assertEquals(IntervalUtils.getSpanningInterval(locs), expectedResult);
    }

    // -------------------------------------------------------------------------------------
    //
    // tests to ensure the quality of the interval cuts of the interval cutting functions
    //
    // -------------------------------------------------------------------------------------

    private class IntervalSlicingTest extends TestDataProvider {
        public int parts;
        public double maxAllowableVariance;

        private IntervalSlicingTest(final int parts, final double maxAllowableVariance) {
            super(IntervalSlicingTest.class);
            this.parts = parts;
            this.maxAllowableVariance = maxAllowableVariance;
        }

        public String toString() {
            return String.format("IntervalSlicingTest parts=%d maxVar=%.2f", parts, maxAllowableVariance);
        }
    }

    @DataProvider(name = "intervalslicingdata")
    public Object[][] createTrees() {
        new IntervalSlicingTest(1, 0);
        new IntervalSlicingTest(2, 1);
        new IntervalSlicingTest(5, 1);
        new IntervalSlicingTest(10, 1);
        new IntervalSlicingTest(67, 1);
        new IntervalSlicingTest(100, 1);

        return IntervalSlicingTest.getTests(IntervalSlicingTest.class);
    }

    @Test(dataProvider = "intervalslicingdata")
    public void testFixedScatterIntervalsAlgorithm(IntervalSlicingTest test) {
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(hg19exomeIntervals, test.parts);

        long totalSize = IntervalUtils.intervalSize(hg19exomeIntervals);
        long idealSplitSize = totalSize / test.parts;

        long sumOfSplitSizes = 0;
        int counter = 0;
        for ( final List<GenomeLoc> split : splits ) {
            long splitSize = IntervalUtils.intervalSize(split);
            double sigma = (splitSize - idealSplitSize) / (1.0 * idealSplitSize);
            //logger.warn(String.format("Split %d size %d ideal %d sigma %.2f", counter, splitSize, idealSplitSize, sigma));
            counter++;
            sumOfSplitSizes += splitSize;
            Assert.assertTrue(Math.abs(sigma) <= test.maxAllowableVariance, String.format("Interval %d (size %d ideal %d) has a variance %.2f outside of the tolerated range %.2f", counter, splitSize, idealSplitSize, sigma, test.maxAllowableVariance));
        }

        Assert.assertEquals(totalSize, sumOfSplitSizes, "Split intervals don't contain the exact number of bases in the origianl intervals");
    }

    // -------------------------------------------------------------------------------------
    //
    // splitLocusIntervals tests
    //
    // -------------------------------------------------------------------------------------

    /** large scale tests for many intervals */
    private class SplitLocusIntervalsTest extends TestDataProvider {
        final List<GenomeLoc> originalIntervals;
        final public int parts;

        private SplitLocusIntervalsTest(final String name, List<GenomeLoc> originalIntervals, final int parts) {
            super(SplitLocusIntervalsTest.class, name);
            this.parts = parts;
            this.originalIntervals = originalIntervals;
        }

        public String toString() {
            return String.format("%s parts=%d", super.toString(), parts);
        }
    }

    @DataProvider(name = "IntervalRepartitionTest")
    public Object[][] createIntervalRepartitionTest() {
        for ( int parts : Arrays.asList(1, 2, 3, 10, 13, 100, 151, 1000, 10000) ) {
            //for ( int parts : Arrays.asList(10) ) {
            new SplitLocusIntervalsTest("hg19RefLocs", hg19ReferenceLocs, parts);
            new SplitLocusIntervalsTest("hg19ExomeLocs", hg19exomeIntervals, parts);
        }

        return SplitLocusIntervalsTest.getTests(SplitLocusIntervalsTest.class);
    }

    @Test(dataProvider = "IntervalRepartitionTest")
    public void testIntervalRepartition(SplitLocusIntervalsTest test) {
        List<List<GenomeLoc>> splitByLocus = IntervalUtils.splitLocusIntervals(test.originalIntervals, test.parts);
        Assert.assertEquals(splitByLocus.size(), test.parts, "SplitLocusIntervals failed to generate correct number of intervals");
        List<GenomeLoc> flat = IntervalUtils.flattenSplitIntervals(splitByLocus);

        // test overall size
        final long originalSize = IntervalUtils.intervalSize(test.originalIntervals);
        final long flatSize = IntervalUtils.intervalSize(flat);
        Assert.assertEquals(flatSize, originalSize, "SplitLocusIntervals locs cover an incorrect number of bases");

        // test size of each split
        final long ideal = (long)Math.floor(originalSize / (1.0 * test.parts));
        final long maxSize = ideal + (originalSize % test.parts) * test.parts; // no more than N * rounding error in size
        for ( final List<GenomeLoc> split : splitByLocus ) {
            final long splitSize = IntervalUtils.intervalSize(split);
            Assert.assertTrue(splitSize >= ideal && splitSize <= maxSize,
                    String.format("SplitLocusIntervals interval (start=%s) has size %d outside of bounds ideal=%d, max=%d",
                            split.get(0), splitSize, ideal, maxSize));
        }

        // test that every base in original is covered once by a base in split by locus intervals
        String diff = IntervalUtils.equateIntervals(test.originalIntervals, flat);
        Assert.assertNull(diff, diff);
    }

    /** small scale tests where the expected cuts are enumerated upfront for testing */
    private class SplitLocusIntervalsSmallTest extends TestDataProvider {
        final List<GenomeLoc> original;
        final public int parts;
        final public int expectedParts;
        final List<GenomeLoc> expected;

        private SplitLocusIntervalsSmallTest(final String name, List<GenomeLoc> originalIntervals, final int parts, List<GenomeLoc> expected) {
            this(name, originalIntervals, parts,  expected, parts);
        }

        private SplitLocusIntervalsSmallTest(final String name, List<GenomeLoc> originalIntervals, final int parts, List<GenomeLoc> expected, int expectedParts) {
            super(SplitLocusIntervalsSmallTest.class, name);
            this.parts = parts;
            this.expectedParts = expectedParts;
            this.original = originalIntervals;
            this.expected = expected;
        }

        public String toString() {
            return String.format("%s parts=%d", super.toString(), parts);
        }
    }

    @DataProvider(name = "SplitLocusIntervalsSmallTest")
    public Object[][] createSplitLocusIntervalsSmallTest() {
        GenomeLoc bp01_10 = TestResources.getGenomeLocParser().createGenomeLoc("1", 1, 10);

        GenomeLoc bp1_5 = TestResources.getGenomeLocParser().createGenomeLoc("1", 1, 5);
        GenomeLoc bp6_10 = TestResources.getGenomeLocParser().createGenomeLoc("1", 6, 10);
        new SplitLocusIntervalsSmallTest("cut into two", Arrays.asList(bp01_10), 2, Arrays.asList(bp1_5, bp6_10));

        GenomeLoc bp20_30 = TestResources.getGenomeLocParser().createGenomeLoc("1", 20, 30);
        new SplitLocusIntervalsSmallTest("two in two", Arrays.asList(bp01_10, bp20_30), 2, Arrays.asList(bp01_10, bp20_30));

        GenomeLoc bp1_7 = TestResources.getGenomeLocParser().createGenomeLoc("1", 1, 7);
        GenomeLoc bp8_10 = TestResources.getGenomeLocParser().createGenomeLoc("1", 8, 10);
        GenomeLoc bp20_23 = TestResources.getGenomeLocParser().createGenomeLoc("1", 20, 23);
        GenomeLoc bp24_30 = TestResources.getGenomeLocParser().createGenomeLoc("1", 24, 30);
        new SplitLocusIntervalsSmallTest("two in three", Arrays.asList(bp01_10, bp20_30), 3,
                Arrays.asList(bp1_7, bp8_10, bp20_23, bp24_30));

        GenomeLoc bp1_2 = TestResources.getGenomeLocParser().createGenomeLoc("1", 1, 2);
        GenomeLoc bp1_1 = TestResources.getGenomeLocParser().createGenomeLoc("1", 1, 1);
        GenomeLoc bp2_2 = TestResources.getGenomeLocParser().createGenomeLoc("1", 2, 2);
        new SplitLocusIntervalsSmallTest("too many pieces", Arrays.asList(bp1_2), 5, Arrays.asList(bp1_1, bp2_2), 2);

        new SplitLocusIntervalsSmallTest("emptyList", Collections.<GenomeLoc>emptyList(), 5, Collections.<GenomeLoc>emptyList(), 0);

        return SplitLocusIntervalsSmallTest.getTests(SplitLocusIntervalsSmallTest.class);
    }

    @Test(dataProvider = "SplitLocusIntervalsSmallTest")
    public void splitLocusIntervalsSmallTest(SplitLocusIntervalsSmallTest test) {
        List<List<GenomeLoc>> splitByLocus = IntervalUtils.splitLocusIntervals(test.original, test.parts);
        Assert.assertEquals(splitByLocus.size(), test.expectedParts, "SplitLocusIntervals failed to generate correct number of intervals");
        List<GenomeLoc> flat = IntervalUtils.flattenSplitIntervals(splitByLocus);

        // test sizes
        final long originalSize = IntervalUtils.intervalSize(test.original);
        final long splitSize = IntervalUtils.intervalSize(flat);
        Assert.assertEquals(splitSize, originalSize, "SplitLocusIntervals locs cover an incorrect number of bases");

        Assert.assertEquals(flat, test.expected, "SplitLocusIntervals locs not expected intervals");
    }

    //
    // Misc. tests
    //

    @Test(expectedExceptions=UserException.class)
    public void testMergeListsBySetOperatorNoOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> listEveryTwoFromOne = new ArrayList<>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(TestResources.getGenomeLocParser().createGenomeLoc("chr1",x,x));
            else
                listEveryTwoFromOne.add(TestResources.getGenomeLocParser().createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 100);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, null);
        Assert.assertEquals(ret.size(), 100);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 0);
    }

    @Test
    public void testMergeListsBySetOperatorAllOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(TestResources.getGenomeLocParser().createGenomeLoc("1",x,x));
            allSites.add(TestResources.getGenomeLocParser().createGenomeLoc("1",x,x));
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 150);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, null);
        Assert.assertEquals(ret.size(), 150);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 50);
    }

    @Test
    public void testMergeListsBySetOperator() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 5 == 0) {
                listEveryTwoFromTwo.add(TestResources.getGenomeLocParser().createGenomeLoc("1",x,x));
                allSites.add(TestResources.getGenomeLocParser().createGenomeLoc("1",x,x));
            }
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, null);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 20);
    }

    @Test
    public void testOverlappingIntervalsFromSameSourceWithIntersection() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> source1 = new ArrayList<>();
        List<GenomeLoc> source2 = new ArrayList<>();

        source1.add(TestResources.getGenomeLocParser().createGenomeLoc("1", 10, 20));
        source1.add(TestResources.getGenomeLocParser().createGenomeLoc("1", 15, 25));

        source2.add(TestResources.getGenomeLocParser().createGenomeLoc("1", 16, 18));
        source2.add(TestResources.getGenomeLocParser().createGenomeLoc("1", 22, 24));

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(source1, source2, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 2);
    }

    @Test
    public void testGetContigLengths() {
        Map<String, Integer> lengths = IntervalUtils.getContigSizes(new File(TestResources.exampleReference));
        Assert.assertEquals((long)lengths.get("1"), 16000);
        Assert.assertEquals((long)lengths.get("2"), 16000);
        Assert.assertEquals((long) lengths.get("3"), 16000);
        Assert.assertEquals((long)lengths.get("4"), 16000);
    }

    @Test
    public void testParseIntervalArguments() {
        Assert.assertEquals(hg19ReferenceLocs.size(), 4);
        Assert.assertEquals(TestResources.intervalStringsToGenomeLocs("1", "2", "3").size(), 3);
        Assert.assertEquals(TestResources.intervalStringsToGenomeLocs("1:1-2", "1:4-5", "2:1-1", "3:2-2").size(), 4);
    }

    @Test
    public void testParseUnmappedIntervalArgument() {
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord("1", 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final GenomeLocParser parser = new GenomeLocParser(dictionary);

        List<GenomeLoc> unmappedLoc = IntervalUtils.parseIntervalArguments(parser, "unmapped");
        Assert.assertTrue(unmappedLoc.size() == 1);
        Assert.assertTrue(unmappedLoc.get(0).isUnmapped());
    }

    @Test
    public void testParseUnmappedIntervalArgumentInIntervalFile() {
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord("1", 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final GenomeLocParser parser = new GenomeLocParser(dictionary);

        List<GenomeLoc> unmappedLoc = IntervalUtils.parseIntervalArguments(parser, TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped.intervals");
        Assert.assertTrue(unmappedLoc.size() == 1);
        Assert.assertTrue(unmappedLoc.get(0).isUnmapped());
    }

    @Test
    public void testIsIntervalFile() {
        Assert.assertTrue(IntervalUtils.isIntervalFile(emptyIntervals));
        Assert.assertTrue(IntervalUtils.isIntervalFile(emptyIntervals, true));

        Assert.assertFalse(IntervalUtils.isIntervalFile(INTERVAL_TEST_DATA + "intervals_from_features_test.vcf"));
        Assert.assertFalse(IntervalUtils.isIntervalFile(INTERVAL_TEST_DATA + "intervals_from_features_test.bed"));

        List<String> extensions = Arrays.asList("interval_list", "intervals", "list", "picard");
        for (String extension: extensions) {
            Assert.assertTrue(IntervalUtils.isIntervalFile("test_intervals." + extension, false), "Tested interval file extension: " + extension);
        }
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testMissingIntervalFile() {
        IntervalUtils.isIntervalFile(TestResources.publicTestDir + "no_such_intervals.list");
    }

    @Test
    public void testParseIntervalWithPeriodInContigName() {
        // Make sure that we don't interpret contigs with periods in their name as files
        final String contigName = "GL000249.1";
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord(contigName, 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final GenomeLocParser parser = new GenomeLocParser(dictionary);

        final List<GenomeLoc> result = IntervalUtils.parseIntervalArguments(parser, contigName);
        Assert.assertEquals(result.size(), 1);
        Assert.assertEquals(result.get(0).getContig(), contigName);
        Assert.assertEquals(result.get(0).getStart(), 1);
        Assert.assertEquals(result.get(0).getEnd(), 100);
    }

    @Test
    public void testFixedScatterIntervalsBasic() {
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3");

        List<File> files = testFiles("basic.", 3, ".intervals");

        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs("1", "2", "3");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterFixedIntervalsLessFiles() {
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3");
        GenomeLoc chr4 = TestResources.getGenomeLocParser().parseGenomeLoc("4");

        List<File> files = testFiles("less.", 3, ".intervals");

        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs("1", "2", "3", "4");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
        Assert.assertEquals(locs3.get(1), chr4);
    }

    @Test(expectedExceptions=CommandLineException.BadArgumentValue.class)
    public void testSplitFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs("1", "2");
        IntervalUtils.splitFixedIntervals(locs, files.size());
    }

    @Test(expectedExceptions=CommandLineException.BadArgumentValue.class)
    public void testScatterFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs("1", "2");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, locs.size()); // locs.size() instead of files.size()
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);
    }
    @Test
    public void testScatterFixedIntervalsStart() {
        List<String> intervals = Arrays.asList("1:1-2", "1:4-5", "2:1-1", "3:2-2");
        GenomeLoc chr1a = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-2");
        GenomeLoc chr1b = TestResources.getGenomeLocParser().parseGenomeLoc("1:4-5");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs2.get(0), chr1b);
        Assert.assertEquals(locs3.get(0), chr2);
        Assert.assertEquals(locs3.get(1), chr3);
    }

    @Test
    public void testScatterFixedIntervalsMiddle() {
        List<String> intervals = Arrays.asList("1:1-1", "2:1-2", "2:4-5", "3:2-2");
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-1");
        GenomeLoc chr2a = TestResources.getGenomeLocParser().parseGenomeLoc("2:1-2");
        GenomeLoc chr2b = TestResources.getGenomeLocParser().parseGenomeLoc("2:4-5");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs3.get(0), chr2b);
        Assert.assertEquals(locs3.get(1), chr3);
    }

    @Test
    public void testScatterFixedIntervalsEnd() {
        List<String> intervals = Arrays.asList("1:1-1", "2:2-2", "3:1-2", "3:4-5");
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2:2-2");
        GenomeLoc chr3a = TestResources.getGenomeLocParser().parseGenomeLoc("3:1-2");
        GenomeLoc chr3b = TestResources.getGenomeLocParser().parseGenomeLoc("3:4-5");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = TestResources.intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3a);
        Assert.assertEquals(locs3.get(0), chr3b);
    }

    @Test
    public void testScatterFixedIntervalsFile() {
        List<File> files = testFiles("sg.", 10, ".intervals");
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(TestResources.hg19MiniIntervalFile));
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());

        int[] counts = {
                33, 29, 68, 58, 64,
                31, 27, 62, 59, 52
        };

        //String splitCounts = "";
        for (int i = 0; i < splits.size(); i++) {
            int splitCount = splits.get(i).size();
            Assert.assertEquals(splitCount, counts[i], "Num intervals in split " + i);
        }
        //System.out.println(splitCounts.substring(2));

        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        int locIndex = 0;
        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), counts[i], "Intervals in " + file);
            for (GenomeLoc parsedLoc: parsedLocs)
                Assert.assertEquals(parsedLoc, locs.get(locIndex), String.format("Genome loc %d from file %d", locIndex++, i));
        }
        Assert.assertEquals(locIndex, locs.size(), "Total number of GenomeLocs");
    }

    @Test
    public void testScatterFixedIntervalsMax() {
        List<File> files = testFiles("sg.", 4, ".intervals");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(hg19ReferenceLocs, files.size());
        IntervalUtils.scatterFixedIntervals(TestResources.getHg19Header(), splits, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    @Test
    public void testScatterContigIntervalsOrder() {
        List<String> intervals = Arrays.asList("2:1-1", "1:1-1", "3:2-2");
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr2);
        Assert.assertEquals(locs2.get(0), chr1);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsBasic() {
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3");

        List<File> files = testFiles("contig_basic.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs("1", "2", "3"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsLessFiles() {
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3");
        GenomeLoc chr4 = TestResources.getGenomeLocParser().parseGenomeLoc("4");

        List<File> files = testFiles("contig_less.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs("1", "2", "3", "4"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3);
        Assert.assertEquals(locs3.get(0), chr4);
    }

    @Test(expectedExceptions=UserException.BadInput.class)
    public void testScatterContigIntervalsMoreFiles() {
        List<File> files = testFiles("contig_more.", 3, ".intervals");
        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs("1", "2"), files);
    }

    @Test
    public void testScatterContigIntervalsStart() {
        List<String> intervals = Arrays.asList("1:1-2", "1:4-5", "2:1-1", "3:2-2");
        GenomeLoc chr1a = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-2");
        GenomeLoc chr1b = TestResources.getGenomeLocParser().parseGenomeLoc("1:4-5");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3:2-2");

        List<File> files = testFiles("contig_split_start.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs1.get(1), chr1b);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsMiddle() {
        List<String> intervals = Arrays.asList("1:1-1", "2:1-2", "2:4-5", "3:2-2");
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-1");
        GenomeLoc chr2a = TestResources.getGenomeLocParser().parseGenomeLoc("2:1-2");
        GenomeLoc chr2b = TestResources.getGenomeLocParser().parseGenomeLoc("2:4-5");
        GenomeLoc chr3 = TestResources.getGenomeLocParser().parseGenomeLoc("3:2-2");

        List<File> files = testFiles("contig_split_middle.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 2);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs2.get(1), chr2b);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsEnd() {
        List<String> intervals = Arrays.asList("1:1-1", "2:2-2", "3:1-2", "3:4-5");
        GenomeLoc chr1 = TestResources.getGenomeLocParser().parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = TestResources.getGenomeLocParser().parseGenomeLoc("2:2-2");
        GenomeLoc chr3a = TestResources.getGenomeLocParser().parseGenomeLoc("3:1-2");
        GenomeLoc chr3b = TestResources.getGenomeLocParser().parseGenomeLoc("3:4-5");

        List<File> files = testFiles("contig_split_end.", 3 ,".intervals");

        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), TestResources.intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3a);
        Assert.assertEquals(locs3.get(1), chr3b);
    }

    @Test
    public void testScatterContigIntervalsMax() {
        List<File> files = testFiles("sg.", 4, ".intervals");
        IntervalUtils.scatterContigIntervals(TestResources.getHg19Header(), hg19ReferenceLocs, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    private List<File> testFiles(String prefix, int count, String suffix) {
        ArrayList<File> files = new ArrayList<>();
        for (int i = 1; i <= count; i++) {
            files.add(createTempFile(prefix + i, suffix));
        }
        return files;
    }

    @DataProvider(name="unmergedIntervals")
    public Object[][] getUnmergedIntervals() {
        return new Object[][] {
                new Object[] {"small_unmerged_picard_intervals.list"},
                new Object[] {"small_unmerged_gatk_intervals.list"}
        };
    }

    @Test(dataProvider="unmergedIntervals")
    public void testUnmergedIntervals(String unmergedIntervals) {
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(TestResources.getGenomeLocParser(), Collections.singletonList(TestResources.publicTestDir + unmergedIntervals));
        Assert.assertEquals(locs.size(), 2);

        List<GenomeLoc> merged;

        merged = IntervalUtils.mergeIntervalLocations(locs, IntervalMergingRule.ALL);
        Assert.assertEquals(merged.size(), 1);

        // Test that null means the same as ALL
        merged = IntervalUtils.mergeIntervalLocations(locs, null);
        Assert.assertEquals(merged.size(), 1);
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class, dataProvider="invalidIntervalTestData")
    public void testInvalidGATKFileIntervalHandling(GenomeLocParser genomeLocParser,
                                                    String contig, int intervalStart, int intervalEnd ) throws Exception {

        File gatkIntervalFile = createTempFile("testInvalidGATKFileIntervalHandling", ".intervals",
                String.format("%s:%d-%d", contig, intervalStart, intervalEnd));

        List<String> intervalArgs = new ArrayList<>(1);
        intervalArgs.add(gatkIntervalFile.getAbsolutePath());

        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    @Test(expectedExceptions=UserException.MalformedFile.class, dataProvider="invalidIntervalTestData")
    public void testLoadIntervalsInvalidPicardIntervalHandling(GenomeLocParser genomeLocParser,
                                                  String contig, int intervalStart, int intervalEnd ) throws Exception {

        SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("1"));
        // Intentionally add an extra contig not in our genomeLocParser's sequence dictionary
        // so that we can add intervals to the Picard interval file for contigs not in our dictionary
        picardFileHeader.addSequence(new SAMSequenceRecord("2", 100000));

        IntervalList picardIntervals = new IntervalList(picardFileHeader);
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname"));

        File picardIntervalFile = createTempFile("testLoadIntervalsInvalidPicardIntervalHandling", ".intervals");
        picardIntervals.write(picardIntervalFile);

        List<String> intervalArgs = new ArrayList<>(1);
        intervalArgs.add(picardIntervalFile.getAbsolutePath());

        // loadIntervals() will validate all intervals against the sequence dictionary in our genomeLocParser,
        // and should throw for all intervals in our invalidIntervalTestData set
        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    // TODO - remove once a corrected version of the exome interval list is released.
    @Test(dataProvider="negativeOneLengthIntervalTestData")
    public void testIntervalFileToListNegativeOneLength(GenomeLocParser genomeLocParser,
                                                  String contig, int intervalStart, int intervalEnd ) throws Exception {

        final SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("1"));

        final IntervalList picardIntervals = new IntervalList(picardFileHeader);
        // Need one good interval or else a UserException.MalformedFile( is thrown if no intervals
        picardIntervals.add(new Interval(contig, 1, 2, true, "dummyname0"));
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname1"));

        final File picardIntervalFile = createTempFile("testIntervalFileToListNegativeOneLength", ".intervals");
        picardIntervals.write(picardIntervalFile);

        IntervalUtils.intervalFileToList(genomeLocParser, picardIntervalFile.getAbsolutePath());
    }

    @Test(expectedExceptions=UserException.CouldNotReadInputFile.class, dataProvider="invalidIntervalTestData")
    public void testIntervalFileToListInvalidPicardIntervalHandling(GenomeLocParser genomeLocParser,
                                       String contig, int intervalStart, int intervalEnd ) throws Exception {

        final SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("1"));
        picardFileHeader.addSequence(new SAMSequenceRecord("2", 100000));

        final IntervalList picardIntervals = new IntervalList(picardFileHeader);
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname"));

        final File picardIntervalFile = createTempFile("testIntervalFileToListInvalidPicardIntervalHandling", ".intervals");
        picardIntervals.write(picardIntervalFile);

        // loadIntervals() will validate all intervals against the sequence dictionary in our genomeLocParser,
        // and should throw for all intervals in our invalidIntervalTestData set
        IntervalUtils.intervalFileToList(genomeLocParser, picardIntervalFile.getAbsolutePath());
    }

    @DataProvider(name="invalidIntervalTestData")
    public Object[][] invalidIntervalDataProvider() throws Exception {
        File fastaFile = new File(TestResources.publicTestDir + "exampleFASTA.fasta");
        GenomeLocParser genomeLocParser = new GenomeLocParser(new IndexedFastaSequenceFile(fastaFile));

        return new Object[][] {
                new Object[] {genomeLocParser, "1", 10000000, 20000000},
                new Object[] {genomeLocParser, "2", 1, 2},
                new Object[] {genomeLocParser, "1", -1, 50}
        };
    }

    // TODO - remove once a corrected version of the exome interval list is released.
    @DataProvider(name="negativeOneLengthIntervalTestData")
    public Object[][] negativeOneLengthIntervalDataProvider() throws Exception {
        File fastaFile = new File(TestResources.publicTestDir + "exampleFASTA.fasta");
        GenomeLocParser genomeLocParser = new GenomeLocParser(new IndexedFastaSequenceFile(fastaFile));

        return new Object[][] {
                new Object[] {genomeLocParser, "1", 2, 1},
                new Object[] {genomeLocParser, "1", 500, 499},
                new Object[] {genomeLocParser, "1", 1000, 999}
        };
    }

    private File createTempFile( String tempFilePrefix, String tempFileExtension, String... lines ) throws Exception {
        File tempFile = BaseTest.createTempFile(tempFilePrefix, tempFileExtension);
        FileUtils.writeLines(tempFile, Arrays.asList(lines));
        return tempFile;
    }

    @DataProvider(name = "sortAndMergeIntervals")
    public Object[][] getSortAndMergeIntervals() {
        return new Object[][] {
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, TestResources.intervalStringsToGenomeLocs("1:1", "1:3", "1:2"), TestResources.intervalStringsToGenomeLocs("1:1", "1:2", "1:3") },
                new Object[] { IntervalMergingRule.ALL, TestResources.intervalStringsToGenomeLocs("1:1", "1:3", "1:2"), TestResources.intervalStringsToGenomeLocs("1:1-3") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, TestResources.intervalStringsToGenomeLocs("1:1", "1:3", "2:2"), TestResources.intervalStringsToGenomeLocs("1:1", "1:3", "2:2") },
                new Object[] { IntervalMergingRule.ALL, TestResources.intervalStringsToGenomeLocs("1:1", "1:3", "2:2"), TestResources.intervalStringsToGenomeLocs("1:1", "1:3", "2:2") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, TestResources.intervalStringsToGenomeLocs("1:1", "1"), TestResources.intervalStringsToGenomeLocs("1") },
                new Object[] { IntervalMergingRule.ALL, TestResources.intervalStringsToGenomeLocs("1:1", "1"), TestResources.intervalStringsToGenomeLocs("1") }
        };
    }

    @Test(dataProvider = "sortAndMergeIntervals")
    public void testSortAndMergeIntervals(IntervalMergingRule merge, List<GenomeLoc> unsorted, List<GenomeLoc> expected) {
        List<GenomeLoc> sorted = IntervalUtils.sortAndMergeIntervals(TestResources.getGenomeLocParser(), unsorted, merge).toList();
        Assert.assertEquals(sorted, expected);
    }

    @DataProvider(name="loadintervals")
    public Object[][] getloadIntervals(){
        final String severalIntervals = "src/test/resources/org/broadinstitute/hellbender/utils/interval/example_intervals.list";

        return new Object[][]{
                new Object[]{Arrays.asList("1:1-2"), IntervalSetRule.UNION, 0, TestResources.intervalStringsToGenomeLocs("1:1-2")},
                new Object[]{Arrays.asList("1:1-2", "2:1-2"), IntervalSetRule.UNION, 0, TestResources.intervalStringsToGenomeLocs("1:1-2", "2:1-2")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), IntervalSetRule.UNION, 0, TestResources.intervalStringsToGenomeLocs("1:1-15")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), IntervalSetRule.INTERSECTION, 0, TestResources.intervalStringsToGenomeLocs("1:5-10")},
                new Object[]{Arrays.asList("1:5-5"), IntervalSetRule.UNION, 5, TestResources.intervalStringsToGenomeLocs("1:1-10")},
                new Object[]{Arrays.asList(severalIntervals), IntervalSetRule.UNION, 0, TestResources.intervalStringsToGenomeLocs("1:100-200", "2:20-30", "4:50")}
        };
    }

    @Test( dataProvider = "loadintervals")
    public void loadIntervalsTest(List<String> intervals, IntervalSetRule setRule, int padding, List<GenomeLoc> results){
        GenomeLocSortedSet loadedIntervals = IntervalUtils.loadIntervals(intervals, setRule, IntervalMergingRule.ALL, padding, TestResources.getGenomeLocParser());
        Assert.assertEquals(loadedIntervals, results);
    }

    @DataProvider(name = "loadIntervalsFromFeatureFileData")
    public Object[][] loadIntervalsFromFeatureFileData() {
        final List<GenomeLoc> expectedIntervals = TestResources.intervalStringsToGenomeLocs("1:100-100", "1:203-206", "1:1000-1003", "2:200-200", "2:548-550", "4:776-779");

        return new Object[][] {
                { new File(INTERVAL_TEST_DATA + "intervals_from_features_test.vcf"), expectedIntervals },
                { new File(INTERVAL_TEST_DATA + "intervals_from_features_test.bed"), expectedIntervals }
        };
    }

    @Test(dataProvider = "loadIntervalsFromFeatureFileData")
    public void testLoadIntervalsFromFeatureFile( final File featureFile, final List<GenomeLoc> expectedIntervals ) {
        final GenomeLocSortedSet actualIntervals = IntervalUtils.loadIntervals(Arrays.asList(featureFile.getAbsolutePath()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, TestResources.getGenomeLocParser());
        Assert.assertEquals(actualIntervals, expectedIntervals, "Wrong intervals loaded from Feature file " + featureFile.getAbsolutePath());
    }

    // Note: because the file does not exist and all characters are allowed in contig names,
    // we will not know that this is supposed to be interpreted as a file.
    // So we'll blow up with MalformedGenomeLoc and not anything related to files
    @Test(expectedExceptions = UserException.MalformedGenomeLoc.class)
    public void testLoadIntervalsFromNonExistentFile() {
        IntervalUtils.loadIntervals(Arrays.asList(BaseTest.getSafeNonExistentFile("non_existent_file.vcf").getAbsolutePath()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, TestResources.getGenomeLocParser());
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testLoadIntervalsFromUnrecognizedFormatFile() {
        IntervalUtils.loadIntervals(Arrays.asList(INTERVAL_TEST_DATA + "unrecognized_format_file.xyz"), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, TestResources.getGenomeLocParser());
    }

    @Test(expectedExceptions = UserException.MalformedFile.class)
    public void loadIntervalsEmptyFile(){
        IntervalUtils.loadIntervals(Arrays.asList(emptyIntervals), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, TestResources.getGenomeLocParser());
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData",expectedExceptions = IllegalArgumentException.class)
    public void testGenomeLocsFromLocatablesNullParser(final GenomeLocParser parser,
                                                       final List<? extends Locatable> locatables) {
        IntervalUtils.genomeLocsFromLocatables(null, locatables);
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData",expectedExceptions = IllegalArgumentException.class)
    public void testGenomeLocsFromLocatablesNullLocatables(final GenomeLocParser parser,
                                                       final List<? extends Locatable> locatables) {
        IntervalUtils.genomeLocsFromLocatables(parser,null);
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData",expectedExceptions = IllegalArgumentException.class)
    public void testGenomeLocsFromLocatablesNullContainingLocatables(final GenomeLocParser parser,
                                                           final List<? extends Locatable> locatables) {
        if (locatables.size() == 0) {
            IntervalUtils.genomeLocsFromLocatables(parser,null);
        } else {
            final List<? extends Locatable> withNull = new ArrayList<>(locatables);
            withNull.set(withNull.size() / 2,null);
            IntervalUtils.genomeLocsFromLocatables(parser,withNull);
        }
    }

    @Test(dataProvider = "genomeLocsFromLocatablesData")
    public void testGenomeLocsFromLocatables(final GenomeLocParser parser, final List<? extends Locatable> locatables) {
        final List<GenomeLoc> result = IntervalUtils.genomeLocsFromLocatables(parser, locatables);
        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), locatables.size());
        for (int i = 0; i < result.size(); i++) {
            final GenomeLoc resultLoc = result.get(i);
            final Locatable inputLoc = locatables.get(i);
            Assert.assertEquals(resultLoc.getContig(),inputLoc.getContig());
            Assert.assertEquals(resultLoc.getStart(),inputLoc.getStart());
            Assert.assertEquals(resultLoc.getStop(),inputLoc.getEnd());
        }
        // is it real only:
        try {
            result.add(parser.createGenomeLoc("1",1,2));
            Assert.fail("the result collection should not allow to call add");
        } catch (final UnsupportedOperationException ex) {
            // ok.
        }
        if (result.size() > 0) {
            try {
                result.set(0,parser.createGenomeLoc("1",1,2));
                Assert.fail("the result collection should not allow to call set");
            } catch (final UnsupportedOperationException ex) {
                // ok.
            }
        }
    }

    @DataProvider(name="genomeLocsFromLocatablesData")
    public Object[][] genomeLocsFromLocatablesData() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
        final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());

        return new Object[][] {
                { genomeLocParser, Collections.emptyList() },
                { genomeLocParser, Collections.singletonList(new SimpleInterval("1",1,2)) },
                { genomeLocParser, Arrays.asList(new SimpleInterval("1", 1, 2), new SimpleInterval("1", 1, 1), new SimpleInterval("1", 5, 9)) },
                { genomeLocParser, Arrays.asList(new SimpleInterval("1",9,10),new SimpleInterval("1",5,9), new SimpleInterval("1",1,9)) }
        };
    }

    @Test(dataProvider = "overlapData")
    public void testOverlap(final SimpleInterval l, final SimpleInterval r, final boolean expected) {
        Assert.assertEquals(IntervalUtils.overlaps(l, r), expected);
        Assert.assertEquals(IntervalUtils.overlaps(
                new SimpleFeature(l.getContig(),l.getStart(),l.getEnd()),
                new SimpleFeature(r.getContig(),r.getStart(),r.getEnd())),expected);
        Assert.assertEquals(IntervalUtils.overlaps(
                new SimpleFeature(l.getContig(),l.getStart(),l.getEnd()),
                r),expected);
        Assert.assertEquals(IntervalUtils.overlaps(
                l,
                new SimpleFeature(r.getContig(),r.getStart(),r.getEnd())),expected);
        Assert.assertFalse(IntervalUtils.overlaps(
                new SimpleFeature(null, r.getStart(), r.getEnd()),
                new SimpleFeature(null, l.getStart(), l.getEnd())));
    }

    @DataProvider(name = "overlapData")
    public Object[][] overlapData() {
        // removing cases where the second interval is null.
        return Stream.of(SimpleIntervalUnitTest.getIntervalOverlapData()).filter(a -> a[1] != null)
                .collect(Collectors.toList()).toArray(new Object[0][]);
    }

    @Test(dataProvider = "lexicographicalOrderComparatorData")
    public void testLexicographicalOrderComparator(final List<Locatable> unsorted) {
        final List<Locatable> sorted = unsorted.stream().sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        for (int i = 0; i < sorted.size() - 1; i++) {
           final Locatable thisLoc = sorted.get(i);
           final Locatable nextLoc = sorted.get(i + 1);
           if (thisLoc.getContig() == nextLoc.getContig()) {
               Assert.assertTrue(thisLoc.getStart() <= nextLoc.getStart());
               if (thisLoc.getStart() == nextLoc.getStart()) {
                   Assert.assertTrue(thisLoc.getEnd() <= nextLoc.getEnd());
               }
           } else if (thisLoc.getContig() == null) {
               Assert.fail("Null contig must go last");
           } else if (nextLoc.getContig() != null) {
               Assert.assertTrue(thisLoc.getContig().compareTo(nextLoc.getContig()) <= 0);
               if (thisLoc.getContig().equals(nextLoc.getContig())) {
                   Assert.assertTrue(thisLoc.getStart() <= nextLoc.getStart());
                   if (thisLoc.getStart() == nextLoc.getStart()) {
                       Assert.assertTrue(thisLoc.getEnd() <= nextLoc.getEnd());
                   }
               }
           }
        }
    }

    @DataProvider(name = "lexicographicalOrderComparatorData")
    public Object[][] lexicographicalOrderComparatorData() {
        final String[] CONTIG_NAMES = new String[] {"A","AA","B","C", null};
        final int[] CONTIG_SIZES = new int[] {200,300,400,500,600};
        final int MIN_INTERVAL_SIZE = 1;
        final int MAX_INTERVAL_SIZE = 100;
        final Random rdn = new Random(1131312131);
        final int CASE_NUMBER = 100;
        final List<Object[]> result = new ArrayList<>(100);
        for (int i = 0; i < CASE_NUMBER; i++) {
            final int locatableCount = rdn.nextInt(100) + 1;
            final List<Locatable> locatables = new ArrayList<>(locatableCount);
            for (int j = 0; j < locatableCount; j++) {
                final int contigIdx = rdn.nextInt(CONTIG_NAMES.length);
                final String contig = CONTIG_NAMES[contigIdx];

                final boolean useSimpleInterval = contig == null ? false : rdn.nextBoolean();
                final int intervalSize = rdn.nextInt(MAX_INTERVAL_SIZE - MIN_INTERVAL_SIZE + 1) + MIN_INTERVAL_SIZE;
                final int start = rdn.nextInt(CONTIG_SIZES[contigIdx] - intervalSize) + 1;
                final int end = start + intervalSize - 1;
                final Locatable loc = useSimpleInterval ? new SimpleInterval(contig,start,end) : new SimpleFeature(contig,start,end);
                locatables.add(loc);
            }
            result.add(new Object[] { locatables });
        }
        return result.toArray(new Object[result.size()][]);
    }


    @DataProvider(name="intervals")
    public Object[][] intervals(){
        ArrayList<SimpleInterval> input = Lists.newArrayList(new SimpleInterval("1", 10, 100));
        ArrayList<SimpleInterval> output = input;
        ArrayList<SimpleInterval> hundred = Lists.newArrayList(new SimpleInterval("1", 1, 100));

        return new Object[][]{
                // input fits all in one shard
                new Object[]{input, 1000, output},
                // input could fit in the shard, but is cut anyways so it's aligned to integer multiples of shard boundary
                new Object[]{input, 90, Lists.newArrayList(new SimpleInterval("1", 10, 90),new SimpleInterval("1", 91, 100))},
                // shard ends just after input
                new Object[]{hundred, 101, hundred},
                // shard ends at input
                new Object[]{hundred, 100, hundred},
                // shard ends just before input does
                new Object[]{hundred, 99, Lists.newArrayList(new SimpleInterval("1", 1, 99),new SimpleInterval("1", 100, 100))},
                // input overlaps three shards
                new Object[]{hundred, 40, Lists.newArrayList(new SimpleInterval("1", 1, 40),new SimpleInterval("1", 41, 80),new SimpleInterval("1", 81, 100))},
                // input covers four shards exactly
                new Object[]{hundred, 25, Lists.newArrayList(
                        new SimpleInterval("1", 1, 25),new SimpleInterval("1", 26, 50),new SimpleInterval("1", 51, 75),new SimpleInterval("1", 76, 100))},
        };
    }

    @Test(dataProvider = "intervals")
    public void testCutToShards(ArrayList<SimpleInterval> input, int shardSize, ArrayList<SimpleInterval> expected) throws Exception {
        List<SimpleInterval> actual = IntervalUtils.cutToShards(input, shardSize);
        Assert.assertEquals(
                actual,
                expected
        );
    }

    @DataProvider(name="shardIndex")
    public Object[][] shardIndex(){
        return new Object[][]{
                // offset, shardSize, expected shardIndex
                new Object[]{1,10,0},
                new Object[]{10,10,0},
                new Object[]{11,10,1},
                new Object[]{20,10,1},
                new Object[]{21,10,2}
        };
    }

    @Test(dataProvider = "shardIndex")
    public void testShardIndex(int offset, int shardSize, int shardIndex) throws Exception {
        Assert.assertEquals(
                IntervalUtils.shardIndex(offset, shardSize),
                shardIndex
        );
    }

    @DataProvider(name="shardBegin")
    public Object[][] shardBegin(){
        return new Object[][]{
                // shard index, shard size, shard begin
                new Object[]{0,10,1},
                new Object[]{1,10,11},
                new Object[]{1,42,43},
                new Object[]{2,10,21},
                new Object[]{2,42,85}
        };
    }

    @Test(dataProvider = "shardBegin")
    public void testBeginOfShard(int index, int shardSize, int shardBegin) throws Exception {
        Assert.assertEquals(
                IntervalUtils.beginOfShard(index, shardSize),
                shardBegin
        );
    }

    @DataProvider(name="shardEnd")
    public Object[][] shardEnd(){
        return new Object[][]{
                // shard index, shard size, shard end
                new Object[]{0,10,10},
                new Object[]{1,10,20},
                new Object[]{2,10,30},
                new Object[]{0,42,42},
                new Object[]{1,42,84},
        };
    }

    @Test(dataProvider = "shardEnd")
    public void testEndOfShard(int index, int shardSize, int shardEnd) throws Exception {
        Assert.assertEquals(
                IntervalUtils.endOfShard(index, shardSize),
                shardEnd
        );
    }

    @Test
    public void testConvertSimpleIntervalToQueryInterval() {
        final SAMSequenceRecord contigRecord = new SAMSequenceRecord("1", 100);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(contigRecord));
        final SimpleInterval originalInterval = new SimpleInterval("1", 5, 10);

        final QueryInterval convertedInterval = IntervalUtils.convertSimpleIntervalToQueryInterval(originalInterval, dictionary);
        Assert.assertEquals(convertedInterval.referenceIndex, 0);
        Assert.assertEquals(convertedInterval.start, 5);
        Assert.assertEquals(convertedInterval.end, 10);
    }

    @DataProvider(name = "TrimIntervalToContigData")
    public Object[][] trimIntervalToContigData() {
        return new Object[][] {
                { "1", 5, 10, 100, new SimpleInterval("1", 5, 10) },
                { "1", 1, 100, 100, new SimpleInterval("1", 1, 100) },
                { "1", -1, 100, 100, new SimpleInterval("1", 1, 100) },
                { "1", -5, 100, 100, new SimpleInterval("1", 1, 100) },
                { "1", -5, 90, 100, new SimpleInterval("1", 1, 90) },
                { "1", 1, 101, 100, new SimpleInterval("1", 1, 100) },
                { "1", 1, 105, 100, new SimpleInterval("1", 1, 100) },
                { "1", 5, 105, 100, new SimpleInterval("1", 5, 100) },
                { "1", -1, 10, 100, new SimpleInterval("1", 1, 10) },
                { "1", -5, 10, 100, new SimpleInterval("1", 1, 10) },
                { "1", 90, 101, 100, new SimpleInterval("1", 90, 100) },
                { "1", 90, 105, 100, new SimpleInterval("1", 90, 100) }
        };
    }

    @Test(dataProvider = "TrimIntervalToContigData")
    public void testTrimIntervalToContig( final String contig, final int start, final int stop, final int contigLength, final SimpleInterval expectedInterval ) {
        Assert.assertEquals(IntervalUtils.trimIntervalToContig(contig, start, stop, contigLength), expectedInterval);
    }

    @DataProvider(name = "TrimIntervalToContigInvalidData")
    public Object[][] trimIntervalToContigInvalidData() {
        return new Object[][] {
                { null, 1, 10, 100 },
                { "1", 1, 10, 0 },
                { "1", 1, 10, -1 }
        };
    }

    @Test(dataProvider = "TrimIntervalToContigInvalidData", expectedExceptions = IllegalArgumentException.class)
    public void testTrimIntervalToContigInvalidArg( final String contig, final int start, final int stop, final int contigLength ) {
        IntervalUtils.trimIntervalToContig(contig, start, stop, contigLength);
    }

    @DataProvider(name = "IntervalIsOnDictionaryContigData")
    public Object[][] intervalIsOnDictionaryContigData() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 1000)));

        return new Object[][] {
                { new SimpleInterval("1", 1, 10), dictionary, true },
                { new SimpleInterval("1", 50, 100), dictionary, true },
                { new SimpleInterval("1", 1, 1000), dictionary, true },
                { new SimpleInterval("2", 1, 10), dictionary, false },
                { new SimpleInterval("1", 1, 1001), dictionary, false },
                { new SimpleInterval("1", 1, 2000), dictionary, false }
        };
    }

    @Test(dataProvider = "IntervalIsOnDictionaryContigData")
    public void testIntervalIsOnDictionaryContig( final SimpleInterval interval, final SAMSequenceDictionary dictionary, final boolean expectedResult ) {
        Assert.assertEquals(IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary), expectedResult);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalIsOnDictionaryContigNullInterval() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 1000)));
        IntervalUtils.intervalIsOnDictionaryContig(null, dictionary);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalIsOnDictionaryContigNullDictionary() {
        IntervalUtils.intervalIsOnDictionaryContig(new SimpleInterval("1", 1, 10), null);
    }

    @Test
    public void testToGatkIntervalString(){
            final SimpleInterval interval = new SimpleInterval("1",1,100);
            Assert.assertEquals(IntervalUtils.locatableToString(interval), "1:1-100");
    }
}
