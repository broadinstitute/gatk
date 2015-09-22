package org.broadinstitute.hellbender.utils.interval;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.SimpleFeature;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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
    public static final String INTERVAL_TEST_DATA = publicTestDir + "org/broadinstitute/hellbender/utils/interval/";
    public static final String emptyIntervals = BaseTest.publicTestDir + "empty_intervals.list";
    private List<GenomeLoc> hg19ReferenceLocs;
    private List<GenomeLoc> hg19exomeIntervals;


    @BeforeClass
    public void init() {
        hg19ReferenceLocs = Collections.unmodifiableList(GenomeLocSortedSet.createSetFromSequenceDictionary(hg19Header.getSequenceDictionary()).toList()) ;
        hg19exomeIntervals = Collections.unmodifiableList(IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(hg19MiniIntervalFile)));
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
        GenomeLoc bp01_10 = hg19GenomeLocParser.createGenomeLoc("1", 1, 10);

        GenomeLoc bp1_5 = hg19GenomeLocParser.createGenomeLoc("1", 1, 5);
        GenomeLoc bp6_10 = hg19GenomeLocParser.createGenomeLoc("1", 6, 10);
        new SplitLocusIntervalsSmallTest("cut into two", Arrays.asList(bp01_10), 2, Arrays.asList(bp1_5, bp6_10));

        GenomeLoc bp20_30 = hg19GenomeLocParser.createGenomeLoc("1", 20, 30);
        new SplitLocusIntervalsSmallTest("two in two", Arrays.asList(bp01_10, bp20_30), 2, Arrays.asList(bp01_10, bp20_30));

        GenomeLoc bp1_7 = hg19GenomeLocParser.createGenomeLoc("1", 1, 7);
        GenomeLoc bp8_10 = hg19GenomeLocParser.createGenomeLoc("1", 8, 10);
        GenomeLoc bp20_23 = hg19GenomeLocParser.createGenomeLoc("1", 20, 23);
        GenomeLoc bp24_30 = hg19GenomeLocParser.createGenomeLoc("1", 24, 30);
        new SplitLocusIntervalsSmallTest("two in three", Arrays.asList(bp01_10, bp20_30), 3,
                Arrays.asList(bp1_7, bp8_10, bp20_23, bp24_30));

        GenomeLoc bp1_2 = hg19GenomeLocParser.createGenomeLoc("1", 1, 2);
        GenomeLoc bp1_1 = hg19GenomeLocParser.createGenomeLoc("1", 1, 1);
        GenomeLoc bp2_2 = hg19GenomeLocParser.createGenomeLoc("1", 2, 2);
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
                listEveryTwoFromTwo.add(hg19GenomeLocParser.createGenomeLoc("chr1",x,x));
            else
                listEveryTwoFromOne.add(hg19GenomeLocParser.createGenomeLoc("chr1",x,x));
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
                listEveryTwoFromTwo.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
            allSites.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
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
                listEveryTwoFromTwo.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
                allSites.add(hg19GenomeLocParser.createGenomeLoc("1",x,x));
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

        source1.add(hg19GenomeLocParser.createGenomeLoc("1", 10, 20));
        source1.add(hg19GenomeLocParser.createGenomeLoc("1", 15, 25));

        source2.add(hg19GenomeLocParser.createGenomeLoc("1", 16, 18));
        source2.add(hg19GenomeLocParser.createGenomeLoc("1", 22, 24));

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(source1, source2, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 2);
    }

    @Test
    public void testGetContigLengths() {
        Map<String, Integer> lengths = IntervalUtils.getContigSizes(new File(BaseTest.exampleReference));
        Assert.assertEquals((long)lengths.get("1"), 16000);
        Assert.assertEquals((long)lengths.get("2"), 16000);
        Assert.assertEquals((long) lengths.get("3"), 16000);
        Assert.assertEquals((long)lengths.get("4"), 16000);
    }

    @Test
    public void testParseIntervalArguments() {
        Assert.assertEquals(hg19ReferenceLocs.size(), 4);
        Assert.assertEquals(intervalStringsToGenomeLocs("1", "2", "3").size(), 3);
        Assert.assertEquals(intervalStringsToGenomeLocs("1:1-2", "1:4-5", "2:1-1", "3:2-2").size(), 4);
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
        IntervalUtils.isIntervalFile(BaseTest.publicTestDir + "no_such_intervals.list");
    }

    @Test
    public void testFixedScatterIntervalsBasic() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");

        List<File> files = testFiles("basic.", 3, ".intervals");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2", "3");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterFixedIntervalsLessFiles() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");
        GenomeLoc chr4 = hg19GenomeLocParser.parseGenomeLoc("4");

        List<File> files = testFiles("less.", 3, ".intervals");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2", "3", "4");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
        Assert.assertEquals(locs3.get(1), chr4);
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testSplitFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2");
        IntervalUtils.splitFixedIntervals(locs, files.size());
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testScatterFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = intervalStringsToGenomeLocs("1", "2");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, locs.size()); // locs.size() instead of files.size()
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);
    }
    @Test
    public void testScatterFixedIntervalsStart() {
        List<String> intervals = Arrays.asList("1:1-2", "1:4-5", "2:1-1", "3:2-2");
        GenomeLoc chr1a = hg19GenomeLocParser.parseGenomeLoc("1:1-2");
        GenomeLoc chr1b = hg19GenomeLocParser.parseGenomeLoc("1:4-5");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2a = hg19GenomeLocParser.parseGenomeLoc("2:1-2");
        GenomeLoc chr2b = hg19GenomeLocParser.parseGenomeLoc("2:4-5");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:2-2");
        GenomeLoc chr3a = hg19GenomeLocParser.parseGenomeLoc("3:1-2");
        GenomeLoc chr3b = hg19GenomeLocParser.parseGenomeLoc("3:4-5");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = intervalStringsToGenomeLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(hg19MiniIntervalFile));
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

        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        int locIndex = 0;
        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
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
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    @Test
    public void testScatterContigIntervalsOrder() {
        List<String> intervals = Arrays.asList("2:1-1", "1:1-1", "3:2-2");
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr2);
        Assert.assertEquals(locs2.get(0), chr1);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsBasic() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");

        List<File> files = testFiles("contig_basic.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs("1", "2", "3"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsLessFiles() {
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3");
        GenomeLoc chr4 = hg19GenomeLocParser.parseGenomeLoc("4");

        List<File> files = testFiles("contig_less.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs("1", "2", "3", "4"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs("1", "2"), files);
    }

    @Test
    public void testScatterContigIntervalsStart() {
        List<String> intervals = Arrays.asList("1:1-2", "1:4-5", "2:1-1", "3:2-2");
        GenomeLoc chr1a = hg19GenomeLocParser.parseGenomeLoc("1:1-2");
        GenomeLoc chr1b = hg19GenomeLocParser.parseGenomeLoc("1:4-5");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:1-1");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("contig_split_start.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2a = hg19GenomeLocParser.parseGenomeLoc("2:1-2");
        GenomeLoc chr2b = hg19GenomeLocParser.parseGenomeLoc("2:4-5");
        GenomeLoc chr3 = hg19GenomeLocParser.parseGenomeLoc("3:2-2");

        List<File> files = testFiles("contig_split_middle.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        GenomeLoc chr1 = hg19GenomeLocParser.parseGenomeLoc("1:1-1");
        GenomeLoc chr2 = hg19GenomeLocParser.parseGenomeLoc("2:2-2");
        GenomeLoc chr3a = hg19GenomeLocParser.parseGenomeLoc("3:1-2");
        GenomeLoc chr3b = hg19GenomeLocParser.parseGenomeLoc("3:4-5");

        List<File> files = testFiles("contig_split_end.", 3 ,".intervals");

        IntervalUtils.scatterContigIntervals(hg19Header, intervalStringsToGenomeLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(files.get(2).toString()));

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
        IntervalUtils.scatterContigIntervals(hg19Header, hg19ReferenceLocs, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
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
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Collections.singletonList(publicTestDir + unmergedIntervals));
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
    public void testInvalidPicardIntervalHandling(GenomeLocParser genomeLocParser,
                                                  String contig, int intervalStart, int intervalEnd ) throws Exception {

        SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("1"));
        // Intentionally add an extra contig not in our genomeLocParser's sequence dictionary
        // so that we can add intervals to the Picard interval file for contigs not in our dictionary
        picardFileHeader.addSequence(new SAMSequenceRecord("2", 100000));

        IntervalList picardIntervals = new IntervalList(picardFileHeader);
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname"));

        File picardIntervalFile = createTempFile("testInvalidPicardIntervalHandling", ".intervals");
        picardIntervals.write(picardIntervalFile);

        List<String> intervalArgs = new ArrayList<>(1);
        intervalArgs.add(picardIntervalFile.getAbsolutePath());

        // loadIntervals() will validate all intervals against the sequence dictionary in our genomeLocParser,
        // and should throw for all intervals in our invalidIntervalTestData set
        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    /*
    Split into tests that can be written to files and tested by writeFlankingIntervals,
    and lists that cannot but are still handled by getFlankingIntervals.
    */
    private static abstract class FlankingIntervalsTestData extends TestDataProvider {
        final public File referenceFile;
        final public GenomeLocParser parser;
        final int basePairs;
        final List<GenomeLoc> original;
        final List<GenomeLoc> expected;

        protected FlankingIntervalsTestData(Class<?> clazz, String name, File referenceFile, GenomeLocParser parser,
                                            int basePairs, List<String> original, List<String> expected) {
            super(clazz, name);
            this.referenceFile = referenceFile;
            this.parser = parser;
            this.basePairs = basePairs;
            this.original = parse(parser, original);
            this.expected = parse(parser, expected);
        }

        private static List<GenomeLoc> parse(GenomeLocParser parser, List<String> locs) {
            List<GenomeLoc> parsed = new ArrayList<>();
            for (String loc: locs)
                parsed.add("unmapped".equals(loc) ? GenomeLoc.UNMAPPED : parser.parseGenomeLoc(loc));
            return parsed;
        }
    }

    private static class FlankingIntervalsFile extends FlankingIntervalsTestData {
        public FlankingIntervalsFile(String name, File referenceFile, GenomeLocParser parser,
                                     int basePairs, List<String> original, List<String> expected) {
            super(FlankingIntervalsFile.class, name, referenceFile, parser, basePairs, original, expected);
        }
    }

    private static class FlankingIntervalsList extends FlankingIntervalsTestData {
        public FlankingIntervalsList(String name, File referenceFile, GenomeLocParser parser,
                                     int basePairs, List<String> original, List<String> expected) {
            super(FlankingIntervalsList.class, name, referenceFile, parser, basePairs, original, expected);
        }
    }

    /* Intervals where the original and the flanks can be written to files. */
    @DataProvider(name = "flankingIntervalsFiles")
    public Object[][] getFlankingIntervalsFiles() {
        File hg19ReferenceFile = new File(BaseTest.exampleReference);
        int hg19Length1 = hg19GenomeLocParser.getContigInfo("1").getSequenceLength();

        new FlankingIntervalsFile("atStartBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:1"),
                Arrays.asList("1:2"));

        new FlankingIntervalsFile("atStartBase50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:1"),
                Arrays.asList("1:2-51"));

        new FlankingIntervalsFile("atStartRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:1-10"),
                Arrays.asList("1:11-60"));

        new FlankingIntervalsFile("atEndBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:" + hg19Length1),
                Arrays.asList("1:" + (hg19Length1 - 1)));

        new FlankingIntervalsFile("atEndBase50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:" + hg19Length1),
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 50, hg19Length1 - 1)));

        new FlankingIntervalsFile("atEndRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 10, hg19Length1)),
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 60, hg19Length1 - 11)));

        new FlankingIntervalsFile("nearStartBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:2"),
                Arrays.asList("1:1", "1:3"));

        new FlankingIntervalsFile("nearStartRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:21-30"),
                Arrays.asList("1:1-20", "1:31-80"));

        new FlankingIntervalsFile("nearEndBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:" + (hg19Length1 - 1)),
                Arrays.asList("1:" + (hg19Length1 - 2), "1:" + hg19Length1));

        new FlankingIntervalsFile("nearEndRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 30, hg19Length1 - 21)),
                Arrays.asList(
                        String.format("1:%d-%d", hg19Length1 - 80, hg19Length1 - 31),
                        String.format("1:%d-%d", hg19Length1 - 20, hg19Length1)));

        new FlankingIntervalsFile("beyondStartBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:3"),
                Arrays.asList("1:2", "1:4"));

        new FlankingIntervalsFile("beyondStartRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200"),
                Arrays.asList("1:51-100", "1:201-250"));

        new FlankingIntervalsFile("beyondEndBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:" + (hg19Length1 - 3)),
                Arrays.asList("1:" + (hg19Length1 - 4), "1:" + (hg19Length1 - 2)));

        new FlankingIntervalsFile("beyondEndRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 200, hg19Length1 - 101)),
                Arrays.asList(
                        String.format("1:%d-%d", hg19Length1 - 250, hg19Length1 - 201),
                        String.format("1:%d-%d", hg19Length1 - 100, hg19Length1 - 51)));

        new FlankingIntervalsFile("betweenFar50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:401-500"),
                Arrays.asList("1:51-100", "1:201-250", "1:351-400", "1:501-550"));

        new FlankingIntervalsFile("betweenSpan50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:301-400"),
                Arrays.asList("1:51-100", "1:201-300", "1:401-450"));

        new FlankingIntervalsFile("betweenOverlap50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:271-400"),
                Arrays.asList("1:51-100", "1:201-270", "1:401-450"));

        new FlankingIntervalsFile("betweenShort50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:221-400"),
                Arrays.asList("1:51-100", "1:201-220", "1:401-450"));

        new FlankingIntervalsFile("betweenNone50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:121-400"),
                Arrays.asList("1:51-100", "1:401-450"));

        new FlankingIntervalsFile("twoContigs", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "2:301-400"),
                Arrays.asList("1:51-100", "1:201-250", "2:251-300", "2:401-450"));

        // Explicit testing a problematic agilent target pair
        new FlankingIntervalsFile("badAgilent", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("2:6257-6411", "2:6487-6628"),
                // wrong!    ("2:6206-6256", "2:6412-6462", "2:6436-6486", "2:6629-6679")
                Arrays.asList("2:6207-6256", "2:6412-6486", "2:6629-6678"));

        return TestDataProvider.getTests(FlankingIntervalsFile.class);
    }

    /* Intervals where either the original and/or the flanks cannot be written to a file. */
    @DataProvider(name = "flankingIntervalsLists")
    public Object[][] getFlankingIntervalsLists() {
        File hg19ReferenceFile = new File(BaseTest.exampleReference);
        List<String> empty = Collections.emptyList();

        new FlankingIntervalsList("empty", hg19ReferenceFile, hg19GenomeLocParser, 50,
                empty,
                empty);

        new FlankingIntervalsList("unmapped", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("unmapped"),
                empty);

        new FlankingIntervalsList("fullContig", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1"),
                empty);

        new FlankingIntervalsList("fullContigs", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1", "2", "3"),
                empty);

        new FlankingIntervalsList("betweenWithUnmapped", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:301-400", "unmapped"),
                Arrays.asList("1:51-100", "1:201-300", "1:401-450"));

        return TestDataProvider.getTests(FlankingIntervalsList.class);
    }

    @Test(dataProvider = "flankingIntervalsFiles")
    public void testWriteFlankingIntervals(FlankingIntervalsTestData data) throws Exception {
        File originalFile = createTempFile("original.", ".intervals");
        File flankingFile = createTempFile("flanking.", ".intervals");
        try {
            List<String> lines = new ArrayList<>();
            for (GenomeLoc loc: data.original)
                lines.add(loc.toString());
            FileUtils.writeLines(originalFile, lines);

            IntervalUtils.writeFlankingIntervals(data.referenceFile, originalFile, flankingFile, data.basePairs);

            List<GenomeLoc> actual = IntervalUtils.intervalFileToList(data.parser, flankingFile.getAbsolutePath());

            String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                    data.toString(), data.original, actual, data.expected);
            Assert.assertEquals(actual, data.expected, description);
        } finally {
            FileUtils.deleteQuietly(originalFile);
            FileUtils.deleteQuietly(flankingFile);
        }
    }

    @Test(dataProvider = "flankingIntervalsLists", expectedExceptions = UserException.class)
    public void testWritingBadFlankingIntervals(FlankingIntervalsTestData data) throws Exception {
        File originalFile = createTempFile("original.", ".intervals");
        File flankingFile = createTempFile("flanking.", ".intervals");
        try {
            List<String> lines = new ArrayList<>();
            for (GenomeLoc loc: data.original)
                lines.add(loc.toString());
            FileUtils.writeLines(originalFile, lines);

            // Should throw a user exception on bad input if either the original
            // intervals are empty or if the flanking intervals are empty
            IntervalUtils.writeFlankingIntervals(data.referenceFile, originalFile, flankingFile, data.basePairs);
        } finally {
            FileUtils.deleteQuietly(originalFile);
            FileUtils.deleteQuietly(flankingFile);
        }
    }

    @Test(dataProvider = "flankingIntervalsLists")
    public void testGetFlankingIntervals(FlankingIntervalsTestData data) {
        List<GenomeLoc> actual = IntervalUtils.getFlankingIntervals(data.parser, data.original, data.basePairs);
        String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                data.toString(), data.original, actual, data.expected);
        Assert.assertEquals(actual, data.expected, description);
    }


    @DataProvider(name="invalidIntervalTestData")
    public Object[][] invalidIntervalDataProvider() throws Exception {
        File fastaFile = new File(publicTestDir + "exampleFASTA.fasta");
        GenomeLocParser genomeLocParser = new GenomeLocParser(new IndexedFastaSequenceFile(fastaFile));

        return new Object[][] {
                new Object[] {genomeLocParser, "1", 10000000, 20000000},
                new Object[] {genomeLocParser, "2", 1, 2},
                new Object[] {genomeLocParser, "1", -1, 50}
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
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, intervalStringsToGenomeLocs("1:1", "1:3", "1:2"), intervalStringsToGenomeLocs("1:1", "1:2", "1:3") },
                new Object[] { IntervalMergingRule.ALL, intervalStringsToGenomeLocs("1:1", "1:3", "1:2"), intervalStringsToGenomeLocs("1:1-3") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, intervalStringsToGenomeLocs("1:1", "1:3", "2:2"), intervalStringsToGenomeLocs("1:1", "1:3", "2:2") },
                new Object[] { IntervalMergingRule.ALL, intervalStringsToGenomeLocs("1:1", "1:3", "2:2"), intervalStringsToGenomeLocs("1:1", "1:3", "2:2") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, intervalStringsToGenomeLocs("1:1", "1"), intervalStringsToGenomeLocs("1") },
                new Object[] { IntervalMergingRule.ALL, intervalStringsToGenomeLocs("1:1", "1"), intervalStringsToGenomeLocs("1") }
        };
    }

    @Test(dataProvider = "sortAndMergeIntervals")
    public void testSortAndMergeIntervals(IntervalMergingRule merge, List<GenomeLoc> unsorted, List<GenomeLoc> expected) {
        List<GenomeLoc> sorted = IntervalUtils.sortAndMergeIntervals(hg19GenomeLocParser, unsorted, merge).toList();
        Assert.assertEquals(sorted, expected);
    }

    @DataProvider(name="loadintervals")
    public Object[][] getloadIntervals(){
        final String severalIntervals = "src/test/resources/org/broadinstitute/hellbender/utils/interval/example_intervals.list";

        return new Object[][]{
                new Object[]{Arrays.asList("1:1-2"), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:1-2")},
                new Object[]{Arrays.asList("1:1-2", "2:1-2"), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:1-2", "2:1-2")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:1-15")},
                new Object[]{Arrays.asList("1:1-10", "1:5-15"), IntervalSetRule.INTERSECTION, 0, intervalStringsToGenomeLocs("1:5-10")},
                new Object[]{Arrays.asList("1:5-5"), IntervalSetRule.UNION, 5, intervalStringsToGenomeLocs("1:1-10")},
                new Object[]{Arrays.asList(severalIntervals), IntervalSetRule.UNION, 0, intervalStringsToGenomeLocs("1:100-200", "2:20-30", "4:50")}
        };
    }

    @Test( dataProvider = "loadintervals")
    public void loadIntervalsTest(List<String> intervals, IntervalSetRule setRule, int padding, List<GenomeLoc> results){
        GenomeLocSortedSet loadedIntervals = IntervalUtils.loadIntervals(intervals, setRule, IntervalMergingRule.ALL, padding, hg19GenomeLocParser);
        Assert.assertEquals(loadedIntervals, results);
    }

    @DataProvider(name = "loadIntervalsFromFeatureFileData")
    public Object[][] loadIntervalsFromFeatureFileData() {
        final List<GenomeLoc> expectedIntervals = intervalStringsToGenomeLocs("1:100-100", "1:203-206", "1:1000-1003", "2:200-200", "2:548-550", "4:776-779");

        return new Object[][] {
                { new File(INTERVAL_TEST_DATA + "intervals_from_features_test.vcf"), expectedIntervals },
                { new File(INTERVAL_TEST_DATA + "intervals_from_features_test.bed"), expectedIntervals }
        };
    }

    @Test(dataProvider = "loadIntervalsFromFeatureFileData")
    public void testLoadIntervalsFromFeatureFile( final File featureFile, final List<GenomeLoc> expectedIntervals ) {
        final GenomeLocSortedSet actualIntervals = IntervalUtils.loadIntervals(Arrays.asList(featureFile.getAbsolutePath()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
        Assert.assertEquals(actualIntervals, expectedIntervals, "Wrong intervals loaded from Feature file " + featureFile.getAbsolutePath());
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testLoadIntervalsFromNonExistentFile() {
        IntervalUtils.loadIntervals(Arrays.asList(publicTestDir + "non_existent_file.vcf"), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testLoadIntervalsFromUnrecognizedFormatFile() {
        IntervalUtils.loadIntervals(Arrays.asList(INTERVAL_TEST_DATA + "unrecognized_format_file.xyz"), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
    }

    @Test(expectedExceptions = UserException.MalformedFile.class)
    public void loadIntervalsEmptyFile(){
        IntervalUtils.loadIntervals(Arrays.asList(emptyIntervals), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, hg19GenomeLocParser);
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
}
