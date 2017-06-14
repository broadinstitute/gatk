package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.ojalgo.netio.BufferedInputStreamReader;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

public class PSScoreUtilsTest extends CommandLineProgramTest {

    final static double MIN_COV = 0.90;
    final static double MIN_IDENT = 0.90;
    Map<String, String> refNameToTax;

    private static String getTestCigar(final int readLength, final int insertionLength, final int deletionLength, final int clipLength) {
        final int numMatchOrMismatch = readLength - insertionLength - clipLength;
        if (numMatchOrMismatch < 0) {
            throw new TestException("Clip length plus insertions plus clipping was greater than read length");
        }
        String cigar = "";
        if (clipLength > 0) {
            cigar = cigar + clipLength + "S";
        }
        if (insertionLength > 0) {
            cigar = cigar + insertionLength + "I";
        }
        if (deletionLength > 0) {
            cigar = cigar + deletionLength + "D";
        }
        if (numMatchOrMismatch > 0) {
            cigar = cigar + numMatchOrMismatch + "M";
        }
        return cigar;
    }

    private static GATKRead generateRead(final int length, final List<Integer> NM, final List<Integer> clip,
                                         final List<Integer> insert, final List<Integer> delete,
                                         final List<String> contig, final String altTag) {
        final GATKRead read = ArtificialReadUtils.createRandomRead(length);
        read.setAttribute("NM", NM.get(0));
        read.setPosition(contig.get(0), 1);
        read.setCigar(getTestCigar(length, insert.get(0), delete.get(0), clip.get(0)));
        if (NM.size() > 1) {
            String tagValue = "";
            for (int i = 1; i < NM.size(); i++) {
                if (altTag.equals("XA")) {
                    tagValue += contig.get(i) + ",+1," + getTestCigar(length, insert.get(i), delete.get(i), clip.get(i)) + "," + NM.get(i) + ";";
                } else if (altTag.equals("SA")) {
                    tagValue += contig.get(i) + ",1,+," + getTestCigar(length, insert.get(i), delete.get(i), clip.get(i)) + ",0," + NM.get(i) + ";";
                } else {
                    throw new TestException("Unknown tag " + altTag);
                }
            }
            read.setAttribute(altTag, tagValue);
        }
        return read;
    }

    private static Iterable<GATKRead> generateUnpairedRead(final int length, final List<Integer> NM, final List<Integer> clip,
                                                           final List<Integer> insert, final List<Integer> delete,
                                                           final List<String> contig, final String altTag) {

        final List<GATKRead> pair = new ArrayList<>(1);
        final GATKRead read1 = generateRead(length, NM, clip, insert, delete, contig, altTag);
        read1.setIsPaired(true);
        pair.add(read1);
        return pair;
    }

    private static Iterable<GATKRead> generateReadPair(final int length, final List<Integer> NM1, final List<Integer> NM2,
                                                       final List<Integer> clip1, final List<Integer> clip2,
                                                       final List<Integer> insert1, final List<Integer> insert2,
                                                       final List<Integer> delete1, final List<Integer> delete2,
                                                       final List<String> contig1, final List<String> contig2,
                                                       final String altTag) {
        final List<GATKRead> pair = new ArrayList<>(2);
        final GATKRead read1 = generateRead(length, NM1, clip1, insert1, delete1, contig1, altTag);
        read1.setIsPaired(true);
        pair.add(read1);

        final GATKRead read2 = generateRead(length, NM2, clip2, insert2, delete2, contig2, altTag);
        read2.setIsPaired(true);
        pair.add(read2);

        return pair;
    }

    @BeforeMethod
    public void before() {
        refNameToTax = new HashMap<>(3);
        refNameToTax.put("recordA", "1");
        refNameToTax.put("recordB", "2");
        refNameToTax.put("recordC", "3");
    }

    @Test(groups = "spark")
    public void testDoHeaderWarnings() {

        final PipelineOptions options = null;
        final File file = createTempFile("header_warnings", "txt");
        final PSTree tree = new PSTree("1");
        final Map<String, String> map = new HashMap<>(10);
        map.put("seq1", "2");
        map.put("seq2", "2");
        map.put("seq3", "3");
        final PSTaxonomyDatabase taxDB = new PSTaxonomyDatabase(tree, map);

        //Header and taxDB with matching records -> no warnings expected
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("seq1", 1000));
        header.addSequence(new SAMSequenceRecord("seq2", 1000));
        header.addSequence(new SAMSequenceRecord("seq3", 1000));
        PSScoreUtils.writeHeaderWarnings(file.getAbsolutePath(), header, taxDB, logger);
        try {
            final BufferedReader inputStream = new BufferedInputStreamReader(new FileInputStream(file));
            Assert.assertNull(inputStream.readLine());
        } catch (IOException e) {
            throw new TestException("Could not open temporary file", e);
        }

        //Add sequence not in taxDB -> expect warning to be written
        header.addSequence(new SAMSequenceRecord("seq4", 1000));
        PSScoreUtils.writeHeaderWarnings(file.getAbsolutePath(), header, taxDB, logger);
        try {
            final BufferedReader inputStream = new BufferedInputStreamReader(new FileInputStream(file));
            Assert.assertNotNull(inputStream.readLine());
        } catch (IOException e) {
            throw new TestException("Could not open temporary file", e);
        }
    }

    @DataProvider(name = "mapPairs")
    public Object[][] getMapPairData() {
        return new Object[][]{

                //Perfect match
                {101, Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("recordA"), Arrays.asList("1")},

                //First read sub-threshold mismatches
                {101, Arrays.asList(10), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("recordA"), Arrays.asList("1")},

                //First read too many mismatches
                {101, Arrays.asList(11), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("recordA"), Arrays.asList()},

                //Second read sub-threshold mismatches
                {101, Arrays.asList(0), Arrays.asList(10), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("recordA"), Arrays.asList("1")},

                //Second read too many mismatches
                {101, Arrays.asList(0), Arrays.asList(11), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("recordA"), Arrays.asList()},

                //Indistinct hits
                {101, Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordA"), Arrays.asList("recordA", "recordA"), Arrays.asList("1")},

                //Chimeric hits
                {101, Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordB"), Arrays.asList("recordA", "recordA"), Arrays.asList("1")},

                //Dual hits
                {101, Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordB"), Arrays.asList("recordA", "recordB"), Arrays.asList("1", "2")},

                //Sequence name truncation
                {101, Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA 1", "recordB 2"), Arrays.asList("recordA 1", "recordB 2"), Arrays.asList("1", "2")},

                //Second hit has too many mismatches
                {101, Arrays.asList(0, 11), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordB"), Arrays.asList("recordA", "recordB"), Arrays.asList("1")},

                //Alternate alignment better than primary alignment
                {101, Arrays.asList(11, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordB"), Arrays.asList("recordA", "recordB"), Arrays.asList("2")}
        };
    }

    @Test(dataProvider = "mapPairs", groups = "spark")
    public void testMapGroupedReadsToTax(final int readLength, final List<Integer> NM1, final List<Integer> NM2,
                                         final List<Integer> clip1, final List<Integer> clip2,
                                         final List<Integer> insert1, final List<Integer> insert2,
                                         final List<Integer> delete1, final List<Integer> delete2,
                                         final List<String> contig1, final List<String> contig2,
                                         final List<String> truthTax) {

        if (!(NM1.size() == clip1.size() && NM1.size() == insert1.size() && NM1.size() == delete1.size() && NM1.size() == contig1.size())
                || !(NM2.size() == clip2.size() && NM2.size() == insert2.size() && NM2.size() == delete2.size() && NM2.size() == contig2.size())) {
            throw new TestException("Input lists for each read must be of uniform length");
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final Broadcast<Map<String, String>> refNameToTaxBroadcast = ctx.broadcast(refNameToTax);

        //Test with alternate alignments assigned to the XA tag
        final List<Iterable<GATKRead>> readListXA = new ArrayList<>();
        readListXA.add(generateReadPair(readLength, NM1, NM2, clip1, clip2, insert1, insert2, delete1, delete2, contig1, contig2, "XA"));
        final JavaRDD<Iterable<GATKRead>> pairsXA = ctx.parallelize(readListXA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultXA = PSScoreUtils.mapGroupedReadsToTax(pairsXA, MIN_COV, MIN_IDENT, refNameToTaxBroadcast);
        final PSPathogenAlignmentHit infoXA = resultXA.first()._2;

        Assert.assertNotNull(infoXA);
        Assert.assertEquals(infoXA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoXA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoXA.numMates, 2);

        //Test SA tag
        final List<Iterable<GATKRead>> readListSA = new ArrayList<>();
        readListSA.add(generateReadPair(readLength, NM1, NM2, clip1, clip2, insert1, insert2, delete1, delete2, contig1, contig2, "SA"));
        final JavaRDD<Iterable<GATKRead>> pairsSA = ctx.parallelize(readListSA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultSA = PSScoreUtils.mapGroupedReadsToTax(pairsSA, MIN_COV, MIN_IDENT, refNameToTaxBroadcast);
        final PSPathogenAlignmentHit infoSA = resultSA.first()._2;

        Assert.assertNotNull(infoSA);
        Assert.assertEquals(infoSA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoSA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoSA.numMates, 2);
    }

    @DataProvider(name = "mapUnpaired")
    public Object[][] getMapPairDataUnpaired() {
        return new Object[][]{

                //Perfect match
                {101, Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Sub-threshold mismatches
                {101, Arrays.asList(10), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Too many mismatches
                {101, Arrays.asList(11), Arrays.asList(0), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Sub-threshold clip
                {101, Arrays.asList(0), Arrays.asList(10), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Too much clip
                {101, Arrays.asList(0), Arrays.asList(11), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Sub-threshold insertions
                {101, Arrays.asList(0), Arrays.asList(0), Arrays.asList(10), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Too many insertions
                {101, Arrays.asList(0), Arrays.asList(0), Arrays.asList(11), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Sub-threshold deletions
                {101, Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(10),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Too many deletions
                {101, Arrays.asList(0), Arrays.asList(0), Arrays.asList(0), Arrays.asList(11),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Sub-threshold errors
                {101, Arrays.asList(3), Arrays.asList(0), Arrays.asList(3), Arrays.asList(4),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Too many errors
                {101, Arrays.asList(4), Arrays.asList(0), Arrays.asList(3), Arrays.asList(4),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Sub-threshold NM with large clip
                {101, Arrays.asList(9), Arrays.asList(10), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Critical NM with large clip
                {101, Arrays.asList(10), Arrays.asList(10), Arrays.asList(0), Arrays.asList(0),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Sub-threshold with large clip
                {104, Arrays.asList(3), Arrays.asList(10), Arrays.asList(3), Arrays.asList(3),
                        Arrays.asList("recordA"), Arrays.asList("1")},

                //Critical with large clip
                {101, Arrays.asList(8), Arrays.asList(10), Arrays.asList(3), Arrays.asList(4),
                        Arrays.asList("recordA"), Arrays.asList()},

                //Distinct hits
                {101, Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordB"), Arrays.asList("1", "2")},

                //Indistinct hits
                {101, Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0), Arrays.asList(0, 0),
                        Arrays.asList("recordA", "recordA"), Arrays.asList("1")},
        };
    }

    @Test(dataProvider = "mapUnpaired", groups = "spark")
    public void testMapGroupedReadsToTaxUnpaired(final int readLength, final List<Integer> NM, final List<Integer> clip,
                                                 final List<Integer> insert, final List<Integer> delete,
                                                 final List<String> contig, final List<String> truthTax) {

        if (!(NM.size() == clip.size() && NM.size() == insert.size() && NM.size() == delete.size() && NM.size() == contig.size())) {
            throw new TestException("Input lists for read must be of uniform length");
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final Broadcast<Map<String, String>> refNameToTaxBroadcast = ctx.broadcast(refNameToTax);

        //Test with alternate alignments assigned to the XA tag
        final List<Iterable<GATKRead>> readListXA = new ArrayList<>();
        readListXA.add(generateUnpairedRead(readLength, NM, clip, insert, delete, contig, "XA"));
        final JavaRDD<Iterable<GATKRead>> pairsXA = ctx.parallelize(readListXA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultXA = PSScoreUtils.mapGroupedReadsToTax(pairsXA, MIN_COV, MIN_IDENT, refNameToTaxBroadcast);
        final PSPathogenAlignmentHit infoXA = resultXA.first()._2;

        Assert.assertNotNull(infoXA);
        Assert.assertEquals(infoXA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoXA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoXA.numMates, 1);

        //Test SA tag
        final List<Iterable<GATKRead>> readListSA = new ArrayList<>();
        readListSA.add(generateUnpairedRead(readLength, NM, clip, insert, delete, contig, "SA"));
        final JavaRDD<Iterable<GATKRead>> pairsSA = ctx.parallelize(readListSA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultSA = PSScoreUtils.mapGroupedReadsToTax(pairsSA, MIN_COV, MIN_IDENT, refNameToTaxBroadcast);
        final PSPathogenAlignmentHit infoSA = resultSA.first()._2;

        Assert.assertNotNull(infoSA);
        Assert.assertEquals(infoSA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoSA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoSA.numMates, 1);
    }

    @Test
    public void testComputeTaxScores() {

        PSTree tree = new PSTree("1");
        List<PSPathogenAlignmentHit> readTaxHits = new ArrayList<>(1);
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("4"), 2));
        try {
            final Map<String, PSPathogenTaxonScore> result = PSScoreUtils.computeTaxScores(readTaxHits, tree);
            Assert.assertTrue(result.isEmpty(), "Result should be empty since the hit does not exist in the tree");
        } catch (Exception e) {
            Assert.fail("Threw and exception when a HitInfo references a tax ID not in the tree, or vice versa", e);
        }

        tree.addNode("2", "n2", "1", 0, "genus");
        tree.addNode("3", "n3", "2", 100, "species");
        readTaxHits.clear();
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("3"), 2));
        Map<String, PSPathogenTaxonScore> result = PSScoreUtils.computeTaxScores(readTaxHits, tree);
        Assert.assertEquals(result.size(), 3);
        Assert.assertEquals(result.get("1").score, result.get("2").score);
        Assert.assertEquals(result.get("2").score, result.get("3").score);
        Assert.assertEquals(result.get("1").score, PSScoreUtils.SCORE_NORMALIZATION_FACTOR * 2.0 / 100);
        Assert.assertEquals(result.get("1").scoreNormalized, 100.0);
        Assert.assertEquals(result.get("2").scoreNormalized, 100.0);
        Assert.assertEquals(result.get("3").scoreNormalized, 100.0);
        Assert.assertEquals(result.get("3").reads, 2);
        Assert.assertEquals(result.get("3").unambiguous, 2);
        Assert.assertEquals(result.get("3").refLength, 100);

        tree.addNode("4", "n4", "1", 0, "genus");
        tree.addNode("5", "n5", "2", 100, "species");
        tree.addNode("6", "n6", "4", 100, "species");
        tree.addNode("7", "n7", "4", 100, "species");
        readTaxHits.clear();
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("4"), 2)); //Invalid hit, ref length 0
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("3"), 2));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("3", "6"), 2));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("5"), 2));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("6"), 1));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList("8"), 2)); //Invalid hit, not in tree
        result = PSScoreUtils.computeTaxScores(readTaxHits, tree);

        final double score3 = PSScoreUtils.SCORE_NORMALIZATION_FACTOR * (0.5 * 2.0 + 2.0) / 100;
        final double score5 = PSScoreUtils.SCORE_NORMALIZATION_FACTOR * 2.0 / 100;
        final double score6 = PSScoreUtils.SCORE_NORMALIZATION_FACTOR * (0.5 * 2.0 + 1.0) / 100;
        final double score2 = score3 + score5;
        final double score4 = score6;
        final double score1 = score2 + score4;

        final int reads1 = 7;
        final int reads2 = 6;
        final int reads3 = 4;
        final int reads4 = 3;
        final int reads5 = 2;
        final int reads6 = 3;

        final int unambiguous1 = 7;
        final int unambiguous2 = 4;
        final int unambiguous3 = 2;
        final int unambiguous4 = 1;
        final int unambiguous5 = 2;
        final int unambiguous6 = 1;

        Assert.assertEquals(result.size(), 6);
        Assert.assertFalse(result.containsKey("7"));
        Assert.assertEquals(result.get("1").score, score1);
        Assert.assertEquals(result.get("2").score, score2);
        Assert.assertEquals(result.get("3").score, score3);
        Assert.assertEquals(result.get("4").score, score4);
        Assert.assertEquals(result.get("5").score, score5);
        Assert.assertEquals(result.get("6").score, score6);

        Assert.assertEquals(result.get("1").scoreNormalized, 100.0 * score1 / score1);
        Assert.assertEquals(result.get("2").scoreNormalized, 100.0 * score2 / score1);
        Assert.assertEquals(result.get("3").scoreNormalized, 100.0 * score3 / score1);
        Assert.assertEquals(result.get("4").scoreNormalized, 100.0 * score4 / score1);
        Assert.assertEquals(result.get("5").scoreNormalized, 100.0 * score5 / score1);
        Assert.assertEquals(result.get("6").scoreNormalized, 100.0 * score6 / score1);

        Assert.assertEquals(result.get("1").reads, reads1);
        Assert.assertEquals(result.get("2").reads, reads2);
        Assert.assertEquals(result.get("3").reads, reads3);
        Assert.assertEquals(result.get("4").reads, reads4);
        Assert.assertEquals(result.get("5").reads, reads5);
        Assert.assertEquals(result.get("6").reads, reads6);

        Assert.assertEquals(result.get("1").unambiguous, unambiguous1);
        Assert.assertEquals(result.get("2").unambiguous, unambiguous2);
        Assert.assertEquals(result.get("3").unambiguous, unambiguous3);
        Assert.assertEquals(result.get("4").unambiguous, unambiguous4);
        Assert.assertEquals(result.get("5").unambiguous, unambiguous5);
        Assert.assertEquals(result.get("6").unambiguous, unambiguous6);

        Assert.assertEquals(result.get("1").refLength, 0);
        Assert.assertEquals(result.get("3").refLength, 100);
        Assert.assertEquals(result.get("4").refLength, 0);

    }

    @DataProvider(name = "tupleSecond")
    public Object[][] getTupleSecondData() {
        return new Object[][]{
                {Arrays.asList()},
                {Arrays.asList(0)},
                {Arrays.asList(83, 19, 43, 0, -5, 73, 2383, 748, 3, 1, 1, 21)}
        };
    }

    @Test(dataProvider = "tupleSecond", groups = "spark")
    public void testCollectTupleSecond(final List<Integer> intsList) {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final int numInts = intsList.size();

        //Create RDD with tuples of (null,Integer)
        final List<Tuple2<Object, Integer>> data = new ArrayList<>(numInts);
        intsList.stream().forEach(val -> data.add(new Tuple2<>(null, val)));
        final JavaRDD<Tuple2<Object, Integer>> dataRDD = ctx.parallelize(data);

        //Test
        final Iterable<Integer> result = PSScoreUtils.collectTupleSecond(dataRDD);

        //Check that the result is equal to the original integer list, allowing for different order
        final List<Integer> resultList = StreamSupport.stream(result.spliterator(), false)
                .collect(Collectors.toList());
        Collections.sort(resultList);
        Collections.sort(intsList);
        Assert.assertEquals(resultList, intsList);
    }

    @DataProvider(name = "firstIterable")
    public Object[][] getFirstIterable() {
        return new Object[][]{
                {Arrays.asList()},
                {Arrays.asList(new Integer[]{})},
                {Arrays.asList(new int[]{8, 21, -5, 0, 11})},
                {Arrays.asList(new int[]{83, 19, 43, 0, -5, 73, 2383, 748, 3, 1, 1, 21}, new int[]{4, 7, 0, 1, -34, 100})}
        };
    }

    @Test(dataProvider = "firstIterable", groups = "spark")
    public void testFlattenTupleFirstIterable(final List<int[]> arrayList) {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final int numIterables = arrayList.size();

        //Create RDD with tuples of (Iterable<Integer>,null)
        final List<Tuple2<Iterable<Integer>, Object>> data = new ArrayList<>(numIterables);
        final List<List<Integer>> iterableList = arrayList.stream()
                .map(array -> IntStream.of(array)
                        .mapToObj(Integer::valueOf)
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        iterableList.stream().forEach(list -> data.add(new Tuple2<>(list, null)));
        final JavaRDD<Tuple2<Iterable<Integer>, Object>> dataRDD = ctx.parallelize(data);

        //Check that the result is equal to the original integer list, allowing for different order
        final List<Integer> resultList = new ArrayList<>(PSScoreUtils.flattenTupleFirstIterable(dataRDD).collect());
        final List<Integer> truthList = iterableList.stream()
                .flatMap(list -> list.stream())
                .collect(Collectors.toList());
        Collections.sort(resultList);
        Collections.sort(truthList);
        Assert.assertEquals(resultList, truthList);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testReadDatabase() throws Exception {
        final File tempFile = createTempFile("test", ".dat");
        final PSTree tree_in = new PSTree("1");
        tree_in.addNode("2", "n2", "1", 100, "rank_2");
        final Map<String, String> map_in = new HashMap<>();
        map_in.put("acc_2a", "2");
        map_in.put("acc_2b", "2");
        final PSTaxonomyDatabase expectedDatabase = new PSTaxonomyDatabase(tree_in, map_in);

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);
        final Output output = new Output(new FileOutputStream(tempFile.getPath()));
        kryo.writeObject(output, expectedDatabase);
        output.close();

        final PSTaxonomyDatabase taxonomyDatabase = PSScoreUtils.readTaxonomyDatabase(tempFile.getPath());
        Assert.assertEquals(expectedDatabase.tree, taxonomyDatabase.tree);
        Assert.assertEquals(expectedDatabase.accessionToTaxId, taxonomyDatabase.accessionToTaxId);

        try {
            PSScoreUtils.readTaxonomyDatabase("/bad_dir");
            Assert.fail("Did not throw UserException for bad path to readTaxonomyDatabase()");
        } catch (UserException e) {
        }
    }

}