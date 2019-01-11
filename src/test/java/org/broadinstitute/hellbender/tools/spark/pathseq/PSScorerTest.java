package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
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

public class PSScorerTest extends CommandLineProgramTest {

    private static final double MIN_IDENT = 0.80;
    private static final double IDENT_MARGIN = 0.05;
    private static final double SCORE_ABSOLUTE_ERROR_TOLERANCE = 1e-6;
    private PSTaxonomyDatabase taxonomyDatabase;
    private Map<String, Integer> refNameToTax;

    @BeforeMethod
    public void before() {
        refNameToTax = new HashMap<>(3);
        refNameToTax.put("recordA", 1);
        refNameToTax.put("recordB", 2);
        refNameToTax.put("recordC", 3);
        taxonomyDatabase = new PSTaxonomyDatabase(new PSTree(PSTaxonomyConstants.ROOT_ID), refNameToTax);
    }

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

    @Test(groups = "spark")
    public void testDoHeaderWarnings() {

        final File file = createTempFile("header_warnings", "txt");
        final PSTree tree = new PSTree(1);
        final Map<String, Integer> map = new HashMap<>(10);
        map.put("seq1", 2);
        map.put("seq2", 2);
        map.put("seq3", 3);
        final PSTaxonomyDatabase taxDB = new PSTaxonomyDatabase(tree, map);

        //Header and taxDB with matching records -> no warnings expected
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("seq1", 1000));
        header.addSequence(new SAMSequenceRecord("seq2", 1000));
        header.addSequence(new SAMSequenceRecord("seq3", 1000));
        PSScorer.writeMissingReferenceAccessions(file.getAbsolutePath(), header, taxDB, logger);
        try {
            final BufferedReader inputStream = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
            Assert.assertNull(inputStream.readLine());
        } catch (IOException e) {
            throw new TestException("Could not open temporary file", e);
        }

        //Add sequence not in taxDB -> expect warning to be written
        header.addSequence(new SAMSequenceRecord("seq4", 1000));
        PSScorer.writeMissingReferenceAccessions(file.getAbsolutePath(), header, taxDB, logger);
        try {
            final BufferedReader inputStream = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
            Assert.assertNotNull(inputStream.readLine());
        } catch (IOException e) {
            throw new TestException("Could not open temporary file", e);
        }
    }

    private static final List<Object> _empty_ = Collections.emptyList();
    private static final List<Integer> _0_ = Collections.singletonList(0);
    private static final List<Integer> _1_ = Collections.singletonList(1);
    private static final List<Integer> _2_ = Collections.singletonList(2);
    private static final List<Integer> _1_2_ = Arrays.asList(1, 2);
    private static final List<Integer> _20_ = Collections.singletonList(20);
    private static final List<Integer> _21_ = Collections.singletonList(21);
    private static final List<Integer> _0_0_ = Arrays.asList(0, 0);
    private static final List<String> _recordA_ = Collections.singletonList("recordA");
    private static final List<String> _recordA_recordA_ = Arrays.asList("recordA", "recordA");
    private static final List<String> _recordA_recordB_ = Arrays.asList("recordA", "recordB");
    @DataProvider(name = "mapPairs")
    public Object[][] getMapPairData() {
        return new Object[][]{

                //Perfect match
                {101, _0_, _0_, _0_, _0_,
                        _0_, _0_, _0_, _0_,
                        _recordA_, _recordA_, _1_},

                //First read too few mismatches
                {101, _20_, _0_, _0_, _0_,
                        _0_, _0_, _0_, _0_,
                        _recordA_, _recordA_, _1_},

                //First read too many mismatches
                {101, _21_, _0_, _0_, _0_,
                        _0_, _0_, _0_, _0_,
                        _recordA_, _recordA_, _empty_},

                //Second read too few mismatches
                {101, _0_, _20_, _0_, _0_,
                        _0_, _0_, _0_, _0_,
                        _recordA_, _recordA_, _1_},

                //Second read too many mismatches
                {101, _0_, _21_, _0_, _0_,
                        _0_, _0_, _0_, _0_,
                        _recordA_, _recordA_, _empty_},

                //Indistinct hits
                {101, _0_0_, _0_0_, _0_0_, _0_0_,
                        _0_0_, _0_0_, _0_0_, _0_0_,
                        _recordA_recordA_, _recordA_recordA_, _1_},

                //Chimeric hits
                {101, _0_0_, _0_0_, _0_0_, _0_0_,
                        _0_0_, _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _recordA_recordA_, _1_},

                //Dual hits
                {101, _0_0_, _0_0_, _0_0_, _0_0_,
                        _0_0_, _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _recordA_recordB_, _1_2_},

                //Sequence name truncation
                {101, _0_0_, _0_0_, _0_0_, _0_0_,
                        _0_0_, _0_0_, _0_0_, _0_0_,
                        Arrays.asList("recordA 1", "recordB 2"), Arrays.asList("recordA 1", "recordB 2"), _1_2_},

                //Second hit has too many mismatches
                {101, Arrays.asList(0, 21), _0_0_, _0_0_, _0_0_,
                        _0_0_, _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _recordA_recordB_, _1_},

                //Alternate alignment better than primary alignment
                {101, Arrays.asList(21, 0), _0_0_, _0_0_, _0_0_,
                        _0_0_, _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _recordA_recordB_, _2_}
        };
    }

    @Test(dataProvider = "mapPairs", groups = "spark")
    public void testMapGroupedReadsToTax(final int readLength, final List<Integer> NM1, final List<Integer> NM2,
                                         final List<Integer> clip1, final List<Integer> clip2,
                                         final List<Integer> insert1, final List<Integer> insert2,
                                         final List<Integer> delete1, final List<Integer> delete2,
                                         final List<String> contig1, final List<String> contig2,
                                         final List<Integer> truthTax) {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final Broadcast<PSTaxonomyDatabase> taxonomyDatabaseBroadcast = ctx.broadcast(taxonomyDatabase);

        //Test with alternate alignments assigned to the XA tag
        final List<Iterable<GATKRead>> readListXA = new ArrayList<>();
        readListXA.add(generateReadPair(readLength, NM1, NM2, clip1, clip2, insert1, insert2, delete1, delete2, contig1, contig2, "XA"));
        final JavaRDD<Iterable<GATKRead>> pairsXA = ctx.parallelize(readListXA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultXA = PSScorer.mapGroupedReadsToTax(pairsXA,
                MIN_IDENT, IDENT_MARGIN, taxonomyDatabaseBroadcast);
        final PSPathogenAlignmentHit infoXA = resultXA.first()._2;

        Assert.assertNotNull(infoXA);
        Assert.assertEquals(infoXA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoXA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoXA.numMates, 2);

        //Test SA tag
        final List<Iterable<GATKRead>> readListSA = new ArrayList<>();
        readListSA.add(generateReadPair(readLength, NM1, NM2, clip1, clip2, insert1, insert2, delete1, delete2, contig1, contig2, "SA"));
        final JavaRDD<Iterable<GATKRead>> pairsSA = ctx.parallelize(readListSA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultSA = PSScorer.mapGroupedReadsToTax(pairsSA,
                MIN_IDENT, IDENT_MARGIN, taxonomyDatabaseBroadcast);
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
               {100, _0_, _0_, _0_, _0_,
                        _recordA_, _1_},

                //Suprathreshold mismatches
                {100, _20_, _0_, _0_, _0_,
                        _recordA_, _1_},

                //Subthreshold mismatches
                {100, _21_, _0_, _0_, _0_,
                        _recordA_, _empty_},

                //Suprathreshold clip
                {100, _0_, _20_, _0_, _0_,
                        _recordA_, _1_},

                //Subthreshold clip + mismatches
                {100, Arrays.asList(1), _20_, _0_, _0_,
                        _recordA_, _empty_},

                //Suprathreshold insertions
                {100, _20_, _0_, _20_, _0_,
                        _recordA_, _1_},

                //Subthreshold insertions + mismatches
                {100, Arrays.asList(21), _0_, _20_, _0_,
                        _recordA_, _empty_},

                //Suprathreshold deletions
                {100, _20_, _0_, _0_, _20_,
                        _recordA_, _1_},

                //Subthreshold deletions
                {100, _21_, _0_, _0_, _21_,
                        _recordA_, _empty_},

                //Two taxa
                {100, _0_0_, _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _1_2_},

                //Two taxa, one with mismatches within the margin
                {100, Arrays.asList(5, 0), _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _1_2_},

                //Two taxa, one with mismatches below the margin
                {100, Arrays.asList(6, 0), _0_0_, _0_0_, _0_0_,
                        _recordA_recordB_, _2_},
        };
    }

    @Test(dataProvider = "mapUnpaired", groups = "spark")
    public void testMapGroupedReadsToTaxUnpaired(final int readLength, final List<Integer> NM, final List<Integer> clip,
                                                 final List<Integer> insert, final List<Integer> delete,
                                                 final List<String> contig, final List<Integer> truthTax) {

        if (!(NM.size() == clip.size() && NM.size() == insert.size() && NM.size() == delete.size() && NM.size() == contig.size())) {
            throw new TestException("Input lists for read must be of uniform length");
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final Broadcast<PSTaxonomyDatabase> taxonomyDatabaseBroadcast = ctx.broadcast(taxonomyDatabase);

        //Test with alternate alignments assigned to the XA tag
        final List<Iterable<GATKRead>> readListXA = new ArrayList<>();
        readListXA.add(generateUnpairedRead(readLength, NM, clip, insert, delete, contig, "XA"));
        final JavaRDD<Iterable<GATKRead>> pairsXA = ctx.parallelize(readListXA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultXA = PSScorer.mapGroupedReadsToTax(pairsXA,
                MIN_IDENT, IDENT_MARGIN, taxonomyDatabaseBroadcast);
        final PSPathogenAlignmentHit infoXA = resultXA.first()._2;

        Assert.assertNotNull(infoXA);
        Assert.assertEquals(infoXA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoXA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoXA.numMates, 1);

        //Test SA tag
        final List<Iterable<GATKRead>> readListSA = new ArrayList<>();
        readListSA.add(generateUnpairedRead(readLength, NM, clip, insert, delete, contig, "SA"));
        final JavaRDD<Iterable<GATKRead>> pairsSA = ctx.parallelize(readListSA);
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> resultSA = PSScorer.mapGroupedReadsToTax(pairsSA,
                MIN_IDENT, IDENT_MARGIN, taxonomyDatabaseBroadcast);
        final PSPathogenAlignmentHit infoSA = resultSA.first()._2;

        Assert.assertNotNull(infoSA);
        Assert.assertEquals(infoSA.taxIDs.size(), truthTax.size());
        Assert.assertTrue(infoSA.taxIDs.containsAll(truthTax));
        Assert.assertEquals(infoSA.numMates, 1);
    }

    private static Map<Integer,PSPathogenTaxonScore> scoreIteratorToMap(final Iterator<Tuple2<Integer, PSPathogenTaxonScore>> iter) {
        final Map<Integer,PSPathogenTaxonScore> map = new HashMap<>();
        iter.forEachRemaining(pair -> map.put(pair._1, pair._2));
        return map;
    }

    @Test
    public void testComputeTaxScores() {

        PSTree tree = new PSTree(1);
        List<PSPathogenAlignmentHit> readTaxHits = new ArrayList<>(1);
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(4), 2));
        boolean divideByGenomeLength = true;
        boolean notNormalizedByKingdom = true;
        final PSTaxonomyDatabase testDatabase = new PSTaxonomyDatabase(tree, null);
        try {
            final Iterator<Tuple2<Integer, PSPathogenTaxonScore>> resultIter = PSScorer.computeTaxScores(readTaxHits.iterator(), testDatabase, divideByGenomeLength);
            Map<Integer,PSPathogenTaxonScore> resultMap = PSScorer.computeNormalizedScores(scoreIteratorToMap(resultIter), testDatabase.tree, notNormalizedByKingdom);
            Assert.assertTrue(resultMap.isEmpty(), "Result should be empty since the hit does not exist in the tree");
        } catch (Exception e) {
            Assert.fail("Threw an exception when a HitInfo references a tax ID not in the tree, or vice versa", e);
        }

        tree.addNode(2, "n2", 1, 0, PSTaxonomyConstants.KINGDOM_RANK_NAME);
        tree.addNode(3, "n3", 2, 100, "species");
        readTaxHits.clear();
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(3), 2));
        Iterator<Tuple2<Integer, PSPathogenTaxonScore>> resultIter = PSScorer.computeTaxScores(readTaxHits.iterator(), testDatabase, divideByGenomeLength);
        Map<Integer,PSPathogenTaxonScore> resultMap = PSScorer.computeNormalizedScores(scoreIteratorToMap(resultIter), testDatabase.tree, notNormalizedByKingdom);
        Assert.assertEquals(resultMap.size(), 3);
        Assert.assertEquals(resultMap.get(1).getSelfScore(), resultMap.get(2).getSelfScore());
        Assert.assertEquals(resultMap.get(1).getDescendentScore(), resultMap.get(2).getDescendentScore());
        Assert.assertEquals(resultMap.get(2).getDescendentScore(), resultMap.get(3).getSelfScore());
        Assert.assertEquals(resultMap.get(3).getSelfScore(), PSScorer.SCORE_GENOME_LENGTH_UNITS * 2.0 / 100);
        Assert.assertEquals(resultMap.get(1).getScoreNormalized(), 100.0);
        Assert.assertEquals(resultMap.get(2).getScoreNormalized(), 100.0);
        Assert.assertEquals(resultMap.get(3).getScoreNormalized(), 100.0);
        Assert.assertEquals(resultMap.get(3).getTotalReads(), 2);
        Assert.assertEquals(resultMap.get(3).getUnambiguousReads(), 2);
        Assert.assertEquals(resultMap.get(3).getReferenceLength(), 100);
        Assert.assertEquals(resultMap.get(3).getKingdomTaxonId(), 1);

        tree.addNode(4, "n4", 1, 0, PSTaxonomyConstants.SUPERKINGDOM_RANK_NAME);
        tree.addNode(5, "n5", 2, 100, "species");
        tree.addNode(6, "n6", 4, 100, "species");
        tree.addNode(7, "n7", 4, 100, "species");
        readTaxHits.clear();
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(4), 2)); //Invalid hit, ref length 0
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(3), 2));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(3, 6), 2));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(5), 2));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(6), 1));
        readTaxHits.add(new PSPathogenAlignmentHit(Arrays.asList(8), 2)); //Invalid hit, not in tree
        resultIter = PSScorer.computeTaxScores(readTaxHits.iterator(), testDatabase, divideByGenomeLength);
        resultMap = PSScorer.computeNormalizedScores(scoreIteratorToMap(resultIter), tree, notNormalizedByKingdom);
        checkComputedScores(resultMap, divideByGenomeLength, notNormalizedByKingdom);

        //Test after switching genome length and kingdom normalization
        divideByGenomeLength = false;
        notNormalizedByKingdom = false;
        resultIter = PSScorer.computeTaxScores(readTaxHits.iterator(), testDatabase, divideByGenomeLength);
        resultMap = PSScorer.computeNormalizedScores(scoreIteratorToMap(resultIter), tree, notNormalizedByKingdom);
        checkComputedScores(resultMap, divideByGenomeLength, notNormalizedByKingdom);
    }

    private static void checkComputedScores(final Map<Integer,PSPathogenTaxonScore> resultMap, final boolean divideByGenomeLength,
                                           final boolean notNormalizeByKingdom) {
        double score3 = 0.5 * 2.0 + 2.0;
        double score5 = 2.0;
        double score6 = 0.5 * 2.0 + 1.0;
        if (divideByGenomeLength) {
            score3 *= PSScorer.SCORE_GENOME_LENGTH_UNITS / 100.;
            score5 *= PSScorer.SCORE_GENOME_LENGTH_UNITS / 100.;
            score6 *= PSScorer.SCORE_GENOME_LENGTH_UNITS / 100.;
        }
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

        Assert.assertEquals(resultMap.size(), 6);
        Assert.assertFalse(resultMap.containsKey(7));
        Assert.assertEquals(resultMap.get(1).getSelfScore() + resultMap.get(1).getDescendentScore(), score1);
        Assert.assertEquals(resultMap.get(2).getSelfScore() + resultMap.get(2).getDescendentScore(), score2);
        Assert.assertEquals(resultMap.get(3).getSelfScore() + resultMap.get(3).getDescendentScore(), score3);
        Assert.assertEquals(resultMap.get(4).getSelfScore() + resultMap.get(4).getDescendentScore(), score4);
        Assert.assertEquals(resultMap.get(5).getSelfScore() + resultMap.get(5).getDescendentScore(), score5);
        Assert.assertEquals(resultMap.get(6).getSelfScore() + resultMap.get(6).getDescendentScore(), score6);

        final double denominator2 = notNormalizeByKingdom ? score1 : score2;
        final double denominator4 = notNormalizeByKingdom ? score1 : score4;
        final int numKingdoms = notNormalizeByKingdom ? 1: 2;
        Assert.assertTrue(PathSeqTestUtils.equalWithinTolerance(resultMap.get(1).getScoreNormalized(), 100.0 * numKingdoms, SCORE_ABSOLUTE_ERROR_TOLERANCE));
        Assert.assertTrue(PathSeqTestUtils.equalWithinTolerance(resultMap.get(2).getScoreNormalized(), 100.0 * score2 / denominator2, SCORE_ABSOLUTE_ERROR_TOLERANCE));
        Assert.assertTrue(PathSeqTestUtils.equalWithinTolerance(resultMap.get(3).getScoreNormalized(), 100.0 * score3 / denominator2, SCORE_ABSOLUTE_ERROR_TOLERANCE));
        Assert.assertTrue(PathSeqTestUtils.equalWithinTolerance(resultMap.get(4).getScoreNormalized(), 100.0 * score4 / denominator4, SCORE_ABSOLUTE_ERROR_TOLERANCE));
        Assert.assertTrue(PathSeqTestUtils.equalWithinTolerance(resultMap.get(5).getScoreNormalized(), 100.0 * score5 / denominator2, SCORE_ABSOLUTE_ERROR_TOLERANCE));
        Assert.assertTrue(PathSeqTestUtils.equalWithinTolerance(resultMap.get(6).getScoreNormalized(), 100.0 * score6 / denominator4, SCORE_ABSOLUTE_ERROR_TOLERANCE));

        Assert.assertEquals(resultMap.get(1).getTotalReads(), reads1);
        Assert.assertEquals(resultMap.get(2).getTotalReads(), reads2);
        Assert.assertEquals(resultMap.get(3).getTotalReads(), reads3);
        Assert.assertEquals(resultMap.get(4).getTotalReads(), reads4);
        Assert.assertEquals(resultMap.get(5).getTotalReads(), reads5);
        Assert.assertEquals(resultMap.get(6).getTotalReads(), reads6);

        Assert.assertEquals(resultMap.get(1).getUnambiguousReads(), unambiguous1);
        Assert.assertEquals(resultMap.get(2).getUnambiguousReads(), unambiguous2);
        Assert.assertEquals(resultMap.get(3).getUnambiguousReads(), unambiguous3);
        Assert.assertEquals(resultMap.get(4).getUnambiguousReads(), unambiguous4);
        Assert.assertEquals(resultMap.get(5).getUnambiguousReads(), unambiguous5);
        Assert.assertEquals(resultMap.get(6).getUnambiguousReads(), unambiguous6);

        Assert.assertEquals(resultMap.get(1).getReferenceLength(), 0);
        Assert.assertEquals(resultMap.get(3).getReferenceLength(), 100);
        Assert.assertEquals(resultMap.get(4).getReferenceLength(), 0);

        if (!notNormalizeByKingdom) {
            Assert.assertEquals(resultMap.get(1).getKingdomTaxonId(), 1);
            Assert.assertEquals(resultMap.get(2).getKingdomTaxonId(), 2);
            Assert.assertEquals(resultMap.get(3).getKingdomTaxonId(), 2);
            Assert.assertEquals(resultMap.get(4).getKingdomTaxonId(), 4);
            Assert.assertEquals(resultMap.get(5).getKingdomTaxonId(), 2);
            Assert.assertEquals(resultMap.get(6).getKingdomTaxonId(), 4);
        } else {
            Assert.assertEquals(resultMap.get(1).getKingdomTaxonId(), 1);
            Assert.assertEquals(resultMap.get(2).getKingdomTaxonId(), 1);
            Assert.assertEquals(resultMap.get(3).getKingdomTaxonId(), 1);
            Assert.assertEquals(resultMap.get(4).getKingdomTaxonId(), 1);
            Assert.assertEquals(resultMap.get(5).getKingdomTaxonId(), 1);
            Assert.assertEquals(resultMap.get(6).getKingdomTaxonId(), 1);
        }
    }

    @DataProvider(name = "tupleSecond")
    public Object[][] getTupleSecondData() {
        return new Object[][]{
                {_empty_},
                {_0_},
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
        final Iterable<Integer> result = PSScorer.collectValues(dataRDD);

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
                {_empty_},
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
        final List<Integer> resultList = new ArrayList<>(PSScorer.flattenIterableKeys(dataRDD).collect());
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
        final PSTree tree_in = new PSTree(1);
        tree_in.addNode(2, "n2", 1, 100, "rank_2");
        final Map<String, Integer> map_in = new HashMap<>();
        map_in.put("acc_2a", 2);
        map_in.put("acc_2b", 2);
        final PSTaxonomyDatabase expectedDatabase = new PSTaxonomyDatabase(tree_in, map_in);

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);
        final Output output = new Output(new FileOutputStream(tempFile.getPath()));
        kryo.writeObject(output, expectedDatabase);
        output.close();

        final PSTaxonomyDatabase taxonomyDatabase = PSScorer.readTaxonomyDatabase(tempFile.getPath());
        Assert.assertEquals(expectedDatabase.tree, taxonomyDatabase.tree);
        Assert.assertEquals(expectedDatabase.accessionToTaxId, taxonomyDatabase.accessionToTaxId);

        try {
            PSScorer.readTaxonomyDatabase("/bad_dir");
            Assert.fail("Did not throw UserException for bad path to readTaxonomyDatabase()");
        } catch (UserException e) {
        }
    }

}