package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.pathseq.PSUtils.addReferenceSequencesToHeader;

public class PSUtilsTest extends BaseTest {

    @Override
    public String getTestedClassName() {
        return PSUtils.class.getSimpleName();
    }

    @DataProvider(name = "maskData")
    public Object[][] getAlignmentData() {
        return new Object[][]{
                {"", 0, new byte[0], Boolean.FALSE},
                {"0", 31, null, Boolean.TRUE},
                {"0,1", 31, new byte[]{0, 1}, Boolean.FALSE},
                {"0,1,", 31, new byte[]{0, 1}, Boolean.FALSE},
                {"0,1,10,8", 31, new byte[]{0, 1, 10, 8}, Boolean.FALSE},
                {"0,-1", 31, null, Boolean.TRUE},
                {"0,31", 31, null, Boolean.TRUE}
        };
    }

    @Test(dataProvider = "maskData")
    public void testParseMask(final String maskArg, final int kSize, final byte[] expected, final Boolean throwsException) {
        if (throwsException) {
            try {
                PSUtils.parseMask(maskArg, kSize);
                Assert.fail("Did not throw error for mask " + maskArg + " and kSize " + kSize);
            } catch (Exception e) {
            }
        } else {
            byte[] result = PSUtils.parseMask(maskArg, kSize);
            Assert.assertNotNull(result, "Parser should return empty array instead of null value");
            Assert.assertEquals(result, expected);
        }
    }

    @Test(groups = "spark")
    public void testPrimaryReads() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<GATKRead> readList = new ArrayList<>(2);

        final GATKRead primaryRead = ArtificialReadUtils.createRandomRead(101);
        readList.add(primaryRead);

        final GATKRead secondaryRead = ArtificialReadUtils.createRandomRead(101);
        secondaryRead.setIsSecondaryAlignment(true);
        readList.add(secondaryRead);

        final GATKRead supplementaryRead = ArtificialReadUtils.createRandomRead(101);
        supplementaryRead.setIsSupplementaryAlignment(true);
        readList.add(supplementaryRead);

        final JavaRDD<GATKRead> readRDD = ctx.parallelize(readList);
        final List<GATKRead> result = PSUtils.primaryReads(readRDD).collect();

        Assert.assertTrue(result.contains(primaryRead));
        Assert.assertTrue(result.size() == 1);

    }

    @Test(groups = "spark")
    public void testGroupReadsByName() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<GATKRead> readList = new ArrayList<>(2);

        final SAMSequenceDictionary seq = new SAMSequenceDictionary();
        seq.addSequence(new SAMSequenceRecord("test_seq", 1000));
        final SAMFileHeader header = new SAMFileHeader(seq);

        final List<GATKRead> readPair1 = ArtificialReadUtils.createPair(header, "paired_1", 101, 1, 102, false, false);
        readList.addAll(readPair1);

        final List<GATKRead> readPair2 = ArtificialReadUtils.createPair(header, "paired_2", 101, 1, 102, false, false);
        readList.addAll(readPair2);

        final GATKRead unpairedRead = ArtificialReadUtils.createRandomRead(101);
        unpairedRead.setName("unpaired_1");
        readList.add(unpairedRead);
        final List<GATKRead> unpairedReadIterable = new ArrayList<>(1);
        unpairedReadIterable.add(unpairedRead);

        final JavaRDD<GATKRead> readRDD = ctx.parallelize(readList);
        final List<Iterable<GATKRead>> result = PSUtils.groupReadPairs(readRDD).collect();

        Assert.assertTrue(result.size() == 3);
        for (final Iterable<GATKRead> group : result) {
            String groupName = null;
            for (final GATKRead read : group) {
                final String readName = read.getName();
                Assert.assertTrue(readName.equals("paired_1") || readName.equals("paired_2") || readName.equals("unpaired_1"),
                        "Unrecognized read name: " + readName);
                if (groupName == null) {
                    groupName = readName;
                } else {
                    Assert.assertEquals(readName, groupName);
                }
            }
            Assert.assertNotNull(groupName, "Returned an empty group");
        }

    }

    @Test(groups = "spark")
    public void testReadPairing() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<GATKRead> readList = new ArrayList<>(2);

        final SAMSequenceDictionary seq = new SAMSequenceDictionary();
        seq.addSequence(new SAMSequenceRecord("test_seq", 1000));
        final SAMFileHeader header = new SAMFileHeader(seq);

        final List<GATKRead> readPair1 = ArtificialReadUtils.createPair(header, "paired_1", 101, 1, 102, false, false);
        readList.addAll(readPair1);

        final List<GATKRead> readPair2 = ArtificialReadUtils.createPair(header, "paired_2", 101, 1, 102, false, false);
        readList.addAll(readPair2);

        final GATKRead unpairedRead = ArtificialReadUtils.createRandomRead(101);
        unpairedRead.setName("unpaired_1");
        readList.add(unpairedRead);

        final List<GATKRead> readPair3 = ArtificialReadUtils.createPair(header, "paired_3", 101, 1, 102, false, false);
        final GATKRead formerlyPairedRead = readPair3.get(0);
        readList.add(formerlyPairedRead);

        final JavaRDD<GATKRead> readRDD = ctx.parallelize(readList);

        final List<GATKRead> resultPaired = PSUtils.pairedReads(readRDD).collect();
        Assert.assertEquals(resultPaired.size(), 4);
        Assert.assertTrue(resultPaired.contains(readPair1.get(0)));
        Assert.assertTrue(resultPaired.contains(readPair1.get(1)));
        Assert.assertTrue(resultPaired.contains(readPair2.get(0)));
        Assert.assertTrue(resultPaired.contains(readPair2.get(1)));

        final List<GATKRead> resultUnpaired = PSUtils.unpairedReads(readRDD).collect();
        Assert.assertEquals(resultUnpaired.size(), 2);
        Assert.assertTrue(resultUnpaired.contains(unpairedRead));
        Assert.assertTrue(resultUnpaired.contains(formerlyPairedRead));
    }

    @Test(groups = "spark")
    public void testUnpairedReads() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<GATKRead> readList = new ArrayList<>(2);

        final SAMSequenceDictionary seq = new SAMSequenceDictionary();
        seq.addSequence(new SAMSequenceRecord("test_seq", 1000));
        final SAMFileHeader header = new SAMFileHeader(seq);

        final List<GATKRead> readPair1 = ArtificialReadUtils.createPair(header, "paired_1", 101, 1, 102, false, false);
        readList.addAll(readPair1);

        final List<GATKRead> readPair2 = ArtificialReadUtils.createPair(header, "paired_2", 101, 1, 102, false, false);
        readList.addAll(readPair2);

        final GATKRead unpairedRead = ArtificialReadUtils.createRandomRead(101);
        unpairedRead.setName("unpaired_1");
        readList.add(unpairedRead);

        final List<GATKRead> readPair3 = ArtificialReadUtils.createPair(header, "paired_3", 101, 1, 102, false, false);
        final GATKRead formerlyPairedRead = readPair3.get(0);
        readList.add(formerlyPairedRead);

        final JavaRDD<GATKRead> readRDD = ctx.parallelize(readList);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testWriteReadDatabase() throws Exception {
        final File tempFile = createTempFile("test", ".dat");
        PSTree tree_in = new PSTree("1");
        tree_in.addNode("2", "n2", "1", 100, "rank_2");
        Map<String, String> map_in = new HashMap<>();
        map_in.put("acc_2a", "2");
        map_in.put("acc_2b", "2");
        PSUtils.writeKryoTwo(tempFile.getPath(), tree_in, map_in);

        PSTaxonomyDatabase db = PSUtils.readTaxonomyDatabase(tempFile.getPath(), null);
        PSTree tree_out = db.tree;
        Map<String, String> map_out = db.contigToTaxIDMap;
        Assert.assertEquals(tree_in, tree_out);
        Assert.assertEquals(map_in, map_out);

        try {
            PSUtils.writeKryoTwo("/bad_dir", tree_in, map_in);
            Assert.fail("Did not throw UserException for bad path to writeKryoTwo()");
        } catch (UserException e) {
        }

        try {
            PSUtils.readTaxonomyDatabase("/bad_dir", null);
            Assert.fail("Did not throw UserException for bad path to readTaxonomyDatabase()");
        } catch (UserException e) {
        }
    }

    @Test
    public void testAddReferenceSequencesToHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        final SerializableFunction<GATKRead, SimpleInterval> windowFunction = ReferenceWindowFunctions.IDENTITY_FUNCTION;
        final String referencePath = hg19MiniReference;
        addReferenceSequencesToHeader(header, referencePath, windowFunction, null);

        final ReferenceMultiSource ref = new ReferenceMultiSource((PipelineOptions)null, referencePath, windowFunction);
        Assert.assertEquals(ref.getReferenceSequenceDictionary(null).size(), header.getSequenceDictionary().size());
        for (final SAMSequenceRecord rec : ref.getReferenceSequenceDictionary(null).getSequences()) {
            final SAMSequenceRecord recTest = header.getSequenceDictionary().getSequence(rec.getSequenceName());
            Assert.assertNotNull(recTest);
            Assert.assertEquals(rec, recTest);
        }
    }

    @Test
    public void testLogItemizedWarning() {
        final Collection<String> items = new ArrayList<>(3);
        PSUtils.logItemizedWarning(logger, items, "Test warning statement");
        items.add("x");
        items.add("y");
        items.add("z");
        PSUtils.logItemizedWarning(logger, items, "Test warning statement");
    }

}