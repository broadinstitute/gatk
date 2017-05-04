package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

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
    public void testWriteTwoKryo() throws Exception {
        final File tempFile = createTempFile("test", ".dat");
        final Integer int_in = 29382;
        final String str_in = "test string";
        PSUtils.writeKryoTwo(tempFile.getPath(), int_in, str_in);

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);
        final Input input = new Input(BucketUtils.openFile(tempFile.getPath(), null));
        final Integer int_out = (Integer) kryo.readClassAndObject(input);
        final String str_out = (String) kryo.readClassAndObject(input);
        input.close();

        Assert.assertEquals(int_in, int_out);
        Assert.assertEquals(str_in, str_out);

        try {
            PSUtils.writeKryoTwo("/bad_dir", int_out, str_out);
            Assert.fail("Did not throw UserException for bad path to writeKryoTwo()");
        } catch (UserException e) {
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