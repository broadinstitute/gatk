package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Random;

public class PathSeqBuildKmersSparkIntegrationTest extends CommandLineProgramTest {

    private final static double BLOOM_FPP = 0.01;
    private final static int NUM_FPP_TRIALS = 10000;

    @Override
    public String getTestedClassName() {
        return PathSeqBuildKmers.class.getSimpleName();
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testHopscotchSetFromFasta() throws Exception {

        final String libraryPath = publicTestDir + PathSeqBuildKmers.class.getPackage().getName().replace(".", "/") + "/hg19mini.hss";
        final File expectedFile = new File(libraryPath);
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(PathSeqBuildKmers.REFERENCE_LONG_NAME, ref);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input inputExpected = new Input(FileUtils.openInputStream(expectedFile));
        final Input inputTest = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        final PSKmerSet expectedKmerLib = kryo.readObject(inputExpected, PSKmerSet.class);
        final PSKmerSet testKmerLib = kryo.readObject(inputTest, PSKmerSet.class);

        Assert.assertEquals(testKmerLib, expectedKmerLib);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testBloomFilterFromFasta() throws Exception {

        final String libraryPath = publicTestDir + PathSeqBuildKmers.class.getPackage().getName().replace(".", "/") + "/hg19mini.hss";
        final File expectedFile = new File(libraryPath);
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".bfi");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(PathSeqBuildKmers.REFERENCE_LONG_NAME, ref);
        args.addArgument(PathSeqBuildKmers.BLOOM_FILTER_FALSE_POSITIVE_P_LONG_NAME, Double.toString(BLOOM_FPP));
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input inputExpected = new Input(FileUtils.openInputStream(expectedFile));
        final Input inputTest = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        final PSKmerSet expectedKmerLib = kryo.readObject(inputExpected, PSKmerSet.class);
        final PSKmerBloomFilter testKmerLib = kryo.readObject(inputTest, PSKmerBloomFilter.class);

        final LongIterator itr = expectedKmerLib.iterator();
        while (itr.hasNext()) {
            Assert.assertTrue(testKmerLib.contains(new SVKmerShort(itr.next())));
        }

        final Random rand = new Random(72939);
        int numFP = 0;
        for (int i = 0; i < NUM_FPP_TRIALS; i++) {
            final long randomValue = rand.nextLong() >>> 2;
            if (testKmerLib.contains(new SVKmerShort(randomValue)) && !expectedKmerLib.contains(new SVKmerShort(randomValue))) {
                numFP++;
            }
        }
        Assert.assertTrue(numFP < 1.2 * NUM_FPP_TRIALS * BLOOM_FPP);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testMaskedHopscotchSetFromFasta() throws Exception {
        final File expectedFile = getTestFile("hg19mini.mask_4_15.hss");
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(PathSeqBuildKmers.REFERENCE_LONG_NAME, ref);
        args.addArgument(PathSeqBuildKmers.KMER_MASK_LONG_NAME, "4,15");
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input inputExpected = new Input(FileUtils.openInputStream(expectedFile));
        final Input inputTest = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        final PSKmerSet expectedKmerLib = kryo.readObject(inputExpected, PSKmerSet.class);
        final PSKmerSet testKmerLib = kryo.readObject(inputTest, PSKmerSet.class);

        Assert.assertEquals(testKmerLib, expectedKmerLib);
    }

    @DataProvider(name = "badArgs")
    public Object[][] getBadArguments() {
        return new Object[][]{
                {PathSeqBuildKmers.KMER_SIZE_LONG_NAME, "0"},
                {PathSeqBuildKmers.KMER_SIZE_LONG_NAME, "2"},
                {PathSeqBuildKmers.KMER_SIZE_LONG_NAME, "32"},
                {PathSeqBuildKmers.KMER_SIZE_LONG_NAME, "33"},
                {PathSeqBuildKmers.BLOOM_FILTER_FALSE_POSITIVE_P_LONG_NAME, "-0.1"},
                {PathSeqBuildKmers.BLOOM_FILTER_FALSE_POSITIVE_P_LONG_NAME, "1"},
                {PathSeqBuildKmers.KMER_MASK_LONG_NAME, "0,32"},
                {PathSeqBuildKmers.KMER_MASK_LONG_NAME, "-1,15"},
                {PathSeqBuildKmers.KMER_SPACING_LONG_NAME, "0"}
        };
    }

    @Test(dataProvider = "badArgs", expectedExceptions = Exception.class)
    public void testBadInputs(final String argName, final String argVal) throws Exception {
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(PathSeqBuildKmers.REFERENCE_LONG_NAME, ref);
        args.addArgument(argName, argVal);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());
    }

}
