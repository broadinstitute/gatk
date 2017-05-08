package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongBloomFilter;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Random;
import java.util.function.BiPredicate;

public class PathSeqKmerSparkIntegrationTest extends CommandLineProgramTest {

    private final static double TEST_DOWNSAMPLE_PROB = 0.5;
    private final static double BLOOM_FPP = 0.01;
    private final static int NUM_FPP_TRIALS = 1000;

    private static boolean checkDownsampleSet(final LargeLongHopscotchSet resultSet, final LargeLongHopscotchSet originalSet) {
        return resultSet.size() < originalSet.size() * TEST_DOWNSAMPLE_PROB * 1.5
                && resultSet.size() > originalSet.size() * TEST_DOWNSAMPLE_PROB * 0.5;
    }

    @Override
    public String getTestedClassName() {
        return PathSeqKmerSpark.class.getSimpleName();
    }

    @DataProvider(name = "argChecks")
    public Object[][] getTestArguments() {
        return new Object[][]{

                {null, null,
                        (BiPredicate<LargeLongHopscotchSet, LargeLongHopscotchSet>) (result, original) -> result.equals(original)},

                {"downsampleProbability", Double.toString(TEST_DOWNSAMPLE_PROB),
                        (BiPredicate<LargeLongHopscotchSet, LargeLongHopscotchSet>) (result, original) -> checkDownsampleSet(result, original)},

        };
    }

    @SuppressWarnings("unchecked")
    @Test(dataProvider = "argChecks", groups = "spark")
    public void testHopscotchSetFromFasta(final String argName, final String argVal,
                                          final BiPredicate<LargeLongHopscotchSet, LargeLongHopscotchSet> testFunction) throws Exception {
        final File expectedFile = getTestFile("kmer.hss");
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("reference", ref);
        if (argName != null && argVal != null) {
            args.addArgument(argName, argVal);
        }
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input input_expected = new Input(FileUtils.openInputStream(expectedFile));
        final Input input_test = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);

        final LargeLongHopscotchSet expectedKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected);
        final LargeLongHopscotchSet testKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_test);

        Assert.assertTrue(testFunction.test(testKmerLib, expectedKmerLib));
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testBloomFilterFromFasta() throws Exception {
        final File expectedFile = getTestFile("kmer.hss");
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".bfi");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("reference", ref);
        args.addArgument("bloomFalsePositiveProbability", new Double(BLOOM_FPP).toString());
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input input_expected = new Input(FileUtils.openInputStream(expectedFile));
        final Input input_test = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);

        final LargeLongHopscotchSet expectedKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected);
        final LargeLongBloomFilter testKmerLib = (LargeLongBloomFilter) kryo.readClassAndObject(input_test);

        final LongIterator itr = expectedKmerLib.iterator();
        while (itr.hasNext()) {
            Assert.assertTrue(testKmerLib.contains(itr.next()));
        }

        final Random rand = new Random(72939);
        int numFP = 0;
        for (int i = 0; i < 100; i++) {
            if (testKmerLib.contains(rand.nextLong() >>> 1)) {
                numFP++;
            }
        }
        Assert.assertTrue(numFP < 10 * NUM_FPP_TRIALS * BLOOM_FPP);
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testMaskedHopscotchSetFromFasta() throws Exception {
        final File expectedFile = getTestFile("kmer.mask_0_15.hss");
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("reference", ref);
        args.addOutput(output);
        args.addArgument("kmerMask", "0,15");
        this.runCommandLine(args.getArgsArray());

        final Input input_expected = new Input(FileUtils.openInputStream(expectedFile));
        final Input input_test = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);

        final LargeLongHopscotchSet expectedKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected);
        final LargeLongHopscotchSet testKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_test);

        Assert.assertTrue(testKmerLib.equals(expectedKmerLib));
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testBloomFilterFromHopscotchSet() throws Exception {
        final File expectedFile = getTestFile("kmer.hss");
        final File output = createTempFile("test", ".bfi");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("kmerFiles", expectedFile);
        args.addArgument("bloomFalsePositiveProbability", new Double(BLOOM_FPP).toString());
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input input_expected = new Input(FileUtils.openInputStream(expectedFile));
        final Input input_test = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);

        final LargeLongHopscotchSet expectedKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected);
        final LargeLongBloomFilter testKmerLib = (LargeLongBloomFilter) kryo.readClassAndObject(input_test);

        final LongIterator itr = expectedKmerLib.iterator();
        while (itr.hasNext()) {
            Assert.assertTrue(testKmerLib.contains(itr.next()));
        }

        final Random rand = new Random(72939);
        int numFP = 0;
        for (int i = 0; i < 100; i++) {
            if (testKmerLib.contains(rand.nextLong() >>> 1)) {
                numFP++;
            }
        }
        Assert.assertTrue(numFP < 10 * NUM_FPP_TRIALS * BLOOM_FPP);
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testHopscotchSetFromMultipleSets() throws Exception {
        final File expectedFile1 = getTestFile("kmer.hss");
        final File expectedFile2 = getTestFile("exampleFASTA.hss");
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("kmerFiles", expectedFile1.getAbsolutePath() + "," + expectedFile2.getAbsolutePath());
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final Input input_expected1 = new Input(FileUtils.openInputStream(expectedFile1));
        final Input input_expected2 = new Input(FileUtils.openInputStream(expectedFile2));
        final Input input_test = new Input(FileUtils.openInputStream(output));

        final Kryo kryo = new Kryo();
        kryo.setReferences(false);

        final LargeLongHopscotchSet expectedKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected1);
        final LargeLongHopscotchSet expectedKmerLib_2 = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected2);
        final LargeLongHopscotchSet testKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_test);

        LongIterator iter = expectedKmerLib_2.iterator();
        while (iter.hasNext()) {
            expectedKmerLib.add(iter.next());
        }

        Assert.assertEquals(testKmerLib.size(), expectedKmerLib.size());

        iter = testKmerLib.iterator();
        while (iter.hasNext()) {
            final Long val = iter.next();
            Assert.assertTrue(expectedKmerLib.contains(val));
        }
    }

    @DataProvider(name = "badArgs")
    public Object[][] getBadArguments() {
        return new Object[][]{
                {"kSize", "0"},
                {"kSize", "2"},
                {"kSize", "32"},
                {"kSize", "33"},
                {"bloomFalsePositiveProbability", "-0.1"},
                {"bloomFalsePositiveProbability", "1"},
                {"kmerMask", "0"},
                {"kmerMask", "0,32"},
                {"kmerMask", "-1,15"},
                {"maxPartitionSize", "0"},
                {"maxPartitionSize", "5000"},
                {"kmerSpacing", "0"},
                {"downsampleProbability", "0"},
                {"downsampleProbability", "1.1"}
        };
    }

    @Test(dataProvider = "badArgs", expectedExceptions = Exception.class)
    public void testBadInputs(final String argName, final String argVal) throws Exception {
        final File ref = new File(hg19MiniReference);
        final File output = createTempFile("test", ".hss");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("reference", ref);
        args.addArgument(argName, argVal);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());
    }

}
