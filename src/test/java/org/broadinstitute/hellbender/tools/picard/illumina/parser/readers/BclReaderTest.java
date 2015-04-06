package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.BclData;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class BclReaderTest {

    public static final File TestDataDir = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/readerTests");
    public static final File PASSING_BCL_FILE = new File(TestDataDir, "bcl_passing.bcl");
    public static final File QUAL_0FAILING_BCL_FILE = new File(TestDataDir, "bcl_failing.bcl");
    public static final File QUAL_1FAILING_BCL_FILE = new File(TestDataDir, "bcl_failing2.bcl");
    public static final File FILE_TOO_LONG = new File(TestDataDir, "bcl_tooLong.bcl");
    public static final File FILE_TOO_SHORT = new File(TestDataDir, "bcl_tooShort.bcl");

    public static final char[] expectedBases = new char[]{
            'C', 'A', 'A', 'A', 'T', 'C', 'T', 'G', 'T', 'A', 'A', 'G', 'C', 'C', 'A', 'A',
            'C', 'A', 'C', 'C', 'A', 'A', 'C', 'G', 'A', 'T', 'A', 'C', 'A', 'A', 'C', 'A',
            'T', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'G', 'C', 'A', 'A', 'G', 'T', 'G', 'C',
            'A', 'C', 'G', 'T', 'A', 'C', 'A', 'A', 'C', 'G', 'C', 'A', 'C', 'A', 'T', 'T',
            'T', 'A', 'A', 'G', 'C', 'G', 'T', 'C', 'A', 'T', 'G', 'A', 'G', 'C', 'T', 'C',
            'T', 'A', 'C', 'G', 'A', 'A', 'C', 'C', 'C', 'A', 'T', 'A', 'T', 'G', 'G', 'G',
            'C', 'T', 'G', 'A', 'A', '.', '.', 'G', 'A', 'C', 'C', 'G', 'T', 'A', 'C', 'A',
            'G', 'T', 'G', 'T', 'A', '.'
    };

    public static final int[] expectedQuals = new int[]{
            18, 29, 8, 17, 27, 25, 28, 27, 9, 29, 8, 20, 25, 24, 27, 27,
            30, 8, 19, 24, 29, 29, 25, 28, 8, 29, 26, 24, 29, 8, 18, 8,
            29, 28, 26, 29, 25, 8, 26, 25, 28, 25, 8, 28, 28, 27, 29, 26,
            25, 26, 27, 25, 8, 18, 8, 26, 24, 29, 25, 8, 24, 8, 25, 27,
            27, 25, 8, 28, 24, 27, 25, 25, 8, 27, 25, 8, 16, 24, 28, 25,
            28, 8, 24, 27, 25, 8, 20, 29, 24, 27, 28, 8, 23, 10, 23, 11,
            15, 11, 10, 12, 12, 2, 2, 31, 24, 8, 4, 36, 12, 17, 21, 4,
            8, 12, 18, 23, 27, 2
    };

    public byte[] qualsAsBytes() {
        final byte[] byteVals = new byte[expectedQuals.length];
        for (int i = 0; i < byteVals.length; i++) {
            byteVals[i] = (byte) expectedQuals[i];
        }
        return byteVals;
    }

    @Test
    public void readValidFile() {
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);
        final BclReader reader = new BclReader(PASSING_BCL_FILE, bclQualityEvaluationStrategy, false);
        final byte[] quals = qualsAsBytes();

        Assert.assertEquals(reader.numClustersPerCycle[0], expectedBases.length);

        int readNum = 0;
        while (readNum < reader.numClustersPerCycle[0]) {
            final BclData bv = reader.next();
            Assert.assertEquals(bv.bases[0][0], expectedBases[readNum], " On num cluster: " + readNum);
            Assert.assertEquals(bv.qualities[0][0], quals[readNum], " On num cluster: " + readNum);
            ++readNum;
        }
        bclQualityEvaluationStrategy.assertMinimumQualities();
        reader.close();
    }

    @DataProvider(name = "failingFiles")
    public Object[][] failingFiles() {
        return new Object[][]{
                {QUAL_0FAILING_BCL_FILE},
                {QUAL_1FAILING_BCL_FILE},
                {new File(TestDataDir, "SomeNoneExistentFile.bcl")},
                {FILE_TOO_LONG},
                {FILE_TOO_SHORT}
        };
    }

    @Test(expectedExceptions = IlluminaReaderException.class, dataProvider = "failingFiles")
    public void failingFileTest(final File failingFile) {
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);
        final BclReader reader = new BclReader(failingFile, bclQualityEvaluationStrategy, false);
        Assert.assertEquals(reader.numClustersPerCycle[0], expectedBases.length);
        while (reader.hasNext()) {
            reader.next();
        }
        reader.close();
        bclQualityEvaluationStrategy.assertMinimumQualities();
    }

    /**
     * Asserts appropriate functionality of a quality-minimum-customized BLC reader, such that (1) if sub-Q2 qualities are found, the BCL
     * reader does not throw an exception, (2) sub-minimum calls are set to quality 1 and (3) sub-minimum calls are counted up properly.
     */
    @Test
    public void lowQualityButPassingTest() throws ExecutionException, InterruptedException {
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(1);

        // Build a list of callables, then submit them and check for errors.
        final Collection<Callable<Void>> callables = new LinkedList<Callable<Void>>();
        for (int i = 0; i < 10; i++) {
            final boolean even_i = i % 2 == 0;
            callables.add(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    final BclReader reader = new BclReader(even_i ? QUAL_1FAILING_BCL_FILE : QUAL_0FAILING_BCL_FILE,
                            bclQualityEvaluationStrategy, false);
                    Assert.assertEquals(reader.numClustersPerCycle[0], expectedBases.length);
                    while (reader.hasNext()) {
                        reader.next();
                    }
                    reader.close();
                    return null;
                }
            });
        }
        final ExecutorService executorService = Executors.newFixedThreadPool(callables.size());
        final Collection<Future<Void>> futures = new LinkedList<Future<Void>>();
        for (final Callable<Void> callable : callables) {
            futures.add(executorService.submit(callable));
        }
        for (final Future<Void> future : futures) {
            future.get();
        }
        bclQualityEvaluationStrategy.assertMinimumQualities();
        Assert.assertEquals((int) bclQualityEvaluationStrategy.getPoorQualityFrequencies().get((byte) 0), 25);
        Assert.assertEquals((int) bclQualityEvaluationStrategy.getPoorQualityFrequencies().get((byte) 1), 25);
    }

    @Test(expectedExceptions = IlluminaReaderException.class)
    public void lowQualityAndFailingTest() throws ExecutionException, InterruptedException {
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);

        // Build a list of callables, then submit them and check for errors.
        final Collection<Callable<Void>> callables = new LinkedList<Callable<Void>>();
        for (int i = 0; i < 10; i++) {
            final boolean even_i = i % 2 == 0;
            callables.add(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    final BclReader reader = new BclReader(even_i ? QUAL_1FAILING_BCL_FILE : QUAL_0FAILING_BCL_FILE,
                            bclQualityEvaluationStrategy, false);
                    Assert.assertEquals(reader.numClustersPerCycle[0], expectedBases.length);
                    while (reader.hasNext()) {
                        reader.next();
                    }
                    reader.close();
                    return null;
                }
            });
        }
        final ExecutorService executorService = Executors.newFixedThreadPool(callables.size());
        final Collection<Future<Void>> futures = new LinkedList<Future<Void>>();
        for (final Callable<Void> callable : callables) {
            futures.add(executorService.submit(callable));
        }
        for (final Future<Void> future : futures) {
            future.get();
        }
        Assert.assertEquals((int) bclQualityEvaluationStrategy.getPoorQualityFrequencies().get((byte) 0), 25);
        Assert.assertEquals((int) bclQualityEvaluationStrategy.getPoorQualityFrequencies().get((byte) 1), 25);
        bclQualityEvaluationStrategy.assertMinimumQualities();
    }
}
