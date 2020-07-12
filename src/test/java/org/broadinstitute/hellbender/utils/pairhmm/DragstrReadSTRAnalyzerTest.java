package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ExponentialDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.PoissonDistribution;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPInputStream;

/**
 * This test uses the matlab routine shared by Illumina to generate use cases.
 *
 * Unfortunately we cannot share the routine code, as is not ditributable under this software license. You still
 * can find the matlab script that invoke such a routine, which was written by the Broad. You may be able to get
 * the routine source code if you contact Illumina directly.
 */
public final class DragstrReadSTRAnalyzerTest extends BaseTest {

    private static final int MAX_PERIOD = 8;
    private static final int MAX_REPEATS = 20;
    private static final int GENERATED_TEST_CASES = 1000;
    private static final double AVERAGE_NUMBER_OF_REPEATS_FOR_100_BPS = 1.0;

    private static final int GENERATED_AVG_SEQ_LENGTH = 300;
    private static final int MIN_LENGTH = 10;
    private static final double GENERATED_SD_SEQ_LENGTH = 1.0;
    private static final double REPEAT_LENGTH_MEAN = 3.0;

    private static final String GENERATED_TEST_CASES_FILENAME = "./src/test/resources/large/dragstr-str-read-period-repeat-test-cases-with-induced-strs.txt.gz";

    @Test(dataProvider = "testCasesData")
    public void testSingleSequence(final String seq, final String periods, final String repeats) {
        final DragstrReadSTRAnalyzer analizer = DragstrUtils.repeatPeriodAndCounts(seq.length(), MAX_PERIOD);
        assertCorrectAnalyzerResults(seq, periods, repeats, analizer);
    }

    private void assertCorrectAnalyzerResults(String seq, String periods, String repeats, DragstrReadSTRAnalyzer analizer) {
        final int[] periodsAsInts = Arrays.stream(periods.split("\\s+"))
                .mapToInt(s -> Integer.parseInt(s, 16)).toArray();
        final int[] repeatsAsInts = Arrays.stream(repeats.split("\\s+"))
                .mapToInt(s -> Integer.parseInt(s, 16)).toArray();
        analizer.load(seq.getBytes());
        assertCorrectAnalyzerResults(periodsAsInts, repeatsAsInts, analizer);
    }

    private void assertCorrectAnalyzerResults(final int[] periodsAsInts, final int[] repeatsAsInts, final DragstrReadSTRAnalyzer analizer) {
        if (MathUtils.arrayMax(periodsAsInts) > MAX_PERIOD || MathUtils.arrayMax(repeatsAsInts) > MAX_REPEATS) {
            throw new SkipException("");
        }
        for (int i = 0; i < periodsAsInts.length; i++) {
            Assert.assertEquals(analizer.mostRepeatedPeriod(i), periodsAsInts[i], "wrong period at position " + i + " where exepected repeat is " + repeatsAsInts[i] + " and/but it is " + analizer.numberOfMostRepeats(i) + "; ");
            Assert.assertEquals(Math.min(MAX_REPEATS, analizer.numberOfMostRepeats(i)), repeatsAsInts[i], "wrong repeat number at position " + i + " where exepected period is " + periodsAsInts[i] + " and/but it is " + analizer.mostRepeatedPeriod(i) + "; ");
        }
    }

    @Test
    public void testMultiSequence() {
        final Object[][] cases = testCasesData();
        final int maxSeqLength = Arrays.stream(cases).map(o -> (String)o[0]).mapToInt(String::length).max().orElse(0);
        final DragstrReadSTRAnalyzer subject = DragstrUtils.repeatPeriodAndCounts(maxSeqLength, MAX_PERIOD);
        for (final Object[] params : cases) {
            final String seq = (String) params[0];
            final String periods = (String) params[1];
            final String repeats = (String) params[2];
            assertCorrectAnalyzerResults(seq, periods, repeats, subject);
        }
    }

    @DataProvider(name = "testCasesData")
    public Object[][] testCasesData() {
        final List<Object[]> result = new ArrayList<>();
        // An empirical test case.
        result.add(new Object[] {
                "ATTTTTTCAATGTTTACACATTTCCTTCCTCCCTCCCTCCTTCCTTTCCTCCCTTCCTCCCTTCCTCCCTTCCTTCCTGTTTGCTTTATTATTGTATTG",
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 4 4 4 1 1 1 1 1 1 1 1 4 4 4 1 1 1 1 1 1 1 8 8 1 1 1 1 8 8 8 8 1 1 1 1 8 8 8 8 1 1 1 1 8 8 8 8 8 1 1 1 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 1 1 1",
                "6 6 6 6 6 6 6 2 2 2 1 3 3 3 3 2 2 2 2 3 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 3 3 3 3 1 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2"
        });
        try (final BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(GENERATED_TEST_CASES_FILENAME))))) {
            String seq;
            while ((seq = reader.readLine()) != null) {
                final String periods = reader.readLine();
                final String repeats = reader.readLine();
                result.add(new Object[] {seq, periods, repeats});
            }

        } catch (final IOException e) {
            throw new UncheckedIOException(e);
        }
        return result.toArray(new Object[result.size()][]);
    }


    public static void main(final String[] args) throws IOException, MathException {
        final File outfile = File.createTempFile("input", ".txt");
        System.err.println("output file is: " + outfile);
        final int seed = args.length == 0 ? 13 : Integer.parseInt(args[0]);
        final int count = args.length < 2 ? GENERATED_TEST_CASES : Integer.parseInt(args[1]);
        final Random rdn = new Random(seed);
        final RandomDNA rdnDNA = new RandomDNA(rdn);
        final ExponentialDistribution repeatLenDist = new ExponentialDistributionImpl(REPEAT_LENGTH_MEAN);
        try (final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile)))) {
            for (int i = 0; i < count; i++) {
                if (i > 0) writer.write('\n');
                final int length = (int) Math.max(MIN_LENGTH, (rdn.nextGaussian() * GENERATED_SD_SEQ_LENGTH + 1) * GENERATED_AVG_SEQ_LENGTH);
                final byte[] seq = rdnDNA.nextBases(length);
                final PoissonDistribution repeatsDist = new PoissonDistributionImpl(AVERAGE_NUMBER_OF_REPEATS_FOR_100_BPS  * seq.length / 100.0);
                final int numRepeats = Math.round(repeatsDist.inverseCumulativeProbability(rdn.nextDouble()));
                for (int j = 0; j < numRepeats; j++) {
                    final int start = rdn.nextInt(seq.length);
                    final int unitLength = rdn.nextInt(MAX_PERIOD) + 1;
                    int repeatLengthInBases = (int) (unitLength * Math.round(repeatLenDist.inverseCumulativeProbability(rdn.nextDouble())));
                    if (repeatLengthInBases + start > seq.length) {
                        repeatLengthInBases = seq.length - start;
                    }
                    for (int k = unitLength; k < repeatLengthInBases; k++) {
                        seq[start + k] = seq[start + k % unitLength];
                    }
                }
                writer.write(new String(seq));
            }
        }
    }
}
