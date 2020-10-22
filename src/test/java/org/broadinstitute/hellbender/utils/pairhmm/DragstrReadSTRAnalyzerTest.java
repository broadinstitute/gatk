package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
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
        final DragstrReadSTRAnalyzer analizer = DragstrReadSTRAnalyzer.of(seq.getBytes(), MAX_PERIOD);
        assertCorrectAnalyzerResults(periods, repeats, analizer);
    }

    private void assertCorrectAnalyzerResults(String periods, String repeats, DragstrReadSTRAnalyzer analizer) {
        final int[] periodsAsInts = Arrays.stream(periods.split("\\s+"))
                .mapToInt(s -> Integer.parseInt(s, 16)).toArray();
        final int[] repeatsAsInts = Arrays.stream(repeats.split("\\s+"))
                .mapToInt(s -> Integer.parseInt(s, 16)).toArray();
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

    @Test(dataProvider = "testSequenceAndMaxPeriodData")
    public void testRepeatPeriodAndCount(final String sequenceStr, final int maxPeriod) {
        final byte[] sequence = sequenceStr.getBytes();
        final DragstrReadSTRAnalyzer rpc = DragstrReadSTRAnalyzer.of(sequence, maxPeriod);
        for (int position = 0; position < sequence.length; position++) {
            for (int period = 1; period <= maxPeriod; period++) {
                final int repeatCount = rpc.numberOfRepeats(position, period);
                final int expected = calculateExpectedRepeatLength(sequence, position, period);
                Assert.assertEquals(repeatCount, expected, new String(sequence) + " " + position + " " + period);
            }
        }
    }

    @Test(dataProvider = "testSequenceAndMaxPeriodData")
    public void testRepeatBestPeriodAndCount(final String sequenceStr, final int maxPeriod) {
        final byte[] sequence = sequenceStr.getBytes();
        final DragstrReadSTRAnalyzer rpc = DragstrReadSTRAnalyzer.of(sequence, maxPeriod);
        for (int position = 0; position < sequence.length; position++) {
            final int[] expected = calculateBestPeriodAndRepeat(sequence, position, maxPeriod);
            final int bestPeriod = rpc.mostRepeatedPeriod(position);
            final int bestRepeat = rpc.numberOfMostRepeats(position);
            Assert.assertEquals(bestPeriod, expected[0], new String(sequence) + " " + position);
            Assert.assertEquals(bestRepeat, expected[1], new String(sequence) + " " + position);
        }
    }

    /**
     * Brute force approach to calcuate the "best" repeat unit length and repeat count. The best is the one
     * that has a maximum number of repeated units. In case of a tie, the smaller unit is considered a better
     * one.
     */
    private static int[] calculateBestPeriodAndRepeat(final byte[] sequence, final int position, final int maxPeriod) {
        final int[] result = new int[2];
        result[0] = 1;
        result[1] = calculateExpectedRepeatLength(sequence, position, 1);
        for (int period = 2; period <= maxPeriod; period++) {
            final int candidate = calculateExpectedRepeatLength(sequence, position, period);
            if (candidate > result[1]) {
                result[0] = period;
                result[1] = candidate;
            }
        }
        return result;
    }

    /**
     * Brute force calculation of the repeat length of a particular period at a position of a sequence.
     */
    public static int calculateExpectedRepeatLength(final byte[] sequence, final int position, final int period) {
        if (period > sequence.length) {
            return 0;
        }
        // The simple story:
        // we must calculate the maxium number of adjacent repeats for any of the overlapping k-mers (k = period)
        // at that position + the k-mer that immediately follows it
        // Details:
        //     The outer loop goes thru the start of each of those k-mers. start = [position - period ... position + 1]
        //        Then two inner loops simply calculate the number of copies downstream (forward) and upstream (backwards)
        //     There is an adjustment a the end due to a quirkiness of DRAGEN's STR analyzer in where a copy is discounted
        //        in the first k-1 position of the sequence.
        int start = Math.max(0, position - period + 1);
        int max = 0;
        for (int end = start + period; start <= position + 1 && end <= sequence.length; start++, end++) {
            final byte[] unit = Arrays.copyOfRange(sequence, start, end);
            int forward = 1; // = 1 since we must count the unit starting at position.
            for (int offset = end; offset <= sequence.length - period; offset += period) {
                final byte[] other = Arrays.copyOfRange(sequence, offset, offset + period);
                if (!Arrays.equals(unit, other)) {
                    break;
                }
                forward++;
            }
            int backward = 0;
            for (int offset = start - period; offset >= 0; offset -= period) {
                final byte[] other = Arrays.copyOfRange(sequence, offset, offset + period);
                if (!Arrays.equals(unit, other)) {
                    break;
                }
                backward++;
            }
            int candidate = forward + backward;
            // Matching strange behavior in DRAGEN where repeat lengths at first position for period larger than one
            // are reduced by one (minimum one. (but no less than 1.
            if (position < period - 1 && start == 0) {
                candidate = Math.max(1, candidate - 1);
            }

            max = Math.max(max, candidate);
        }
        return max;
    }

    @DataProvider
    public static Object[][] testSequenceAndMaxPeriodData() {
        final List<Object[]> result = new ArrayList<>();
        final String[] fixSequences = {
                "TGATTTGCTCTGTCTGCTGCTGCTGCCTTCAGTAGGGTTGCACGCCTGGGCACGCCTGGAAT",
                "AGTATACTGAT",
                "GTCTATATATATTTTAATTAATTAATTAATTAAATATATTTTCTGCTGCCTTTTGGAT",
                "AAAAA",
                "A",
                "ACGTAGATCTGTAGCACTATCGAGC"};
        final Random rdn = new Random(131);
        final RandomDNA rdnDNA = new RandomDNA(rdn);

        final String[] randomSequences = new String[100];
        for (int i = 0; i < randomSequences.length; i++) {
            final int length = rdn.nextInt(195)  + 5;
            final byte[] bases = rdnDNA.nextBases(length);
            for (int position = 0; position < bases.length; position++) {
                if (rdn.nextDouble() < 1.0 / length) {
                    final int unitLength = rdn.nextInt(10) + 1;
                    final int intendedRepeats = rdn.nextInt(10) + 2;
                    final byte[] unit = rdnDNA.nextBases(unitLength);
                    for (int repeat = 0, offset = position; repeat < intendedRepeats && offset < length - unitLength; offset += unitLength, repeat++) {
                        System.arraycopy(unit, 0, bases, offset, Math.min(unitLength, bases.length - offset));
                    }
                }
            }
            randomSequences[i] = new String(bases);
        }
        for (final String fixSequence : fixSequences) {
            result.add(new Object[] { fixSequence, Math.max(5, fixSequence.length() / 4)});
        }
        for (final String randomSequence : randomSequences) {
            result.add(new Object[] { randomSequence, Math.max(5, randomSequence.length() / 4 )});
        }
        return result.toArray(new Object[result.size()][]);
    }

    /**
     * App main program to generate the input file to be used with DRAGEN matlab script to generate the truth.
     */
    public static void main(final String[] args) throws IOException {
        final File outfile = File.createTempFile("input", ".txt");
        System.err.println("output file is: " + outfile);
        final int seed = args.length == 0 ? 13 : Integer.parseInt(args[0]);
        final int count = args.length < 2 ? GENERATED_TEST_CASES : Integer.parseInt(args[1]);
        final Random rdn = new Random(seed);
        final RandomDNA rdnDNA = new RandomDNA(rdn);
        final ExponentialDistribution repeatLenDist = new ExponentialDistribution(REPEAT_LENGTH_MEAN);
        try (final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile)))) {
            for (int i = 0; i < count; i++) {
                if (i > 0) writer.write('\n');
                final int length = (int) Math.max(MIN_LENGTH, (rdn.nextGaussian() * GENERATED_SD_SEQ_LENGTH + 1) * GENERATED_AVG_SEQ_LENGTH);
                final byte[] seq = rdnDNA.nextBases(length);
                final PoissonDistribution repeatsDist = new PoissonDistribution(AVERAGE_NUMBER_OF_REPEATS_FOR_100_BPS  * seq.length / 100.0);
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
