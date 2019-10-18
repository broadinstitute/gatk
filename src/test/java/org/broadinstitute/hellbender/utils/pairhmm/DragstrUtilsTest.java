package org.broadinstitute.hellbender.utils.pairhmm;

import breeze.stats.distributions.Rand;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;


public class DragstrUtilsTest {

    @Test(dataProvider = "testSequenceAndMaxPeriodData")
    public void testRepeatPeriodAndCount(final String sequenceStr, final int maxPeriod) {
        final byte[] sequence = sequenceStr.getBytes();
        final DragstrUtils.STRSequenceAnalyzer rpc = DragstrUtils.repeatPeriodAndCounts(sequence.length, maxPeriod);
        rpc.load(sequence);
        final Random rdn = new Random(Arrays.hashCode(sequence) * 31 + maxPeriod);
        final int[] positions = new int[sequence.length];
        for (int i = 0; i < positions.length; i++) {
            positions[i] = i;
        }
        ArrayUtils.shuffle(positions, rdn);
        for (int position : positions) {
            final int[] periods = IntStream.range(1, maxPeriod + 1).toArray();
            ArrayUtils.shuffle(periods, rdn);
            for (final int period : periods) {
                final int repeatCount = rpc.numberOfRepeats(position, period);
                final int expected = calculate(sequence, position, period);
                Assert.assertEquals(repeatCount, expected, new String(sequence) + " " + position + " " + period);
            }
        }
    }

    @Test(dataProvider = "testSequenceAndMaxPeriodData")
    public void testRepeatBestPeriodAndCount(final String sequenceStr, final int maxPeriod) {
        final byte[] sequence = sequenceStr.getBytes();
        final DragstrUtils.STRSequenceAnalyzer rpc = DragstrUtils.repeatPeriodAndCounts(sequence.length, maxPeriod);
        rpc.load(sequence);
        final Random rdn = new Random(Arrays.hashCode(sequence) * 31 + maxPeriod);
        final int[] positions = new int[sequence.length];
        for (int i = 0; i < positions.length; i++) {
            positions[i] = i;
        }
        ArrayUtils.shuffle(positions, rdn);
        for (int position : positions) {
            final int[] periods = IntStream.range(1, maxPeriod + 1).toArray();
            ArrayUtils.shuffle(periods, rdn);
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
     * @param sequence
     * @param position
     * @param maxPeriod
     * @return
     */
    public static int[] calculateBestPeriodAndRepeat(final byte[] sequence, final int position, final int maxPeriod) {
        final int[] result = new int[2];
        result[0] = 1;
        result[1] = calculate(sequence, position, 1);
        for (int period = 2; period <= maxPeriod; period++) {
            final int candidate = calculate(sequence, position, period);
            if (candidate > result[1]) {
                result[0] = period;
                result[1] = candidate;
            }
        }
        return result;
    }


    public static int calculate(final byte[] sequence, final int position, final int period) {
        final int start = Math.max(0, position - period + 1);
        final int windows = position - period + 1 >= 0 ? period : position + 1;
        int max = 0;
        for (int i = 0; i < windows; i++) {
            final byte[] unit = Arrays.copyOfRange(sequence, start + i, start + i + period);
            int forward = 0;
            for (int offset = start + i + period; offset <= sequence.length - period; offset += period) {
                final byte[] other = Arrays.copyOfRange(sequence, offset, offset + period);
                if (!Arrays.equals(unit, other)) {
                    break;
                }
                forward++;
            }
            int backward = 0;
            for (int offset = start + i - period; offset >= 0; offset -= period) {
                final byte[] other = Arrays.copyOfRange(sequence, offset, offset + period);
                if (!Arrays.equals(unit, other)) {
                    break;
                }
                backward++;
            }
            final int candidate = forward + backward + 1;
            if (candidate > max) {
                max = candidate;
            }
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
            result.add(new Object[] { fixSequence, Math.max(5, fixSequence.length() / 4) });
        }
        for (final String randomSequence : randomSequences) {
            result.add(new Object[] { randomSequence, Math.max(5, randomSequence.length() / 4 )});
        }
        return result.toArray(new Object[result.size()][]);
    }
}
