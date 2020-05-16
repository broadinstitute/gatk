package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;


public class DragstrReferenceSTRsTest {

    @Test(dataProvider = "testSequenceAndMaxPeriodData")
    public void testRepeatPeriodAndCountFullSequence(final String sequenceStr, final int maxPeriod) {
        final byte[] sequence = sequenceStr.getBytes();
        final DragstrReferenceSTRs subject = DragstrReferenceSTRs.of(sequence, 0, sequence.length, maxPeriod);
        assertCorrectness(subject, sequenceStr, maxPeriod, 0, sequence.length);
    }

    private void assertCorrectness(final DragstrReferenceSTRs subject, final String sequence, final int maxPeriod, final int start, final int end) {
        for (int position = start; position < end; position++) {
            final int[] expectedBestPeriodAndRepeat = calculateBestPeriodAndRepeat(sequence, position, maxPeriod);
            Assert.assertEquals(subject.period(position), expectedBestPeriodAndRepeat[0], "" + position);
            Assert.assertEquals(subject.repeatLength(position), expectedBestPeriodAndRepeat[1], "" + position);
            final String expectedUnit = sequence.substring(position, position + expectedBestPeriodAndRepeat[0]);
            Assert.assertEquals(subject.repeatUnitAsString(position), expectedUnit);
            Assert.assertEquals(new String(subject.repeatUnit(position)), expectedUnit);
        }
    }

    @Test(dataProvider = "testSequenceAndMaxPeriodDataWithStartEnd")
    public void testRepeatBestPeriodAndCountPartialSequence(final int start, final int end, final int maxPeriod, final String sequenceStr) {
        for (int position = start; position < end; position++) {
            final byte[] sequence = sequenceStr.getBytes();
            final DragstrReferenceSTRs subject = DragstrReferenceSTRs.of(sequence, 0, sequence.length, maxPeriod);
            assertCorrectness(subject, sequenceStr, maxPeriod, start, end);
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
    public static int[] calculateBestPeriodAndRepeat(final String sequence, final int position, final int maxPeriod) {
        final int[] result = new int[2];
        result[0] = 1;
        result[1] = calculateRepeats(sequence, position, 1, true);
        for (int period = 2; period <= maxPeriod; period++) {
            final int candidate = calculateRepeats(sequence, position, period, true);
            if (candidate > result[1]) {
                result[0] = period;
                result[1] = candidate;
            }
        }
        result[1] = calculateRepeats(sequence, position, result[0], false);
        return result;
    }

    public static int calculateRepeats(String sequence, final int position, final int period, final boolean onlyForward) {
        if (period + position > sequence.length()) {
            return 0;
        }
        int result = 1; // the first match is a given.
        // forward counting.
        outter:
        for (int offset = period + position; offset <= sequence.length() - period; offset += period) {
            for (int i = 0; i < period; i++) {
                if (sequence.charAt(offset + i) != sequence.charAt(position + i)) {
                    break outter;
                }
            }
            result++;
        }

        if (onlyForward) {
            return result;
        }
        // backward counting.
        outter:
        for (int offset = position - period; offset >= 0; offset -= period) {
            for (int i = 0; i < period; i++) {
                if (sequence.charAt(offset + i) != sequence.charAt(position + i)) {
                    break outter;
                }
            }
            result++;
        }
        return result;
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
                "",
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
            result.add(new Object[] { randomSequence, Math.max(5, randomSequence.length() / 4 ) });
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider
    public static Object[][] testSequenceAndMaxPeriodDataWithStartEnd() {
        final Object[][] baseData = testSequenceAndMaxPeriodData();
        final Random rdn = new Random(313111241);
        final List<Object[]>  all = Arrays.stream(baseData)
                .flatMap(in -> {
                    final String seq = (String) in[0];
                    final int maxPeriod = (Integer) in[1];
                    return IntStream.range(0, seq.length())
                            .boxed()
                            .flatMap(start -> IntStream.range(start, seq.length() + 1).mapToObj(end -> new int[] {start, end}))
                            .map(startEnd -> new Object[] {startEnd[0], startEnd[1], maxPeriod, seq});
                }).collect(Collectors.toList());
        Collections.shuffle(all, rdn);
        return all.subList(0, 10000).toArray(new Object[10000][]);
    }
}
