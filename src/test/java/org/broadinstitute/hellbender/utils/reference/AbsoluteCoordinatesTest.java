package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

public final class AbsoluteCoordinatesTest extends BaseTest {

    @Test(dataProvider = "randomDictionaries")
    public void testAll(final SAMSequenceDictionary dictionary) {
        final AbsoluteCoordinates subject = AbsoluteCoordinates.of(dictionary);
        final Random rdn = new Random(31 * dictionary.hashCode());
        for (int i = 0; i < 10000; i++) {
            final SAMSequenceRecord record = dictionary.getSequence(rdn.nextInt(dictionary.getSequences().size()));
            final int pos = 1 + rdn.nextInt(record.getSequenceLength());
            final long abs = subject.toAbsolute(record.getSequenceIndex(), pos);
            final long expectedAbs = slowToAbsolute(dictionary, record.getSequenceIndex(), pos);
            Assert.assertEquals(abs, expectedAbs, "" + record);
        }
        for (int i = 0; i < dictionary.getSequences().size(); i++) {
            final SAMSequenceRecord seq = dictionary.getSequence(i);
            final long lastAbs = subject.toAbsolute(i, seq.getSequenceLength());
            final long firstAbs = subject.toAbsolute(i, 1);
            Assert.assertEquals(firstAbs, slowToAbsolute(dictionary, i, 1));
            Assert.assertEquals(lastAbs, slowToAbsolute(dictionary, i, seq.getSequenceLength()));
        }

        for (int i = 0; i < 100000; i++) {
            final long abs = ThreadLocalRandom.current().nextLong(dictionary.getReferenceLength())  +1;
            final SimpleInterval relative = subject.toSimpleInterval(abs, 1);
            final SimpleInterval expected = slowToRelative(dictionary, abs);
            Assert.assertEquals(relative, expected, "" + i);
        }

        for (int i = 0; i < dictionary.getSequences().size(); i++) {
            final long absStart = 1 + IntStream.range(0, i).map(n -> dictionary.getSequence(n).getSequenceLength()).sum();
            final long absEnd = absStart + dictionary.getSequence(i).getSequenceLength() - 1;
            final SimpleInterval relEnd = subject.toSimpleInterval(absEnd, 1);
            final SimpleInterval relStart = subject.toSimpleInterval(absStart, 1);
            Assert.assertEquals(relStart, new SimpleInterval(dictionary.getSequence(i).getSequenceName(), 1, 1));
            final int seqLength =  dictionary.getSequence(i).getSequenceLength();
            Assert.assertEquals(relEnd, new SimpleInterval(dictionary.getSequence(i).getSequenceName(),seqLength, seqLength));
        }
    }

    @Test(dataProvider = "randomDictionaries", enabled = false)
    public void testContinuos(final SAMSequenceDictionary dictionary) {
        final AbsoluteCoordinates subject = AbsoluteCoordinates.of(dictionary);
        final Random rdn = new Random(31 * dictionary.hashCode());
        final int firstCtgLength = dictionary.getSequence(0).getSequenceLength();
        for (int i = 1; i <= firstCtgLength; i++) {
            final SimpleInterval si = subject.toSimpleInterval(i, 1);
            Assert.assertEquals(si, new SimpleInterval(dictionary.getSequence(0).getSequenceName(), i, i));
        }
    }

    @Test(dataProvider = "randomDictionaries", enabled = false)
    public void testRandom(final SAMSequenceDictionary dictionary) {
        final AbsoluteCoordinates subject = AbsoluteCoordinates.of(dictionary);
        final int firstCtgLength = dictionary.getSequence(0).getSequenceLength();
        final ThreadLocalRandom rdn = ThreadLocalRandom.current();
        final long size = dictionary.getReferenceLength();
        for (int i = 1; i <= firstCtgLength; i++) {
            final long p = rdn.nextLong(size) + 1;
            subject.toSimpleInterval(p, 1);
        }
    }

    @Test(dataProvider = "randomDictionaries", enabled = false)
    public void benckMarkRelative(final SAMSequenceDictionary dictionary) {
        final AbsoluteCoordinates subject = AbsoluteCoordinates.of(dictionary);
        final long total = dictionary.getReferenceLength();
        for (long i =1; i <= total; i++) {
            subject.toSimpleInterval(i, 1);
        }
    }

    public long slowToAbsolute(final SAMSequenceDictionary dict, final int ctg, final int pos) {
        long accu = 0;
        for (final SAMSequenceRecord seq : dict.getSequences()) {
            if (seq.getSequenceIndex() == ctg) {
                break;
            }
            accu += seq.getSequenceLength();
        }
        return accu + pos;
    }

    public SimpleInterval slowToRelative(final SAMSequenceDictionary dict, final long abs) {
        long remaining = abs;
        for (final SAMSequenceRecord record : dict.getSequences()) {
            if (remaining <= record.getSequenceLength()) {
                final int pos = (int) remaining;
                return new SimpleInterval(record.getSequenceName(), pos, pos);
            } else {
                remaining -= record.getSequenceLength();
            }
        }
        throw new IllegalArgumentException("abs to large " + abs + " > " + dict.getReferenceLength());
    }

    @DataProvider
    public Object[][] randomDictionaries() {
        final List<Object[]> result = new ArrayList<>(100);
        final Random rdn = new Random(13);
        final IntegerDistribution seqLengthDistr = new BoundedDicreteParetoDistribution(1000, 100_000_000, 0.1);
        final SAMSequenceDictionary oneSeqDictionary =randomDictionary(rdn, 1, seqLengthDistr, true);
        result.add(new Object[] {oneSeqDictionary});
        final IntegerDistribution numberSeqDist = new PoissonDistribution(100);
        while (result.size() < 100) {
            final int numberSeq = 1 + numberSeqDist.inverseCumulativeProbability(rdn.nextDouble());
            final SAMSequenceDictionary randomDictionary = randomDictionary(rdn, numberSeq, seqLengthDistr, true);
            result.add(new Object[]{randomDictionary});
        }
        return result.stream().toArray(Object[][]::new);
    }

    public SAMSequenceDictionary randomDictionary(final Random rdn, final int seqNum, final IntegerDistribution lengthDistribution, final boolean sortBySize) {
        final List<SAMSequenceRecord> contigs = new ArrayList<>(seqNum);
        for (int i = 0; i < seqNum; i++) {
            final SAMSequenceRecord contig = new SAMSequenceRecord("seq" + (i + 1), lengthDistribution.inverseCumulativeProbability(rdn.nextDouble()));
            contigs.add(contig);
        }
        if (sortBySize) {
            Collections.sort(contigs, Comparator.comparingInt(SAMSequenceRecord::getSequenceLength).reversed());
        }
//        for (final SAMSequenceRecord record : contigs) {
//            System.err.println(record.getSequenceLength());
//        }
        return new SAMSequenceDictionary(contigs);
    }


    private static class BoundedDicreteParetoDistribution extends AbstractIntegerDistribution {

        private static final long serialVersionUID = -1;

        private final int L, H;
        private final double alpha;

        private BoundedDicreteParetoDistribution(final int L, final int H, final double alpha) {
            super(new JDKRandomGenerator());
            this.L = L;
            this.H = H;
            this.alpha = alpha;
        }

        @Override
        public double probability(int x) {
            return (alpha * Math.pow(L, alpha) * Math.pow(x, -alpha -1))
                    / (1 - Math.pow((L/(double)H), alpha));
        }

        @Override
        public double cumulativeProbability(int x) {
            return (1- Math.pow(L, alpha)*Math.pow(x, -alpha)) / (1 - Math.pow(L/(double)H, alpha));
        }

        @Override
        public double cumulativeProbability(int x0, int x1) {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public int inverseCumulativeProbability(final double p) {
            final double result = L / Math.pow((1 - p * (1 - Math.pow(L/(double)H,alpha))), 1/alpha);

            //final int result =  (int) Math.round(Math.pow(- (p * Math.pow(H, alpha) - p * Math.pow(L, alpha) - Math.pow(H * L, alpha))
            //        / (Math.pow(H * L, alpha)), - 1.0 / alpha));
            if (result < L || result > H) {
                throw new IllegalStateException("bad output " + result);
            }
//            System.err.println((int) result);
            return (int) Math.round(result);
        }

        @Override
        public double getNumericalMean() {
            return alpha != 1.0
                    ? ((Math.pow(L, alpha) * alpha) / ((1 - Math.pow(L / (double) H, alpha)) * (alpha - 1))) * ( Math.pow(L, - alpha + 1) - Math.pow(H, -alpha + 1) )
                    : ((H * L) / ((double) H - L)) * (Math.log(H) - Math.log(L));
        }

        @Override
        public double getNumericalVariance() {
            return alpha != 2.0
                    ? ((Math.pow(L, alpha) * alpha) / ((1 - Math.pow(L / (double) H, alpha)) * (alpha - 2))) * ( Math.pow(L, - alpha + 2) - Math.pow(H, -alpha + 2) )
                    : ((2 * H * H * L * L) / ((double) H * H - L * L)) * (Math.log(H) - Math.log(L));
        }

        @Override
        public int getSupportLowerBound() {
            return L;
        }

        @Override
        public int getSupportUpperBound() {
            return H;
        }

        @Override
        public boolean isSupportConnected() {
            return true;
        }

    }


}
