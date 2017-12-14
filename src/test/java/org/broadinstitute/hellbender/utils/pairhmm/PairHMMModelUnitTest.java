package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;

/**
 * Unit tests for {@link PairHMMModel}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class PairHMMModelUnitTest extends GATKBaseTest {

    final double TOLERANCE = 1E-9;

    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbs(final int insQual, final int delQual, final int gcp, final double[] expected) {
        final double[] actual = PairHMMModel.qualToTransProbs((byte)insQual,(byte)delQual,(byte)gcp);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, PairHMMModel.TRANS_PROB_ARRAY_LENGTH);
        assertEqualsDoubleArray(actual,expected,TOLERANCE);
        Assert.assertEquals(actual.length, PairHMMModel.TRANS_PROB_ARRAY_LENGTH);
    }

    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbsLog10(final int insQuals, final int delQual, final int gcp, final double[] expected) {
        final double[] logExpected = new double[expected.length];
        for (int i = 0; i < logExpected.length; i++)
            logExpected[i] = Math.log10(expected[i]);
        final double[] actual = PairHMMModel.qualToTransProbsLog10((byte) insQuals, (byte) delQual, (byte) gcp);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, PairHMMModel.TRANS_PROB_ARRAY_LENGTH);
        assertEqualsDoubleArray(actual,logExpected,TOLERANCE);
    }

    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbsFill(final int insQual, final int delQual, final int gcp, final double[] expected) {
        final double[] actual = new double[PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
        PairHMMModel.qualToTransProbs(actual, (byte) insQual, (byte) delQual, (byte) gcp);
        assertEqualsDoubleArray(actual,expected,TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbs(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.qualToTransProbs(insQuals,delQuals,gapQuals);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, expected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length, expected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],expected[i],TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbsLog10(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.qualToTransProbsLog10(insQuals, delQuals, gapQuals);
        final double[][] logExpected = new double[expected.length][expected[0].length];
        for (int i = 1; i < expected.length; i++)
            for (int j = 0; j < expected[0].length; j++)
                logExpected[i][j] = Math.log10(expected[i][j]);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, logExpected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length, logExpected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],logExpected[i],TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbsLog10Fill(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.createTransitionMatrix(insQuals.length);
        PairHMMModel.qualToTransProbsLog10(actual, insQuals, delQuals, gapQuals);
        final double[][] logExpected = new double[expected.length][expected[0].length];
        for (int i = 1; i < expected.length; i++)
            for (int j = 0; j < expected[0].length; j++)
                logExpected[i][j] = Math.log10(expected[i][j]);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, logExpected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length, logExpected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],logExpected[i],TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbsFill(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.createTransitionMatrix(insQuals.length);
        PairHMMModel.qualToTransProbs(actual,insQuals,delQuals,gapQuals);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, expected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length, expected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],expected[i],TOLERANCE);
    }
    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbsLog10Fill(final int insQuals, final int delQual, final int gcp, final double[] expected) {
        final double[] logExpected = new double[expected.length];
        for (int i = 0; i < logExpected.length; i++)
            logExpected[i] = Math.log10(expected[i]);
        final double[] actual = new double[PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
        PairHMMModel.qualToTransProbsLog10(actual, (byte) insQuals, (byte) delQual, (byte) gcp);
        assertEqualsDoubleArray(actual,logExpected,TOLERANCE);
    }


    @DataProvider(name="qualToTransDataProvider")
    public Iterator<Object[]> qualToTransDataProvider() {
        return new Iterator<Object[]>() {

            private final Iterator<Integer> readLengthIterator = readLengthIterator();
            private Iterator<int[]> qualsIterator = qualIterator();

            @Override
            public boolean hasNext() {
                return readLengthIterator.hasNext();
            }

            @Override
            public Object[] next() {
                final int readLength = readLengthIterator.next();
                double[][] matrix = new double[readLength+1][PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
                final byte[] insQuals = new byte[readLength];
                final byte[] delQuals = new byte[readLength];
                final byte[] gapQuals = new byte[readLength];
                for (int i = 0; i < readLength; i++) {
                    if (!qualsIterator.hasNext())
                        qualsIterator = qualIterator();
                    final int[] quals = qualsIterator.next();
                    final int insQual = quals[0];
                    final int delQual = quals[1];
                    final int gapQual = quals[2];
                    final double[] trans = qualsToProbs(insQual, delQual, gapQual);
                    matrix[i+1] = trans;
                    insQuals[i] = (byte)insQual;
                    delQuals[i] = (byte)delQual;
                    gapQuals[i] = (byte)gapQual;
                }

                return new Object[] { insQuals, delQuals, gapQuals, matrix };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }


    @DataProvider(name="qualToProbsDataProvider")
    public Iterator<Object[]> qualToProbsDataProvider() {
        return new Iterator<Object[]>() {
            private final Iterator<int[]> qualsIterator = qualIterator();

            @Override
            public boolean hasNext() {
                return qualsIterator.hasNext();
            }

            @Override
            public Object[] next() {
                final int[] quals = qualsIterator.next();
                final int insQual = quals[0];
                final int delQual = quals[1];
                final int gapQual = quals[2];

                final double[] trans = qualsToProbs(insQual, delQual, gapQual);


                return new Object[] { insQual, delQual, gapQual, trans };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private double[] qualsToProbs(final int insQual, final int delQual, final int gapQual) {
        final double[] trans = new double[PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
        final double matchToMatch = PairHMMModel.matchToMatchProb(insQual, delQual);
        final double matchToInsert = QualityUtils.qualToErrorProb(insQual);
        final double matchToDeletion = QualityUtils.qualToErrorProb(delQual);
        final double indelToMatch = QualityUtils.qualToProb(gapQual);
        final double indelToIndel = QualityUtils.qualToErrorProb(gapQual);

        trans[PairHMMModel.matchToMatch] = matchToMatch;
        trans[PairHMMModel.matchToInsertion] = matchToInsert;
        trans[PairHMMModel.matchToDeletion] = matchToDeletion;
        trans[PairHMMModel.indelToMatch] = indelToMatch;
        trans[PairHMMModel.deletionToDeletion] = trans[PairHMMModel.insertionToInsertion] = indelToIndel;
        return trans;
    }

    private Iterator<Integer> readLengthIterator() {
        return Arrays.asList(READ_LENGTHS).iterator();
    }

    private Iterator<int[]> qualIterator() {
        final int totalCount = INS_QUALS.length * DEL_QUALS.length * GAP_QUALS.length;

        return new Iterator<int[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < totalCount;
            }

            @Override
            public int[] next() {
                final int gap = i % GAP_QUALS.length;
                final int indelGroup = i / GAP_QUALS.length;
                final int del = indelGroup % DEL_QUALS.length;
                final int ins = indelGroup % DEL_QUALS.length;
                i++;
                return new int[] { INS_QUALS[ins], DEL_QUALS[del], GAP_QUALS[gap]};
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }



    @Test(dataProvider = "dualTestDataProvider")
    public void testDoubleQualToProb(final int insQual, final int delQual, final double log10Expected, final double expected) {
        Assert.assertEquals(PairHMMModel.matchToMatchProb(insQual, delQual), expected, TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProb(delQual, insQual), expected, TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProbLog10(insQual, delQual), log10Expected, TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProbLog10(delQual, insQual), log10Expected, TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProb((byte) insQual, (byte) delQual), expected, TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProbLog10((byte) insQual, (byte) delQual), log10Expected, TOLERANCE);
    }

    @DataProvider(name = "dualTestDataProvider")
    private Iterator<Object[]> dualTestDataProvider() {
        final int[] testQuals = new int[] { 0, 1, 2, 5, 10, 13, 17, 20, 23, 27, 30, 43, 57, 70, 100, 200, 254};

        return new Iterator<Object[]>() {
            private int i = 0;
            private int j = 0;

            @Override
            public Object[] next() {

                final int qual1 =  testQuals[i];
                final int qual2 =  testQuals[j];

                final double errorProb1 = Math.pow(10, -0.1 * qual1);
                final double errorProb2 = Math.pow(10, -0.1 * qual2);

                //Note: matchToMatchProb assumes that insertion and deletion never coincide, so the probability that neither occurs
                // is 1 - errorProb1 - errorProb2, without the correction + errorProb1*errorProb2
                final double expected = Math.max(0, (1 - (errorProb1 + errorProb2)));
                final Object[] result = new Object[] { qual1, qual2, Math.log10(Math.min(1, expected)), Math.min(1, expected)};

                if (++j >= testQuals.length) {
                    i++;
                    j = i;
                }
                return result;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }

            @Override
            public boolean hasNext() {
                return i < testQuals.length;
            }
        };
    }


    private static int[] INS_QUALS = {30, 45, 20, 10, 5, 60, 123 };

    private static int[] DEL_QUALS = {30, 45, 20, 10, 5, 60, 123 };

    private static int[] GAP_QUALS = {10, 20, 5};

    private static Integer[] READ_LENGTHS = { 0, 1, 5, 20, 100, 250};
}
