package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

public final class RandomDNAUnitTest {

    private int[] counts(final byte[] bytes){
        final int[] b= new int[4];
        for(int i=0; i < bytes.length; i++){
            switch (bytes[i]){
                case 'A': b[0]++; break;
                case 'C': b[1]++; break;
                case 'G': b[2]++; break;
                case 'T': b[3]++; break;
                default: throw new IllegalStateException("illegal base:" + bytes[i]);
            }
        }
        return b;
    }
    @Test
    public void testBases1(){
        int[] results = new int[4];

        final int n = 1000;
        final int m = 13;
        for (int i= 0; i < n; i++) {
            final byte[] b = new RandomDNA().nextBases(m);
            final int[] b0 = counts(b);
            results = pairwiseAdd(results, b0);
        }

        checkResults(results, n, m);
    }

    @Test
    public void testBases(){
        int[] results = new int[4];

        final int n = 1000;
        final int m = 13;
        for (int i= 0; i < n; i++) {
            final byte[] b = new byte[m];
            new RandomDNA().nextBases(b);
            final int[] b0 = counts(b);
            results = pairwiseAdd(results, b0);
        }

        checkResults(results, n, m);
    }

    public void checkResults(final int[] results, final int n, final int m) {
        final double[] dresults = MathUtils.promote(results);
        final double mean = MathUtils.mean(dresults, 0, dresults.length);
        final double std = MathUtils.stddev(dresults, 0, dresults.length);
        final double expectedMean = (n*m)/4.0;
        final double s = std; // not really because it's the population not the sample dtd but it'll do
        Assert.assertTrue(mean < expectedMean + 2 * s / Math.sqrt(n * m), "unexpected mean:" + mean);
        Assert.assertTrue(mean > expectedMean-2*s/Math.sqrt(n*m), "unexpected mean:" +mean);
    }

    private int[] pairwiseAdd(int[] a, int[] b) {
        if (a.length != b.length){
            throw  new IllegalArgumentException();
        }
        final int[] results = new int[a.length];
        for (int i = 0; i < a.length; i++) {
            results[i] = a[i] + b[i];
        }
        return results;
    }
}
