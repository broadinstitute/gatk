package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class KmerUnitTest extends BaseTest {
    @DataProvider(name = "KMerCreationData")
    public Object[][] makeKMerCreationData() {
        List<Object[]> tests = new ArrayList<>();

        final String bases = "ACGTAACCGGTTAAACCCGGGTTT";
        for ( int start = 0; start < bases.length(); start++ ) {
            for ( int length = 1; start + length < bases.length(); length++ ) {
                final String myBases = bases.substring(start, start+length);
                tests.add(new Object[]{bases.getBytes(), start, length, myBases});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "KMerCreationData")
    public void testFullConstructor(final byte[] allBases, final int start, final int length, final String expected) {
        testKmerCreation(new Kmer(allBases, start, length), start, length, expected);
    }

    @Test(dataProvider = "KMerCreationData")
    public void testByteConstructor(final byte[] allBases, final int start, final int length, final String expected) {
        testKmerCreation(new Kmer(Arrays.copyOfRange(allBases, start, start + length)), 0, length, expected);
    }

    @Test(dataProvider = "KMerCreationData")
    public void testStringConstructor(final byte[] allBases, final int start, final int length, final String expected) {
        testKmerCreation(new Kmer(new String(Arrays.copyOfRange(allBases, start, start + length))), 0, length, expected);
    }

    private void testKmerCreation(final Kmer kmer, final int start, final int length, final String expected) {
        Assert.assertEquals(kmer.start(), start);
        Assert.assertEquals(kmer.length(), length);
        Assert.assertEquals(new String(kmer.bases()), expected);

        // check that the caching is working by calling again
        Assert.assertEquals(kmer.start(), 0);
        Assert.assertEquals(kmer.length(), length);
        Assert.assertEquals(new String(kmer.bases()), expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeStart() throws Exception {
        final byte[] bases = "ACGTACGT".getBytes();
        new Kmer(bases, -2, 3);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeLength() throws Exception {
        final byte[] bases = "ACGTACGT".getBytes();
        new Kmer(bases, 2, -3);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEndTooFar() throws Exception {
        final byte[] bases = "ACGTACGT".getBytes();
        new Kmer(bases, 2, 10);
    }

    @Test
    public void testEquals() {
        final byte[] bases = "ACGTACGT".getBytes();
        final Kmer eq1 = new Kmer(bases, 0, 3);
        final Kmer eq2 = new Kmer(bases, 4, 3);
        final Kmer eq4 = new Kmer(new Kmer(bases, 4, 3).bases());
        final Kmer neq = new Kmer(bases, 1, 3);

        Assert.assertNotEquals(eq1, "FRED");

//        for ( final Kmer eq : Arrays.asList(eq1, eq2) ) { // TODO -- deal with me
        for ( final Kmer eq : Arrays.asList(eq1, eq2, eq4) ) {
            Assert.assertEquals(eq1, eq, "Should have been equal but wasn't: " + eq1.hashCode() + " vs " + eq.hashCode()); // , "should be equals " + eq1 + " with " + eq);
            Assert.assertEquals(eq1.hashCode(), eq.hashCode());
            Assert.assertNotEquals(eq, neq, "incorrectly equals " + eq + " with " + neq);
        }
    }

    @Test
    public void testSubkmer() {
        final String bases = "ACGT";
        final Kmer one = new Kmer(bases.getBytes());

        for ( int start = 0; start < bases.length(); start++ ) {
            for ( int length = 0; start + length < bases.length(); length++ ) {
                Assert.assertEquals(new String(one.subKmer(start, length).bases()), bases.substring(start, start + length));
            }
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void getDifferingPositions() {
        final byte[] bases = "ACGTACGT".getBytes();
        final Kmer eq1 = new Kmer(bases, 0, 3);
        final byte[] newBases = bases.clone();
        final int[] differingIndices = new int[newBases.length];
        final byte[] differingBases = new byte[newBases.length];

        eq1.getDifferingPositions(eq1, -2, differingIndices, differingBases);
    }

    @Test
    public void testStringsTooFar() throws Exception {
        final byte[] bases  = "ACGTACGT".getBytes();
        final byte[] bases2 = "TGCATGCA".getBytes();
        final Kmer eq1 = new Kmer(bases, 0, 3);
        final Kmer eq2 = new Kmer(bases2, 0, 3);
        final int[] differingIndices = new int[bases2.length];
        final byte[] differingBases = new byte[bases2.length];

        Assert.assertEquals(eq1.getDifferingPositions(eq2, 1, differingIndices, differingBases), -1);
    }

    @Test
    public void testDifferingPositions() {
        final String bases = "ACGTCAGACGTACGTTTGACGTCAGACGTACGT";
        final Kmer baseKmer = new Kmer(bases.getBytes());


        final int NUM_TEST_CASES = 30;

        for (int test = 0; test < NUM_TEST_CASES; test++) {

            final int numBasesToChange =  test % bases.length();

            // changes numBasesToChange bases - spread regularly through read string
            final int step = (numBasesToChange > 0? Math.min(bases.length() / numBasesToChange, 1) : 1);

            final byte[] newBases = bases.getBytes().clone();
            int actualChangedBases =0; // could be different from numBasesToChange due to roundoff
            for (int idx=0; idx < numBasesToChange; idx+=step) {
                // now change given positions
                newBases[idx] = (newBases[idx] == (byte)'A'? (byte)'T':(byte)'A');
                actualChangedBases++;
            }

            // compute changed positions
            final int[] differingIndices = new int[newBases.length];
            final byte[] differingBases = new byte[newBases.length];
            final int numDiffs = baseKmer.getDifferingPositions(new Kmer(newBases),newBases.length,differingIndices,differingBases);
            Assert.assertEquals(numDiffs, actualChangedBases);
            for (int k=0; k < numDiffs; k++) {
                final int idx = differingIndices[k];
                Assert.assertTrue(newBases[idx] != bases.getBytes()[idx]);
                Assert.assertEquals(differingBases[idx], newBases[idx]);
            }
        }
    }
}
