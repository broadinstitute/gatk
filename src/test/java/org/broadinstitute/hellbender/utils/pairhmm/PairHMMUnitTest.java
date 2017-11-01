package org.broadinstitute.hellbender.utils.pairhmm;

import com.google.common.base.Strings;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class PairHMMUnitTest extends GATKBaseTest {
    private final static boolean ALLOW_READS_LONGER_THAN_HAPLOTYPE = true;
    private final static boolean DEBUG = false;
    final static boolean EXTENSIVE_TESTING = true;
    final N2MemoryPairHMM exactHMM = new Log10PairHMM(true); // the log truth implementation
    final N2MemoryPairHMM originalHMM = new Log10PairHMM(false); // the reference implementation
    final N2MemoryPairHMM loglessHMM = new LoglessPairHMM();

    @BeforeClass
    public void initialize() {
        exactHMM.doNotUseTristateCorrection();
        originalHMM.doNotUseTristateCorrection();
        loglessHMM.doNotUseTristateCorrection();
    }

    private List<N2MemoryPairHMM> getHMMs() {
        return Arrays.asList(exactHMM, originalHMM, loglessHMM);
    }

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class BasicLikelihoodTestProvider {
        final String ref, nextRef, read;
        final byte[] refBasesWithContext, nextRefBasesWithContext, readBasesWithContext;
        final int baseQual, insQual, delQual, gcp;
        final int expectedQual;
        final boolean left, right;
        final static String CONTEXT = "ACGTAATGACGATTGCA";
        final static String LEFT_FLANK = "GATTTATCATCGAGTCTGC";
        final static String RIGHT_FLANK = "CATGGATCGTTATCAGCTATCTCGAGGGATTCACTTAACAGTTTTA";

        public BasicLikelihoodTestProvider(final String ref, final String nextRef, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp ) {
            this(ref, nextRef, read, baseQual, insQual, delQual, expectedQual, gcp, false, false);
        }

        public BasicLikelihoodTestProvider(final String ref, final String nextRef, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp, final boolean left, final boolean right) {
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
            this.gcp = gcp;
            this.read = read;
            this.ref = ref;
            this.nextRef = nextRef;
            this.expectedQual = expectedQual;
            this.left = left;
            this.right = right;

            refBasesWithContext = asBytes(ref, left, right);
            nextRefBasesWithContext = asBytes(nextRef, left, right);
            readBasesWithContext = asBytes(read, false, false);
        }

        @Override
        public String toString() {
            return String.format("ref=%s nextRef=%s read=%s b/i/d/c quals = %d/%d/%d/%d l/r flank = %b/%b e[qual]=%d", ref, nextRef, read, baseQual, insQual, delQual, gcp, left, right, expectedQual);
        }

        public double expectedLogLikelihood() {
            final double log10MagicConstant = 0.03;
            return (expectedQual / -10.0) + 0.03 + Math.log10(1.0 /
                    refBasesWithContext.length);
        }

        public double getTolerance(final PairHMM hmm) {
            if ( hmm instanceof Log10PairHMM) {
                return ((Log10PairHMM)hmm).isDoingExactLog10Calculations() ? toleranceFromExact() : toleranceFromReference();
            } else
                return toleranceFromTheoretical();
        }

        // NOTE: natural log units rescale log10 units by a factor of log(10), so we rescale tolerances by this factor
        public double toleranceFromTheoretical() {
            return 0.2;
        }

        public double toleranceFromReference() {
            return 1E-3; // has to be very tolerant -- this approximation is quite approximate
        }

        public double toleranceFromExact() {
            return 1E-9;
        }

        public double calcLog10Likelihood(final PairHMM pairHMM, boolean anchorIndel) {
            pairHMM.initialize(readBasesWithContext.length, refBasesWithContext.length);
            return pairHMM.computeReadLikelihoodGivenHaplotypeLog10(
                    refBasesWithContext, readBasesWithContext,
                    qualAsBytes(baseQual, false, anchorIndel), qualAsBytes(insQual, true, anchorIndel), qualAsBytes(delQual, true, anchorIndel),
                    qualAsBytes(gcp, false, anchorIndel), true, nextRefBasesWithContext);
        }

        private byte[] asBytes(final String bases, final boolean left, final boolean right) {
            if(bases == null)
                return null;
            else
                return ( (left ? LEFT_FLANK : "") + CONTEXT + bases + CONTEXT + (right ? RIGHT_FLANK : "")).getBytes();
        }

        private byte[] qualAsBytes(final int phredQual, final boolean doGOP, final boolean anchorIndel) {
            final byte phredQuals[] = new byte[readBasesWithContext.length];

            if( anchorIndel ) {
                // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
                Arrays.fill(phredQuals, (byte) 100);

                // update just the bases corresponding to the provided micro read with the quality scores
                if( doGOP ) {
                    phredQuals[CONTEXT.length()] = (byte)phredQual;
                } else {
                    for ( int i = 0; i < read.length(); i++)
                        phredQuals[i + CONTEXT.length()] = (byte)phredQual;
                }
            } else {
                Arrays.fill(phredQuals, (byte) phredQual);
            }

            return phredQuals;
        }
    }

    @DataProvider(name = "BasicLikelihoodTestProvider")
    public Object[][] makeBasicLikelihoodTests() {
        // context on either side is ACGTTGCA REF ACGTTGCA
        // test all combinations
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30, 40, 50) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(20, 30, 40, 50) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(8, 10, 20) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(2, 3, 4, 5, 7, 8, 9, 10, 20, 30, 35) : Arrays.asList(2);

        final List<Object[]> tests = new ArrayList<>();

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {

                    // test substitutions
                    for ( final byte refBase : BaseUtils.BASES ) {
                        for ( final byte readBase : BaseUtils.BASES ) {
                            final String ref  = new String(new byte[]{refBase});
                            final String read = new String(new byte[]{readBase});
                            final int expected = refBase == readBase ? 0 : baseQual;
                            // runBasicLikelihoodTests uses calcLogLikelihood(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                            tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp)});
                        }
                    }

                    // test insertions and deletions
                    for ( final int size : sizes ) {
                        for ( final byte base : BaseUtils.BASES ) {
                            final int expected = indelQual + (size - 2) * gcp;

                            for ( boolean insertionP : Arrays.asList(true, false)) {
                                final String small = Utils.dupChar((char) base, 1);
                                final String big = Utils.dupChar((char) base, size);

                                final String ref = insertionP ? small : big;
                                final String read = insertionP ? big : small;

                                // runBasicLikelihoodTests uses calcLogLikelihood(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp)});
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp, true, false)});
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp, false, true)});
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp, true, true)});
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "OptimizedLikelihoodTestProvider")
    public Object[][] makeOptimizedLikelihoodTests() {
        Utils.resetRandomGenerator();
        final Random random = Utils.getRandomGenerator();
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 30, 40, 60) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(20, 40, 60) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(3, 20, 50, 90, 160) : Arrays.asList(2);

        final List<Object[]> tests = new ArrayList<>();

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {
                    for ( final int refSize : sizes ) {
                        for ( final int readSize : sizes ) {
                            String ref = "";
                            String read = "";
                            for( int iii = 0; iii < refSize; iii++) {
                                ref += (char) BaseUtils.BASES[random.nextInt(4)];
                            }
                            for( int iii = 0; iii < readSize; iii++) {
                                read += (char) BaseUtils.BASES[random.nextInt(4)];
                            }

                            for ( final boolean leftFlank : Arrays.asList(true, false) )
                                for ( final boolean rightFlank : Arrays.asList(true, false) )
                                    // runOptimizedLikelihoodTests uses calcLogLikelihood(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                                    tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, -0, gcp, leftFlank, rightFlank)});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "BasicLikelihoodTestProvider")
    public void testBasicLikelihoods(BasicLikelihoodTestProvider cfg) {
        if ( ALLOW_READS_LONGER_THAN_HAPLOTYPE || cfg.read.length() <= cfg.ref.length() ) {
            final double exactLogL = cfg.calcLog10Likelihood(exactHMM, true);

            for ( final PairHMM hmm : getHMMs() ) {
                final double actualLogL = cfg.calcLog10Likelihood(hmm, true);
                final double expectedLogL = cfg.expectedLogLikelihood();

                // compare to our theoretical expectation with appropriate tolerance
                Assert.assertEquals(actualLogL, expectedLogL, cfg.toleranceFromTheoretical(), "Failed with hmm " + hmm);

                // compare to the exact reference implementation with appropriate tolerance
                Assert.assertEquals(actualLogL, exactLogL, cfg.getTolerance(hmm), "Failed with hmm " + hmm);
                Assert.assertTrue(MathUtils.goodLog10Probability(actualLogL), "Bad log likelihood " + actualLogL);
            }
        }
    }

    @Test(enabled = !DEBUG, dataProvider = "OptimizedLikelihoodTestProvider")
    public void testOptimizedLikelihoods(BasicLikelihoodTestProvider cfg) {
        if ( ALLOW_READS_LONGER_THAN_HAPLOTYPE || cfg.read.length() <= cfg.ref.length() ) {
            final double exactLogL = cfg.calcLog10Likelihood(exactHMM, false);

            for ( final PairHMM hmm : getHMMs() ) {
                final double calculatedLogL = cfg.calcLog10Likelihood(hmm, false);
                // compare to the exact reference implementation with appropriate tolerance
                Assert.assertEquals(calculatedLogL, exactLogL, cfg.getTolerance(hmm), String.format("Test: logL calc=%.2f expected=%.2f for %s with hmm %s", calculatedLogL, exactLogL, cfg.toString(), hmm));
                Assert.assertTrue(MathUtils.goodLog10Probability(calculatedLogL), "Bad log10 likelihood " + calculatedLogL);
            }
        }
    }

    @Test(enabled = !DEBUG)
    public void testMismatchInEveryPositionInTheReadWithCenteredHaplotype() {
        final byte[] haplotype1 = "TTCTCTTCTGTTGTGGCTGGTT".getBytes();
        final byte matchQual = 90;
        final byte mismatchQual = 20;
        final byte indelQual = 80;
        final int offset = 2;
        final byte[] gop = new byte[haplotype1.length - 2 * offset];
        Arrays.fill(gop, indelQual);
        final byte[] gcp = new byte[haplotype1.length - 2 * offset];
        Arrays.fill(gcp, indelQual);
        loglessHMM.initialize(gop.length, haplotype1.length);

        for( int k = 0; k < haplotype1.length - 2 * offset; k++ ) {
            final byte[] quals = new byte[haplotype1.length - 2 * offset];
            Arrays.fill(quals, matchQual);
            // one base mismatches the haplotype
            quals[k] = mismatchQual;
            final byte[] mread = Arrays.copyOfRange(haplotype1,offset,haplotype1.length-offset);
            // change single base at position k to C. If it's a C, change to T
            mread[k] = ( mread[k] == (byte)'C' ? (byte)'T' : (byte)'C');
            final double res1 = loglessHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotype1, mread, quals, gop, gop, gcp, true, null);
            final double expected = Math.log10(1.0 / haplotype1.length * Math.pow(QualityUtils.qualToProb(matchQual), mread.length - 1) * QualityUtils.qualToErrorProb(mismatchQual));
            Assert.assertEquals(res1, expected, 1e-2);
        }
    }


    @Test(enabled = ! DEBUG)
    public void testMismatchInEveryPositionInTheRead() {
        final byte[] haplotype1 = "TTCTCTTCTGTTGTGGCTGGTT".getBytes();
        final byte matchQual = 90;
        final byte mismatchQual = 20;
        final byte indelQual = 80;

        final int offset = 2;
        final byte[] gop = new byte[haplotype1.length - offset];
        Arrays.fill(gop, indelQual);
        final byte[] gcp = new byte[haplotype1.length - offset];
        Arrays.fill(gcp, indelQual);
        loglessHMM.initialize(gop.length, haplotype1.length);

        for( int k = 0; k < haplotype1.length - offset; k++ ) {
            final byte[] quals = new byte[haplotype1.length - offset];
            Arrays.fill(quals, matchQual);
            // one base mismatches the haplotype with low qual
            quals[k] = mismatchQual;
            final byte[] mread = Arrays.copyOfRange(haplotype1,offset,haplotype1.length);
            // change single base at position k to C. If it's a C, change to T
            mread[k] = ( mread[k] == (byte)'C' ? (byte)'T' : (byte)'C');
            final double res1 = loglessHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotype1, mread, quals, gop, gop, gcp, true, null);
            final double expected = Math.log10(1.0 / haplotype1.length * Math.pow(QualityUtils.qualToProb(matchQual), mread.length - 1) * QualityUtils.qualToErrorProb(mismatchQual));
            Assert.assertEquals(res1, expected, 1e-2);
        }
    }

    @DataProvider(name = "HMMProvider")
    public Object[][] makeHMMProvider() {
        List<Object[]> tests = new ArrayList<>();

        for ( final int readSize : Arrays.asList(1, 2, 5, 10) ) {
            for ( final int refSize : Arrays.asList(1, 2, 5, 10) ) {
                if ( refSize > readSize ) {
                    for ( final PairHMM hmm : getHMMs() )
                        tests.add(new Object[]{hmm, readSize, refSize});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "HMMProvider")
    void testMultipleReadMatchesInHaplotype(final PairHMM hmm, final int readSize, final int refSize) {
        final byte[] readBases =  Utils.dupBytes((byte)'A', readSize);
        final byte[] refBases = ("CC" + new String(Utils.dupBytes((byte)'A', refSize)) + "GGA").getBytes();
        final byte baseQual = 20;
        final byte insQual = 37;
        final byte delQual = 37;
        final byte gcp = 10;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        final double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
        Assert.assertTrue(d <= 0.0, "Likelihoods should be <= 0 but got " + d);
    }

    @DataProvider(name = "HMMProviderSimple")
    public Object[][] HMMProviderSimple() {
        List<Object[]> tests = new ArrayList<>();

        for ( final int readSize : Arrays.asList(1, 2, 5, 10) ) {
            for ( final PairHMM hmm : getHMMs() )
                tests.add(new Object[]{hmm, readSize});
        }

        return tests.toArray(new Object[][]{});
    }
    @Test(enabled = !DEBUG, dataProvider = "HMMProviderSimple")
    void testReadSameAsHaplotype(final PairHMM hmm, final int readSize) {
        final byte[] readBases =  Utils.dupBytes((byte)'A', readSize);
        final byte[] refBases = readBases;
        final byte baseQual = 20;
        final byte insQual = 37;
        final byte delQual = 37;
        final byte gcp = 10;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        final double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
        Assert.assertTrue(d <= 0.0, "Likelihoods should be <= 0 but got " + d);
    }

    @Test(enabled = !DEBUG, dataProvider = "HMMProvider")
    void testAllMatchingRead(final PairHMM hmm, final int readSize, final int refSize) {
        final byte[] readBases =  Utils.dupBytes((byte)'A', readSize);
        final byte[] refBases = Utils.dupBytes((byte) 'A', refSize);
        final byte baseQual = 20;
        final byte insQual = 100;
        final byte delQual = 100;
        final byte gcp = 100;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
        double expected = getExpectedMatchingLogLikelihood(readBases, refBases, baseQual, insQual);
        Assert.assertEquals(d, expected, 1e-3, "Likelihoods should sum to just the error prob of the read " + String.format("readSize=%d refSize=%d", readSize, refSize));
    }

    private static double getExpectedMatchingLogLikelihood(final byte[] readBases, final byte[] refBases, final byte baseQual, final byte insQual) {
        double expected =  0;
        final double initialCondition = ((double) Math.abs(refBases.length - readBases.length + 1))/refBases.length;
        if (readBases.length < refBases.length) {
            expected = Math.log10(initialCondition * Math.pow(QualityUtils.qualToProb(baseQual), readBases.length));
        } else if (readBases.length > refBases.length) {
            expected = Math.log10(initialCondition * Math.pow(QualityUtils.qualToProb(baseQual), refBases.length) * Math.pow(QualityUtils.qualToErrorProb(insQual), readBases.length - refBases.length));
        }
        return expected;
    }

    @DataProvider(name = "HMMProviderWithBigReads")
    public Object[][] makeBigReadHMMProvider() {
        List<Object[]> tests = new ArrayList<>();

        final String read1 = "ACCAAGTAGTCACCGT";
        final String ref1  = "ACCAAGTAGTCACCGTAACG";

        for ( final int nReadCopies : Arrays.asList(1, 2, 10, 20, 50) ) {
            for ( final int nRefCopies : Arrays.asList(1, 2, 10, 20, 100) ) {
                if ( nRefCopies > nReadCopies ) {
                    for ( final PairHMM hmm : getHMMs() ) {
                        final String read = Strings.repeat(read1, nReadCopies);
                        final String ref  = Strings.repeat(ref1, nRefCopies);
                        tests.add(new Object[]{hmm, read, ref});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "HMMProviderWithBigReads")
    void testReallyBigReads(final PairHMM hmm, final String read, final String ref) {
        final byte[] readBases =  read.getBytes();
        final byte[] refBases = ref.getBytes();
        final byte baseQual = 30;
        final byte insQual = 40;
        final byte delQual = 40;
        final byte gcp = 10;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
    }

    @Test(enabled = !DEBUG)
    void testPreviousBadValue() {
        final byte[] readBases = "A".getBytes();
        final byte[] refBases =  "AT".getBytes();
        final byte baseQual = 30;
        final byte insQual = 40;
        final byte delQual = 40;
        final byte gcp = 10;
        // running HMMs with no haplotype caching. Should therefore pass null in place of nextRef bases
        exactHMM.initialize(readBases.length, refBases.length);
        exactHMM.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);

        loglessHMM.initialize(readBases.length, refBases.length);
        loglessHMM.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
    }

    @DataProvider(name = "JustHMMProvider")
    public Object[][] makeJustHMMProvider() {
        List<Object[]> tests = new ArrayList<>();

        for ( final PairHMM hmm : getHMMs() ) {
            tests.add(new Object[]{hmm});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "JustHMMProvider")
    void testMaxLengthsBiggerThanProvidedRead(final PairHMM hmm) {
        final byte[] readBases = "CTATCTTAGTAAGCCCCCATACCTGCAAATTTCAGGATGTCTCCTCCAAAAATCAACA".getBytes();
        final byte[] refBases =  "CTATCTTAGTAAGCCCCCATACCTGCAAATTTCAGGATGTCTCCTCCAAAAATCAAAACTTCTGAGAAAAAAAAAAAAAATTAAATCAAACCCTGATTCCTTAAAGGTAGTAAAAAAACATCATTCTTTCTTAGTGGAATAGAAACTAGGTCAAAAGAACAGTGATTC".getBytes();

        final byte[] quals = new byte[]{35,34,31,32,35,34,32,31,36,30,31,32,36,34,33,32,32,32,33,32,30,35,33,35,36,36,33,33,33,32,32,32,37,33,36,35,33,32,34,31,36,35,35,35,35,33,34,31,31,30,28,27,26,29,26,25,29,29};
        final byte[] insQual = new byte[]{46,46,46,46,46,47,45,46,45,48,47,44,45,48,46,43,43,42,48,48,45,47,47,48,48,47,48,45,38,47,45,39,47,48,47,47,48,46,49,48,49,48,46,47,48,44,44,43,39,32,34,36,46,48,46,44,45,45};
        final byte[] delQual = new byte[]{44,44,44,43,45,44,43,42,45,46,45,43,44,47,45,40,40,40,45,46,43,45,45,44,46,46,46,43,35,44,43,36,44,45,46,46,44,44,47,43,47,45,45,45,46,45,45,46,44,35,35,35,45,47,45,44,44,43};
        final byte[] gcp = Utils.dupBytes((byte) 10, delQual.length);
        hmm.initialize(readBases.length + 100, refBases.length + 100);
        for ( int nExtraMaxSize = 0; nExtraMaxSize < 100; nExtraMaxSize++ ) {
            // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
            hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases, quals, insQual, delQual, gcp, true, null);
        }
    }

    @Test(dataProvider = "JustHMMProvider")
    public void testLikelihoodsFromHaplotypes(final PairHMM hmm){
        final int readSize= 10;
        final int refSize = 20;
        final byte[] readBases =  Utils.dupBytes((byte)'A', readSize);
        final byte[] refBases = Utils.dupBytes((byte) 'A', refSize);
        final byte baseQual = 20;
        final byte insQual = 100;
        final byte gcp = 100;

        final Haplotype refH= new Haplotype(refBases, true);

        final byte[] readQuals= Utils.dupBytes(baseQual, readBases.length);
        final List<GATKRead> reads = Arrays.asList(ArtificialReadUtils.createArtificialRead(readBases, readQuals, readBases.length + "M"));
        final Map<GATKRead, byte[]> gpcs = buildGapContinuationPenalties(reads, gcp);

        hmm.computeLog10Likelihoods(matrix(Arrays.asList(refH)), Collections.emptyList(), gpcs);
        Assert.assertEquals(hmm.getLogLikelihoodArray(), null);

        hmm.computeLog10Likelihoods(matrix(Arrays.asList(refH)), reads, gpcs);
        final double expected = getExpectedMatchingLogLikelihood(readBases, refBases, baseQual, insQual);
        final double[] la = hmm.getLogLikelihoodArray();

        Assert.assertEquals(la.length, 1);
        Assert.assertEquals(la[0], expected, 1e-3, "Likelihoods should sum to just the error prob of the read " + String.format("readSize=%d refSize=%d", readSize, refSize));

    }

    private LikelihoodMatrix<Haplotype> matrix(final List<Haplotype> haplotypes) {
        return new LikelihoodMatrix<Haplotype>() {
            @Override
            public List<GATKRead> reads() {
                throw new UnsupportedOperationException();
            }

            @Override
            public List<Haplotype> alleles() {
                return haplotypes;
            }

            @Override
            public void set(int alleleIndex, int readIndex, double value) {
//                throw new UnsupportedOperationException();
            }

            @Override
            public double get(int alleleIndex, int readIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int indexOfAllele(Haplotype allele) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int indexOfRead(GATKRead read) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int numberOfAlleles() {
                throw new UnsupportedOperationException();
            }

            @Override
            public int numberOfReads() {
                throw new UnsupportedOperationException();
            }

            @Override
            public Haplotype getAllele(int alleleIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public GATKRead getRead(int readIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void copyAlleleLikelihoods(int alleleIndex, double[] dest, int offset) {
                throw new UnsupportedOperationException();
            }
        };
    }

    private static Map<GATKRead, byte[]> buildGapContinuationPenalties(final List<GATKRead> processedReads, final byte gcp) {
        final Map<GATKRead,byte[]> result = new LinkedHashMap<>(processedReads.size());
        for (final GATKRead read : processedReads) {
            final byte[] readGcpArray = new byte[read.getLength()];
            Arrays.fill(readGcpArray,gcp);
            result.put(read,readGcpArray);
        }
        return result;
    }
    
    @DataProvider(name = "HaplotypeIndexingProvider")
    public Object[][] makeHaplotypeIndexingProvider() {
        List<Object[]> tests = new ArrayList<>();

        // First difference (root2, root3) is the base position immediately following first difference (root1, root2)
        final String root1    = "ACGTGTCAAACCGGGTT";
        final String root2    = "ACGTGTCACACTGGGTT"; // differs in two locations from root1
        final String root3    = "ACGTGTCACTCCGCGTT"; // differs in two locations from root2

        final String read1    = "ACGTGTCACACTGGATT"; // 1 diff from 2, 2 diff from root1, 2 diff from root3
        final String read2    = root1; // same as root1
        final String read3    = root2; // same as root2
        final String read4    = "ACGTGTCACACTGGATTCGAT";
        final String read5    = "CCAGTAACGTGTCACACTGGATTCGAT";

//        for ( final String read : Arrays.asList(read2) ) {
        for ( final String read : Arrays.asList(read1, read2, read3, read4, read5) ) {
            for ( final PairHMM hmm : getHMMs() ) {
//                int readLength = read.length(); {
                for ( int readLength = 10; readLength < read.length(); readLength++ ) {
                    final String myRead = read.substring(0, readLength);
                    tests.add(new Object[]{hmm, root1, root2, root3, myRead});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "HaplotypeIndexingProvider")
    void testHaplotypeIndexing(final PairHMM hmm, final String root1, final String root2, final String root3, final String read) {
        final double TOLERANCE = 1e-9;
        final String prefix   = "AACCGGTTTTTGGGCCCAAACGTACGTACAGTTGGTCAACATCGATCAGGTTCCGGAGTAC";

        final int maxReadLength = read.length();
        final int maxHaplotypeLength = prefix.length() + root1.length();

        // the initialization occurs once, at the start of the evalution of reads
        hmm.initialize(maxReadLength, maxHaplotypeLength);

        for ( int prefixStart = prefix.length(); prefixStart >= 0; prefixStart-- ) {
            final String myPrefix = prefix.substring(prefixStart, prefix.length());
            final String hap1 = myPrefix + root1;
            final String hap2 = myPrefix + root2;
            final String hap3 = myPrefix + root3;

            final int hapStart = PairHMM.findFirstPositionWhereHaplotypesDiffer(hap1.getBytes(), hap2.getBytes());

            // Run the HMM on the first haplotype, peaking ahead the second, to set up caching
            // Then run on the second haplotype in both cached and uncached mode, and verify that results are the same
            // When evaluating actual2, it is important that we both apply old caching from hap1 and set up new caching for hap3, to ensure read/write operations do not cause conflicts
            final double actual1 = testHaplotypeIndexingCalc(hmm, hap1, hap2, read, 0, true);
            final double actual2 = testHaplotypeIndexingCalc(hmm, hap2, hap3, read, hapStart, false);
            final double expected2 = testHaplotypeIndexingCalc(hmm, hap2, null, read, 0, true);
            Assert.assertEquals(actual2, expected2, TOLERANCE, "HMM " + hmm.getClass() + " Caching calculation failed for read " + read + " against haplotype with prefix '" + myPrefix
                    + "' expected " + expected2 + " but got " + actual2 + " with hapStart of " + hapStart);
        }
    }

    private double testHaplotypeIndexingCalc(final PairHMM hmm, final String hap, final String nextHap, final String read, final int hapStart, final boolean recache) {
        final byte[] readBases = read.getBytes();
        // if not peaking ahead to capture info for a future cache run, the next haplotype will be null, and this should be passed to HMM
        final byte[] nextHapBases = nextHap == null ? null : nextHap.getBytes();
        final byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);
        final byte[] insQuals = Utils.dupBytes((byte)45, readBases.length);
        final byte[] delQuals = Utils.dupBytes((byte)40, readBases.length);
        final byte[] gcp = Utils.dupBytes((byte) 10, readBases.length);
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(hap.getBytes(), readBases, baseQuals, insQuals, delQuals, gcp, recache, nextHapBases);
        Assert.assertTrue(MathUtils.goodLog10Probability(d), "Likelihoods = " + d + " was bad for read " + read + " and ref " + hap + " with hapStart " + hapStart);
        return d;
    }

    @Test(enabled = !DEBUG)
    public void testFindFirstPositionWhereHaplotypesDiffer() {
        for ( int haplotypeSize1 = 10; haplotypeSize1 < 30; haplotypeSize1++ ) {
            for ( int haplotypeSize2 = 10; haplotypeSize2 < 50; haplotypeSize2++ ) {
                final int maxLength = Math.max(haplotypeSize1, haplotypeSize2);
                final int minLength = Math.min(haplotypeSize1, haplotypeSize2);
                for ( int differingSite = 0; differingSite < maxLength + 1; differingSite++) {
                    for ( final boolean oneIsDiff : Arrays.asList(true, false) ) {
                        final byte[] hap1 = Utils.dupBytes((byte)'A', haplotypeSize1);
                        final byte[] hap2 = Utils.dupBytes((byte)'A', haplotypeSize2);
                        final int expected = oneIsDiff
                                ? makeDiff(hap1, differingSite, minLength)
                                : makeDiff(hap2, differingSite, minLength);
                        final int actual = PairHMM.findFirstPositionWhereHaplotypesDiffer(hap1, hap2);
                        Assert.assertEquals(actual, expected, "Bad differing site for " + new String(hap1) + " vs. " + new String(hap2));
                    }
                }
            }
        }
    }

    private int makeDiff(final byte[] bytes, final int site, final int minSize) {
        if ( site < bytes.length ) {
            bytes[site] = 'C';
            return Math.min(site, minSize);
        } else
            return minSize;
    }

    @DataProvider(name = "UninitializedHMMs")
    public Object[][] makeUninitializedHMMs() {
        List<Object[]> tests = new ArrayList<>();

        final Log10PairHMM myLog10PairHMM = new Log10PairHMM(true);
        myLog10PairHMM.doNotUseTristateCorrection();
        tests.add(new Object[]{myLog10PairHMM});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, expectedExceptions = IllegalStateException.class, dataProvider = "UninitializedHMMs")
    public void testNoInitializeCall(final PairHMM hmm) {
        byte[] readBases = "A".getBytes();
        byte[] refBases =  "AT".getBytes();
        byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);

        // didn't call initialize => should exception out
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                baseQuals, baseQuals, baseQuals, baseQuals, true, null);
    }

    @Test(enabled = true, expectedExceptions = IllegalArgumentException.class, dataProvider = "JustHMMProvider")
    public void testHapTooLong(final PairHMM hmm) {
        byte[] readBases = "AAA".getBytes();
        byte[] refBases =  "AAAT".getBytes();
        byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);

        hmm.initialize(3, 3);
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                baseQuals, baseQuals, baseQuals, baseQuals, true, null);
    }

    @Test(enabled = true, expectedExceptions = IllegalArgumentException.class, dataProvider = "JustHMMProvider")
    public void testReadTooLong(final PairHMM hmm) {
        byte[] readBases = "AAA".getBytes();
        byte[] refBases =  "AAAT".getBytes();
        byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);

        hmm.initialize(2, 3);
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(refBases, readBases,
                baseQuals, baseQuals, baseQuals, baseQuals, true, null);
    }

    @Test
    public void dumpMatrices(){
        //doesn't test anything other than not-blowing up
        getHMMs().forEach(hmm -> hmm.initialize(3, 3));
        getHMMs().forEach(hmm -> hmm.dumpMatrices());
    }
}