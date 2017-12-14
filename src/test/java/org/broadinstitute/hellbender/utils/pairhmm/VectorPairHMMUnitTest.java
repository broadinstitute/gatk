package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import picard.util.BasicInputParser;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.ImmutablePair;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.*;

public final class VectorPairHMMUnitTest extends GATKBaseTest {

    private static final String pairHMMTestData = publicTestDir + "pairhmm-testdata.txt";

    // Return a list of supported VectorLoglessPairHMM implementations, skip the test if none are supported
    private List<Pair<PairHMM, Boolean> > getHMMs() {
        List<Pair<PairHMM, Boolean> > list = new ArrayList<>();
        PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.useDoublePrecision = false;
        args.maxNumberOfThreads = 1;

        for (VectorLoglessPairHMM.Implementation imp : VectorLoglessPairHMM.Implementation.values()) {
            boolean loaded = true;
            PairHMM avxPairHMM = null;
            try {
                avxPairHMM = new VectorLoglessPairHMM(imp, args);
                //avxPairHMM.doNotUseTristateCorrection();
            }
            catch (UserException.HardwareFeatureException e ) {
                loaded = false;
            }

            final Pair<PairHMM, Boolean> hmm_load = new ImmutablePair<PairHMM, Boolean>(avxPairHMM, new Boolean(loaded));
            list.add(hmm_load);
        }

        return list;
    }

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "JustHMMProvider")
    public Object[][] makeJustHMMProvider() {
        List<Object[]> tests = new ArrayList<>();

        for ( final Pair<PairHMM, Boolean> hmm_load : getHMMs() ) {
            tests.add(new Object[]{hmm_load.getLeft(), hmm_load.getRight()});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "JustHMMProvider")
    public void testLikelihoodsFromHaplotypes(final PairHMM hmm, Boolean loaded){

        // skip if not loaded
        if(!loaded.booleanValue()) {
            throw new SkipException("AVX PairHMM is not supported on this system or the library is not available");
        }

        BasicInputParser parser = null;
        try {
            parser = new BasicInputParser(true, new FileInputStream(pairHMMTestData));
        } catch (FileNotFoundException e) {
            Assert.fail("PairHMM test data not found : " + pairHMMTestData);
        }

        while (parser.hasNext()) {
            String tokens[] = parser.next();

            final Haplotype hap = new Haplotype(tokens[0].getBytes(), true);

            final byte[] bases = tokens[1].getBytes();
            final byte[] baseQuals = normalize(tokens[2].getBytes(), 6);
            final byte[] insertionQuals = normalize(tokens[3].getBytes());
            final byte[] deletionQuals = normalize(tokens[4].getBytes());
            final byte[] gcp = normalize(tokens[5].getBytes());
            final double expectedResult = Double.parseDouble(tokens[6]);

            final int readLength = bases.length;
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, baseQuals, readLength + "M");
            ReadUtils.setInsertionBaseQualities(read, insertionQuals);
            ReadUtils.setDeletionBaseQualities(read, deletionQuals);

            final Map<GATKRead,byte[]> gpcs = new LinkedHashMap<>(readLength);
            gpcs.put(read, gcp);

            hmm.initialize(Arrays.asList(hap), null, 0, 0);
            hmm.computeLog10Likelihoods(matrix(Arrays.asList(hap)), Arrays.asList(read), gpcs);

            final double[] la = hmm.getLogLikelihoodArray();

            Assert.assertEquals(la[0], expectedResult, 1e-5, "Likelihood not in expected range.");
        }

        hmm.close();
    }

    private static byte[] normalize(byte[] scores) {
        return normalize(scores, 0);
    }

    private static byte[] normalize(byte[] scores, int min) {
        for (int i = 0; i < scores.length; i++) {
            scores[i] -= 33;
            scores[i] = scores[i] < min ? (byte)min : scores[i];
        }
        return scores;
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
                return haplotypes.size();
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

}
