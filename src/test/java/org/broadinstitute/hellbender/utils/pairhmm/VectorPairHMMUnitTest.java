package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.BasicInputParser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class VectorPairHMMUnitTest extends GATKBaseTest {

    private static final String pairHMMTestData = publicTestDir + "pairhmm-testdata.txt";

    // This method originally used a DataProvider to individually test each available implementation,
    // but was refactored to avoid resultant intermittent failures
    // (possibly caused by concurrency issues when loading libraries;
    // see https://github.com/broadinstitute/gatk/pull/5026#issuecomment-607596205).
    @Test
    public void testLikelihoodsFromHaplotypesForAvailableImplementations() {
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.useDoublePrecision = false;
        args.maxNumberOfThreads = 1;

        for (final VectorLoglessPairHMM.Implementation imp : VectorLoglessPairHMM.Implementation.values()) {
            PairHMM hmm;
            try {
                logger.warn(String.format("Trying PairHMM implementation %s...", imp.name()));
                hmm = new VectorLoglessPairHMM(imp, args);
                //hmm.doNotUseTristateCorrection();
            } catch (final UserException.HardwareFeatureException e ) {
                logger.warn(String.format("PairHMM implementation %s not available, skipping test...", imp.name()));
                continue;
            }

            BasicInputParser parser = null;
            try {
                parser = new BasicInputParser(true, new FileInputStream(pairHMMTestData));
            } catch (final FileNotFoundException e) {
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

                final Map<GATKRead, byte[]> gpcs = new LinkedHashMap<>(readLength);
                gpcs.put(read, gcp);

                hmm.initialize(Arrays.asList(hap), null, 0, 0);
                hmm.computeLog10Likelihoods(matrix(Arrays.asList(hap)), Arrays.asList(read), gpcs);

                final double[] la = hmm.getLogLikelihoodArray();

                Assert.assertEquals(la[0], expectedResult, 1e-5,
                        String.format("Likelihood not in expected range for PairHMM implementation: %s.", imp.name()));
            }

            hmm.close();
        }
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

    private LikelihoodMatrix<GATKRead, Haplotype> matrix(final List<Haplotype> haplotypes) {
        return new LikelihoodMatrix<GATKRead, Haplotype>() {
            @Override
            public List<GATKRead> evidence() {
                throw new UnsupportedOperationException();
            }

            @Override
            public List<Haplotype> alleles() {
                return haplotypes;
            }

            @Override
            public void set(int alleleIndex, int evidenceIndex, double value) {
//                throw new UnsupportedOperationException();
            }

            @Override
            public double get(int alleleIndex, int evidenceIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int indexOfAllele(Haplotype allele) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int indexOfEvidence(GATKRead evidence) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int numberOfAlleles() {
                return haplotypes.size();
            }

            @Override
            public int evidenceCount() {
                throw new UnsupportedOperationException();
            }

            @Override
            public Haplotype getAllele(int alleleIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public GATKRead getEvidence(int evidenceIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void copyAlleleLikelihoods(int alleleIndex, double[] dest, int offset) {
                throw new UnsupportedOperationException();
            }
        };
    }

}
