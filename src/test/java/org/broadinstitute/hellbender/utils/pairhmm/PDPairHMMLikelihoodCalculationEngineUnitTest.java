package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PDPairHMMNativeArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Checks the Java and native partially determined PairHMMs against likelihoods recorded from a DRAGEN-mode
 * HaplotypeCaller run, and against each other.
 */
public class PDPairHMMLikelihoodCalculationEngineUnitTest extends GATKBaseTest  {

    static final String DRAGEN_GATK_TEST_ASSERT_FILE = largeFileTestDir + "expected.PDHMM.hmmresults.txt";

    /** Largest allowed difference in log10 likelihood between the native and Java PDPairHMMs, which both compute in double precision. */
    static final double NATIVE_TOLERANCE = 1e-10;

    /** Largest allowed difference from a recorded log10 likelihood, set by the precision they were recorded at. */
    static final double RECORDED_TOLERANCE = 1e-5;

    private static final String SAMPLE = "sample";

    /**
     * A recorded read and the partially determined haplotypes it was scored against, scored in a single PairHMM call
     * as the HaplotypeCaller does.
     */
    private static final class RecordedRead {
        final GATKRead read;
        final List<PartiallyDeterminedHaplotype> haplotypes;

        /** Per-base gap continuation penalties, which DRAGEN mode derives from the read's STR context. */
        final byte[] gapContinuationPenalties;

        /** Recorded log10 likelihood of the read against each haplotype. */
        final double[] likelihoods;

        RecordedRead(final GATKRead read, final List<PartiallyDeterminedHaplotype> haplotypes, final byte[] gapContinuationPenalties,
                     final double[] likelihoods) {
            this.read = read;
            this.haplotypes = haplotypes;
            this.gapContinuationPenalties = gapContinuationPenalties;
            this.likelihoods = likelihoods;
        }

        LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> computeLikelihoods(final PDPairHMM hmm) {
            final List<GATKRead> reads = Collections.singletonList(read);
            final AlleleLikelihoods<GATKRead, PartiallyDeterminedHaplotype> result = new AlleleLikelihoods<>(
                    SampleList.singletonSampleList(SAMPLE), new IndexedAlleleList<>(haplotypes), Collections.singletonMap(SAMPLE, reads));
            hmm.initialize(read.getLength(), haplotypes.stream().mapToInt(Haplotype::length).max().getAsInt());
            hmm.computeLog10Likelihoods(result.sampleMatrix(0), reads, this::scores, -1);
            return result.sampleMatrix(0);
        }

        private PairHMMInputScoreImputation scores(final GATKRead read) {
            return new PairHMMInputScoreImputation() {
                @Override
                public byte[] delOpenPenalties() {
                    return ReadUtils.getBaseDeletionQualities(read);
                }

                @Override
                public byte[] insOpenPenalties() {
                    return ReadUtils.getBaseInsertionQualities(read);
                }

                @Override
                public byte[] gapContinuationPenalties() {
                    return gapContinuationPenalties;
                }
            };
        }
    }

    @Test
    public void testJavaPDPairHMMMatchesRecordedLikelihoods() throws IOException {
        try (final PDPairHMM hmm = makeJavaHmm()) {
            for (final RecordedRead recorded : readRecordedReads()) {
                assertMatchesRecorded(recorded.computeLikelihoods(hmm), recorded, "Java PDPairHMM");
            }
        }
    }

    @Test
    public void testNativePDPairHMMMatchesRecordedLikelihoods() throws IOException {
        try (final PDPairHMM hmm = makeNativeHmm()) {
            for (final RecordedRead recorded : readRecordedReads()) {
                assertMatchesRecorded(recorded.computeLikelihoods(hmm), recorded, "native PDPairHMM");
            }
        }
    }

    @Test
    public void testNativePDPairHMMMatchesJava() throws IOException {
        try (final PDPairHMM nativeHmm = makeNativeHmm();
             final PDPairHMM javaHmm = makeJavaHmm()) {
            for (final RecordedRead recorded : readRecordedReads()) {
                final LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> nativeLikelihoods = recorded.computeLikelihoods(nativeHmm);
                final LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> javaLikelihoods = recorded.computeLikelihoods(javaHmm);
                for (int haplotype = 0; haplotype < recorded.haplotypes.size(); haplotype++) {
                    Assert.assertEquals(nativeLikelihoods.get(haplotype, 0), javaLikelihoods.get(haplotype, 0), NATIVE_TOLERANCE,
                                        "native PDPairHMM differs from the Java PDPairHMM for " + describe(recorded, haplotype));
                }
            }
        }
    }

    private static PDPairHMM makeJavaHmm() {
        return PDPairHMM.Implementation.LOGLESS_CACHING.makeNewHMM(PDPairHMMNativeArgumentCollection.getDefaultPDPairHMMArgs());
    }

    /**
     * Creates the native PDPairHMM, or skips the calling test when it cannot be loaded on this machine.
     */
    private static PDPairHMM makeNativeHmm() {
        try {
            return PDPairHMM.Implementation.AVX_LOGLESS_CACHING.makeNewHMM(PDPairHMMNativeArgumentCollection.getDefaultPDPairHMMArgs());
        } catch (final UserException.HardwareFeatureException e) {
            throw new SkipException("The native PDPairHMM is not available on this machine");
        }
    }

    private static void assertMatchesRecorded(final LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> likelihoods,
                                              final RecordedRead recorded, final String hmmDescription) {
        for (int haplotype = 0; haplotype < recorded.haplotypes.size(); haplotype++) {
            // A recorded -Infinity marks a pair whose likelihood the caller skipped rather than computed.
            if (recorded.likelihoods[haplotype] != Double.NEGATIVE_INFINITY) {
                Assert.assertEquals(likelihoods.get(haplotype, 0), recorded.likelihoods[haplotype], RECORDED_TOLERANCE,
                                    hmmDescription + " differs from the recorded likelihood for " + describe(recorded, haplotype));
            }
        }
    }

    private static String describe(final RecordedRead recorded, final int haplotypeIndex) {
        final PartiallyDeterminedHaplotype haplotype = recorded.haplotypes.get(haplotypeIndex);
        return "read " + recorded.read.getBasesString() + " vs haplotype " + haplotype.getBaseString() + " with PD bases " + Arrays.toString(haplotype.getAlternateBases());
    }

    /**
     * Reads the recorded file, one tab-separated line per pair: haplotype bases, its PD bases as a bracketed list,
     * read bases, base, insertion and deletion qualities and gap continuation penalties in FASTQ encoding, and the
     * log10 likelihood. Consecutive lines scoring the same read against distinct haplotypes form one
     * {@link RecordedRead}.
     */
    private static List<RecordedRead> readRecordedReads() throws IOException {
        final List<RecordedRead> recordedReads = new ArrayList<>();
        List<PartiallyDeterminedHaplotype> haplotypes = new ArrayList<>();
        List<Double> likelihoods = new ArrayList<>();
        GATKRead read = null;
        String readKey = null;
        byte[] gapContinuationPenalties = null;
        for (final String line : Files.readAllLines(Paths.get(DRAGEN_GATK_TEST_ASSERT_FILE))) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            final String[] tokens = line.split("\t");
            final PartiallyDeterminedHaplotype haplotype = recordedHaplotype(tokens[0], tokens[1]);
            final String key = String.join("\t", Arrays.copyOfRange(tokens, 2, 7));
            if (!key.equals(readKey) || haplotypes.contains(haplotype)) {
                if (read != null) {
                    recordedReads.add(new RecordedRead(read, haplotypes, gapContinuationPenalties, toArray(likelihoods)));
                }
                haplotypes = new ArrayList<>();
                likelihoods = new ArrayList<>();
                readKey = key;
                read = recordedRead(tokens);
                gapContinuationPenalties = SAMUtils.fastqToPhred(tokens[6]);
            }
            haplotypes.add(haplotype);
            likelihoods.add(Double.parseDouble(tokens[7]));
        }
        recordedReads.add(new RecordedRead(read, haplotypes, gapContinuationPenalties, toArray(likelihoods)));
        return recordedReads;
    }

    /**
     * Builds a partially determined haplotype from its bases and PD bases. Only those two arrays reach the PairHMMs,
     * so the events and position are placeholders.
     */
    private static PartiallyDeterminedHaplotype recordedHaplotype(final String bases, final String bracketedPDBases) {
        final byte[] haplotypeBases = bases.getBytes();
        final String[] pdTokens = bracketedPDBases.substring(1, bracketedPDBases.length() - 1).split(",");
        final byte[] pdBases = new byte[pdTokens.length];
        for (int i = 0; i < pdTokens.length; i++) {
            pdBases[i] = Byte.parseByte(pdTokens[i].trim());
        }
        final Haplotype base = new Haplotype(haplotypeBases, new SimpleInterval("1", 1, haplotypeBases.length));
        final Cigar cigar = new Cigar(Collections.singletonList(new CigarElement(haplotypeBases.length, CigarOperator.M)));
        return new PartiallyDeterminedHaplotype(base, pdBases, Collections.emptyList(), Collections.emptySet(), cigar, 1, Collections.emptyList(), 0);
    }

    private static GATKRead recordedRead(final String[] tokens) {
        final byte[] bases = tokens[2].getBytes();
        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, SAMUtils.fastqToPhred(tokens[3]), bases.length + "M");
        ReadUtils.setInsertionBaseQualities(read, SAMUtils.fastqToPhred(tokens[4]));
        ReadUtils.setDeletionBaseQualities(read, SAMUtils.fastqToPhred(tokens[5]));
        return read;
    }

    private static double[] toArray(final List<Double> values) {
        return values.stream().mapToDouble(Double::doubleValue).toArray();
    }
}
