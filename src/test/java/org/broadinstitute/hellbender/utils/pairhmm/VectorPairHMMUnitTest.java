package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.RandomSequenceTestUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.StandardPairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Checks that the native PairHMM implementations agree with the Java {@link LoglessPairHMM}, and that both agree
 * with previously recorded likelihoods, on recorded read/haplotype pairs and on randomly generated ones.
 */
public final class VectorPairHMMUnitTest extends GATKBaseTest {

    /** Largest allowed difference in log10 likelihood between the Java PairHMM and a single-precision native one. */
    public static final double FLOAT_TOLERANCE = 1e-5;

    /** Largest allowed difference in log10 likelihood between the Java PairHMM and a double-precision native one. */
    public static final double DOUBLE_TOLERANCE = 1e-5;

    /** Largest allowed difference from a recorded log10 likelihood, set by the precision they were recorded at. */
    public static final double RECORDED_TOLERANCE = 1e-5;

    private static final String SAMPLE = "sample";
    private static final int RANDOM_BATCH_COUNT = 50;

    /** The read/haplotype sets the PairHMMs are compared on. */
    private enum Dataset {
        /** Pairs and likelihoods recorded when the native PairHMM was first added to GATK. */
        PAIR_HMM_TEST_DATA,
        /** Pairs and likelihoods recorded from a HaplotypeCaller run with the Java PairHMM. */
        HAPLOTYPE_CALLER_RESULTS,
        /** Random haplotype sets with reads drawn from them; nothing recorded. */
        RANDOM;

        List<Batch> batches() throws IOException {
            switch (this) {
                case PAIR_HMM_TEST_DATA:       return readRecordedBatches(publicTestDir + "pairhmm-testdata.txt");
                case HAPLOTYPE_CALLER_RESULTS: return readRecordedBatches(toolsTestDir + "haplotypecaller/expected.Java.hmmresults.txt");
                default:                       return randomBatches();
            }
        }
    }

    /**
     * Haplotypes and reads whose likelihoods are computed in a single PairHMM call, as the HaplotypeCaller does for
     * each active region.
     */
    private static final class Batch {
        final List<Haplotype> haplotypes;
        final List<GATKRead> reads;
        final byte gapContinuationPenalty;

        /** Recorded log10 likelihoods indexed by read then haplotype, or null when none were recorded. */
        final double[] recorded;

        Batch(final List<Haplotype> haplotypes, final List<GATKRead> reads, final byte gapContinuationPenalty, final double[] recorded) {
            this.haplotypes = haplotypes;
            this.reads = reads;
            this.gapContinuationPenalty = gapContinuationPenalty;
            this.recorded = recorded;
        }

        LikelihoodMatrix<GATKRead, Haplotype> computeLikelihoods(final PairHMM hmm) {
            final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(
                    SampleList.singletonSampleList(SAMPLE), new IndexedAlleleList<>(haplotypes), Collections.singletonMap(SAMPLE, reads));
            final int maxReadLength = reads.stream().mapToInt(GATKRead::getLength).max().getAsInt();
            final int maxHaplotypeLength = haplotypes.stream().mapToInt(Haplotype::length).max().getAsInt();
            hmm.initialize(haplotypes, Collections.singletonMap(SAMPLE, reads), maxReadLength, maxHaplotypeLength);
            hmm.computeLog10Likelihoods(likelihoods.sampleMatrix(0), reads, StandardPairHMMInputScoreImputator.newInstance(gapContinuationPenalty));
            return likelihoods.sampleMatrix(0);
        }

        double recorded(final int readIndex, final int haplotypeIndex) {
            return recorded[readIndex * haplotypes.size() + haplotypeIndex];
        }
    }

    @DataProvider
    public Object[][] recordedDatasets() {
        return new Object[][]{
                {Dataset.PAIR_HMM_TEST_DATA},
                {Dataset.HAPLOTYPE_CALLER_RESULTS},
        };
    }

    @DataProvider
    public Object[][] nativeImplementationsAndDatasets() {
        final List<Object[]> cases = new ArrayList<>();
        for (final VectorLoglessPairHMM.Implementation implementation : VectorLoglessPairHMM.Implementation.values()) {
            for (final boolean doublePrecision : new boolean[]{false, true}) {
                for (final Dataset dataset : Dataset.values()) {
                    cases.add(new Object[]{implementation, doublePrecision, dataset});
                }
            }
        }
        return cases.toArray(new Object[0][]);
    }

    @Test(dataProvider = "recordedDatasets")
    public void testJavaPairHMMMatchesRecordedLikelihoods(final Dataset dataset) throws IOException {
        try (final PairHMM hmm = new LoglessPairHMM()) {
            for (final Batch batch : dataset.batches()) {
                assertMatchesRecorded(batch.computeLikelihoods(hmm), batch, RECORDED_TOLERANCE, "Java PairHMM on " + dataset);
            }
        }
    }

    @Test(dataProvider = "nativeImplementationsAndDatasets")
    public void testNativePairHMMMatchesJava(final VectorLoglessPairHMM.Implementation implementation, final boolean doublePrecision,
                                             final Dataset dataset) throws IOException {
        final double tolerance = doublePrecision ? DOUBLE_TOLERANCE : FLOAT_TOLERANCE;
        final String hmmDescription = implementation + (doublePrecision ? " double" : " float") + " PairHMM on " + dataset;
        try (final PairHMM nativeHmm = makeNativeHmm(implementation, doublePrecision);
             final PairHMM javaHmm = new LoglessPairHMM()) {
            for (final Batch batch : dataset.batches()) {
                final LikelihoodMatrix<GATKRead, Haplotype> nativeLikelihoods = batch.computeLikelihoods(nativeHmm);
                final LikelihoodMatrix<GATKRead, Haplotype> javaLikelihoods = batch.computeLikelihoods(javaHmm);
                for (int read = 0; read < batch.reads.size(); read++) {
                    for (int haplotype = 0; haplotype < batch.haplotypes.size(); haplotype++) {
                        Assert.assertEquals(nativeLikelihoods.get(haplotype, read), javaLikelihoods.get(haplotype, read), tolerance,
                                            hmmDescription + " differs from the Java PairHMM for " + describe(batch, read, haplotype));
                    }
                }
                if (batch.recorded != null) {
                    assertMatchesRecorded(nativeLikelihoods, batch, RECORDED_TOLERANCE, hmmDescription);
                }
            }
        }
    }

    /**
     * Creates the requested native PairHMM, or skips the calling test when it cannot be loaded on this machine.
     */
    private static PairHMM makeNativeHmm(final VectorLoglessPairHMM.Implementation implementation, final boolean doublePrecision) {
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.useDoublePrecision = doublePrecision;
        args.maxNumberOfThreads = 1;
        try {
            return new VectorLoglessPairHMM(implementation, args);
        } catch (final UserException.HardwareFeatureException e) {
            throw new SkipException("Native PairHMM implementation " + implementation + " is not available on this machine");
        }
    }

    private static void assertMatchesRecorded(final LikelihoodMatrix<GATKRead, Haplotype> likelihoods, final Batch batch,
                                              final double tolerance, final String hmmDescription) {
        for (int read = 0; read < batch.reads.size(); read++) {
            for (int haplotype = 0; haplotype < batch.haplotypes.size(); haplotype++) {
                Assert.assertEquals(likelihoods.get(haplotype, read), batch.recorded(read, haplotype), tolerance,
                                    hmmDescription + " differs from the recorded likelihood for " + describe(batch, read, haplotype));
            }
        }
    }

    private static String describe(final Batch batch, final int readIndex, final int haplotypeIndex) {
        return "read " + batch.reads.get(readIndex).getBasesString() + " vs haplotype " + batch.haplotypes.get(haplotypeIndex).getBaseString();
    }

    /**
     * Reads recorded pairs, one per line as "haplotype read baseQuals insertionQuals deletionQuals gcp log10Likelihood"
     * with FASTQ-encoded qualities. Consecutive lines that score the same read against distinct haplotypes become one
     * batch, so that the recorded files also exercise the multi-haplotype path.
     */
    private static List<Batch> readRecordedBatches(final String path) throws IOException {
        final List<Batch> batches = new ArrayList<>();
        List<Haplotype> haplotypes = new ArrayList<>();
        List<Double> likelihoods = new ArrayList<>();
        GATKRead read = null;
        String readKey = null;
        byte gapContinuationPenalty = 0;
        for (final String line : Files.readAllLines(Paths.get(path))) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            final String[] tokens = line.trim().split("\\s+");
            final Haplotype haplotype = new Haplotype(tokens[0].getBytes());
            final String key = String.join(" ", Arrays.copyOfRange(tokens, 1, 6));
            if (!key.equals(readKey) || haplotypes.contains(haplotype)) {
                if (read != null) {
                    batches.add(new Batch(haplotypes, Collections.singletonList(read), gapContinuationPenalty, toArray(likelihoods)));
                }
                haplotypes = new ArrayList<>();
                likelihoods = new ArrayList<>();
                readKey = key;
                read = recordedRead(tokens);
                gapContinuationPenalty = constantQuality(SAMUtils.fastqToPhred(tokens[5]));
            }
            haplotypes.add(haplotype);
            likelihoods.add(Double.parseDouble(tokens[6]));
        }
        batches.add(new Batch(haplotypes, Collections.singletonList(read), gapContinuationPenalty, toArray(likelihoods)));
        return batches;
    }

    /**
     * Builds the read on a recorded line. Base qualities are floored at {@link QualityUtils#MIN_USABLE_Q_SCORE} as the
     * likelihood engines do before calling the PairHMM, which is how the recorded likelihoods were computed.
     */
    private static GATKRead recordedRead(final String[] tokens) {
        final byte[] bases = tokens[1].getBytes();
        final byte[] baseQualities = SAMUtils.fastqToPhred(tokens[2]);
        for (int i = 0; i < baseQualities.length; i++) {
            baseQualities[i] = (byte) Math.max(baseQualities[i], QualityUtils.MIN_USABLE_Q_SCORE);
        }
        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, baseQualities, bases.length + "M");
        ReadUtils.setInsertionBaseQualities(read, SAMUtils.fastqToPhred(tokens[3]));
        ReadUtils.setDeletionBaseQualities(read, SAMUtils.fastqToPhred(tokens[4]));
        return read;
    }

    private static byte constantQuality(final byte[] qualities) {
        for (final byte quality : qualities) {
            Utils.validate(quality == qualities[0], "gap continuation penalties must be constant within a read");
        }
        return qualities[0];
    }

    private static double[] toArray(final List<Double> values) {
        return values.stream().mapToDouble(Double::doubleValue).toArray();
    }

    /**
     * Random haplotype sets with reads drawn from them, spanning the read lengths, haplotype counts and quality ranges
     * the HaplotypeCaller produces.
     */
    private static List<Batch> randomBatches() {
        Utils.resetRandomGenerator();
        final Random rng = Utils.getRandomGenerator();
        final List<Batch> batches = new ArrayList<>();
        for (int i = 0; i < RANDOM_BATCH_COUNT; i++) {
            final byte[] reference = RandomSequenceTestUtils.randomBases(rng, 150 + rng.nextInt(250));
            final List<Haplotype> haplotypes = new ArrayList<>();
            haplotypes.add(new Haplotype(reference, true));
            final int haplotypeCount = 1 + rng.nextInt(8);
            while (haplotypes.size() < haplotypeCount) {
                final Haplotype haplotype = new Haplotype(RandomSequenceTestUtils.withRandomEdits(rng, reference, 1 + rng.nextInt(3), 8));
                if (!haplotypes.contains(haplotype)) {
                    haplotypes.add(haplotype);
                }
            }
            final int minHaplotypeLength = haplotypes.stream().mapToInt(Haplotype::length).min().getAsInt();
            final List<GATKRead> reads = new ArrayList<>();
            final int readCount = 1 + rng.nextInt(20);
            for (int r = 0; r < readCount; r++) {
                reads.add(randomRead(rng, haplotypes.get(rng.nextInt(haplotypes.size())), minHaplotypeLength - 50));
            }
            batches.add(new Batch(haplotypes, reads, (byte) (5 + rng.nextInt(11)), null));
        }
        return batches;
    }

    /**
     * A read copied from a random slice of {@code source}, never longer than {@code maxLength} plus the inserted
     * bases. Four reads in five carry up to two sequencing errors. The rest carry up to twelve, which spreads the
     * likelihoods down through the range where single precision underflows.
     */
    private static GATKRead randomRead(final Random rng, final Haplotype source, final int maxLength) {
        final int editCount = rng.nextInt(5) == 0 ? rng.nextInt(13) : rng.nextInt(3);
        // Heavily edited reads start longer so that deletions cannot consume them
        final int minLength = editCount > 2 ? 60 : 20;
        final int length = Math.min(minLength + rng.nextInt(251 - minLength), maxLength);
        final int start = rng.nextInt(source.length() - length + 1);
        final byte[] slice = Arrays.copyOfRange(source.getBases(), start, start + length);
        final byte[] bases = RandomSequenceTestUtils.withRandomEdits(rng, slice, editCount, 4);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, randomQualities(rng, bases.length, 18, 45, 6, 17), bases.length + "M");
        ReadUtils.setInsertionBaseQualities(read, randomQualities(rng, bases.length, 45, 45, 6, 50));
        ReadUtils.setDeletionBaseQualities(read, randomQualities(rng, bases.length, 45, 45, 6, 50));
        return read;
    }

    /**
     * Qualities drawn uniformly from [{@code usualMin}, {@code usualMax}], except that one in ten is drawn from
     * [{@code rareMin}, {@code rareMax}] instead.
     */
    private static byte[] randomQualities(final Random rng, final int length, final int usualMin, final int usualMax,
                                          final int rareMin, final int rareMax) {
        final byte[] qualities = new byte[length];
        for (int i = 0; i < length; i++) {
            final boolean rare = rng.nextInt(10) == 0;
            final int min = rare ? rareMin : usualMin;
            final int max = rare ? rareMax : usualMax;
            qualities[i] = (byte) (min + rng.nextInt(max - min + 1));
        }
        return qualities;
    }
}
