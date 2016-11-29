package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/*
 * Classic likelihood computation: full pair-hmm all haplotypes vs all reads.
 */
public final class PairHMMLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {

    private static final Logger logger = LogManager.getLogger(PairHMMLikelihoodCalculationEngine.class);

    private static final int MAX_STR_UNIT_LENGTH = 8;
    private static final int MAX_REPEAT_LENGTH   = 20;
    private static final int MIN_ADJUSTED_QSCORE = 10;

    @VisibleForTesting
    static final double INITIAL_QSCORE = 40.0;

    private final byte constantGCP;

    private final double log10globalReadMismappingRate;

    private final PairHMM pairHMM;

    @VisibleForTesting
    static boolean writeLikelihoodsToFile = false;

    public static final String LIKELIHOODS_FILENAME = "likelihoods.txt";
    private final PrintStream likelihoodsStream;

    public enum PCRErrorModel {
        /** no specialized PCR error model will be applied; if base insertion/deletion qualities are present they will be used */
        NONE(0.0),
        /** a most aggressive model will be applied that sacrifices true positives in order to remove more false positives */
        HOSTILE(1.0),
        /** a more aggressive model will be applied that sacrifices true positives in order to remove more false positives */
        AGGRESSIVE(2.0),
        /** a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives */
        CONSERVATIVE(3.0);

        private final double rateFactor;

        /** rate factor is applied to the PCR error model.  Can be 0.0 to imply no correction */
        PCRErrorModel(final double rateFactor) {
            this.rateFactor = rateFactor;
        }
        public double getRateFactor() { return rateFactor; }
        public boolean hasRateFactor() { return rateFactor != 0.0; }
    }

    private final PCRErrorModel pcrErrorModel;
    
    private final byte baseQualityScoreThreshold;

    /**
     * The expected rate of random sequencing errors for a read originating from its true haplotype.
     *
     * For example, if this is 0.01, then we'd expect 1 error per 100 bp.
     */
    private static final double EXPECTED_ERROR_RATE_PER_BASE = 0.02;
    
    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP the gap continuation penalty to use with the PairHMM
     * @param hmmType the type of the HMM to use
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param pcrErrorModel model to correct for PCR indel artifacts
     */
    public PairHMMLikelihoodCalculationEngine(final byte constantGCP,
                                              final PairHMM.Implementation hmmType,
                                              final double log10globalReadMismappingRate,
                                              final PCRErrorModel pcrErrorModel) {
        this( constantGCP, hmmType, log10globalReadMismappingRate, pcrErrorModel, PairHMM.BASE_QUALITY_SCORE_THRESHOLD );
    }

    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP the gap continuation penalty to use with the PairHMM
     * @param hmmType the type of the HMM to use
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param pcrErrorModel model to correct for PCR indel artifacts
     * @param baseQualityScoreThreshold Base qualities below this threshold will be reduced to the minimum usable base
     *                                  quality.
     */
    public PairHMMLikelihoodCalculationEngine(final byte constantGCP,
                                              final PairHMM.Implementation hmmType,
                                              final double log10globalReadMismappingRate,
                                              final PCRErrorModel pcrErrorModel,
                                              final byte baseQualityScoreThreshold) {
        Utils.nonNull(hmmType, "hmmType is null");
        Utils.nonNull(pcrErrorModel, "pcrErrorModel is null");
        if (constantGCP < 0){
            throw new IllegalArgumentException("gap continuation penalty must be non-negative");
        }
        if (log10globalReadMismappingRate > 0){
            throw new IllegalArgumentException("log10globalReadMismappingRate must be negative");
        }
        this.constantGCP = constantGCP;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.pcrErrorModel = pcrErrorModel;
        pairHMM = hmmType.makeNewHMM();

        initializePCRErrorModel();

        likelihoodsStream = makeLikelihoodStream();

        if (baseQualityScoreThreshold < QualityUtils.MIN_USABLE_Q_SCORE) {
            throw new IllegalArgumentException("baseQualityScoreThreshold must be greater than or equal to " + QualityUtils.MIN_USABLE_Q_SCORE + " (QualityUtils.MIN_USABLE_Q_SCORE)");
        }
        this.baseQualityScoreThreshold = baseQualityScoreThreshold;
    }

    private PrintStream makeLikelihoodStream() {
        try {
            return writeLikelihoodsToFile ? new PrintStream(new FileOutputStream(new File(LIKELIHOODS_FILENAME))) : null;
        } catch ( final FileNotFoundException e ) {
            throw new GATKException("can't open a file to write likelihoods to", e);
        }
    }

    @Override
    public void close() {
        if ( likelihoodsStream != null ) {
            likelihoodsStream.close();
        }
        pairHMM.close();
    }

    @Override
    public ReadLikelihoods<Haplotype> computeReadLikelihoods( final AssemblyResultSet assemblyResultSet, final SampleList samples, final Map<String, List<GATKRead>> perSampleReadList ) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        initializePairHMM(haplotypeList, perSampleReadList);

        // Add likelihoods for each sample's reads to our result
        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i));
        }

        result.normalizeLikelihoods(false, log10globalReadMismappingRate);
        result.filterPoorlyModeledReads(EXPECTED_ERROR_RATE_PER_BASE);
        return result;
    }

    /**
     * Creates a new GATKRead with the source read's header, read group and mate
     * information, but with the following fields set to user-supplied values:
     *  - Read Bases
     *  - Base Qualities
     *  - Base Insertion Qualities
     *  - Base Deletion Qualities
     *
     *  Cigar string is empty (not-null)
     *
     * Use this method if you want to create a new GATKRead based on
     * another GATKRead, but with modified bases and qualities
     *
     * @param read a read to copy the header from
     * @param readBases an array containing the new bases you wish use in place of the originals
     * @param baseQualities an array containing the new base qualities you wish use in place of the originals
     * @param baseInsertionQualities an array containing the new base insertion qaulities
     * @param baseDeletionQualities an array containing the new base deletion qualities
     * @return a read with modified bases and qualities, safe for the GATK
     */
    private static GATKRead createQualityModifiedRead(final GATKRead read,
                                                      final byte[] readBases,
                                                      final byte[] baseQualities,
                                                      final byte[] baseInsertionQualities,
                                                      final byte[] baseDeletionQualities) {
        Utils.validateArg( baseQualities.length == readBases.length && baseInsertionQualities.length == readBases.length && baseDeletionQualities.length == readBases.length,
                () -> String.format("Read bases and read quality arrays aren't the same size: Bases: %d vs Base Q's: %d vs Insert Q's: %d vs Delete Q's: %d.",
                        readBases.length, baseQualities.length, baseInsertionQualities.length, baseDeletionQualities.length));

        final GATKRead processedRead = ReadUtils.emptyRead(read);
        processedRead.setBases(readBases);
        processedRead.setBaseQualities(baseQualities);
        ReadUtils.setInsertionBaseQualities(processedRead, baseInsertionQualities);
        ReadUtils.setDeletionBaseQualities(processedRead, baseDeletionQualities);
        return processedRead;
    }

    /**
     * Initialize our pairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     *
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *
     * @param haplotypes a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    private void initializePairHMM(final List<Haplotype> haplotypes, final Map<String, List<GATKRead>> perSampleReadList) {
        final int readMaxLength = perSampleReadList.entrySet().stream().flatMap(e -> e.getValue().stream()).mapToInt(read -> read.getLength()).max().orElse(0);
        final int haplotypeMaxLength = haplotypes.stream().mapToInt(h -> h.getBases().length).max().orElse(0);

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        pairHMM.initialize(haplotypes, perSampleReadList, readMaxLength, haplotypeMaxLength);
    }

    private void computeReadLikelihoods(final LikelihoodMatrix<Haplotype> likelihoods) {
        // Modify the read qualities by applying the PCR error model and capping the minimum base,insertion,deletion qualities
        final List<GATKRead> processedReads = modifyReadQualities(likelihoods.reads());

        final Map<GATKRead, byte[]> gapContinuationPenalties = buildGapContinuationPenalties(processedReads, constantGCP);

        // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
        pairHMM.computeLog10Likelihoods(likelihoods, processedReads, gapContinuationPenalties);

        writeDebugLikelihoods(likelihoods);
    }

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */
    private List<GATKRead> modifyReadQualities(final List<GATKRead> reads) {
        final List<GATKRead> result = new ArrayList<>(reads.size());

        for (final GATKRead read : reads) {
            final byte[] readBases = read.getBases();

            // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
            //Using close here is justified - it's an array of primitives.
            final byte[] readQuals = read.getBaseQualities().clone();
            final byte[] readInsQuals = ReadUtils.getBaseInsertionQualities(read).clone();
            final byte[] readDelQuals = ReadUtils.getBaseDeletionQualities(read).clone();

            applyPCRErrorModel(readBases, readInsQuals, readDelQuals);
            capMinimumReadQualities(read, readQuals, readInsQuals, readDelQuals, baseQualityScoreThreshold);

            // Create a new copy of the read and sets its base qualities to the modified versions.
            result.add(createQualityModifiedRead(read, readBases, readQuals, readInsQuals, readDelQuals));
        }
        return result;
    }

    private static void capMinimumReadQualities(final GATKRead read, final byte[] readQuals, final byte[] readInsQuals, final byte[] readDelQuals, final byte baseQualityScoreThreshold) {
        for( int i = 0; i < readQuals.length; i++ ) {
            readQuals[i] = (byte) Math.min(0xff & readQuals[i], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
            readQuals[i] =    setToFixedValueIfTooLow( readQuals[i],    baseQualityScoreThreshold,             QualityUtils.MIN_USABLE_Q_SCORE );
            readInsQuals[i] = setToFixedValueIfTooLow( readInsQuals[i], QualityUtils.MIN_USABLE_Q_SCORE,       QualityUtils.MIN_USABLE_Q_SCORE );
            readDelQuals[i] = setToFixedValueIfTooLow( readDelQuals[i], QualityUtils.MIN_USABLE_Q_SCORE,       QualityUtils.MIN_USABLE_Q_SCORE );
        }
    }

    private static byte setToFixedValueIfTooLow(final byte currentVal, final byte minQual, final byte fixedQual){
        return currentVal < minQual ? fixedQual : currentVal;
    }

    private static Map<GATKRead, byte[]> buildGapContinuationPenalties(final List<GATKRead> reads, final byte gapPenalty) {
        return reads.stream().collect(Collectors.toMap(read -> read,
                                                       read -> Utils.dupBytes(gapPenalty, read.getLength())));
    }

    private void writeDebugLikelihoods(final LikelihoodMatrix<Haplotype> likelihoods) {
        if (!writeLikelihoodsToFile || likelihoodsStream == null) {
            return;
        }

        final List<GATKRead> reads = likelihoods.reads();
        final List<Haplotype> haplotypes = likelihoods.alleles();
        for (int i = 0; i < reads.size(); i++) {
            for (int j = 0; j < haplotypes.size(); j++) {
                writeDebugLikelihoods(reads.get(i), haplotypes.get(j), likelihoods.get(j, i));
            }
        }
        likelihoodsStream.flush();
    }

    private void writeDebugLikelihoods(final GATKRead processedRead, final Haplotype haplotype, final double log10l){
        // Note: the precision of log10l in the debug output is only ~6 digits (ie., not all digits are necessarily printed)
        likelihoodsStream.printf("%s %s %s %s %s %s %f%n",
                haplotype.getBaseString(),
                new String(processedRead.getBases()),
                SAMUtils.phredToFastq(processedRead.getBaseQualities()),
                SAMUtils.phredToFastq(ReadUtils.getBaseInsertionQualities(processedRead)),
                SAMUtils.phredToFastq(ReadUtils.getBaseDeletionQualities(processedRead)),
                SAMUtils.phredToFastq(constantGCP),
                log10l);
    }

    /* --------------------------------------------------------------------------------
    *
    * Experimental attempts at PCR error rate modeling
    *
    -------------------------------------------------------------------------------- */

    private byte[] pcrIndelErrorModelCache;

    private void initializePCRErrorModel() {
        if ( !pcrErrorModel.hasRateFactor() ) {
            return;
        }

        pcrIndelErrorModelCache = new byte[MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel.getRateFactor();

        for( int i = 0; i <= MAX_REPEAT_LENGTH; i++ ) {
            pcrIndelErrorModelCache[i] = getErrorModelAdjustedQual(i, rateFactor);
        }
    }

    static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor) {
        return (byte) Math.max(MIN_ADJUSTED_QSCORE, MathUtils.fastRound(INITIAL_QSCORE - Math.exp(repeatLength / (rateFactor * Math.PI)) + 1.0));
    }

    @VisibleForTesting
    void applyPCRErrorModel( final byte[] readBases, final byte[] readInsQuals, final byte[] readDelQuals ) {
        if ( pcrErrorModel == PCRErrorModel.NONE ) {
            return;
        }

        for ( int i = 1; i < readBases.length; i++ ) {
            final int repeatLength = findTandemRepeatUnits(readBases, i-1).getRight();
            readInsQuals[i-1] = (byte) Math.min(0xff & readInsQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
            readDelQuals[i-1] = (byte) Math.min(0xff & readDelQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        }
    }

    @VisibleForTesting
    static Pair<byte[], Integer> findTandemRepeatUnits(final byte[] readBases, final int offset) {
        int maxBW = 0;
        byte[] bestBWRepeatUnit = {readBases[offset]};
        for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
            // fix repeat unit length
            //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
            if (offset+1-str < 0) {
                break;
            }

            // get backward repeat unit and # repeats
            maxBW = GATKVariantContextUtils.findNumberOfRepetitions(readBases, offset - str + 1,  str , readBases, 0, offset + 1, false);
            if (maxBW > 1) {
                bestBWRepeatUnit = Arrays.copyOfRange(readBases, offset - str + 1, offset + 1);
                break;
            }
        }
        byte[] bestRepeatUnit = bestBWRepeatUnit;
        int maxRL = maxBW;

        if (offset < readBases.length-1) {
            byte[] bestFWRepeatUnit = {readBases[offset+1]};
            int maxFW = 0;

            for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
                // fix repeat unit length
                //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
                if (offset+str+1 > readBases.length) {
                    break;
                }

                // get forward repeat unit and # repeats
                maxFW = GATKVariantContextUtils.findNumberOfRepetitions(readBases, offset + 1, str, readBases, offset + 1, readBases.length-offset -1, true);
                if (maxFW > 1) {
                    bestFWRepeatUnit = Arrays.copyOfRange(readBases, offset + 1, offset+str+1);
                    break;
                }
            }
            // if FW repeat unit = BW repeat unit it means we're in the middle of a tandem repeat - add FW and BW components
            if (Arrays.equals(bestFWRepeatUnit, bestBWRepeatUnit)) {
                maxRL = maxBW + maxFW;
                bestRepeatUnit = bestFWRepeatUnit; // arbitrary
            } else {
                // tandem repeat starting forward from current offset.
                // It could be the case that best BW unit was different from FW unit, but that BW still contains FW unit.
                // For example, TTCTT(C) CCC - at (C) place, best BW unit is (TTC)2, best FW unit is (C)3.
                // but correct representation at that place might be (C)4.
                // Hence, if the FW and BW units don't match, check if BW unit can still be a part of FW unit and add
                // representations to total
                final byte[] testString = Arrays.copyOfRange(readBases, 0, offset + 1);
                maxBW = GATKVariantContextUtils.findNumberOfRepetitions(bestFWRepeatUnit, testString, false);
                maxRL = maxFW + maxBW;
                bestRepeatUnit = bestFWRepeatUnit;
            }
        }

        if(maxRL > MAX_REPEAT_LENGTH) {
            maxRL = MAX_REPEAT_LENGTH;
        }
        return Pair.of(bestRepeatUnit, maxRL);
    }
}
