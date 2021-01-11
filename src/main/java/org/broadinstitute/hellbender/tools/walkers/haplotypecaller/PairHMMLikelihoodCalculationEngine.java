package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.ToDoubleFunction;

/*
 * Classic likelihood computation: full pair-hmm all haplotypes vs all reads.
 */
public final class PairHMMLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {

    static final double DEFAULT_DYNAMIC_DISQUALIFICATION_SCALE_FACTOR = 1.0;
    private static final Logger logger = LogManager.getLogger(PairHMMLikelihoodCalculationEngine.class);

    private static final int MAX_STR_UNIT_LENGTH = 8;
    private static final int MAX_REPEAT_LENGTH   = 20;
    private static final int MIN_ADJUSTED_QSCORE = 10;

    @VisibleForTesting
    static final double INITIAL_QSCORE = 40.0;
    public static final String HMM_BASE_QUALITIES_TAG = "HMMQuals";

    private final byte constantGCP;

    private final double log10globalReadMismappingRate;

    private final PairHMM pairHMM;

    // DRAGEN-GATK related parameters
    private final DragstrParams dragstrParams;
    private final boolean dynamicDisqualification;
    private final double readDisqualificationScale;
    private final double expectedErrorRatePerBase;
    private final boolean disableCapReadQualitiesToMapQ;
    private final boolean symmetricallyNormalizeAllelesToReference;
    private final boolean modifySoftclippedBases;

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
    public static final double DEFAULT_EXPECTED_ERROR_RATE_PER_BASE = 0.02;
    
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
                                              final DragstrParams dragstrParams,
                                              final PairHMMNativeArguments arguments,
                                              final PairHMM.Implementation hmmType,
                                              final double log10globalReadMismappingRate,
                                              final PCRErrorModel pcrErrorModel) {
        this( constantGCP, dragstrParams, arguments, hmmType, log10globalReadMismappingRate, pcrErrorModel, PairHMM.BASE_QUALITY_SCORE_THRESHOLD, false, DEFAULT_DYNAMIC_DISQUALIFICATION_SCALE_FACTOR, DEFAULT_EXPECTED_ERROR_RATE_PER_BASE, true, false, true);
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
                                              final DragstrParams dragstrParams,
                                              final PairHMMNativeArguments arguments,
                                              final PairHMM.Implementation hmmType,
                                              final double log10globalReadMismappingRate,
                                              final PCRErrorModel pcrErrorModel,
                                              final byte baseQualityScoreThreshold,
                                              final boolean dynamicReadDisqualificaiton,
                                              final double readDisqualificationScale,
                                              final double expectedErrorRatePerBase,
                                              final boolean symmetricallyNormalizeAllelesToReference,
                                              final boolean disableCapReadQualitiesToMapQ,
                                              final boolean modifySoftclippedBases) {
        Utils.nonNull(hmmType, "hmmType is null");
        Utils.nonNull(pcrErrorModel, "pcrErrorModel is null");
        if (constantGCP < 0){
            throw new IllegalArgumentException("gap continuation penalty must be non-negative");
        }
        if (log10globalReadMismappingRate > 0){
            throw new IllegalArgumentException("log10globalReadMismappingRate must be negative");
        }
        this.dragstrParams = dragstrParams;
        this.constantGCP = constantGCP;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.pcrErrorModel = this.dragstrParams == null ? pcrErrorModel : PCRErrorModel.NONE;
        this.pairHMM = hmmType.makeNewHMM(arguments);
        this.dynamicDisqualification = dynamicReadDisqualificaiton;
        this.readDisqualificationScale = readDisqualificationScale;
        this.symmetricallyNormalizeAllelesToReference = symmetricallyNormalizeAllelesToReference;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;
        this.disableCapReadQualitiesToMapQ = disableCapReadQualitiesToMapQ;
        this.modifySoftclippedBases = modifySoftclippedBases;

        initializePCRErrorModel();

        if (baseQualityScoreThreshold < QualityUtils.MIN_USABLE_Q_SCORE) {
            throw new IllegalArgumentException("baseQualityScoreThreshold must be greater than or equal to " + QualityUtils.MIN_USABLE_Q_SCORE + " (QualityUtils.MIN_USABLE_Q_SCORE)");
        }
        this.baseQualityScoreThreshold = baseQualityScoreThreshold;
    }

    @Override
    public void close() {
        pairHMM.close();
    }

    @Override
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods( final AssemblyResultSet assemblyResultSet, final SampleList samples, final Map<String, List<GATKRead>> perSampleReadList) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        initializePairHMM(haplotypeList, perSampleReadList);

        // Add likelihoods for each sample's reads to our result
        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
                computeReadLikelihoods(result.sampleMatrix(i));
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate, symmetricallyNormalizeAllelesToReference);

        // sato: disable filtering---for now.
        if (true){
            return result;
        }

        if (dynamicDisqualification) {
            result.filterPoorlyModeledEvidence(daynamicLog10MinLiklihoodModel(readDisqualificationScale, log10MinTrueLikelihood(expectedErrorRatePerBase, false)));
        } else {
            result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase, true)); // sato: 0.02 by default
        }
        return result;
    }

    private ToDoubleFunction<GATKRead> daynamicLog10MinLiklihoodModel(final double dynamicRadQualConstant, final ToDoubleFunction<GATKRead> log10MinTrueLikelihood) {
        return read -> {
            final double dynamicThreshold = calculateLog10DynamicReadQualThreshold(read, dynamicRadQualConstant);
            final double log10MaxLikelihoodForTrueAllele = log10MinTrueLikelihood.applyAsDouble(read);
            if (dynamicThreshold < log10MaxLikelihoodForTrueAllele ) {
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("For read "+ read.getName() + " replacing old threshold ("+log10MaxLikelihoodForTrueAllele+") with new threshold: "+dynamicThreshold);
                }
                return dynamicThreshold;
            } else {
                return log10MaxLikelihoodForTrueAllele;
            }
        };
    }

    private static double calculateLog10DynamicReadQualThreshold(final GATKRead read, final double dynamicReadQualConstant) {
        double sumMean = 0;
        double sumVariance = 0;

        final byte[] baseQualities = read.getOptionalTransientAttribute(HMM_BASE_QUALITIES_TAG, byte[].class)
                                         .orElseGet(read::getBaseQualities);

        for (final int qualByte : baseQualities) {
            final int bq = 0xFF & qualByte; // making sure that larger BQ are not casted into negatives.
            // bound the base qualities for lookup between 1 and 40
            final int entryIndex = bq <= 1 ? 0 : Math.min(MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ, bq) - 1;
            final int meanOffset = entryIndex * DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH + DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET;
            final int varOffset = meanOffset + 1;
            sumMean +=      dynamicReadQualThreshLookupTable[meanOffset];
            sumVariance +=  dynamicReadQualThreshLookupTable[varOffset];
        }

        final double threshold = sumMean + dynamicReadQualConstant * Math.sqrt(sumVariance);
        return QualityUtils.qualToErrorProbLog10(threshold); // = threshold * -.1;
    }

    // TODO i don't like having a lookup table be static like this, i would prefer this be computed at initialization (with the default values being saved as a test)
    // table used for disqualifying reads for genotyping
    // Format for each row of table: baseQ, mean, variance
    // Actual threshold is calculated over the length of the read as:
    // sum(means) + K * sqrt(sum(variances))
    private static double dynamicReadQualThreshLookupTable[] = {
            //baseQ,mean,variance
            1,  5.996842844, 0.196616587, 2,  5.870018422, 1.388545569, 3,  5.401558531, 5.641990128,
            4,  4.818940919, 10.33176216, 5,  4.218758304, 14.25799688, 6,  3.646319832, 17.02880749,
            7,  3.122346753, 18.64537883, 8,  2.654731979, 19.27521677, 9,  2.244479156, 19.13584613,
            10, 1.88893867,  18.43922003, 11, 1.583645342, 17.36842261, 12, 1.3233807, 16.07088712,
            13, 1.102785365, 14.65952563, 14, 0.916703025, 13.21718577, 15, 0.760361881, 11.80207947,
            16, 0.629457387, 10.45304833, 17, 0.520175654, 9.194183767, 18, 0.42918208,  8.038657241,
            19, 0.353590663, 6.991779595, 20, 0.290923699, 6.053379213, 21, 0.23906788,  5.219610436,
            22, 0.196230431, 4.484302033, 23, 0.160897421, 3.839943445, 24, 0.131795374, 3.27839108,
            25, 0.1078567,   2.791361596, 26, 0.088189063, 2.370765375, 27, 0.072048567, 2.008921719,
            28, 0.058816518, 1.698687797, 29, 0.047979438, 1.433525748, 30, 0.039111985, 1.207526336,
            31, 0.031862437, 1.015402928, 32, 0.025940415, 0.852465956, 33, 0.021106532, 0.714585285,
            34, 0.017163711, 0.598145851, 35, 0.013949904, 0.500000349, 36, 0.011332027, 0.41742159,
            37, 0.009200898, 0.348056286, 38, 0.007467036, 0.289881373, 39, 0.006057179, 0.241163527,
            40, 0.004911394, 0.200422214};

    private static final int MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ = 40;
    private static final int DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH = 3;
    private static final int DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET = 1;

    private ToDoubleFunction<GATKRead> log10MinTrueLikelihood(final double maximumErrorPerBase, final boolean capLikelihoods) {
        return read -> {
            // TODO this might be replaced by an explicit calculation
            final int qualifiedReadLength = read.getTransientAttribute(HMM_BASE_QUALITIES_TAG) != null ? ((byte[])read.getTransientAttribute(HMM_BASE_QUALITIES_TAG)).length : read.getLength();
            final double maxErrorsForRead = capLikelihoods ? Math.min(2.0, Math.ceil(qualifiedReadLength * maximumErrorPerBase)) : Math.ceil(qualifiedReadLength * maximumErrorPerBase);
            final double log10QualPerBase = -4.0;
            return maxErrorsForRead * log10QualPerBase;
        };
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
        final int readMaxLength = perSampleReadList.entrySet().stream().flatMap(e -> e.getValue().stream()).mapToInt(GATKRead::getLength).max().orElse(0);
        final int haplotypeMaxLength = haplotypes.stream().mapToInt(h -> h.getBases().length).max().orElse(0);

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        pairHMM.initialize(haplotypes, perSampleReadList, readMaxLength, haplotypeMaxLength);
    }

    private void computeReadLikelihoods(final LikelihoodMatrix<GATKRead, Haplotype> likelihoods) {
        // Modify the read qualities by applying the PCR error model and capping the minimum base,insertion,deletion qualities
        final List<GATKRead> processedReads = modifyReadQualities(likelihoods.evidence());

        for(int counter = 0; counter < processedReads.size(); counter++) {
            GATKRead read = processedReads.get(counter);
            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println("read "+counter +": "+read.getName()+" cigar: "+read.getCigar()+" mapQ: "+read.getMappingQuality()+" loc: ["+read.getStart() +"-"+ read.getEnd()+"] unclippedloc: ["+read.getUnclippedStart()+"-"+read.getUnclippedEnd()+"]");
                HaplotypeCallerGenotypingDebugger.println(Arrays.toString(read.getBaseQualitiesNoCopy()));
            }
        }
        // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
        pairHMM.computeLog10Likelihoods(likelihoods, processedReads, inputScoreImputator);
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
            final GATKRead maybeUnclipped = modifySoftclippedBases ? read : ReadClipper.hardClipSoftClippedBases(read);    //Clip the bases here to remove hap
            final byte[] readBases = maybeUnclipped.getBases();

            // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
            //Using close here is justified - it's an array of primitives.
            final byte[] readQuals = maybeUnclipped.getBaseQualities().clone();
            final byte[] readInsQuals = ReadUtils.getBaseInsertionQualities(maybeUnclipped).clone();
            final byte[] readDelQuals = ReadUtils.getBaseDeletionQualities(maybeUnclipped).clone();

            applyPCRErrorModel(readBases, readInsQuals, readDelQuals);
            capMinimumReadQualities(maybeUnclipped, readQuals, readInsQuals, readDelQuals, baseQualityScoreThreshold, disableCapReadQualitiesToMapQ);

            // Store the actual qualities
            read.setTransientAttribute(HMM_BASE_QUALITIES_TAG, readQuals);
            // Create a new copy of the read and sets its base qualities to the modified versions.
            result.add(createQualityModifiedRead(maybeUnclipped, readBases, readQuals, readInsQuals, readDelQuals));
        }
        return result;
    }

    private static void capMinimumReadQualities(final GATKRead read, final byte[] readQuals, final byte[] readInsQuals, final byte[] readDelQuals, final byte baseQualityScoreThreshold, final boolean disableCapReadQualitiesToMapQ) {
        for( int i = 0; i < readQuals.length; i++ ) {
            if (!disableCapReadQualitiesToMapQ) {
                readQuals[i] = (byte) Math.min(0xff & readQuals[i], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
            }
            readQuals[i] =    setToFixedValueIfTooLow( readQuals[i],    baseQualityScoreThreshold,             QualityUtils.MIN_USABLE_Q_SCORE );
            readInsQuals[i] = setToFixedValueIfTooLow( readInsQuals[i], QualityUtils.MIN_USABLE_Q_SCORE,       QualityUtils.MIN_USABLE_Q_SCORE );
            readDelQuals[i] = setToFixedValueIfTooLow( readDelQuals[i], QualityUtils.MIN_USABLE_Q_SCORE,       QualityUtils.MIN_USABLE_Q_SCORE );
        }
    }

    private static byte setToFixedValueIfTooLow(final byte currentVal, final byte minQual, final byte fixedQual){
        return currentVal < minQual ? fixedQual : currentVal;
    }

    private static Map<GATKRead, byte[]> buildGapContinuationPenalties(final List<GATKRead> reads, final byte gapPenalty) {
        final Map<GATKRead, byte[]> result = new HashMap<>(reads.size());
        reads.stream().forEach(read -> result.put(read, Utils.dupBytes(gapPenalty, read.getLength())));
        return result;
    }

    /* --------------------------------------------------------------------------------
    *
    * Experimental attempts at PCR error rate modeling
    *
    -------------------------------------------------------------------------------- */

    private byte[] pcrIndelErrorModelCache;
    private PairHMMInputScoreImputator inputScoreImputator;

    private void initializePCRErrorModel() {

        inputScoreImputator = dragstrParams == null
                ? StandardPairHMMInputScoreImputator.newInstance(constantGCP)
                : DragstrPairHMMInputScoreImputator.of(dragstrParams) ;

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
