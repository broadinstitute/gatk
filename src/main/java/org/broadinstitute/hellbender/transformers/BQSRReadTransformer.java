package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException.MalformedRead;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CovariateKeyCache;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.PerReadCovariateMatrix;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ReadGroupCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.utils.MathUtils.fastRound;
import static org.broadinstitute.hellbender.utils.QualityUtils.boundQual;
import static org.broadinstitute.hellbender.utils.recalibration.RecalDatum.MAX_RECALIBRATED_Q_SCORE;

public final class BQSRReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;
    // tsato: remove?
    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation
    private final SAMFileHeader header;

    private final int preserveQLessThan;
    private final double constantQualityScorePrior;
    private final boolean emitOriginalQuals;

    //These fields are created to avoid redoing these calculations for every read
    private final int totalCovariateCount;
    private final int specialCovariateCount;

    private static final int BASE_SUBSTITUTION_INDEX = EventType.BASE_SUBSTITUTION.ordinal();

    //Note: varargs allocates a new array every time. We'll pre-allocate one array and reuse it to avoid object allocation for every base of every read.
    private final RecalDatum[] recalDatumsForSpecialCovariates;
    private final boolean useOriginalBaseQualities;

    // TODO: this should be merged with the quantized quals
    private byte[] staticQuantizedMapping;
    private final CovariateKeyCache keyCache;

    // tsato: new flag, needs to be moved to apply BQSR argument
    private boolean allowMissingReadGroup;
    private static final int READ_GROUP_MISSING_IN_RECAL_TABLE_CODE = -1;

    private List<Byte> quantizedQuals;

    /**
     * Constructor using a GATK Report file
     *
     * @param header header for the reads
     * @param bqsrRecalFile a GATK Report file containing the recalibration information
     * @param args ApplyBQSR args
     */
    public BQSRReadTransformer(final SAMFileHeader header, final File bqsrRecalFile, final ApplyBQSRArgumentCollection args) {
        this(header, new RecalibrationReport(bqsrRecalFile), args);
    }

    /**
     * Constructor using separate RecalibrationTables, QuantizationInfo, and StandardCovariateList
     *
     * @param header header for the reads
     * @param recalibrationTables recalibration tables output from BQSR
     * @param quantizationInfo quantization info
     * @param covariates standard covariate set
     * @param args ApplyBQSR arguments
     */
    private BQSRReadTransformer(final SAMFileHeader header, final RecalibrationTables recalibrationTables, final QuantizationInfo quantizationInfo, final StandardCovariateList covariates, final ApplyBQSRArgumentCollection args) {
        this.header = header;
        this.recalibrationTables = recalibrationTables;
        this.covariates = covariates;
        this.quantizationInfo = quantizationInfo;

        if (args.quantizationLevels == 0) { // quantizationLevels == 0 means no quantization, preserve the quality scores
            quantizationInfo.noQuantization();
        } else if (args.quantizationLevels > 0 && args.quantizationLevels != quantizationInfo.getQuantizationLevels()) { // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
            quantizationInfo.quantizeQualityScores(args.quantizationLevels);
        }

        this.preserveQLessThan = args.PRESERVE_QSCORES_LESS_THAN;
        this.constantQualityScorePrior = args.globalQScorePrior;
        this.emitOriginalQuals = args.emitOriginalQuals;
        this.useOriginalBaseQualities = args.useOriginalBaseQualities;

        // staticQuantizedQuals is entirely separate from the dynamic binning that quantizationLevels, and
        // staticQuantizedQuals does not make use of quantizationInfo
        if(args.staticQuantizationQuals != null && !args.staticQuantizationQuals.isEmpty()) {
            staticQuantizedMapping = constructStaticQuantizedMapping(args.staticQuantizationQuals, args.roundDown);
        }

        this.quantizedQuals = quantizationInfo.getQuantizedQuals();
        this.totalCovariateCount = covariates.size();
        this.specialCovariateCount = covariates.numberOfSpecialCovariates();

        //Note: We pre-create the varargs arrays that will be used in the calls. Otherwise we're spending a lot of time allocating those int[] objects
        this.recalDatumsForSpecialCovariates = new RecalDatum[specialCovariateCount];
        this.keyCache = new CovariateKeyCache();//one cache per transformer
        this.allowMissingReadGroup = args.allowMissingReadGroups;
    }

    /**
     * Constructor using a RecalibrationReport
     *
     * @param header header for the reads
     * @param recalInfo the output of BaseRecalibration, containing the recalibration information
     * @param args a set of arguments to control how bqsr is applied
     */
    public BQSRReadTransformer(final SAMFileHeader header, final RecalibrationReport recalInfo, final ApplyBQSRArgumentCollection args) {
        this(header, recalInfo.getRecalibrationTables(), recalInfo.getQuantizationInfo(), recalInfo.getCovariates(), args);
    }

    /**
     * Recalibrates the base qualities of a read
     * <p>
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     * <p>
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     * <p>
     * Given the full recalibration table, we perform the following preprocessing steps:
     * <p>
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     * <p>
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param originalRead the read to recalibrate
     */
    @Override
    public GATKRead apply(final GATKRead originalRead) {
        final GATKRead read = useOriginalBaseQualities ? ReadUtils.resetOriginalBaseQualities(originalRead) : originalRead;
        final SAMReadGroupRecord readGroup = header.getReadGroup(read.getReadGroup());

        // If requested, emit the base qualities of the input bam as the OQ (original qualities) of the output bam.
        // If OQ is already set in the input bam, we do not write over it.
        if (emitOriginalQuals && ! read.hasAttribute(SAMTag.OQ.name())) {
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
            } catch (final IllegalArgumentException e) {
                throw new MalformedRead(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        final PerReadCovariateMatrix perReadCovariateMatrix = RecalUtils.computeCovariates(read, header, covariates, false, keyCache);

        // clear indel qualities // tsato: why? Is this ok? Do we update them?
        read.clearAttribute(ReadUtils.BQSR_BASE_INSERTION_QUALITIES);
        read.clearAttribute(ReadUtils.BQSR_BASE_DELETION_QUALITIES);

        // this array has shape ( read length ) x ( num covariates ). The covariates are encoded as integers.
        final int[][] covariatesForRead = perReadCovariateMatrix.getKeySet(EventType.BASE_SUBSTITUTION);

        final int rgKey = covariates.getReadGroupCovariate().keyFromValue(ReadGroupCovariate.getReadGroupIdentifier(readGroup));

        // the rg key is constant over the whole read, the global deltaQ is too
        final int anyOffset = 0;
        final int rgKey2 = covariatesForRead[anyOffset][StandardCovariateList.READ_GROUP_COVARIATE_DEFAULT_INDEX];
        assert rgKey == rgKey2; // tsato: remove this before merging

        final byte[] recalibratedQuals = read.getBaseQualities(); // recall this returns a new copy of the array
        final byte[] preUpdateQuals = read.getBaseQualities();

        // tsato: if the read group encountered is not in the recalibration table...
        if (rgKey == READ_GROUP_MISSING_IN_RECAL_TABLE_CODE){
            if (allowMissingReadGroup) {
                // Given the way the recalibration code is implemented below, we cannot recalibrate a read with a
                // read group that's not in the recal table. TODO: change the implementation below so we can collapse the read groups i.e. marginalize over it
                for (int i = 0; i < recalibratedQuals.length; i++){
                    recalibratedQuals[i] = staticQuantizedMapping != null ? staticQuantizedMapping[preUpdateQuals[i]] : quantizedQuals.get(preUpdateQuals[i]);
                }
                read.setBaseQualities(recalibratedQuals);

                return read;
            } else {
                throw new GATKException("Read group " + read.getReadGroup() + " not found in the recalibration table." +
                        "Set the " + ApplyBQSRArgumentCollection.ALLOW_MISSING_READ_GROUPS_LONG_NAME + " command line argument to avoid this error and" +
                        "simply copy the original qualities as recalibrated qualities.");
            }
        }

        // Datum for the single covariate, read group, with other covariates marginalized.
        final RecalDatum readGroupDatum = recalibrationTables.getReadGroupTable().get2Keys(rgKey, BASE_SUBSTITUTION_INDEX);

        if (readGroupDatum == null) {
            throw new GATKException("readGroupDatum for " + ReadGroupCovariate.getReadGroupIdentifier(readGroup) + " is null");
        }

        final int readLength = recalibratedQuals.length;

        //Note: this loop is under very heavy use in applyBQSR. Keep it slim.
        for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read

            // only recalibrate usable qualities (default: >= 6) (the original quality will come from the instrument -- reported quality)
            if (recalibratedQuals[offset] < preserveQLessThan) {
                continue;
            }

            // clear and reuse this array to save space
            Arrays.fill(recalDatumsForSpecialCovariates, null);
            final int[] covariatesAtOffset = covariatesForRead[offset];
            final int reportedBaseQualityAtOffset = covariatesAtOffset[StandardCovariateList.BASE_QUALITY_COVARIATE_DEFAULT_INDEX];

            // Datum for the tuple (read group, reported quality score).
            final RecalDatum qualityScoreDatum = recalibrationTables.getQualityScoreTable()
                    .get3Keys(rgKey, reportedBaseQualityAtOffset, BASE_SUBSTITUTION_INDEX);

            for (int j = StandardCovariateList.NUM_REQUIRED_COVARITES; j < totalCovariateCount; j++) {
                // If the covariate is -1 (e.g. the first base in each read should have -1 for the context covariate),
                // we simply leave the corresponding Datum to be null, which will subsequently be ignored when it comes time to recalibrate.
                if (covariatesAtOffset[j] >= 0) {
                    recalDatumsForSpecialCovariates[j - StandardCovariateList.NUM_REQUIRED_COVARITES] = recalibrationTables.getTable(j)
                            .get4Keys(rgKey, reportedBaseQualityAtOffset, covariatesAtOffset[j], BASE_SUBSTITUTION_INDEX);
                }
            }

            // Use the reported quality score of the read group as the prior, which can be non-integer because of collapsing.
            // TODO: Avoid using this "estimated" reported quality, and remove it.
            final double priorQualityScore = constantQualityScorePrior > 0.0 ? constantQualityScorePrior : readGroupDatum.getReportedQuality();
            final double rawRecalibratedQualityScore = hierarchicalBayesianQualityEstimate(priorQualityScore, readGroupDatum, qualityScoreDatum, recalDatumsForSpecialCovariates);
            final byte recalibratedQualityScore = quantizedQuals.get(getBoundedIntegerQual(rawRecalibratedQualityScore));

            // Bin to static quals if requested.
            // TODO: as written the code quantizes *twice* if the static binning is enabled (first time to the dynamic bin), which should be fixed
            // recalibratedQuals[offset] = staticQuantizedMapping == null ? recalibratedQualityScore : staticQuantizedMapping[getBoundedIntegerQual(rawRecalibratedQualityScore)];
            recalibratedQuals[offset] = staticQuantizedMapping == null ? recalibratedQualityScore : staticQuantizedMapping[recalibratedQualityScore];
        }

        read.setBaseQualities(recalibratedQuals);
        return read;
    }

    // recalibrated quality is bound between 1 and MAX_QUAL
    private byte getBoundedIntegerQual(final double recalibratedQualDouble) {
        return boundQual(fastRound(recalibratedQualDouble), MAX_RECALIBRATED_Q_SCORE);
    }


    // PROOF READING START HERE (10/13/24)

    /**
     * Quality score recalibration algorithm works as follows:
     * - Start with the (approximate, or "estimated") reported quality score. (Approximation occurs when marginalizing/collapsing
     * over the reported quality for reach read group).
     * - Compute (retrieve?) the empirical quality score for the read group, usgetBoundedIntegerQualing the reported score as delta (this redundant?)
     * -
     *
     *
     *
     * @param priorReportedQuality the reported quality score (i.e. in log space) for the read group. This value has type
     *                                    double because of the "combine" (or collapse) operation that converts the quality score table to
     *                                    the read group table.
     * @param readGroupDatum the RecalDatum object for a particular read group at hand. May be null. (tsato: can it be null?)
     * @param qualityScoreDatum the RecalDatum object for a particular (read group, reported quality) tuple at hand. May be null.
     * @param specialCovariates the array of RecalDatum objects for the non-required covariates (cycle and context covariates by default).
     * @return
     */
    public static double hierarchicalBayesianQualityEstimate(final double priorReportedQuality,
                                                             final RecalDatum readGroupDatum,
                                                             final RecalDatum qualityScoreDatum,
                                                             final RecalDatum... specialCovariates) {
        // tsato: global in this context means, "across reported qualities"
        // tsato: delta refers to (empirical quality) - (reported quality)
        final double empiricalQualityForReadGroup = readGroupDatum == null ? priorReportedQuality : readGroupDatum.getEmpiricalQuality(priorReportedQuality);
        // tsato: feeding an argument to "getEmpiricalQuality()" here is misleading because the argument will be ignored anyway
        final double posteriorEmpiricalQualityForReportedQuality = qualityScoreDatum == null ? empiricalQualityForReadGroup
                : qualityScoreDatum.getEmpiricalQuality(empiricalQualityForReadGroup); // tsato: this is problematic if the prior is ignored, which I believe it is now.

        double deltaQCovariates = 0.0;
        // tsato: at this point we stop being 'iterative' i.e. the special covariates are treated differently.
        for( final RecalDatum empiricalQualCov : specialCovariates ) {
            if (empiricalQualCov != null) {
                // tsato: again, if we ignore the prior...
                deltaQCovariates += empiricalQualCov.getEmpiricalQuality(posteriorEmpiricalQualityForReportedQuality) - posteriorEmpiricalQualityForReportedQuality;
            }
        }

        return posteriorEmpiricalQualityForReportedQuality + deltaQCovariates;
    }

    /**
     * Constructs an array that maps particular quantized values to a rounded value in staticQuantizedQuals
     *
     * Rounding is done in probability space.
     * For instance, Q34 (error probability 10^(-3.4) = 0.00039) is rounded down to Q40 (error probability 10^-4 = 0.0001).
     *
     * When roundDown is true, we simply round down to the nearest
     * available qual in staticQuantizedQuals
     *
     * @param staticQuantizedQuals the list of qualities to round to
     * @param roundDown round down if true, round to nearest (in probability space) otherwise
     * @return  Array where index representing the quality score to be mapped and the value is the rounded quality score
     */
    public static byte[] constructStaticQuantizedMapping(final List<Integer> staticQuantizedQuals, final boolean roundDown) {
        if (staticQuantizedQuals == null || staticQuantizedQuals.isEmpty()){
            return createIdentityMatrix(QualityUtils.MAX_QUAL);
        }
        Utils.nonNull(staticQuantizedQuals);
        // Create array mapping that maps quals to their rounded value.
        final byte[] mapping = new byte[QualityUtils.MAX_QUAL];

        Collections.sort(staticQuantizedQuals);

        // Fill mapping with one-to-one mappings for values between 0 and MIN_USABLE_Q_SCORE
        // This ensures that quals used as special codes will be preserved
        for(int i = 0 ; i < QualityUtils.MIN_USABLE_Q_SCORE ; i++) {
            mapping[i] = (byte) i;
        }

        // If only one staticQuantizedQual is given, fill mappings larger than QualityUtils.MAX_QUAL with that value
        if(staticQuantizedQuals.size() == 1) {
            final int onlyQual = staticQuantizedQuals.get(0);
            for(int i = QualityUtils.MIN_USABLE_Q_SCORE ; i < QualityUtils.MAX_QUAL ; i++) {
                mapping[i] = (byte) onlyQual;
            }
            return mapping;
        }

        // tsato: this implementation appears to be correct, but highly unreadable and prone to error.
        // better to compare to every quantized score, and take the minimum distance....isn't it.
        final int firstQual = QualityUtils.MIN_USABLE_Q_SCORE;
        int previousQual = firstQual;
        double previousProb = QualityUtils.qualToProb(previousQual); // tsato: note this is 1 - (error probability)
        for (final int nextQual : staticQuantizedQuals) {
            final double nextProb = QualityUtils.qualToProb(nextQual);

            for (int i = previousQual; i < nextQual; i++) {
                if (roundDown) {
                    mapping[i] = (byte) previousQual;
                } else {
                    final double iProb = QualityUtils.qualToProb(i);
                    double diff1 = iProb - previousProb;
                    double diff2 = nextProb - iProb;

                    if (iProb - previousProb > nextProb - iProb) { // tsato: don't we need abs...
                        mapping[i] = (byte) nextQual;
                    } else {
                        mapping[i] = (byte) previousQual;
                    }
                }
            }
            previousQual = nextQual;
            previousProb = nextProb;
        }
        // Round all quals larger than the largest static qual down to the largest static qual
        for(int i = previousQual; i < QualityUtils.MAX_QUAL ; i++) {
            mapping[i] = (byte) previousQual;
        }
        return mapping;
    }

    private static byte[] createIdentityMatrix(int maxQual) {
        final byte[] bytes = new byte[maxQual];
        for (int i = 0; i < bytes.length; i++) {
            bytes[i] = (byte) i;
        }
        return bytes;
    }
}
