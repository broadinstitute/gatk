package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.exceptions.UserException.MalformedRead;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CovariateKeyCache;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ReadCovariates;
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

    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation
    private final SAMFileHeader header;
    
    private final int preserveQLessThan;
    private final double globalQScorePrior;
    private final boolean emitOriginalQuals;

    //These fields are created to avoid redoing these calculations for every read
    private final int totalCovariateCount;
    private final int specialCovariateCount;

    private static final int BASE_SUBSTITUTION_INDEX = EventType.BASE_SUBSTITUTION.ordinal();

    //Note: varargs allocates a new array every time. We'll pre-allocate one array and reuse it to avoid object allocation for every base of every read.
    private final RecalDatum[] empiricalQualCovsArgs;
    private final boolean useOriginalBaseQualities;

    private byte[] staticQuantizedMapping;
    private final CovariateKeyCache keyCache;

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
        this.globalQScorePrior = args.globalQScorePrior;
        this.emitOriginalQuals = args.emitOriginalQuals;
        this.useOriginalBaseQualities = args.useOriginalBaseQualities;

        // staticQuantizedQuals is entirely separate from the dynamic binning that quantizationLevels, and
        // staticQuantizedQuals does not make use of quantizationInfo
        if(args.staticQuantizationQuals != null && !args.staticQuantizationQuals.isEmpty()) {
            staticQuantizedMapping = constructStaticQuantizedMapping(args.staticQuantizationQuals, args.roundDown);
        }

        totalCovariateCount = covariates.size();
        specialCovariateCount = covariates.numberOfSpecialCovariates();

        //Note: We pre-create the varargs arrays that will be used in the calls. Otherwise we're spending a lot of time allocating those int[] objects
        empiricalQualCovsArgs = new RecalDatum[totalCovariateCount - specialCovariateCount];
        keyCache = new CovariateKeyCache();//one cache per transformer
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

        if (emitOriginalQuals && ! read.hasAttribute(SAMTag.OQ.name())) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
            } catch (final IllegalArgumentException e) {
                throw new MalformedRead(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        final ReadCovariates readCovariates = RecalUtils.computeCovariates(read, header, covariates, false, keyCache);

        //clear indel qualities
        read.clearAttribute(ReadUtils.BQSR_BASE_INSERTION_QUALITIES);
        read.clearAttribute(ReadUtils.BQSR_BASE_DELETION_QUALITIES);

        // get the keyset for this base using the error model
        final int[][] fullReadKeySet = readCovariates.getKeySet(EventType.BASE_SUBSTITUTION);

        // the rg key is constant over the whole read, the global deltaQ is too
        final int rgKey = fullReadKeySet[0][0];

        final RecalDatum empiricalQualRG = recalibrationTables.getReadGroupTable().get2Keys(rgKey, BASE_SUBSTITUTION_INDEX);

        if (empiricalQualRG == null) {
            return read;
        }
        final byte[] quals = read.getBaseQualities();

        final int readLength = quals.length;
        final double epsilon = globalQScorePrior > 0.0 ? globalQScorePrior : empiricalQualRG.getEstimatedQReported();

        final NestedIntegerArray<RecalDatum> qualityScoreTable = recalibrationTables.getQualityScoreTable();
        final List<Byte> quantizedQuals = quantizationInfo.getQuantizedQuals();

        //Note: this loop is under very heavy use in applyBQSR. Keep it slim.
        for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read

            // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
            if (quals[offset] < preserveQLessThan) {
                continue;
            }
            Arrays.fill(empiricalQualCovsArgs, null);  //clear the array
            final int[] keySet = fullReadKeySet[offset];


            final RecalDatum empiricalQualQS = qualityScoreTable.get3Keys(keySet[0], keySet[1], BASE_SUBSTITUTION_INDEX);

            for (int i = specialCovariateCount; i < totalCovariateCount; i++) {
                if (keySet[i] >= 0) {
                    empiricalQualCovsArgs[i - specialCovariateCount] = recalibrationTables.getTable(i).get4Keys(keySet[0], keySet[1], keySet[i], BASE_SUBSTITUTION_INDEX);
                }
            }
            final double recalibratedQualDouble = hierarchicalBayesianQualityEstimate(epsilon, empiricalQualRG, empiricalQualQS, empiricalQualCovsArgs);

            final byte recalibratedQualityScore = quantizedQuals.get(getRecalibratedQual(recalibratedQualDouble));

            // Bin to static quals
            quals[offset] = staticQuantizedMapping == null ? recalibratedQualityScore : staticQuantizedMapping[recalibratedQualityScore];
        }
        read.setBaseQualities(quals);
        return read;
    }

    // recalibrated quality is bound between 1 and MAX_QUAL
    private byte getRecalibratedQual(final double recalibratedQualDouble) {
        return boundQual(fastRound(recalibratedQualDouble), MAX_RECALIBRATED_Q_SCORE);
    }

    public static double hierarchicalBayesianQualityEstimate( final double epsilon,
                                                              final RecalDatum empiricalQualRG,
                                                              final RecalDatum empiricalQualQS,
                                                              final RecalDatum... empiricalQualCovs ) {
        final double globalDeltaQ = empiricalQualRG == null ? 0.0 : empiricalQualRG.getEmpiricalQuality(epsilon) - epsilon;
        final double deltaQReported = empiricalQualQS == null ? 0.0 : empiricalQualQS.getEmpiricalQuality(globalDeltaQ + epsilon) - (globalDeltaQ + epsilon);

        double deltaQCovariates = 0.0;
        final double conditionalPrior2 = deltaQReported + globalDeltaQ + epsilon;
        for( final RecalDatum empiricalQualCov : empiricalQualCovs ) {
            if (empiricalQualCov != null) {
                deltaQCovariates += empiricalQualCov.getEmpiricalQuality(conditionalPrior2) - conditionalPrior2;
            }
        }

        return conditionalPrior2 + deltaQCovariates;
    }

    /**
     * Constructs an array that maps particular quantized values to a rounded value in staticQuantizedQuals
     *
     * Rounding is done in probability space.  When roundDown is true, we simply round down to the nearest
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

        final int firstQual = QualityUtils.MIN_USABLE_Q_SCORE;
        int previousQual = firstQual;
        double previousProb = QualityUtils.qualToProb(previousQual);
        for (final int nextQual : staticQuantizedQuals) {
            final double nextProb = QualityUtils.qualToProb(nextQual);

            for (int i = previousQual; i < nextQual; i++) {
                if (roundDown) {
                    mapping[i] = (byte) previousQual;
                } else {
                    final double iProb = QualityUtils.qualToProb(i);
                    if (iProb - previousProb > nextProb - iProb) {
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
