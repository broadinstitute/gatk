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
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CovariateKeyCache;
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

    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation
    private final SAMFileHeader header;
    private final List<SAMReadGroupRecord> readGroups;
    
    private final int preserveQLessThan;
    private final double globalQScorePrior;
    private final boolean emitOriginalQuals;

    //These fields are created to avoid redoing these calculations for every read
    private final int totalCovariateCount;
    private final int specialCovariateCount;
    // tsato: don't like this constant --- just use EventType.BASE_SUBSTITUTION.ordinal().
    private static final int BASE_SUBSTITUTION_INDEX = EventType.BASE_SUBSTITUTION.ordinal();

    //Note: varargs allocates a new array every time. We'll pre-allocate one array and reuse it to avoid object allocation for every base of every read.
    private final RecalDatum[] recalDatumForSpecialCovariates;
    private final boolean useOriginalBaseQualities;

    private byte[] staticQuantizedMapping;
    private final CovariateKeyCache keyCache;

    // tsato: new flag, needs to be moved to apply BQSR argument
    private boolean allowUnknownReadGroup = true;
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
        this.globalQScorePrior = args.globalQScorePrior;
        this.emitOriginalQuals = args.emitOriginalQuals;
        this.useOriginalBaseQualities = args.useOriginalBaseQualities;

        // staticQuantizedQuals is entirely separate from the dynamic binning that quantizationLevels, and
        // staticQuantizedQuals does not make use of quantizationInfo
        if(args.staticQuantizationQuals != null && !args.staticQuantizationQuals.isEmpty()) {
            staticQuantizedMapping = constructStaticQuantizedMapping(args.staticQuantizationQuals, args.roundDown);
        }

        this.quantizedQuals = quantizationInfo.getQuantizedQuals();


        totalCovariateCount = covariates.size();
        specialCovariateCount = covariates.numberOfSpecialCovariates();

        //Note: We pre-create the varargs arrays that will be used in the calls. Otherwise we're spending a lot of time allocating those int[] objects
        recalDatumForSpecialCovariates = new RecalDatum[totalCovariateCount - specialCovariateCount]; // tsato:
        keyCache = new CovariateKeyCache();//one cache per transformer

        readGroups = this.header.getReadGroups();
        allowUnknownReadGroup = args.allowReadGroupsNotInRecalTable;
    }

    /**
     * Constructor using a RecalibrationReport
     *
     * @param header header for the reads
     * @param recalInfo the output of BaseRecalibration, containing the recalibration information
     * @param args a set of arguments to control how bqsr is applied
     */
    public BQSRReadTransformer(final SAMFileHeader header, final RecalibrationReport recalInfo, final ApplyBQSRArgumentCollection args) {
        this(header, recalInfo.getRecalibrationTables(), recalInfo.getQuantizationInfo(), recalInfo.getCovariates(), args); // I see, applyBQSR gets covariate instance form report
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
        if (emitOriginalQuals && ! read.hasAttribute(SAMTag.OQ.name())) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities())); // tsato: phredToFastq?
            } catch (final IllegalArgumentException e) {
                throw new MalformedRead(read, "illegal base quality encountered; " + e.getMessage());
            }
        }
        // tsato: this part needs to be updated if covariates are to be dynamic
        final PerReadCovariateMatrix perReadCovariateMatrix = RecalUtils.computeCovariates(read, header, covariates, false, keyCache);

        // clear indel qualities
        read.clearAttribute(ReadUtils.BQSR_BASE_INSERTION_QUALITIES); // tsato: I forget wtf are indel qualities...how are they computed....
        read.clearAttribute(ReadUtils.BQSR_BASE_DELETION_QUALITIES);

        // this array has dimensions { read length } x { num covariates }. The covariates are encoded as integers.
        final int[][] covariatesForRead = perReadCovariateMatrix.getKeySet(EventType.BASE_SUBSTITUTION); // rename to: cycleByCovariates or something

        final int rgKeyAnotherWay = covariates.getReadGroupCovariate().keyFromValue(ReadGroupCovariate.getReadGroupIdentifier(readGroup)); // tsato: this is ok for now and is better than below
        // the rg key is constant over the whole read, the global deltaQ is too
        final int anyOffset = 0;
        final int rgKey = covariatesForRead[anyOffset][StandardCovariateList.READ_GROUP_COVARIATE_DEFAULT_INDEX]; // tsato: just look up the // tsato: the second 0 identifies the read group.
        // tsato: what is this 2 key thing...


        if (rgKey == READ_GROUP_MISSING_IN_RECAL_TABLE_CODE){
            if (allowUnknownReadGroup) {
                // We cannot hope to recalibrate in this case, since the recalibration code below won't work for such a read.
                final byte[] quals = read.getBaseQualities();
                final byte[] newQuantizedQuals = quals.clone();
                for (int i = 0; i < quals.length; i++){
                    newQuantizedQuals[i] = quantizedQuals.get(quals[i]);
                }
                read.setBaseQualities(newQuantizedQuals);
                // tsato: OQ updated above (as requested)
                return read;
            } else {
                throw new GATKException("copy the same error message as before");
            }
        }
        // tsato: I see. this recalibrationTables is the "cumulative" table...how many tables are there? // The counts for this read group bin
        final RecalDatum readGroupDatum = recalibrationTables.getReadGroupTable().get2Keys(rgKey, BASE_SUBSTITUTION_INDEX);
        // tsato: when can it be null?
        if (readGroupDatum == null) {
            return read; // tsato: ?? so it returns the read as is? no OQ?
        }
        final byte[] quals = read.getBaseQualities();

        final int readLength = quals.length; // tsato: epsilon here is a bad name.
        // tsato: what is this?
        final NestedIntegerArray<RecalDatum> qualityScoreTable = recalibrationTables.getQualityScoreTable();

        //Note: this loop is under very heavy use in applyBQSR. Keep it slim.
        for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read

            // only recalibrate usable qualities (default: >= 6) (the original quality will come from the instrument -- reported quality)
            if (quals[offset] < preserveQLessThan) {
                continue;
            }
            Arrays.fill(recalDatumForSpecialCovariates, null);  //clear the array (tsato: document this array better)
            final int[] covariatesAtOffset = covariatesForRead[offset]; // tsato: rename keySet = covariates

            // tsato: index 0 is read group? 1 is ... base quality?

            final RecalDatum qualityScoreDatum = qualityScoreTable.get3Keys(covariatesAtOffset[StandardCovariateList.READ_GROUP_COVARIATE_DEFAULT_INDEX], covariatesAtOffset[StandardCovariateList.BASE_QUALITY_COVARIATE_DEFAULT_INDEX], BASE_SUBSTITUTION_INDEX);
            // tsato: empiricalQualQS? --- These RecalDatums need to be renamed.
            for (int i = StandardCovariateList.NUM_REQUIRED_COVARITES; i < totalCovariateCount; i++) {
                if (covariatesAtOffset[i] >= 0) {
                    recalDatumForSpecialCovariates[i - StandardCovariateList.NUM_REQUIRED_COVARITES] = recalibrationTables.getTable(i).get4Keys(covariatesAtOffset[0], covariatesAtOffset[1], covariatesAtOffset[i], BASE_SUBSTITUTION_INDEX);
                }
            } // tsato: rename epsilon
            final double epsilon = globalQScorePrior > 0.0 ? globalQScorePrior : readGroupDatum.getEstimatedQReported(); // tsato: wtf are these?
            final double recalibratedQualDouble = hierarchicalBayesianQualityEstimate(epsilon, readGroupDatum, qualityScoreDatum, recalDatumForSpecialCovariates);

            final byte recalibratedQualityScore = quantizedQuals.get(getBoundedIntegerQual(recalibratedQualDouble));

            // Bin to static quals (tsato: in warp pipeline, we use 10, 20, 30, 40, etc.
            quals[offset] = staticQuantizedMapping == null ? recalibratedQualityScore : staticQuantizedMapping[recalibratedQualityScore];
        }
        read.setBaseQualities(quals); // tsato: is this needed? by ref?
        return read;
    }

    // recalibrated quality is bound between 1 and MAX_QUAL // tsato: the return type can be integer, boundQual may be changed.
    private byte getBoundedIntegerQual(final double recalibratedQualDouble) {
        return boundQual(fastRound(recalibratedQualDouble), MAX_RECALIBRATED_Q_SCORE);
    }

    /**
     * Quality score recalibration algorithm works as follows:
     * - Start with the (approximate, or "estimated") reported quality score. (Approximation occurs when marginalizing/collapsing
     * over the reported quality for reach read group).
     * - Compute (retrieve?) the empirical quality score for the read group, usgetBoundedIntegerQualing the reported score as delta (this redundant?)
     * -
     *
     *
     *
     * @param reportedQualityForReadGroup the reported quality score (i.e. in log space) for the read group. This value has type
     *                                    double because of the "combine" (or collapse) operation that converts the quality score table to
     *                                    the read group table.
     * @param readGroupDatum the RecalDatum object for a particular read group at hand. May be null.
     * @param qualityScoreDatum the RecalDatum object for a particular (read group, reported quality) tuple at hand. May be null.
     * @param specialCovariates the array of RecalDatum objects for the non-required covariates.
     * @return
     */
    public static double hierarchicalBayesianQualityEstimate(final double reportedQualityForReadGroup, // tsato: epsilon is by default readGroupDatum.getEstimatedQReported
                                                             final RecalDatum readGroupDatum,
                                                             final RecalDatum qualityScoreDatum,
                                                             final RecalDatum... specialCovariates) {
        // tsato: global in this context means, "across reported qualities"
        // tsato: delta refers to (empirical quality) - (reported quality)
        final double empiricalQualityForReadGroup = readGroupDatum == null ? reportedQualityForReadGroup : readGroupDatum.getEmpiricalQuality(reportedQualityForReadGroup);
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
