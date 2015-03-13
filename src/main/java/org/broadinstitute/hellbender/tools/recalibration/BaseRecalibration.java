package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Utility methods to facilitate on-the-fly base quality score recalibration.
 */

public final class BaseRecalibration {
    private static Logger logger = LogManager.getLogger(BaseRecalibration.class);

    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation

    private final boolean disableIndelQuals;
    private final int preserveQLessThan;
    private final double globalQScorePrior;
    private final boolean emitOriginalQuals;

    /**
     * Constructor using a GATK Report file
     *
     * @param RECAL_FILE         a GATK Report file containing the recalibration information
     * @param quantizationLevels number of bins to quantize the quality scores
     * @param disableIndelQuals  if true, do not emit base indel qualities
     * @param preserveQLessThan  preserve quality scores less than this value
     */
    public BaseRecalibration(final File RECAL_FILE, final int quantizationLevels, final boolean disableIndelQuals, final int preserveQLessThan, final boolean emitOriginalQuals, final double globalQScorePrior) {
        RecalibrationReport recalibrationReport = new RecalibrationReport(RECAL_FILE);

        recalibrationTables = recalibrationReport.getRecalibrationTables();
        covariates = recalibrationReport.getCovariates();
        quantizationInfo = recalibrationReport.getQuantizationInfo();
        if (quantizationLevels == 0) // quantizationLevels == 0 means no quantization, preserve the quality scores
            quantizationInfo.noQuantization();
        else if (quantizationLevels > 0 && quantizationLevels != quantizationInfo.getQuantizationLevels()) // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
            quantizationInfo.quantizeQualityScores(quantizationLevels);

        this.disableIndelQuals = disableIndelQuals;
        this.preserveQLessThan = preserveQLessThan;
        this.globalQScorePrior = globalQScorePrior;
        this.emitOriginalQuals = emitOriginalQuals;
    }

    /**
     * Recalibrates the base qualities of a read
     *
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     *
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     *
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param read the read to recalibrate
     */
    public void recalibrateRead(final SAMRecord read) {
        if (emitOriginalQuals && read.getAttribute(SAMTag.OQ.name()) == null) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
            } catch (IllegalArgumentException e) {
                throw new UserException.MalformedBAM(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        final ReadCovariates readCovariates = RecalUtils.computeCovariates(read, covariates);
        final int readLength = read.getReadLength();

        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
            if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
                ReadUtils.setBaseQualities(read, null, errorModel);
                continue;
            }

            final byte[] quals = ReadUtils.getBaseQualities(read, errorModel);

            // get the keyset for this base using the error model
            final int[][] fullReadKeySet = readCovariates.getKeySet(errorModel);

            // the rg key is constant over the whole read, the global deltaQ is too
            final int rgKey = fullReadKeySet[0][0];
            final RecalDatum empiricalQualRG = recalibrationTables.getReadGroupTable().get(rgKey, errorModel.ordinal());

            if( empiricalQualRG != null ) {
                final double epsilon = ( globalQScorePrior > 0.0 && errorModel.equals(EventType.BASE_SUBSTITUTION) ? globalQScorePrior : empiricalQualRG.getEstimatedQReported() );

                for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read
                    final byte origQual = quals[offset];

                    // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                    if ( origQual >= preserveQLessThan ) {
                        // get the keyset for this base using the error model
                        final int[] keySet = fullReadKeySet[offset];
                        final RecalDatum empiricalQualQS = recalibrationTables.getQualityScoreTable().get(keySet[0], keySet[1], errorModel.ordinal());
                        final List<RecalDatum> empiricalQualCovs = new ArrayList<>();
                        for (int i = 2; i < covariates.size(); i++) {  //XXX the 2 is hard-wired here as the number of special covariates
                            if (keySet[i] < 0) {
                                continue;
                            }
                            empiricalQualCovs.add(recalibrationTables.getTable(i).get(keySet[0], keySet[1], keySet[i], errorModel.ordinal()));
                        }

                        double recalibratedQualDouble = hierarchicalBayesianQualityEstimate( epsilon, empiricalQualRG, empiricalQualQS, empiricalQualCovs );

                        // recalibrated quality is bound between 1 and MAX_QUAL
                        final byte recalibratedQual = QualityUtils.boundQual(MathUtils.fastRound(recalibratedQualDouble), RecalDatum.MAX_RECALIBRATED_Q_SCORE);

                        // return the quantized version of the recalibrated quality
                        final byte recalibratedQualityScore = quantizationInfo.getQuantizedQuals().get(recalibratedQual);

                        quals[offset] = recalibratedQualityScore;
                    }
                }
            }

            // finally update the base qualities in the read
            ReadUtils.setBaseQualities(read, quals, errorModel);
        }
    }

    protected static double hierarchicalBayesianQualityEstimate( final double epsilon, final RecalDatum empiricalQualRG, final RecalDatum empiricalQualQS, final List<RecalDatum> empiricalQualCovs ) {
        final double globalDeltaQ = ( empiricalQualRG == null ? 0.0 : empiricalQualRG.getEmpiricalQuality(epsilon) - epsilon );
        final double deltaQReported = ( empiricalQualQS == null ? 0.0 : empiricalQualQS.getEmpiricalQuality(globalDeltaQ + epsilon) - (globalDeltaQ + epsilon) );
        double deltaQCovariates = 0.0;
        for( final RecalDatum empiricalQualCov : empiricalQualCovs ) {
            deltaQCovariates += ( empiricalQualCov == null ? 0.0 : empiricalQualCov.getEmpiricalQuality(deltaQReported + globalDeltaQ + epsilon) - (deltaQReported + globalDeltaQ + epsilon) );
        }

        return epsilon + globalDeltaQ + deltaQReported + deltaQCovariates;
    }
}
