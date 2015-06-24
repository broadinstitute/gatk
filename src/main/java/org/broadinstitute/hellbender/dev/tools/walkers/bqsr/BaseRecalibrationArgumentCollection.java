package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.baq.BAQ;

import java.io.Serializable;
import java.lang.annotation.Annotation;

/**
 * All the command line arguments for BQSR and its covariates.
 */
public final class BaseRecalibrationArgumentCollection implements Serializable, ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    public BaseRecalibrationArgumentCollection() {}

    /**
     * (Almost) all the command line arguments for BQSR and its covariates.
     */
    @ArgumentCollection(doc="all the command line arguments for BQSR and its covariates")
    public final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    @ArgumentCollection
    public final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();


    @ArgumentCollection
    public final ReferenceInputArgumentCollection referenceArguments = new OptionalReferenceInputArgumentCollection();


    @Argument(fullName = "bqsrBAQGapOpenPenalty", shortName = "bqsrBAQGOP", doc = "BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets", optional = true)
    public double BAQGOP = BAQ.DEFAULT_GOP;

    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out unmodified in the recalibrated BAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    @Argument(fullName = "preserve_qscores_less_than", shortName = "preserveQ", doc = "Don't recalibrate bases with quality scores less than this threshold (with -BQSR)", optional = true)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;


    // --------------------------------------------------------------------------------------------------------------
    //
    // quality encoding checking arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * This flag tells GATK to use the original base qualities (that were in the data before BQSR/recalibration) which
     * are stored in the OQ tag, if they are present, rather than use the post-recalibration quality scores. If no OQ
     * tag is present for a read, the standard qual score will be used.
     */
    @Argument(fullName = "useOriginalQualities", shortName = "OQ", doc = "Use the base quality scores from the OQ tag", optional = true)
    public Boolean useOriginalBaseQualities = false;

    /**
     * If reads are missing some or all base quality scores, this value will be used for all base quality scores.
     * By default this is set to -1 to disable default base quality assignment.
     */
    //TODO: minValue = 0, maxValue = Byte.MAX_VALUE)
    @Argument(fullName = "defaultBaseQualities", shortName = "DBQ", doc = "Assign a default base quality", optional = true)
    public byte defaultBaseQualities = -1;

}
