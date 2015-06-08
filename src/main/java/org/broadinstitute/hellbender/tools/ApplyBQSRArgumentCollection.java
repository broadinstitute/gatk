package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.Serializable;
import java.lang.annotation.Annotation;

public class ApplyBQSRArgumentCollection implements ArgumentCollection, Serializable {
    private static final long serialVersionUID = 1L;

    @Override
    public String doc() {
        return "Options related to applying quality score recalibration";
    }

    /**
     * Turns on the base quantization module. It requires a recalibration report.
     *
     * A value of 0 here means "do not quantize".
     * Any value greater than zero will be used to recalculate the quantization using that many levels.
     * Negative values mean that we should quantize using the recalibration report's quantization level.
     */
    @Argument(fullName="quantize_quals", shortName = "qq", doc = "Quantize quality scores to a given number of levels", optional=true)
    public int quantizationLevels = 0;

    /**
     * Turns off printing of the base insertion and base deletion tags.
     * Only the base substitution qualities will be produced.
     */
    @Argument(fullName="disable_indel_quals", shortName = "DIQ", doc = "Disable printing of base insertion and deletion tags", optional=true)
    public boolean disableIndelQuals = false;

    /**
     * By default, the OQ tag in not emitted. Use this flag to include OQ tags in the output BAM file.
     * Note that this may results in significant file size increase.
     */
    @Argument(fullName="emit_original_quals", shortName = "EOQ", doc = "Emit the OQ tag with the original base qualities", optional=true)
    public boolean emitOriginalQuals = false;

    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out unmodified in the recalibrated BAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    //TODO: add those when https://github.com/broadinstitute/hellbender/issues/143 is fixed
    //TODO: minValue = 0, minRecommendedValue = QualityUtils.MIN_USABLE_Q_SCORE
    @Argument(fullName = "preserve_qscores_less_than", shortName = "preserveQ", doc = "Don't recalibrate bases with quality scores less than this threshold", optional = true)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;

    /**
     * If specified, this value will be used as the prior for all mismatch quality scores instead of the actual reported quality score.
     */
    @Argument(fullName = "globalQScorePrior", shortName = "globalQScorePrior", doc = "Global Qscore Bayesian prior to use for BQSR", optional = true)
    public double globalQScorePrior = -1.0;

    @Override
    public Class<? extends Annotation> annotationType() {
        return ApplyBQSRArgumentCollection.class;
    }
}
