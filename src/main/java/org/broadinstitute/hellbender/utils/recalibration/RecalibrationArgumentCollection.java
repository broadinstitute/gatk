package org.broadinstitute.hellbender.utils.recalibration;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.File;
import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * A collection of the arguments that are used for BQSR. Used to be common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 */

public final class RecalibrationArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    // We always use the same covariates. The field is retained for compatibility with GATK3 reports.
    public static final boolean DO_NOT_USE_STANDARD_COVARIATES = false;

    //We don't support SOLID. The field is retained for compatibility with GATK3 reports.
    public static final String SOLID_RECAL_MODE = "SET_Q_ZERO";
    public static final String SOLID_NOCALL_STRATEGY = "THROW_EXCEPTION";

    //It makes no sense to run BQSR without sites. so we remove this option.
    public static final boolean RUN_WITHOUT_DBSNP = false;

    /**
     * The context covariate will use a context of this size to calculate its covariate value for base mismatches. Must be between 1 and 13 (inclusive). Note that higher values will increase runtime and required java heap size.
     */
    @Argument(fullName = "mismatches-context-size", shortName = "mcs", doc = "Size of the k-mer context to be used for base mismatches", optional = true)
    public int MISMATCHES_CONTEXT_SIZE = 2;

    /**
     * The context covariate will use a context of this size to calculate its covariate value for base insertions and deletions. Must be between 1 and 13 (inclusive). Note that higher values will increase runtime and required java heap size.
     */
    @Argument(fullName = "indels-context-size", shortName = "ics", doc = "Size of the k-mer context to be used for base insertions and deletions", optional = true)
    public int INDELS_CONTEXT_SIZE = 3;

    /**
     * The cycle covariate will generate an error if it encounters a cycle greater than this value.
     * This argument is ignored if the Cycle covariate is not used.
     */
    @Argument(fullName = "maximum-cycle-value", shortName = "max-cycle", doc = "The maximum cycle value permitted for the Cycle covariate", optional = true)
    public int MAXIMUM_CYCLE_VALUE = 500;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off. [default is off]
     */
    @Argument(fullName = "mismatches-default-quality", doc = "default quality for the base mismatches covariate", optional = true)
    public byte MISMATCHES_DEFAULT_QUALITY = -1;

    /**
     * A default base qualities to use as a prior (reported quality) in the insertion covariate model. This parameter is used for all reads without insertion quality scores for each base. [default is on]
     */
    @Argument(fullName = "insertions-default-quality", doc = "default quality for the base insertions covariate", optional = true)
    public byte INSERTIONS_DEFAULT_QUALITY = 45;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off. [default is on]
     */
    @Argument(fullName = "deletions-default-quality", doc = "default quality for the base deletions covariate", optional = true)
    public byte DELETIONS_DEFAULT_QUALITY = 45;

    /**
     * Reads with low quality bases on either tail (beginning or end) will not be considered in the context. This parameter defines the quality below which (inclusive) a tail is considered low quality
     */
    @Argument(fullName = "low-quality-tail", doc = "minimum quality for the bases in the tail of the reads to be considered", optional = true)
    public byte LOW_QUAL_TAIL = 2;

    /**
     * BQSR generates a quantization table for quick quantization later by subsequent tools. BQSR does not quantize the base qualities, this is done by the engine with the -qq or -bqsr options.
     * This parameter tells BQSR the number of levels of quantization to use to build the quantization table.
     */
    @Argument(fullName = "quantizing-levels", optional = true, doc = "number of distinct quality scores in the quantized output")
    public int QUANTIZING_LEVELS = 16;

    /**
     * The tag name for the binary tag covariate (if using it)
     */
    @Argument(fullName = "binary-tag-name", optional = true, doc = "the binary tag covariate name if using it")
    public String BINARY_TAG_NAME = null;

    @Argument(fullName = "bqsr-baq-gap-open-penalty", doc="BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets", optional = true)
    public double BAQGOP = BAQ.DEFAULT_GOP;

    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out unmodified in the recalibrated BAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    @Argument(fullName = "preserve-qscores-less-than", doc = "Don't recalibrate bases with quality scores less than this threshold (with -" + StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME + ")", optional = true)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;

    @Hidden
    @Argument(fullName = "enable-baq", doc = "do BAQ correction")
    public boolean enableBAQ = false;

    @Hidden
    @Argument(fullName = "compute-indel-bqsr-tables", shortName = "indels", doc = "compute indel BQSR tables")
    public boolean computeIndelBQSRTables = false;


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
    @Argument(fullName="use-original-qualities", shortName = "OQ", doc = "Use the base quality scores from the OQ tag", optional = true)
    public Boolean useOriginalBaseQualities = false;

    /**
     * If reads are missing some or all base quality scores, this value will be used for all base quality scores.
     * By default this is set to -1 to disable default base quality assignment.
     */
    //TODO: minValue = 0, maxValue = Byte.MAX_VALUE)
    @Argument(fullName="default-base-qualities", doc = "Assign a default base quality", optional = true)
    public byte defaultBaseQualities = -1;


    /////////////////////////////
    // Debugging-only Arguments
    /////////////////////////////

    @Hidden
    @Argument(fullName = "default-platform", optional = true, doc = "If a read has no platform then default to the provided String. Valid options are illumina, 454, and solid.")
    public String DEFAULT_PLATFORM = null;

    @Hidden
    @Argument(fullName = "force-platform", optional = true, doc = "If provided, the platform of EVERY read will be forced to be the provided String. Valid options are illumina, 454, and solid.")
    public String FORCE_PLATFORM = null;

    public File existingRecalibrationReport = null;

    public GATKReportTable generateReportTable(final String covariateNames) {
        GATKReportTable argumentsTable;
        argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run", 2, GATKReportTable.Sorting.SORT_BY_COLUMN);
        argumentsTable.addColumn("Argument", "%s");
        argumentsTable.addColumn(RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, "");
        argumentsTable.addRowID("covariate", true);
        argumentsTable.set("covariate", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, covariateNames);
        argumentsTable.addRowID("no_standard_covs", true);
        argumentsTable.set("no_standard_covs", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DO_NOT_USE_STANDARD_COVARIATES);
        argumentsTable.addRowID("run_without_dbsnp", true);
        argumentsTable.set("run_without_dbsnp", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, RUN_WITHOUT_DBSNP);
        argumentsTable.addRowID("solid_recal_mode", true);
        argumentsTable.set("solid_recal_mode", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, SOLID_RECAL_MODE);
        argumentsTable.addRowID("solid_nocall_strategy", true);
        argumentsTable.set("solid_nocall_strategy", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, SOLID_NOCALL_STRATEGY);
        argumentsTable.addRowID("mismatches_context_size", true);
        argumentsTable.set("mismatches_context_size", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MISMATCHES_CONTEXT_SIZE);
        argumentsTable.addRowID("indels_context_size", true);
        argumentsTable.set("indels_context_size", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INDELS_CONTEXT_SIZE);
        argumentsTable.addRowID("mismatches_default_quality", true);
        argumentsTable.set("mismatches_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MISMATCHES_DEFAULT_QUALITY);
        argumentsTable.addRowID("deletions_default_quality", true);
        argumentsTable.set("deletions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DELETIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("insertions_default_quality", true);
        argumentsTable.set("insertions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INSERTIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("maximum_cycle_value", true);
        argumentsTable.set("maximum_cycle_value", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MAXIMUM_CYCLE_VALUE);
        argumentsTable.addRowID("low_quality_tail", true);
        argumentsTable.set("low_quality_tail", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, LOW_QUAL_TAIL);
        argumentsTable.addRowID("default_platform", true);
        argumentsTable.set("default_platform", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DEFAULT_PLATFORM);
        argumentsTable.addRowID("force_platform", true);
        argumentsTable.set("force_platform", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, FORCE_PLATFORM);
        argumentsTable.addRowID("quantizing_levels", true);
        argumentsTable.set("quantizing_levels", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, QUANTIZING_LEVELS);
        argumentsTable.addRowID("recalibration_report", true);
        argumentsTable.set("recalibration_report", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, existingRecalibrationReport == null ? "null" : existingRecalibrationReport.getAbsolutePath());
        argumentsTable.addRowID("binary_tag_name", true);
        argumentsTable.set("binary_tag_name", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, BINARY_TAG_NAME == null ? "null" : BINARY_TAG_NAME);
        return argumentsTable;
    }

    /**
     * Returns a map with the arguments that differ between this an
     * another {@link RecalibrationArgumentCollection} instance.
     * <p/>
     * The key is the name of that argument in the report file. The value is a message
     * that explains the difference to the end user.
     * <p/>
     * Thus, a empty map indicates that there is no differences between both argument collection that
     * is relevant to report comparison.
     * <p/>
     * This method should not throw any exception.
     *
     * @param other the argument-collection to compare against.
     * @param thisRole the name used to refer to this RAC report that makes sense to the end user.
     * @param otherRole the name used to refer to the other RAC report that makes sense to the end user.
     *
     * @return never <code>null</code>, but a zero-size collection if there are no differences.
     */
    public Map<String,? extends CharSequence> compareReportArguments(final RecalibrationArgumentCollection other,final String thisRole, final String otherRole) {
        final Map<String,String> result = new LinkedHashMap<>(15);
        compareSimpleReportArgument(result,"no_standard_covs", DO_NOT_USE_STANDARD_COVARIATES, DO_NOT_USE_STANDARD_COVARIATES, thisRole, otherRole);
        compareSimpleReportArgument(result,"run_without_dbsnp",RUN_WITHOUT_DBSNP, RUN_WITHOUT_DBSNP,thisRole,otherRole);
        compareSimpleReportArgument(result,"solid_recal_mode", SOLID_RECAL_MODE, SOLID_RECAL_MODE,thisRole,otherRole);
        compareSimpleReportArgument(result,"solid_nocall_strategy", SOLID_NOCALL_STRATEGY, SOLID_NOCALL_STRATEGY,thisRole,otherRole);
        compareSimpleReportArgument(result,"mismatches_context_size", MISMATCHES_CONTEXT_SIZE,other.MISMATCHES_CONTEXT_SIZE,thisRole,otherRole);
        compareSimpleReportArgument(result,"mismatches_default_quality", MISMATCHES_DEFAULT_QUALITY, other.MISMATCHES_DEFAULT_QUALITY,thisRole,otherRole);
        compareSimpleReportArgument(result,"deletions_default_quality", DELETIONS_DEFAULT_QUALITY, other.DELETIONS_DEFAULT_QUALITY,thisRole,otherRole);
        compareSimpleReportArgument(result,"insertions_default_quality", INSERTIONS_DEFAULT_QUALITY, other.INSERTIONS_DEFAULT_QUALITY,thisRole,otherRole);
        compareSimpleReportArgument(result,"maximum_cycle_value", MAXIMUM_CYCLE_VALUE, other.MAXIMUM_CYCLE_VALUE,thisRole,otherRole);
        compareSimpleReportArgument(result,"low_quality_tail", LOW_QUAL_TAIL, other.LOW_QUAL_TAIL,thisRole,otherRole);
        compareSimpleReportArgument(result,"default_platform", DEFAULT_PLATFORM, other.DEFAULT_PLATFORM,thisRole,otherRole);
        compareSimpleReportArgument(result,"force_platform", FORCE_PLATFORM, other.FORCE_PLATFORM,thisRole,otherRole);
        compareSimpleReportArgument(result,"quantizing_levels", QUANTIZING_LEVELS, other.QUANTIZING_LEVELS,thisRole,otherRole);
        compareSimpleReportArgument(result,"binary_tag_name", BINARY_TAG_NAME, other.BINARY_TAG_NAME,thisRole,otherRole);
        return result;
    }


    /**
     * Annotates a map with any difference encountered in a simple value report argument that differs between this an
     * another {@link RecalibrationArgumentCollection} instance.
     * <p/>
     * The key of the new entry would be the name of that argument in the report file. The value is a message
     * that explains the difference to the end user.
     * <p/>
     *
     * <p/>
     * This method should not return any exception.
     *
     * @param diffs where to annotate the differences.
     * @param name the name of the report argument to compare.
     * @param thisValue this argument collection value for that argument.
     * @param otherValue the other collection value for that argument.
     * @param thisRole the name used to refer to this RAC report that makes sense to the end user.
     * @param otherRole the name used to refer to the other RAC report that makes sense to the end user.
     *
     * @type T the argument Object value type.
     *
     * @return <code>true</code> if a difference has been spotted, thus <code>diff</code> has been modified.
     */
    private <T> boolean compareSimpleReportArgument(final Map<String,String> diffs,
            final String name, final T thisValue, final T otherValue, final String thisRole, final String otherRole) {
        if (thisValue == null && otherValue == null) {
            return false;
        } else if (thisValue != null && thisValue.equals(otherValue)) {
            return false;
        } else {
            diffs.put(name,
                    String.format("differences between '%s' {%s} and '%s' {%s}.",
                            thisRole, thisValue == null ? "" : thisValue,
                            otherRole, otherValue == null ? "" : otherValue));
            return true;
        }

    }
}
