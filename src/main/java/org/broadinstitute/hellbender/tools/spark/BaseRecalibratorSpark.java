package org.broadinstitute.hellbender.tools.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.spark.JoinReadsWithVariants;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.PrintStream;
import java.util.List;

/**
 * Spark version of the first pass of the base quality score recalibration.
 * Generates a recalibration table based on various covariates.
 * The default covariates are read group, reported quality score, machine cycle, and nucleotide context.
 *
 * <p>
 * This walker generates tables based on specified covariates.
 * It does a by-locus traversal operating only at sites that are not in the known-sites resource.
 * ExAc, gnomAD, or dbSNP resources can be used as known sites of variation.
 * We assume that all reference mismatches we see are therefore errors and indicative of poor base quality.
 * Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).
 * </p>
 *

 *
 * <h3>Input</h3>
 * <ol>
 *  <li>The input read data whose base quality scores need to be assessed.</li>
 *  <li>A database of known polymorphic sites to skip over.</li>
 * </ol>
 *
 * <h3>Output</h3>
 * <p>
 * A GATK Report file with many tables:
 * <ol>
 *     <li>The list of arguments</li>
 *     <li>The quantized qualities table</li>
 *     <li>The recalibration table by read group</li>
 *     <li>The recalibration table by quality score</li>
 *     <li>The recalibration table for all the optional covariates</li>
 * </ol>
 *
 * The GATK Report is intended to be easy to read by humans or computers. Check out the documentation of the GATKReport to learn how to manipulate this table.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * gatk BaseRecalibratorSpark \
 *   -I gs://my-gcs-bucket/my_reads.bam \
 *   -R gs://my-gcs-bucket/reference.fasta \
 *   --known-sites gs://my-gcs-bucket/sites_of_variation.vcf \
 *   --known-sites gs://my-gcs-bucket/another/optional/setOfSitesToMask.vcf \
 *   -O gs://my-gcs-bucket/recal_data.table \
 *   -- \
 *   --sparkRunner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 */

@CommandLineProgramProperties(
        summary = BaseRecalibratorSpark.USAGE_SUMMARY,
        oneLineSummary = BaseRecalibratorSpark.USAGE_ONE_LINE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class BaseRecalibratorSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    static final String USAGE_ONE_LINE_SUMMARY = "Generate recalibration table for Base Quality Score Recalibration (BQSR) on Spark";
    static final String USAGE_SUMMARY = "First pass of the Base Quality Score Recalibration (BQSR) on Spark." +
            " Generate a recalibration table based on various user-specified covariates " +
            "(such as read group, reported quality score, machine cycle, and nucleotide context).";

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return BaseRecalibrator.getStandardBQSRReadFilterList();
    }

    @Argument(doc = "the known variants", fullName = BaseRecalibrator.KNOWN_SITES_ARG_FULL_NAME, optional = false)
    private List<String> knownVariants;

    @Argument(doc = "Path to save the final recalibration tables to.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputTablesPath = null;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    @Override
    protected void runTool( JavaSparkContext ctx ) {
        String referenceFileName = addReferenceFilesForSpark(ctx, referenceArguments.getReferenceFileName());
        List<String> localKnownSitesFilePaths = addVCFsForSpark(ctx, knownVariants);

        JavaPairRDD<GATKRead, Iterable<GATKVariant>> readsWithVariants = JoinReadsWithVariants.join(getReads(), localKnownSitesFilePaths);

        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(readsWithVariants, getHeaderForReads(), referenceFileName, bqsrArgs);

        try ( final PrintStream reportStream = new PrintStream(BucketUtils.createFile(outputTablesPath)) ) {
            RecalUtils.outputRecalibrationReport(reportStream, bqsrArgs, bqsrReport.getQuantizationInfo(), bqsrReport.getRecalibrationTables(), bqsrReport.getCovariates());
        }
    }
}
