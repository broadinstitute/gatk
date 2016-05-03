package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.Feature;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadFilterArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.*;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various covariates
 * (such as read group, reported quality score, machine cycle, and nucleotide context).
 *
 * <p>
 * This walker is designed to work as the first pass in a two-pass processing step. It does a by-locus traversal operating
 * only at sites that are not in dbSNP. We assume that all reference mismatches we see are therefore errors and indicative
 * of poor base quality. This walker generates tables based on various user-specified covariates (such as read group,
 * reported quality score, cycle, and context). Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).
 * <p>
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified.
 *
 * <p>
 *
 * <h3>Input</h3>
 * <p>
 * The input read data whose base quality scores need to be assessed.
 * <p>
 * A database of known polymorphic sites to skip over.
 * </p>
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
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T BaseRecalibrator \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -o recal_data.table
 * </pre>
 */

@CommandLineProgramProperties(
        summary = "First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "Generates recalibration table for BQSR",
        programGroup = ReadProgramGroup.class
)
public final class BaseRecalibrator extends ReadWalker {
    protected static final Logger logger = LogManager.getLogger(BaseRecalibrator.class);

    /**
     * All the command line arguments for BQSR and its covariates.
     */
    @ArgumentCollection(doc="all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection recalArgs = new RecalibrationArgumentCollection();

    /**
     * This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference,
     * so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of
     * Feature-containing files (VCF, BCF, BED, etc.) for use as this database. For users wishing to exclude an interval list of known variation simply
     * use -XL my.interval.list to skip over processing those sites. Please note however that the statistics reported by the tool will not accurately
     * reflected those sites skipped by the -XL argument.
     */
    @Argument(fullName = "knownSites", shortName = "knownSites", doc = "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.", optional = false)
    private List<FeatureInput<Feature>> knownSites;

    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.   Use '/dev/stdout' to print to standard out.
     */
    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, doc = "The output recalibration table file to create", optional = false)
    private File recalTableFile = null;

    private BaseRecalibrationEngine recalibrationEngine;

    private ReferenceDataSource referenceDataSource; // datasource for the reference. We're using a different one from the engine itself to avoid messing with its caches.

    /**
     * an object that keeps track of the information necessary for quality score quantization
     */
    private QuantizationInfo quantizationInfo = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    @Override
    public void onTraversalStart() {
        if (recalArgs.FORCE_PLATFORM != null) {
            recalArgs.DEFAULT_PLATFORM = recalArgs.FORCE_PLATFORM;
        }

        Utils.warnOnNonIlluminaReadGroups(getHeaderForReads(), logger);

        recalibrationEngine = new BaseRecalibrationEngine(recalArgs, getHeaderForReads());
        recalibrationEngine.logCovariatesUsed();
        referenceDataSource = ReferenceDataSource.of(referenceArguments.getReferenceFile());
    }

    @Override
    public List<ReadFilterArgumentCollection.CommandLineReadFilter> getDefaultReadFilters() {
        //Note: the order is deliberate - we first check the cheap conditions that do not require decoding the read
        final List<ReadFilterArgumentCollection.CommandLineReadFilter> bqsrFilters = makeBQSRSpecificReadFilters();
        bqsrFilters.add(ReadFilterArgumentCollection.CommandLineReadFilter.WELLFORMED);
        bqsrFilters.add(ReadFilterArgumentCollection.CommandLineReadFilter.ALLOW_ALL);
        return bqsrFilters;
    }

    // used in BaseRecalibratorSpark and BaseRecalibratorSparkSharded as a ReadFilter
    public static ReadFilter getStandardBQSRReadFilter(final SAMFileHeader header ) {
        //Note: the order is deliberate - we first check the cheap conditions that do not require decoding the read
        return makeBQSRSpecificReadFilters()
                .stream()
                .map(f -> f.getFilterInstance(header, null))
                .reduce(ReadFilterLibrary.ALLOW_ALL_READS, (f1, f2) -> f1.and(f2))
                .and(new WellformedReadFilter(header));
    }

    public static List<ReadFilterArgumentCollection.CommandLineReadFilter> makeBQSRSpecificReadFilters() {
        List<ReadFilterArgumentCollection.CommandLineReadFilter> filters = new ArrayList<>(6);
        filters.add(ReadFilterArgumentCollection.CommandLineReadFilter.MAPPING_QUALITY_NOT_ZERO);
        filters.add(ReadFilterArgumentCollection.CommandLineReadFilter.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED);
        filters.add(ReadFilterArgumentCollection.CommandLineReadFilter.PRIMARY_ALIGNMENT);
        filters.add(ReadFilterArgumentCollection.CommandLineReadFilter.NOT_DUPLICATE);
        filters.add(ReadFilterArgumentCollection.CommandLineReadFilter.PASSES_VENDOR_QUALITY_CHECK);
        return filters;
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    @Override
    public void apply( GATKRead read, ReferenceContext ref, FeatureContext featureContext ) {
        recalibrationEngine.processRead(read, referenceDataSource, featureContext.getValues(knownSites));
    }

    @Override
    public Object onTraversalSuccess() {
        recalibrationEngine.finalizeData();

        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();

        logger.info("Writing recalibration report...");
        generateReport();
        logger.info("...done!");

        //logger.info("BaseRecalibrator was able to recalibrate " + result + " reads");
        return recalibrationEngine.getNumReadsProcessed();
    }

    /**
     * go through the quality score table and use the # observations and the empirical quality score
     * to build a quality score histogram for quantization. Then use the QuantizeQual algorithm to
     * generate a quantization map (recalibrated_qual -> quantized_qual)
     */
    private void quantizeQualityScores() {
        quantizationInfo = new QuantizationInfo(recalibrationEngine.getFinalRecalibrationTables(), recalArgs.QUANTIZING_LEVELS);
    }

    private void generateReport() {
        try ( PrintStream recalTableStream = new PrintStream(recalTableFile) ) {
            RecalUtils.outputRecalibrationReport(recalTableStream, recalArgs, quantizationInfo, recalibrationEngine.getFinalRecalibrationTables(), recalibrationEngine.getCovariates());
        }
        catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(recalTableFile, e);
        }
    }
}