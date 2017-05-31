package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn;
import org.broadinstitute.hellbender.tools.exome.pulldown.BalancedHeterozygousPileupPriorModel;
import org.broadinstitute.hellbender.tools.exome.pulldown.BayesianHetPulldownCalculator;
import org.broadinstitute.hellbender.tools.exome.pulldown.HeterogeneousHeterozygousPileupPriorModel;
import org.broadinstitute.hellbender.tools.exome.pulldown.Pulldown;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;

/**
 * Calls heterozygous SNPs, counts ref and alt alleles and gives genotype likelihoods.
 *
 * <p>The tool has three modes. The mode is automatically selected based on the provided arguments as follows. </p>
 *
 * <dl>
 *      <dt>MATCHED_NORMAL_TUMOR</dt><dd>Mode when given both matched normal and tumor BAM files.
 *      Estimates the likelihood of heterozygous genotypes on the normal reads using a balanced prior
 *      (see {@link BalancedHeterozygousPileupPriorModel}). Counts bases for the same sites for the
 *      tumor reads.</dd>
 *
 *      <dt>TUMOR_ONLY</dt><dd>For instances lacking a matched normal.
 *      Mode when given just the tumor BAM file with the --tumor argument. Estimates the likelihood of heterozygous genotypes
 *      using the heterogeneous prior (see {@link HeterogeneousHeterozygousPileupPriorModel}).
 *      This mode has only been tested for low purity tumors, i.e. with high normal contamination.
 *      For high purity tumors, both sensitivity and precision may decrease.</dd>
 *
 *     <dt>NORMAL_ONLY</dt><dd>Extended use-case for germline heterozygous site coverage.
 *     Mode when given just the normal BAM file with the --normal argument.
 *     Estimates the likelihood of heterozygous genotypes using the balanced prior
 *     (see {@link BalancedHeterozygousPileupPriorModel}).</dd>
 * </dl>
 *
 * <p>
 *     The tool writes a table with the following columns:
 *
 *     <pre>
 *     CONTIG,
 *     POSITION,
 *     REF_COUNT,
 *     ALT_COUNT,
 *     REF_NUCLEOTIDE,
 *     ALT_NUCLEOTIDE,
 *     READ_DEPTH,
 *     HET_LOG_ODDS
 *     READ_DEPTH = total pileup depth at the SNP site
 *     HET_LOG_ODDS = log odds of the SNP site being heterozygous
 *     </pre>
 * </p>
 *
 * <p>
 *     The required --snpIntervals file gives the genomic positions of common population Single Nucleotide Polymorphisms (SNPs)
 *     in Picard-style interval list format. See https://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList for details.
 *     Here is an example file:
 * </p>
 *
 * <pre>
 *     @HD	VN:1.4	SO:unsorted
 *     @SQ	SN:1	LN:16000	M5:8c0c38e352d8f3309eabe4845456f274
 *     @SQ	SN:2	LN:16000	M5:5f8388fe3fb34aa38375ae6cf5e45b89
 *     1	10736	10736	+	normal
 *     1	11522	11522	+	normal
 *     1	12098	12098	+	normal
 *     1	12444	12444	+	normal
 *      1	13059	13059	+	normal
 *      1	14630	14630	+	normal
 *     1	15204	15204	+	normal
 *     2	14689	14689	+	normal
 *      2	14982	14982	+	normal
 *     2	15110	15110	+	normal
 *     2	15629	15629	+	normal
 * </pre>
 *
 * <h3>Examples</h3>
 *
 * <p>For MATCHED_NORMAL_TUMOR mode:</p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar GetBayesianHetCoverage \
 *   --reference ref_fasta.fa \
 *   --snpIntervals common_sites.interval_list \
 *   --tumor tumor.bam \
 *   --normal normal.bam \
 *   --tumorHets tumor.hets.tsv
 * </pre>
 *
 * <p>
 *     Detects het sites in the normal, and gives ref and alt counts for corresponding sites in the tumor sample.
 *     The tool will not include the log odds of heterozygosity for the tumor
 *     since the inference is made from normal sample's sites.
 * </p>
 *
 * <p>For TUMOR_ONLY mode:</p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar GetBayesianHetCoverage \
 *   --reference ref_fasta.fa \
 *   --snpIntervals common_sites.interval_list \
 *   --tumor tumor.bam \
 *   --tumorHets tumor.hets.tsv
 * </pre>
 *
 **
 * @author Mehrtash Babadi &lt;mehrtash@brgoadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calls sites of heterozygous single nucleotide polymorphisms (SNPs) and counts reads with reference " +
                "and alternate bases at the called sites. A Picard-style interval_list of expected or common population " +
                "variant sites is required. In the first mode, given both matched normal and tumor BAMs, the tool infers heterozygous " +
                "sites from the normal and produces counts and genotype likelihoods for the corresponding sites in the " +
                "tumor. The tool will not include the log odds of heterozygosity for the tumor since the inference is " +
                "made from the normal reads. \n" +
                "Two other modes are available. Given only the tumor BAM, the tool estimates the likelihood of " +
                "heterozygous genotypes using the heterogeneous prior. In contrast, given only a normal BAM, the tool " +
                "estimates the likelihood of heterozygous genotypes using the balanced prior.",
        oneLineSummary = "Calls heterozygous SNPs, counts ref and alt alleles and gives genotype likelihoods",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class GetBayesianHetCoverage extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(GetBayesianHetCoverage.class);

    protected static final String READ_DEPTH_THRESHOLD_FULL_NAME = "readDepthThreshold";
    protected static final String READ_DEPTH_THRESHOLD_SHORT_NAME = "readThresh";

    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";
    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";

    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";
    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";

    protected static final String HET_CALLING_STRINGENCY_FULL_NAME = "hetCallingStringency";
    protected static final String HET_CALLING_STRINGENCY_SHORT_NAME = "hetS";

    protected static final String QUADRATURE_ORDER_FULL_NAME = "quadratureOrder";
    protected static final String QUADRATURE_ORDER_SHORT_NAME = "quad";

    protected static final String MINIMUM_ABNORMAL_FRACTION_FULL_NAME = "minimumAbnormalFraction";
    protected static final String MINIMUM_ABNORMAL_FRACTION_SHORT_NAME = "minAF";

    protected static final String MAXIMUM_ABNORMAL_FRACTION_FULL_NAME = "maximumAbnormalFraction";
    protected static final String MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME = "maxAF";

    protected static final String MAXIMUM_COPY_NUMBER_FULL_NAME = "maximumCopyNumber";
    protected static final String MAXIMUM_COPY_NUMBER_SHORT_NAME = "maxCN";

    protected static final String ERROR_PROBABILITY_ADJUSTMENT_FACTOR_FULL_NAME = "errorAdjustmentFactor";
    protected static final String ERROR_PROBABILITY_ADJUSTMENT_FACTOR_SHORT_NAME = "errFact";

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REFERENCE_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    @Argument(
            doc = "Normal BAM file (for normal-only and paired normal-tumor cases).",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME,
            optional = true
    )
    protected File normalBamFile;

    @Argument(
            doc = "Tumor BAM file (for tumor-only and paired normal-tumor cases).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME,
            optional = true
    )
    protected File tumorBamFile;

    @Argument(
            doc = "IntervalList file of common SNPs.",
            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
            doc = "Output file for base counts at called heterozygous SNP sites (normal).",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    protected File normalHetOutputFile;

    @Argument(
            doc = "Output file for base counts at called heterozygous SNP sites (tumor).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "Minimum read depth; SNP sites with lower read depths will be ignored.",
            fullName = READ_DEPTH_THRESHOLD_FULL_NAME,
            shortName = READ_DEPTH_THRESHOLD_SHORT_NAME,
            optional = true
    )
    protected int readDepthThreshold = 15;

    @Argument(
            doc = "Minimum mapping quality; reads with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            fullName  = MINIMUM_MAPPING_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumMappingQuality = 30;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_BASE_QUALITY_SHORT_NAME,
            fullName = MINIMUM_BASE_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumBaseQuality = 20;

    @Argument(
            doc = "Validation stringency for all BAM files read by this program.  Setting stringency to SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true)
    protected ValidationStringency VALIDATION_STRINGENCY = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(
            doc = "Heterozygosous SNP calling stringency.",
            fullName = HET_CALLING_STRINGENCY_FULL_NAME,
            shortName = HET_CALLING_STRINGENCY_SHORT_NAME,
            optional = true
    )
    protected double hetCallingStringency = 30;

    @Argument(
            doc = "Estimated minimum abnormal cell fraction (for building the allele fraction prior; to be used for " +
                    "tumor-only cases)",
            fullName = MINIMUM_ABNORMAL_FRACTION_FULL_NAME,
            shortName = MINIMUM_ABNORMAL_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double minimumAbnormalFraction = 0.5;

    @Argument(
            doc = "Estimated maximum abnormal cell fraction (for building the allele fraction prior; to be used for " +
                    "tumor-only cases)",
            fullName = MAXIMUM_ABNORMAL_FRACTION_FULL_NAME,
            shortName = MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double maximumAbnormalFraction = 0.8;

    @Argument(
            doc = "Estimated maximum copy number for either homologs in the abnormal portion of the sample (for " +
                    "building the allele fraction prior; to be used for tumor-only cases)",
            fullName = MAXIMUM_COPY_NUMBER_FULL_NAME,
            shortName = MAXIMUM_COPY_NUMBER_SHORT_NAME,
            optional = true
    )
    protected int maximumCopyNumber = 2;

    @Argument(
            doc = "Integration quadrature order; a good choice is the typical read depth (to be used for tumor-only cases)",
            fullName = QUADRATURE_ORDER_FULL_NAME,
            shortName = QUADRATURE_ORDER_SHORT_NAME,
            optional = true
    )
    protected int quadratureOrder = 200;

    @Argument(
            doc = "(Experimental) error probability adjustment factor (could be any positive value)",
            fullName = ERROR_PROBABILITY_ADJUSTMENT_FACTOR_FULL_NAME,
            shortName = ERROR_PROBABILITY_ADJUSTMENT_FACTOR_SHORT_NAME,
            optional = true
    )
    protected double errorProbabilityAdjustmentFactor = 1.0;

    /**
     * The normal-only workflow
     */
    private void runNormalOnly() {
        final BayesianHetPulldownCalculator normalHetPulldownCalculator;
        final Pulldown normalHetPulldown;

        normalHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor,
                new BalancedHeterozygousPileupPriorModel());

        logger.info("Calculating the Het pulldown from the normal BAM file using the BALANCED prior...");
        normalHetPulldown = normalHetPulldownCalculator.getHetPulldown(normalBamFile, hetCallingStringency);

        logger.info("Writing Het pulldown from normal reads to " + normalHetOutputFile.toString());
        normalHetPulldown.write(normalHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.FULL);
    }

    /**
     * The tumor-only workflow
     */
    private void runTumorOnly() {
        final BayesianHetPulldownCalculator tumorHetPulldownCalculator;
        final Pulldown tumorHetPulldown;

        tumorHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor,
                new HeterogeneousHeterozygousPileupPriorModel(minimumAbnormalFraction, maximumAbnormalFraction,
                        maximumCopyNumber, quadratureOrder));

        logger.info("Calculating the Het pulldown from the tumor BAM file using the HETEROGENEOUS prior...");
        tumorHetPulldown = tumorHetPulldownCalculator.getHetPulldown(tumorBamFile, hetCallingStringency);

        logger.info("Writing Het pulldown from tumor reads to " + tumorHetOutputFile.toString());
        tumorHetPulldown.write(tumorHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.FULL);
    }

    /**
     * The matched norrmal-tumor workflow
     */
    private void runMatchedNormalTumor() {
        final BayesianHetPulldownCalculator normalHetPulldownCalculator, tumorHetPulldownCalculator;
        final Pulldown normalHetPulldown, tumorHetPulldown;

        normalHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor,
                new BalancedHeterozygousPileupPriorModel());

        logger.info("Calculating the Het pulldown from the normal BAM file using the BALANCED prior...");
        normalHetPulldown = normalHetPulldownCalculator.getHetPulldown(normalBamFile, hetCallingStringency);

        logger.info("Writing Het pulldown from normal reads to " + normalHetOutputFile.toString());
        normalHetPulldown.write(normalHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.FULL);

        tumorHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                normalHetPulldown.getIntervals(), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor,
                new BalancedHeterozygousPileupPriorModel());

        logger.info("Calculating the Het pulldown from the tumor BAM file on Hets detected in the normal BAM file...");
        tumorHetPulldown = tumorHetPulldownCalculator.getTumorHetPulldownFromNormalPulldown(tumorBamFile,
                normalHetPulldown);

        logger.info("Writing Het pulldown from tumor reads to " + tumorHetOutputFile.toString());
        tumorHetPulldown.write(tumorHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.INTERMEDIATE);
    }

    @Override
    protected Object doWork() {
        /* check for illegal arguments */
        if (normalBamFile == null && tumorBamFile == null) {
            throw new CommandLineException.BadArgumentValue(String.format("You must provide at least one normal (argument %s)" +
                            " or tumor (argument %s) alignment file.",
                    ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_LONG_NAME,
                    ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_LONG_NAME));
        }
        if (normalBamFile != null && normalHetOutputFile == null) {
            throw new CommandLineException(String.format("If the normal BAM file is provided, the normal" +
                    " Het pulldown output file (argument %s) must also be provided.",
                    ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME));
        }
        if (tumorBamFile != null && tumorHetOutputFile == null) {
            throw new CommandLineException(String.format("If the tumor BAM file is provided, the tumor" +
                    " Het pulldown output file (argument %s) must also be provided.",
                    ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME));
        }

        /* determine the job type */
        if (normalBamFile != null && tumorBamFile != null) { /* both normal and tumor BAM files are provided */
            logger.info("MATCHED_NORMAL_TUMOR mode selected.");
            runMatchedNormalTumor();
        } else if (normalBamFile != null) { /* only the normal BAM file is provided */
            logger.info("NORMAL_ONLY mode selected.");
            runNormalOnly();
        } else { /* only the tumor BAM file is provided */
            logger.info("TUMOR_ONLY mode selected.");
            runTumorOnly();
        }

        return "SUCCESS";
    }
}
