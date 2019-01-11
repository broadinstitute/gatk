package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.denoising.GCBiasCorrector;
import org.broadinstitute.hellbender.tools.copynumber.denoising.HDF5SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;

/**
 * Creates a panel of normals (PoN) for read-count denoising given the read counts for samples in the panel.
 * The resulting PoN can be used with {@link DenoiseReadCounts} to denoise other samples.
 *
 * <p>
 *     The input read counts are first transformed to log2 fractional coverages and preprocessed
 *     according to specified filtering and imputation parameters.  Singular value decomposition (SVD)
 *     is then performed to find the first {@code number-of-eigensamples} principal components,
 *     which are stored in the PoN.  Some or all of these principal components can then be used for
 *     denoising case samples with {@link DenoiseReadCounts}; it is assumed that the principal components used
 *     represent systematic sequencing biases (rather than statistical noise).  Examining the singular values,
 *     which are also stored in the PoN, may be useful in determining the appropriate number
 *     of principal components to use for denoising.
 * </p>
 *
 * <p>
 *     If annotated intervals are provided, explicit GC-bias correction will be performed by {@link GCBiasCorrector}
 *     before filtering and SVD.  GC-content information for the intervals will be stored in the PoN
 *     and used to perform explicit GC-bias correction identically in {@link DenoiseReadCounts}.
 *     Note that if annotated intervals are not provided, it is still likely that GC-bias correction is
 *     implicitly performed by the SVD denoising process (i.e., some of the principal components arise from GC bias).
 * </p>
 *
 * <p>
 *     Note that such SVD denoising cannot distinguish between variance due to systematic sequencing biases and that
 *     due to true common germline CNVs present in the panel; signal from the latter may thus be inadvertently denoised
 *     away.  Furthermore, variance arising from coverage on the sex chromosomes may also significantly contribute
 *     to the principal components if the panel contains samples of mixed sex.  Therefore, if sex chromosomes
 *     are not excluded from coverage collection, it is strongly recommended that users avoid creating panels of
 *     mixed sex and take care to denoise case samples only with panels containing only individuals of the same sex
 *     as the case samples.  (See {@link GermlineCNVCaller}, which avoids these issues by simultaneously learning
 *     a probabilistic model for systematic bias and calling rare and common germline CNVs for samples in the panel.)
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Counts files (TSV or HDF5 output of {@link CollectReadCounts}).
 *     </li>
 *     <li>
 *         (Optional) GC-content annotated-intervals file from {@link AnnotateIntervals}.
 *         Explicit GC-bias correction will be performed on the panel samples and identically for subsequent case samples.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Panel-of-normals file.
 *         This is an HDF5 file containing the panel data in the paths defined in {@link HDF5SVDReadCountPanelOfNormals}.
 *         HDF5 files may be viewed using <a href="https://support.hdfgroup.org/products/java/hdfview/">hdfview</a>
 *         or loaded in python using <a href="http://www.pytables.org/">PyTables</a> or <a href="http://www.h5py.org/">h5py</a>.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk CreateReadCountPanelOfNormals \
 *          -I sample_1.counts.hdf5 \
 *          -I sample_2.counts.hdf5 \
 *          ... \
 *          -O cnv.pon.hdf5
 * </pre>
 *
 * <pre>
 *     gatk CreateReadCountPanelOfNormals \
 *          -I sample_1.counts.hdf5 \
 *          -I sample_2.counts.tsv \
 *          ... \
 *          --annotated-intervals annotated_intervals.tsv \
 *          -O cnv.pon.hdf5
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates a panel of normals for read-count denoising given the read counts for samples in the panel",
        oneLineSummary = "Creates a panel of normals for read-count denoising",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class CreateReadCountPanelOfNormals extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //parameter names
    public static final String MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME = "minimum-interval-median-percentile";
    public static final String MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME = "maximum-zeros-in-sample-percentage";
    public static final String MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME = "maximum-zeros-in-interval-percentage";
    public static final String EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME = "extreme-sample-median-percentile";
    public static final String IMPUTE_ZEROS_LONG_NAME = "do-impute-zeros";
    public static final String EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME = "extreme-outlier-truncation-percentile";
    public static final String MAXIMUM_CHUNK_SIZE = "maximum-chunk-size";

    //default values for filtering
    private static final double DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE = 10.0;
    private static final double DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE = 5.0;
    private static final double DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE = 5.0;
    private static final double DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE = 2.5;
    private static final boolean DEFAULT_DO_IMPUTE_ZEROS = true;
    private static final double DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE = 0.1;

    private static final int DEFAULT_NUMBER_OF_EIGENSAMPLES = 20;
    private static final int DEFAULT_CHUNK_DIVISOR = 16;
    private static final int DEFAULT_MAXIMUM_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / DEFAULT_CHUNK_DIVISOR;

    @Argument(
            doc = "Input TSV or HDF5 files containing integer read counts in genomic intervals for all samples in the panel of normals (output of CollectReadCounts).  " +
                    "Intervals must be identical and in the same order for all samples.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Input file containing annotations for GC content in genomic intervals (output of AnnotateIntervals).  " +
                    "If provided, explicit GC correction will be performed before performing SVD.  " +
                    "Intervals must be identical to and in the same order as those in the input read-counts files.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAnnotatedIntervalsFile = null;

    @Argument(
            doc = "Output file for the panel of normals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputPanelOfNormalsFile;

    @Argument(
            doc = "Genomic intervals with a median (across samples) of fractional coverage (optionally corrected for GC bias) " +
                    "less than or equal to this percentile are filtered out.  " +
                    "(This is the first filter applied.)",
            fullName = MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double minimumIntervalMedianPercentile = DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE;

    @Argument(
            doc = "Samples with a fraction of zero-coverage genomic intervals above this percentage are filtered out.  " +
                    "(This is the second filter applied.)",
            fullName = MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double maximumZerosInSamplePercentage = DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE;

    @Argument(
            doc = "Genomic intervals with a fraction of zero-coverage samples above this percentage are filtered out.  " +
                    "(This is the third filter applied.)",
            fullName = MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double maximumZerosInIntervalPercentage = DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE;

    @Argument(
            doc = "Samples with a median (across genomic intervals) of fractional coverage normalized by genomic-interval medians  " +
                    "below this percentile or above the complementary percentile are filtered out.  " +
                    "(This is the fourth filter applied.)",
            fullName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME,
            minValue = 0.,
            maxValue = 50.,
            optional = true
    )
    private double extremeSampleMedianPercentile = DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE;

    @Argument(
            doc = "If true, impute zero-coverage values as the median of the non-zero values in the corresponding interval.  " +
                    "(This is applied after all filters.)",
            fullName = IMPUTE_ZEROS_LONG_NAME,
            optional = true
    )
    private boolean doImputeZeros = DEFAULT_DO_IMPUTE_ZEROS;

    @Argument(
            doc = "Fractional coverages normalized by genomic-interval medians that are " +
                    "below this percentile or above the complementary percentile are set to the corresponding percentile value.  " +
                    "(This is applied after all filters and imputation.)",
            fullName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME,
            minValue = 0.,
            maxValue = 50.,
            optional = true
    )
    private double extremeOutlierTruncationPercentile = DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE;

    @Argument(
            doc = "Number of eigensamples to use for truncated SVD and to store in the panel of normals.  " +
                    "The number of samples retained after filtering will be used instead if it is smaller than this.",
            fullName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int numEigensamplesRequested = DEFAULT_NUMBER_OF_EIGENSAMPLES;

    @Advanced
    @Argument(
            doc = "Maximum HDF5 matrix chunk size.  Large matrices written to HDF5 are chunked into equally sized " +
                    "subsets of rows (plus a subset containing the remainder, if necessary) to avoid a hard limit in " +
                    "Java HDF5 on the number of elements in a matrix.  However, since a single row is not allowed to " +
                    "be split across multiple chunks, the number of columns must be less than the maximum number of " +
                    "values in each chunk.  Decreasing this number will reduce heap usage when writing chunks.",
            fullName = MAXIMUM_CHUNK_SIZE,
            minValue = 1,
            maxValue = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX,
            optional = true
    )
    private int maximumChunkSize = DEFAULT_MAXIMUM_CHUNK_SIZE;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        //validate parameters and parse optional parameters
        validateArguments();

        //get sample filenames
        final List<String> sampleFilenames = inputReadCountFiles.stream().map(File::getAbsolutePath).collect(Collectors.toList());

        //get sequence dictionary and intervals from the first read-counts file to use to validate remaining files
        //(this first file is read again below, which is slightly inefficient but is probably not worth the extra code)
        final File firstReadCountFile = inputReadCountFiles.get(0);
        logger.info(String.format("Retrieving intervals from first read-counts file (%s)...", firstReadCountFile));
        final SimpleCountCollection firstReadCounts = SimpleCountCollection.read(firstReadCountFile);
        final SAMSequenceDictionary sequenceDictionary = firstReadCounts.getMetadata().getSequenceDictionary();
        final List<SimpleInterval> intervals = firstReadCounts.getIntervals();
        Utils.validateArg(firstReadCounts.size() <= maximumChunkSize,
                String.format("The number of intervals (%d) in each read-counts file cannot exceed the maximum chunk size (%d).",
                        firstReadCounts.size(), maximumChunkSize));

        //get GC content (null if not provided)
        final AnnotatedIntervalCollection annotatedIntervals = CopyNumberArgumentValidationUtils.validateAnnotatedIntervals(
                inputAnnotatedIntervalsFile, firstReadCounts, logger);
        final double[] intervalGCContent = annotatedIntervals == null
                ? null
                : annotatedIntervals.getRecords().stream()
                    .mapToDouble(i -> i.getAnnotationMap().getValue(CopyNumberAnnotations.GC_CONTENT))
                    .toArray();

        //validate input read-counts files (i.e., check intervals and that only integer counts are contained)
        //and aggregate as a RealMatrix with dimensions numIntervals x numSamples
        final RealMatrix readCountMatrix = constructReadCountMatrix(logger, inputReadCountFiles, sequenceDictionary, intervals);

        //create the PoN
        logger.info("Creating the panel of normals...");
        HDF5SVDReadCountPanelOfNormals.create(outputPanelOfNormalsFile, getCommandLine(),
                sequenceDictionary, readCountMatrix, sampleFilenames, intervals, intervalGCContent,
                minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                extremeSampleMedianPercentile, doImputeZeros, extremeOutlierTruncationPercentile, numEigensamplesRequested,
                maximumChunkSize, ctx);

        logger.info("Panel of normals successfully created.");
    }

    private void validateArguments() {
        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-counts files cannot contain duplicates.");
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        if (numEigensamplesRequested > inputReadCountFiles.size()) {
            logger.warn(String.format("Number of eigensamples (%d) is greater than the number of input samples (%d); " +
                            "the number of samples retained after filtering will be used instead.",
                    numEigensamplesRequested, inputReadCountFiles.size()));
        }
    }

    private static RealMatrix constructReadCountMatrix(final Logger logger,
                                                       final List<File> inputReadCountFiles,
                                                       final SAMSequenceDictionary sequenceDictionary,
                                                       final List<SimpleInterval> intervals) {
        logger.info("Validating and aggregating input read-counts files...");
        final int numSamples = inputReadCountFiles.size();
        final int numIntervals = intervals.size();
        final RealMatrix readCountMatrix = new Array2DRowRealMatrix(numSamples, numIntervals);
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-counts file %s (%d / %d)", inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);
            if (!CopyNumberArgumentValidationUtils.isSameDictionary(readCounts.getMetadata().getSequenceDictionary(), sequenceDictionary)) {
                logger.warn(String.format("Sequence dictionary for read-counts file %s does not match those in other read-counts files.", inputReadCountFile));
            }
            Utils.validateArg(readCounts.getIntervals().equals(intervals),
                    String.format("Intervals for read-counts file %s do not match those in other read-counts files.", inputReadCountFile));
            readCountMatrix.setRow(sampleIndex, readCounts.getCounts());
        }
        return readCountMatrix;
    }
}
