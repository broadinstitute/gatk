package org.broadinstitute.hellbender.tools.copynumber.denoising;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.SingularValueDecomposition;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.CreateReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.io.File;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents the SVD panel of normals to be created by {@link CreateReadCountPanelOfNormals}.
 *
 * <p>
 *     Data is stored in the following HDF5 paths:
 * </p>
 * <ul>
 *     <li>
 *         version: /version/value
 *     </li>
 *     <li>
 *         command line: /command_line/value
 *     </li>
 *     <li>
 *         sequence dictionary: /sequence_dictionary/value
 *     </li>
 *     <li>
 *         original read counts: /original_data/read_counts_samples_by_intervals
 *     </li>
 *     <li>
 *         original sample filenames: /original_data/sample_filenames
 *     </li>
 *     <li>
 *         original intervals: /original_data/intervals
 *     </li>
 *     <li>
 *         interval GC content (optional): /original_data/interval_gc_content
 *     </li>
 *     <li>
 *         panel sample filenames (after filtering): /panel/sample_filenames
 *     </li>
 *     <li>
 *         panel intervals (after filtering): /panel/intervals
 *     </li>
 *     <li>
 *         panel interval fractional medians: /panel/interval_fractional_medians
 *     </li>
 *     <li>
 *         panel singular values: /panel/singular_values
 *     </li>
 *     <li>
 *         panel interval fractional medians: /panel/interval_fractional_medians
 *     </li>
 *     <li>
 *         panel eigensamples: /panel/transposed_eigensamples_samples_by_intervals
 *     </li>
 * </ul>
 * <p>
 *     Most attributes are stored as wide matrices (i.e., more columns than rows) when possible.
 *     This avoids a very slow write time in HDF5, since HDF5 writes wide matrices much faster than tall matrices.
 *     See {@link HDF5Utils#writeIntervals} for details on the representation of intervals.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5SVDReadCountPanelOfNormals implements SVDReadCountPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(HDF5SVDReadCountPanelOfNormals.class);

    private static final int CHUNK_DIVISOR = 16;    //limits number of intervals to 16777215
    private static final int NUM_SLICES_FOR_SPARK_MATRIX_CONVERSION = 100;
    private static final double EPSILON = 1E-9;

    /**
     * The version number is a double where the integer part is the
     * major version and the decimal part is the minor version.
     * The minor version should be only a single digit.
     */
    private static final double CURRENT_PON_VERSION = 7.0;
    private static final String PON_VERSION_STRING_FORMAT = "%.1f";

    private static final String VERSION_PATH = "/version/value";    //note that full path names must include a top-level group name ("version" here)
    private static final String COMMAND_LINE_PATH = "/command_line/value";
    private static final String SEQUENCE_DICTIONARY_PATH = "/sequence_dictionary/value";
    private static final String ORIGINAL_DATA_GROUP_NAME = "/original_data";
    private static final String ORIGINAL_READ_COUNTS_PATH = ORIGINAL_DATA_GROUP_NAME + "/read_counts_samples_by_intervals";
    private static final String ORIGINAL_SAMPLE_FILENAMES_PATH = ORIGINAL_DATA_GROUP_NAME + "/sample_filenames";
    private static final String ORIGINAL_INTERVALS_PATH = ORIGINAL_DATA_GROUP_NAME + "/intervals";
    private static final String ORIGINAL_INTERVAL_GC_CONTENT_PATH = ORIGINAL_DATA_GROUP_NAME + "/interval_gc_content";

    private static final String PANEL_GROUP_NAME = "/panel";
    private static final String PANEL_SAMPLE_FILENAMES_PATH = PANEL_GROUP_NAME + "/sample_filenames";
    private static final String PANEL_INTERVALS_PATH = PANEL_GROUP_NAME + "/intervals";
    private static final String PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH = PANEL_GROUP_NAME + "/interval_fractional_medians";
    private static final String PANEL_SINGULAR_VALUES_PATH = PANEL_GROUP_NAME + "/singular_values";
    private static final String PANEL_EIGENSAMPLE_VECTORS_PATH = PANEL_GROUP_NAME + "/transposed_eigensamples_samples_by_intervals";
    private static final String PANEL_NUM_EIGENSAMPLES_PATH = PANEL_EIGENSAMPLE_VECTORS_PATH + HDF5Utils.NUMBER_OF_ROWS_SUB_PATH;

    private final HDF5File file;

    private final Lazy<SAMSequenceDictionary> sequenceDictionary;
    private final Lazy<List<SimpleInterval>> originalIntervals;
    private final Lazy<List<SimpleInterval>> panelIntervals;

    /**
     * DEV NOTE: If you are adding attributes that are neither RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * Otherwise, some operations will hang.
     */
    private HDF5SVDReadCountPanelOfNormals(final HDF5File file) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
        this.file = file;
        sequenceDictionary = new Lazy<>(() -> {
            final String sequenceDictionaryString = file.readStringArray(SEQUENCE_DICTIONARY_PATH)[0];
            return new SAMTextHeaderCodec()
                    .decode(BufferedLineReader.fromString(sequenceDictionaryString), file.getFile().getAbsolutePath())
                    .getSequenceDictionary();
        });
        originalIntervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, ORIGINAL_INTERVALS_PATH));
        panelIntervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, PANEL_INTERVALS_PATH));
    }

    @Override
    public double getVersion() {
        if (!file.isPresent(VERSION_PATH)) {    //the version path may be different in older PoNs
            throw new UserException.BadInput(String.format("The panel of normals is out of date and incompatible.  " +
                    "Please use a panel of normals that was created by CreateReadCountPanelOfNormals and is version " +
                    PON_VERSION_STRING_FORMAT + ".", CURRENT_PON_VERSION));
        }
        return file.readDouble(VERSION_PATH);
    }

    @Override
    public int getNumEigensamples() {
        try {
            return (int) file.readDouble(PANEL_NUM_EIGENSAMPLES_PATH);
        } catch (final HDF5LibException e) {
            //if the panel only contains a single sample or zero eigensamples were requested,
            //we do not store singular values or eigenvectors in the panel;
            //in this case, no value is written to PANEL_NUM_EIGENSAMPLES_PATH and an exception
            //is thrown when we attempt to read it
            return 0;
        }
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary.get();
    }

    @Override
    public double[][] getOriginalReadCounts() {
        return HDF5Utils.readChunkedDoubleMatrix(file, ORIGINAL_READ_COUNTS_PATH);
    }

    @Override
    public List<SimpleInterval> getOriginalIntervals() {
        return originalIntervals.get();
    }

    @Override
    public double[] getOriginalIntervalGCContent() {
        if (!file.isPresent(ORIGINAL_INTERVAL_GC_CONTENT_PATH)) {
            return null;
        }
        return file.readDoubleArray(ORIGINAL_INTERVAL_GC_CONTENT_PATH);
    }

    @Override
    public List<SimpleInterval> getPanelIntervals() {
        return panelIntervals.get();
    }

    @Override
    public double[] getPanelIntervalFractionalMedians() {
        return file.readDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH);
    }

    @Override
    public double[] getSingularValues() {
        if (getNumEigensamples() == 0) {
            throw new UnsupportedOperationException("No singular values were available.  " +
                    "This is because the panel only contains a single sample or no eigensamples were requested upon panel creation.");
        }
        return file.readDoubleArray(PANEL_SINGULAR_VALUES_PATH);
    }

    @Override
    public double[][] getEigensampleVectors() {
        if (getNumEigensamples() == 0) {
            throw new UnsupportedOperationException("No eigensample vectors were available.  " +
                    "This is because the panel only contains a single sample or no eigensamples were requested upon panel creation.");
        }
        return new Array2DRowRealMatrix(
                HDF5Utils.readChunkedDoubleMatrix(file, PANEL_EIGENSAMPLE_VECTORS_PATH), false)
                .transpose().getData();
    }

    /**
     * Create an interface to an HDF5 file.  A version check is performed and a warning message logged if the
     * version number is not up to date.
     */
    public static HDF5SVDReadCountPanelOfNormals read(final HDF5File file) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
        final HDF5SVDReadCountPanelOfNormals pon = new HDF5SVDReadCountPanelOfNormals(file);
        if (pon.getVersion() < CURRENT_PON_VERSION) {
            throw new UserException.BadInput(String.format("The version of the specified panel of normals (%f) is older than the current version (%f).",
                    pon.getVersion(), CURRENT_PON_VERSION));
        }
        return pon;
    }

    /**
     * Create the panel of normals and write it to an HDF5 file.  All inputs are assumed to be valid.
     * The dimensions of {@code originalReadCounts} should be samples x intervals.
     * To reduce memory footprint, {@code originalReadCounts} is modified in place.
     * If {@code intervalGCContent} is null, GC-bias correction will not be performed.
     */
    public static void create(final File outFile,
                              final String commandLine,
                              final SAMSequenceDictionary sequenceDictionary,
                              final RealMatrix originalReadCounts,
                              final List<String> originalSampleFilenames,
                              final List<SimpleInterval> originalIntervals,
                              final double[] intervalGCContent,
                              final double minimumIntervalMedianPercentile,
                              final double maximumZerosInSamplePercentage,
                              final double maximumZerosInIntervalPercentage,
                              final double extremeSampleMedianPercentile,
                              final boolean doImputeZeros,
                              final double extremeOutlierTruncationPercentile,
                              final int numEigensamplesRequested,
                              final int maximumChunkSize,
                              final JavaSparkContext ctx) {
        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            logger.info("Creating " + outFile.getAbsolutePath() + "...");
            final HDF5SVDReadCountPanelOfNormals pon = new HDF5SVDReadCountPanelOfNormals(file);

            logger.info(String.format("Writing version number (" + PON_VERSION_STRING_FORMAT + ")...", CURRENT_PON_VERSION));
            pon.writeVersion(CURRENT_PON_VERSION);

            logger.info("Writing command line...");
            pon.writeCommandLine(commandLine);

            logger.info("Writing sequence dictionary...");
            pon.writeSequenceDictionary(sequenceDictionary);

            logger.info(String.format("Writing original read counts (%d x %d)...",
                    originalReadCounts.getColumnDimension(), originalReadCounts.getRowDimension()));
            pon.writeOriginalReadCountsPath(originalReadCounts, maximumChunkSize);

            logger.info(String.format("Writing original sample filenames (%d)...", originalSampleFilenames.size()));
            pon.writeOriginalSampleFilenames(originalSampleFilenames);

            logger.info(String.format("Writing original intervals (%d)...", originalIntervals.size()));
            pon.writeOriginalIntervals(originalIntervals);

            if (intervalGCContent != null) {
                logger.info(String.format("Writing GC-content annotations for original intervals (%d)...", intervalGCContent.length));
                pon.writeOriginalIntervalGCContent(intervalGCContent);
            }

            //preprocess and standardize read counts and determine filters
            //(originalReadCounts is modified in place and a filtered submatrix is returned)
            logger.info("Preprocessing and standardizing read counts...");
            final SVDDenoisingUtils.PreprocessedStandardizedResult preprocessedStandardizedResult =
                    SVDDenoisingUtils.preprocessAndStandardizePanel(originalReadCounts, intervalGCContent,
                            minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                            extremeSampleMedianPercentile, doImputeZeros, extremeOutlierTruncationPercentile);

            //filter samples and intervals
            final List<String> panelSampleFilenames = IntStream.range(0, originalSampleFilenames.size())
                    .filter(sampleIndex -> !preprocessedStandardizedResult.filterSamples[sampleIndex])
                    .mapToObj(originalSampleFilenames::get).collect(Collectors.toList());
            final List<SimpleInterval> panelIntervals = IntStream.range(0, originalIntervals.size())
                    .filter(intervalIndex -> !preprocessedStandardizedResult.filterIntervals[intervalIndex])
                    .mapToObj(originalIntervals::get).collect(Collectors.toList());

            logger.info(String.format("Writing panel sample filenames (%d)...", panelSampleFilenames.size()));
            pon.writePanelSampleFilenames(panelSampleFilenames);

            logger.info(String.format("Writing panel intervals (%d)...", panelIntervals.size()));
            pon.writePanelIntervals(panelIntervals);

            //get panel interval fractional medians (calculated as an intermediate result during preprocessing)
            final double[] panelIntervalFractionalMedians = preprocessedStandardizedResult.panelIntervalFractionalMedians;

            logger.info(String.format("Writing panel interval fractional medians (%d)...", panelIntervalFractionalMedians.length));
            pon.writePanelIntervalFractionalMedians(panelIntervalFractionalMedians);

            final int numPanelSamples = preprocessedStandardizedResult.preprocessedStandardizedValues.getRowDimension();
            final int numPanelIntervals = preprocessedStandardizedResult.preprocessedStandardizedValues.getColumnDimension();

            //perform SVD, handling number of eigensamples requested vs. that available in filtered panel vs. that available from actual decomposition
            final int numEigensamples = Math.min(numEigensamplesRequested, numPanelSamples);
            if (numEigensamples < numEigensamplesRequested) {
                logger.warn(String.format("%d eigensamples were requested but only %d are available in the panel of normals...",
                        numEigensamplesRequested, numEigensamples));
            }
            logger.info(String.format("Performing SVD (truncated at %d eigensamples) of standardized counts (transposed to %d x %d)...",
                    numEigensamples, numPanelIntervals, numPanelSamples));
            if (numPanelSamples > 1 && numEigensamples > 0) {
                final SingularValueDecomposition<RowMatrix, Matrix> svd = SparkConverter.convertRealMatrixToSparkRowMatrix(
                        ctx, preprocessedStandardizedResult.preprocessedStandardizedValues.transpose(), NUM_SLICES_FOR_SPARK_MATRIX_CONVERSION)
                        .computeSVD(numEigensamples, true, EPSILON);
                final double[] singularValues = svd.s().toArray();    //should be in decreasing order (with corresponding matrices below)
                if (singularValues.length == 0 || Arrays.stream(singularValues).noneMatch(s -> s > EPSILON)) {
                    //if the panel contains more than one sample, we require that at least one non-negligible singular value is found
                    throw new UserException(String.format("No non-zero singular values were found.  It may be necessary to use stricter parameters for filtering.  " +
                            "For example, use a larger value of %s.", CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME));
                }
                if (singularValues.length < numEigensamples) {
                    logger.warn(String.format("Attempted to truncate at %d eigensamples, but only %d non-zero singular values were found...",
                            numEigensamples, singularValues.length));
                }
                final double[][] eigensampleVectors = SparkConverter.convertSparkRowMatrixToRealMatrix(svd.U(), numPanelIntervals).getData();

                logger.info(String.format("Writing singular values (%d)...", singularValues.length));
                pon.writeSingularValues(singularValues);

                logger.info(String.format("Writing eigensample vectors (transposed to %d x %d)...", eigensampleVectors[0].length, eigensampleVectors.length));
                pon.writeEigensampleVectors(eigensampleVectors, maximumChunkSize);
            } else {
                //if the panel only contains a single sample or zero eigensamples were requested,
                //we do not store singular values or eigenvectors in the panel
                logger.info("No eigensamples could be computed because only a single sample was provided or no eigensamples were requested.");
            }
        } catch (final RuntimeException exception) {
            //if any exceptions encountered, delete partial output and rethrow
            logger.warn(String.format("Exception encountered during creation of panel of normals (%s).  Attempting to delete partial output in %s...",
                    exception, outFile.getAbsolutePath()));
            IOUtils.tryDelete(outFile);
            throw new GATKException(String.format("Could not create panel of normals.  It may be necessary to use stricter parameters for filtering.  " +
                    "For example, use a larger value of %s.", CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME),  exception);
        }
        logger.info(String.format("Read-count panel of normals written to %s.", outFile));
    }

    //PRIVATE WRITERS (write values to HDF5 file)
    //these are private to prevent fields from being written individually, which could leave the file in a bad state

    private void writeVersion(final double version) {
        file.makeDouble(VERSION_PATH, version);
    }

    private void writeCommandLine(final String commandLine) {
        file.makeStringArray(COMMAND_LINE_PATH, commandLine);
    }

    private void writeSequenceDictionary(final SAMSequenceDictionary sequenceDictionary) {
        final StringWriter stringWriter = new StringWriter();
        new SAMTextHeaderCodec().encode(stringWriter, new SAMFileHeader(sequenceDictionary));
        file.makeStringArray(SEQUENCE_DICTIONARY_PATH, stringWriter.toString());
    }

    private void writeOriginalReadCountsPath(final RealMatrix originalReadCounts,
                                             final int maximumChunkSize) {
        HDF5Utils.writeChunkedDoubleMatrix(file, ORIGINAL_READ_COUNTS_PATH, originalReadCounts.getData(), maximumChunkSize);
    }

    private void writeOriginalSampleFilenames(final List<String> originalSampleFilenames) {
        file.makeStringArray(ORIGINAL_SAMPLE_FILENAMES_PATH, originalSampleFilenames.toArray(new String[originalSampleFilenames.size()]));
    }

    private void writeOriginalIntervals(final List<SimpleInterval> originalIntervals) {
        HDF5Utils.writeIntervals(file, ORIGINAL_INTERVALS_PATH, originalIntervals);
    }

    private void writeOriginalIntervalGCContent(final double[] originalIntervalGCContent) {
        file.makeDoubleArray(ORIGINAL_INTERVAL_GC_CONTENT_PATH, originalIntervalGCContent);
    }

    private void writePanelSampleFilenames(final List<String> panelSampleFilenames) {
        file.makeStringArray(PANEL_SAMPLE_FILENAMES_PATH, panelSampleFilenames.toArray(new String[panelSampleFilenames.size()]));
    }

    private void writePanelIntervals(final List<SimpleInterval> panelIntervals) {
        HDF5Utils.writeIntervals(file, PANEL_INTERVALS_PATH, panelIntervals);
    }

    private void writePanelIntervalFractionalMedians(final double[] panelIntervalFractionalMedians) {
        file.makeDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH, panelIntervalFractionalMedians);
    }

    private void writeSingularValues(final double[] singularValues) {
        file.makeDoubleArray(PANEL_SINGULAR_VALUES_PATH, singularValues);
    }

    private void writeEigensampleVectors(final double[][] eigensampleVectors,
                                         final int maximumChunkSize) {
        HDF5Utils.writeChunkedDoubleMatrix(file, PANEL_EIGENSAMPLE_VECTORS_PATH,
                new Array2DRowRealMatrix(eigensampleVectors, false).transpose().getData(), maximumChunkSize);
    }
}