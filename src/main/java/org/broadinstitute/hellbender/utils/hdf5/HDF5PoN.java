package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.util.*;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;

/**
 * HDF5 File backed Panel of Normals data structure.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5PoN implements PoN {

    private final static String TARGET_FACTORS_GROUP_NAME = "/target_factors";

    private final static String NORMALIZED_PCOV_GROUP_NAME = "/fnt_control_matrix";

    private final static String MAXIMUM_RATIO_CUTOFF_GROUP_NAME = "/max_ratio_cutoff";

    private final static String VERSION_GROUP_NAME = "/version";

    private final static String LOG_NORMALS_GROUP_NAME = "/log_normals";

    private final static String LOG_NORMALS_PINV_GROUP_NAME = "/log_normals_pinv";

    private final static String REDUCED_PON_GROUP_NAME = "/reduced_pon";

    private final static String REDUCED_PON_PINV_GROUP_NAME = "/reduced_pon_pinv";

    private final static String MAXIMUM_RATIO_CUTOFF_PATH = MAXIMUM_RATIO_CUTOFF_GROUP_NAME + "/values";

    private final static String TARGET_NAMES_PATH = TARGET_FACTORS_GROUP_NAME + "/index";

    private final static String SAMPLE_NAMES_PATH = NORMALIZED_PCOV_GROUP_NAME + "/axis0";

    private final static String TARGET_FACTORS_PATH = TARGET_FACTORS_GROUP_NAME + "/values";

    private final static String NORMALIZED_PCOV_PATH = NORMALIZED_PCOV_GROUP_NAME + "/block0_values";

    private final static String LOG_NORMALS_PATH = LOG_NORMALS_GROUP_NAME + "/block0_values";

    private final static String LOG_NORMALS_SAMPLE_NAMES_PATH = LOG_NORMALS_GROUP_NAME + "/block0_items";

    private final static String LOG_NORMALS_PINV_PATH = LOG_NORMALS_PINV_GROUP_NAME + "/block0_values";

    private final static String REDUCED_PON_PATH = REDUCED_PON_GROUP_NAME + "/block0_values";

    private final static String REDUCED_PON_PINV_PATH = REDUCED_PON_PINV_GROUP_NAME + "/block0_values";

    private final static String VERSION_PATH = VERSION_GROUP_NAME + "/values";

    private final HDF5Reader reader;

    private final List<String> targetNames;

    private final List<String> sampleNames;

    private List<String> logNormalSampleNames;

    public HDF5PoN(final HDF5Reader reader) {
        if (reader == null) {
            throw new IllegalArgumentException("reader is null");
        }
        targetNames = readTargetNames(reader);
        sampleNames = readSampleNames(reader);
        this.reader = reader;
    }

    /**
     * Checks that all sample names in {@code logNormalSampleNames} are present in {@code sampleNames}.
     *
     * <p>
     *     If this is not the case a {@link GATKException} is thrown.
     * </p>
     *
     * @param file the original file, used for the exception message if applies.
     * @param sampleNames the full sample name list.
     * @param logNormalSampleNames the log-normal sample name sub-set.
     * @throws GATKException if elements in {@code logNormalSampleNames} are not present in {@code sampleNames}.
     */
    private static void checkSampleNames(final File file, final List<String> sampleNames, final List<String> logNormalSampleNames) {
        final Set<String> sampleNamesSet = new HashSet<>(sampleNames);
        if (logNormalSampleNames.stream().anyMatch(n -> ! sampleNamesSet.contains(n))) {
            throw new GATKException(
                    String.format("the log-normal sample subset contains samples that are not present in the full sample set for file '%s'; e.g. %s",file,
                            logNormalSampleNames.stream().filter(n -> ! sampleNamesSet.contains(n))
                                    .limit(10).collect(Collectors.joining(", "))));
        }
    }

    /**
     * Reads the log-normal sample names sub-set.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readLogNormalSampleNames(final HDF5Reader reader) {
        final String[] values = reader.readStringArray(LOG_NORMALS_SAMPLE_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    /**
     * Reads the target names.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readTargetNames(final HDF5Reader reader) {
        final String[] values = reader.readStringArray(TARGET_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    /**
     * Reads the sample names.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private List<String> readSampleNames(final HDF5Reader reader) {
        final String[] values = reader.readStringArray(SAMPLE_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    @Override
    public List<String> targetNames() {
        return targetNames;
    }

    @Override
    public List<String> sampleNames() {
        return sampleNames;
    }

    @Override
    public List<String> logNormalSampleNames() {
        if (logNormalSampleNames == null) {
            logNormalSampleNames = readLogNormalSampleNames(reader);
            checkSampleNames(reader.getFile(),sampleNames,logNormalSampleNames);
        }
        return logNormalSampleNames;
    }

    @Override
    public RealMatrix targetFactors() {
        final double[] values = reader.readDoubleArray(TARGET_FACTORS_PATH);
        if (values.length != targetNames.size()) {
            throw new GATKException(String.format("wrong number of elements in the target factors recovered from file '%s': %d != %d",reader.getFile(),values.length,targetNames.size()));
        }
        return new Array2DRowRealMatrix(values);
    }

    @Override
    public RealMatrix normalizedPercentCoverage() {
        return readMatrixAndCheckDimensions(NORMALIZED_PCOV_PATH,targetNames.size(),
                sampleNames.size());
    }

    @Override
    public RealMatrix logNormals() {
        return readMatrixAndCheckDimensions(LOG_NORMALS_PATH,targetNames.size(),
                sampleNames.size());
    }

    @Override
    public RealMatrix logNormalsPseudoInverse() {
        return readMatrixAndCheckDimensions(LOG_NORMALS_PINV_PATH,sampleNames.size(),
                targetNames.size());
    }

    @Override
    public RealMatrix reducedPoN() {
        return readMatrixAndCheckDimensions(REDUCED_PON_PATH,
                r -> r == targetNames.size(),
                c -> c <= logNormalSampleNames().size());
    }

    @Override
    public RealMatrix reducedPoNPseudoInverse() {
        return readMatrixAndCheckDimensions(REDUCED_PON_PINV_PATH,
                r -> r <= logNormalSampleNames().size(),
                c -> c == targetNames.size());
    }

    @Override
    public double maximumRatioCutoff() {
        return reader.readDouble(MAXIMUM_RATIO_CUTOFF_PATH);
    }

    @Override
    public String version() {
        return String.format("%.1f", reader.readDouble(VERSION_PATH));
    }

    /**
     * Reads a matrix from the underlying PoN file.
     * @param fullPath the full path to the matrix data-set within the HDF5 file.
     * @return never {@code null}.
     * @throws GATKException if the matrix does not exist or any other HDF5 level error occurred.
     */
    private RealMatrix readMatrix(final String fullPath) {
        final double[][] values = reader.readDoubleMatrix(fullPath);
        return new Array2DRowRealMatrix(values,false);
    }

    /**
     * Reads a matrix from the underlying PoN file and check its dimensions.
     * @param fullPath the target data-set full path within the HDF5 file.
     * @param expectedRowCount the expected number of rows.
     * @param expectedColumnCount the expected number of columns.
     * @return GATKException if the result matrix dimensions do not match the expectations or
     *  any other cause as described in {@link #readMatrix(String)}.
     */
    private RealMatrix readMatrixAndCheckDimensions(final String fullPath, final int expectedRowCount, final int expectedColumnCount) {
        return readMatrixAndCheckDimensions(fullPath, r -> r == expectedRowCount, c -> c == expectedColumnCount);
    }

    /**
     * Reads a matrix from the underlying PoN file and check its dimensions.
     * @param fullPath the target data-set full path within the HDF5 file.
     * @param expectedRowCount a predicate that returns true iff its argument is an expected number of rows.
     * @param expectedColumnCount a predicate that returns true iff its argument is an expected number of columns.
     * @return GATKException if the result matrix dimensions do not match the expectations or
     *  any other cause as described in {@link #readMatrix(String)}.
     */
    private RealMatrix readMatrixAndCheckDimensions(final String fullPath, final IntPredicate expectedRowCount, final IntPredicate expectedColumnCount) {
        final RealMatrix result = readMatrix(fullPath);
        if (!expectedRowCount.test(result.getRowDimension())) {
            throw new GATKException(String.format("wrong number of rows in '%s' matrix from file '%s': %d",
                    fullPath, reader.getFile(), result.getRowDimension()));
        }
        if (!expectedColumnCount.test(result.getColumnDimension())) {
            throw new GATKException(String.format("wrong number of columns in '%s' from file '%s': %d",
                    fullPath, reader.getFile(), result.getColumnDimension()));
        }
        return result;
    }
}
