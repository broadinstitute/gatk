package org.broadinstitute.hellbender.utils.hdf5;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableColumns;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.IntPredicate;

/**
 * HDF5 File backed Panel of Normals data structure.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5PoN implements PoN {

    public final static String TARGET_FACTORS_GROUP_NAME = "/target_factors";

    private final static String NORMALIZED_PCOV_GROUP_NAME = "/fnt_control_matrix";

    private final static String VERSION_GROUP_NAME = "/version";

    private final static String LOG_NORMALS_GROUP_NAME = "/log_normals";

    private final static String LOG_NORMALS_PINV_GROUP_NAME = "/log_normals_pinv";

    private final static String REDUCED_PON_GROUP_NAME = "/reduced_pon";

    public final static String REDUCED_PON_PINV_GROUP_NAME = "/reduced_pon_pinv";

    public final static String TARGET_NAMES_PATH = TARGET_FACTORS_GROUP_NAME + "/index";

    public final static String SAMPLE_NAMES_PATH = NORMALIZED_PCOV_GROUP_NAME + "/axis0";

    public final static String TARGET_FACTORS_PATH = TARGET_FACTORS_GROUP_NAME + "/values";

    private final static String NORMALIZED_PCOV_PATH = NORMALIZED_PCOV_GROUP_NAME + "/block0_values";

    public final static String LOG_NORMALS_PATH = LOG_NORMALS_GROUP_NAME + "/block0_values";

    private final static String LOG_NORMALS_SAMPLE_NAMES_PATH = LOG_NORMALS_GROUP_NAME + "/block0_items";

    public final static String LOG_NORMALS_PINV_PATH = LOG_NORMALS_PINV_GROUP_NAME + "/block0_values";

    public final static String REDUCED_PON_PATH = REDUCED_PON_GROUP_NAME + "/block0_values";

    private final static String REDUCED_PON_TARGET_PATH = REDUCED_PON_GROUP_NAME + "/axis1";

    public final static String REDUCED_PON_PINV_PATH = REDUCED_PON_PINV_GROUP_NAME + "/block0_values";

    private final static String VERSION_PATH = VERSION_GROUP_NAME + "/values";

    private final static String TARGETS_GROUP_NAME = "/targets";

    private final static String TARGETS_PATH = TARGETS_GROUP_NAME + "/block0_values";

    private final static String RAW_TARGETS_GROUP_NAME = "/raw_targets";

    private final static String RAW_TARGETS_PATH = RAW_TARGETS_GROUP_NAME + "/block0_values";

    public final static String RAW_TARGET_NAMES_PATH = RAW_TARGETS_GROUP_NAME + "/index";

    private final static String LOG_NORMALS_TARGETS_GROUP_NAME = "/log_normals_targets";

    private final static String LOG_NORMALS_TARGETS_PATH = LOG_NORMALS_TARGETS_GROUP_NAME + "/block0_values";

    private final static String TARGET_VARIANCES_GROUP_NAME = "/target_variances";
    private final static String TARGET_VARIANCES_PATH = TARGET_VARIANCES_GROUP_NAME + "/block0_values";

    private final HDF5File file;

    private final Lazy<List<String>> targetNames;

    private final Lazy<List<String>> rawTargetNames;
    private final Lazy<List<Target>> targets;
    private final Lazy<List<Target>> rawTargets;
    private final Lazy<List<Target>> panelTargets;

    private final Lazy<List<String>> reducedTargetNames;

    private final Lazy<List<String>> sampleNames;

    private final Lazy<List<String>> logNormalSampleNames;

    private final static Map<String, Integer> targetColumnToPoNIndex = ImmutableMap.of(
            TargetTableColumns.CONTIG.toString(), 0, TargetTableColumns.START.toString(), 1,
            TargetTableColumns.END.toString(), 2
    );

    private final static int NUM_TARGETS_COLS = targetColumnToPoNIndex.size();

    private final static String NUM_TARGETS_COLS_GROUP_NAME = "/num_target_cols";

    private final static String NUM_TARGETS_COLS_PATH = NUM_TARGETS_COLS_GROUP_NAME + "/values";


    /**
     * Create a new PoN interface to a HDF5 file.
     *
     * <p>DEV NOTE:  If you are adding attributes that are not RealMatrix nor a primitive, you must follow the pattern in the constructor (i.e. the Lazy loading pattern).  See the targetNames private attribute.  Otherwise, some operations will hang.</p>
     * @param file the underlying HDF5 file.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     */
    public HDF5PoN(final HDF5File file) {
        Utils.nonNull(file, "the input file must not be null");
        this.file = file;
        targetNames = new Lazy<>(() -> readTargetNames(file));
        sampleNames = new Lazy<>(() -> readSampleNames(file));
        logNormalSampleNames = new Lazy<>(() -> readLogNormalizedSampleNames(file));
        reducedTargetNames = new Lazy<>(() -> Collections.unmodifiableList(Arrays.asList(file.readStringArray(REDUCED_PON_TARGET_PATH))));
        targets  = new Lazy<>(() -> readTargets(file));
        rawTargets  = new Lazy<>(() -> readRawTargets(file));
        rawTargetNames  = new Lazy<>(() -> readRawTargetNames(file));
        panelTargets = new Lazy<>(() -> readPanelTargets(file));
    }

    /**
     * Reads the log-normalized sample names sub-set.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readLogNormalizedSampleNames(final HDF5File reader) {
        final String[] values = reader.readStringArray(LOG_NORMALS_SAMPLE_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    /**
     * Reads the target names.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readTargetNames(final HDF5File reader) {
        final String[] values = reader.readStringArray(TARGET_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }
    /**
     * Reads the raw target names.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readRawTargetNames(final HDF5File reader) {
        final String[] values = reader.readStringArray(RAW_TARGET_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    /**
     * Reads the sample names.
     * @param reader the source HDF5 reader.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readSampleNames(final HDF5File reader) {
        final String[] values = reader.readStringArray(SAMPLE_NAMES_PATH);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    @Override
    public List<String> getTargetNames() {
        return targetNames.get();
    }

    @Override
    public List<String> getRawTargetNames() {
        return rawTargetNames.get();
    }

    @Override
    public List<String> getPanelTargetNames() {
        return reducedTargetNames.get();
    }

    @Override
    public List<String> getSampleNames() {
        return sampleNames.get();
    }

    @Override
    public List<String> getPanelSampleNames() {
        return logNormalSampleNames.get();

    }

    @Override
    public RealMatrix getTargetFactors() {
        final double[] values = file.readDoubleArray(TARGET_FACTORS_PATH);
        if (values.length != targetNames.get().size()) {
            throw new GATKException(String.format("wrong number of elements in the target factors recovered from file '%s': %d != %d", file.getFile(), values.length, targetNames.get().size()));
        }
        return new Array2DRowRealMatrix(values);
    }

    @Override
    public double[] getTargetVariances() {
        final double[] values = file.readDoubleArray(TARGET_VARIANCES_PATH);
        if (values.length != panelTargets.get().size()) {
            throw new GATKException(String.format("wrong number of elements in the target variances recovered from file '%s': %d != %d", file.getFile(), values.length, panelTargets.get().size()));
        }
        return values;
    }

    @Override
    public List<Target> getTargets() {
        return targets.get();
    }

    @Override
    public List<Target> getRawTargets() {
        return rawTargets.get();
    }

    @Override
    public List<Target> getPanelTargets() {
        return panelTargets.get();
    }

    private static List<Target> readTargets(final HDF5File reader) {
        final String[][] values = reader.readStringMatrix(TARGETS_PATH, NUM_TARGETS_COLS_PATH);

        final List<String> targetNamesToRender = readTargetNames(reader);
        return renderPoNTargets(values, targetNamesToRender, reader);
    }

    private static List<Target> readPanelTargets(final HDF5File reader) {
        final String[][] values = reader.readStringMatrix(LOG_NORMALS_TARGETS_PATH, NUM_TARGETS_COLS_PATH);

        final List<String> targetNamesToRender = Arrays.asList(reader.readStringArray(REDUCED_PON_TARGET_PATH));
        return renderPoNTargets(values, targetNamesToRender, reader);
    }

    private static List<Target> readRawTargets(final HDF5File reader) {
        final String[][] values = reader.readStringMatrix(RAW_TARGETS_PATH, NUM_TARGETS_COLS_PATH);

        final List<String> targetNamesToRender = Arrays.asList(reader.readStringArray(RAW_TARGET_NAMES_PATH));
        return renderPoNTargets(values, targetNamesToRender, reader);
    }

    private static List<Target> renderPoNTargets(final String[][] values, final List<String> targetNamesToRender, final HDF5File reader) {
        if (values.length != targetNamesToRender.size()) {
            throw new GATKException(String.format("wrong number of elements in the targets recovered from file '%s': %d != %d", reader.getFile(), values.length, targetNamesToRender.size()));
        }

        final int numTargetCols = (int) reader.readDouble(NUM_TARGETS_COLS_PATH);

        final List<Target> result = new ArrayList<>(values.length);
        for (int i= 0; i < values.length; i ++) {
            if (values[i].length != numTargetCols) {
                throw new GATKException(String.format("wrong number of column elements in the targets recovered from file '%s': %d != %d", reader.getFile(), values[i].length, numTargetCols));
            }
            result.add(new Target(targetNamesToRender.get(i), new SimpleInterval(values[i][0], Integer.parseInt(values[i][1]), Integer.parseInt(values[i][2]))));
        }
        return result;
    }



    @Override
    public void setTargetFactors(final RealMatrix targetFactors) {
        Utils.nonNull(targetFactors);
        if (targetFactors.getColumnDimension() != 1) {
            throw new IllegalArgumentException("the number of columns in the target factors matrix must be 1 but it is: " + targetFactors.getColumnDimension());
        }
        file.makeGroup(TARGET_FACTORS_GROUP_NAME);
        file.makeDoubleArray(TARGET_FACTORS_PATH, targetFactors.getColumn(0));
    }

    @Override
    public RealMatrix getNormalizedCounts() {
        return readMatrixAndCheckDimensions(NORMALIZED_PCOV_PATH, targetNames.get().size(),
                sampleNames.get().size());
    }

    @Override
    public RealMatrix getLogNormalizedCounts() {
        return readMatrixAndCheckDimensions(LOG_NORMALS_PATH, getPanelTargetNames().size(),
                getPanelSampleNames().size());
    }

    @Override
    public RealMatrix getLogNormalizedPInverseCounts() {
        return readMatrixAndCheckDimensions(LOG_NORMALS_PINV_PATH, getPanelSampleNames().size(),
                getPanelTargetNames().size());
    }

    @Override
    public RealMatrix getReducedPanelCounts() {
        return readMatrixAndCheckDimensions(REDUCED_PON_PATH,
                r -> r == reducedTargetNames.get().size(),
                c -> c <= getPanelSampleNames().size());
    }

    @Override
    public RealMatrix getReducedPanelPInverseCounts() {
        return readMatrixAndCheckDimensions(REDUCED_PON_PINV_PATH,
                r -> r <= getPanelSampleNames().size(),
                c -> c == reducedTargetNames.get().size());
    }

    @Override
    public double getVersion() {
        return file.readDouble(VERSION_PATH);
    }

    /**
     * Reads a matrix from the underlying PoN file.
     * @param fullPath the full path to the matrix data-set within the HDF5 file.
     * @return never {@code null}.
     * @throws GATKException if the matrix does not exist or any other HDF5 level error occurred.
     */
    private RealMatrix readMatrix(final String fullPath) {
        final double[][] values = file.readDoubleMatrix(fullPath);
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
        if (expectedRowCount.test(result.getRowDimension())
                && expectedColumnCount.test(result.getColumnDimension())) {
            return result;
        }
        final RealMatrix transpose = result.transpose();
        if (!expectedRowCount.test(transpose.getRowDimension())) {
            throw new GATKException(String.format("wrong number of rows in '%s' matrix from file '%s': %d",
                    fullPath, file.getFile(), result.getRowDimension()));
        }
        if (!expectedColumnCount.test(transpose.getColumnDimension())) {
            throw new GATKException(String.format("wrong number of columns in '%s' from file '%s': %d",
                    fullPath, file.getFile(), result.getColumnDimension()));
        }
        return transpose;
    }

    /// Write interface:

    /**
     * The version number.
     *
     * The version number is a float point number (double) where the integer part is the
     * major and the decimal part is the minor.
     *
     * Note: the logic choice is to use a free form string but, the Python version was using
     * a HDF5 double, so we are keeping the tradition here.
     *
     * @param version the new version value.
     * @throws UnsupportedOperationException if this PoN instance is read-only access.
     */
    public void setVersion(final double version) {
        file.makeDouble(VERSION_PATH, version);
    }

    /**
     * Changes the sample names in the PoN.
     *
     * @param names the new list of sample names.
     * @throws IllegalArgumentException if {@code name} is {@code null} or contains
     *  any {@code null}.
     */
    public void setSampleNames(final List<String> names) {
        checkNameList(names);
        file.makeStringArray(SAMPLE_NAMES_PATH, names.toArray(new String[names.size()]));
    }

    /**
     * Changes the target names in the PoN.
     */
    public void setTargetNames(final List<String> names) {
        checkNameList(names);
        file.makeStringArray(TARGET_NAMES_PATH, names.toArray(new String[names.size()]));
    }

    /**
     * Changes the raw target names in the PoN.
     */
    public void setRawTargetNames(final List<String> names) {
        checkNameList(names);
        file.makeStringArray(RAW_TARGET_NAMES_PATH, names.toArray(new String[names.size()]));
    }

    /**
     * Set the normalized read counts.
     * @param normalizedCounts the normalized read counts.
     */
    public void setNormalCounts(final RealMatrix normalizedCounts) {
        file.makeDoubleMatrix(NORMALIZED_PCOV_PATH, normalizedCounts.getData());
    }

    public void setReducedPanelCounts(final RealMatrix counts) {
        Utils.nonNull(counts);
        file.makeDoubleMatrix(REDUCED_PON_PATH, counts.getData());
    }

    public void setLogNormalPInverseCounts(final RealMatrix counts) {
        Utils.nonNull(counts);
        file.makeDoubleMatrix(LOG_NORMALS_PINV_PATH, counts.getData());
    }

    public void setLogNormalCounts(final RealMatrix counts) {
        Utils.nonNull(counts);
        file.makeDoubleMatrix(LOG_NORMALS_PATH, counts.getData());
    }

    public void setReducedPanelPInverseCounts(final RealMatrix counts) {
        Utils.nonNull(counts);
        file.makeDoubleMatrix(REDUCED_PON_PINV_PATH, counts.getData());
    }

    public void setPanelSampleNames(final List<String> names) {
        checkNameList(names);
        file.makeStringArray(LOG_NORMALS_SAMPLE_NAMES_PATH, names.toArray(new String[names.size()]));
    }

    public void setPanelTargetNames(final List<String> names) {
        checkNameList(names);
        file.makeStringArray(REDUCED_PON_TARGET_PATH, names.toArray(new String[names.size()]));
    }

    private void checkNameList(final List<String> names) {
        Utils.nonNull(names, "the input names cannot be null");
        if (names.contains(null)) {
            throw new IllegalArgumentException("the input names list cannot contain a null");
        }
    }

    protected void setTargetValues(final List<Target> targets, final String fullPath) {
        Utils.nonNull(targets, "Cannot create targets with null input");

        final String[][] targetValues = new String[targets.size()][NUM_TARGETS_COLS];
        for (int i = 0; i < targets.size(); i++) {
            targetValues[i][targetColumnToPoNIndex.get(TargetTableColumns.CONTIG.toString())] = targets.get(i).getContig();
            targetValues[i][targetColumnToPoNIndex.get(TargetTableColumns.START.toString())] = String.valueOf(targets.get(i).getStart());
            targetValues[i][targetColumnToPoNIndex.get(TargetTableColumns.END.toString())] = String.valueOf(targets.get(i).getEnd());
        }
        file.makeStringMatrix(fullPath, targetValues, NUM_TARGETS_COLS_PATH);
    }

    public void setTargets(final List<Target> targets) {
        setTargetValues(targets, TARGETS_PATH);
    }

    public void setPanelTargets(final List<Target> targets) {
        setTargetValues(targets, LOG_NORMALS_TARGETS_PATH);
    }

    public void setRawTargets(final List<Target> targets) {
        setTargetValues(targets, RAW_TARGETS_PATH);
    }

    public void setTargetVariances(final double[] targetVariances) {
        if (targetVariances.length != panelTargets.get().size()) {
            throw new GATKException(String.format("Writing wrong number of elements in the target variances attempted to file '%s': %d != %d", file.getFile(), targetVariances.length, panelTargets.get().size()));
        }
        file.makeDoubleArray(TARGET_VARIANCES_PATH, targetVariances);
    }

}
