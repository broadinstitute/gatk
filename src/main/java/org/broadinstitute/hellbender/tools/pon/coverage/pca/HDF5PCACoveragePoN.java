package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableColumn;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;

/**
 * HDF5 File backed coverage panel of normals data structure.
 *
 * Several attributes are stored transposed (in other words, the rows and columns are interchanged).
 * This dodges a very slow write time in HDF5, since we usually have many more rows (targets) than columns (samples),
 * and HDF5 writes matrices with few rows and many columns much faster than matrices with many rows and few columns.
 *
 * The following are stored as transposes:
 *
 * <ul>
 *  <li>Normalized Counts</li>
 *  <li>Log-Normalized Counts</li>
 *  <li>Reduced Panel Counts</li>
 *</ul>
 *
 * In these cases, the samples are the rows and the targets are the columns.  No transposing is performed for
 * pseudoinverses since they already have dimensions of samples x targets.
 *
 * This is only for storage.  When saving/loading the above attributes, the transposing is handled transparently.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5PCACoveragePoN implements PCACoveragePoN {
    private final HDF5File file;

    private static final Logger logger = LogManager.getLogger(HDF5PCACoveragePoN.class);

    /**
     * The version number.
     *
     * The version number is a double where the integer part is the
     * major and the decimal part is the minor.
     *
     * Note: the logical choice is to use a free form string but, the Python version was using
     * a HDF5 double, so we are keeping the tradition here.
     */
    public static final double CURRENT_PON_VERSION = 6.0;

    private static final String VERSION_GROUP_NAME = "/version";
    private static final String VERSION_PATH = VERSION_GROUP_NAME + "/values";

    //targets
    private static final String TARGETS_GROUP_NAME = "/targets";
    private static final String TARGETS_PATH = TARGETS_GROUP_NAME + "/block0_values";

    private static final String RAW_TARGETS_GROUP_NAME = "/raw_targets";
    private static final String RAW_TARGETS_PATH = RAW_TARGETS_GROUP_NAME + "/block0_values";

    private static final String PANEL_TARGETS_GROUP_NAME = "/log_normals_targets";
    private static final String PANEL_TARGETS_PATH = PANEL_TARGETS_GROUP_NAME + "/block0_values";

    //target factors and variances
    private static final String TARGET_FACTORS_GROUP_NAME = "/target_factors";
    private static final String TARGET_FACTORS_PATH = TARGET_FACTORS_GROUP_NAME + "/values";

    private static final String TARGET_VARIANCES_GROUP_NAME = "/target_variances";
    private static final String TARGET_VARIANCES_PATH = TARGET_VARIANCES_GROUP_NAME + "/block0_values";

    //normalized/log-normalized counts and log-normalized pseudoinverse for full panel
    private static final String NORMALIZED_COUNTS_GROUP_NAME = "/fnt_control_matrix";
    private static final String NORMALIZED_COUNTS_PATH = NORMALIZED_COUNTS_GROUP_NAME + "/block0_values";

    private static final String LOG_NORMALIZED_COUNTS_GROUP_NAME = "/log_normals";
    private static final String LOG_NORMALIZED_COUNTS_PATH = LOG_NORMALIZED_COUNTS_GROUP_NAME + "/block0_values";

    private static final String LOG_NORMALIZED_PINV_GROUP_NAME = "/log_normals_pinv";
    private static final String LOG_NORMALIZED_PINV_PATH = LOG_NORMALIZED_PINV_GROUP_NAME + "/block0_values";

    //normalized/log-normalized counts and log-normalized pseudoinverse for reduced panel
    private static final String REDUCED_PON_GROUP_NAME = "/reduced_pon";
    private static final String REDUCED_PANEL_COUNTS_PATH = REDUCED_PON_GROUP_NAME + "/block0_values";

    private static final String REDUCED_PANEL_PINV_GROUP_NAME = "/reduced_pon_pinv";
    private static final String REDUCED_PANEL_PINV_PATH = REDUCED_PANEL_PINV_GROUP_NAME + "/block0_values";

    //target names
    private static final String TARGET_NAMES_PATH = TARGET_FACTORS_GROUP_NAME + "/index";
    private static final String RAW_TARGET_NAMES_PATH = RAW_TARGETS_GROUP_NAME + "/index";
    private static final String PANEL_TARGET_NAMES_PATH = REDUCED_PON_GROUP_NAME + "/axis1";

    //sample names
    private static final String SAMPLE_NAMES_PATH = NORMALIZED_COUNTS_GROUP_NAME + "/axis0";
    private static final String PANEL_SAMPLE_NAMES_PATH = LOG_NORMALIZED_COUNTS_GROUP_NAME + "/block0_items";

    private static final String NUM_TARGET_COLUMNS_GROUP_NAME = "/num_target_cols";
    private static final String NUM_TARGET_COLUMNS_PATH = NUM_TARGET_COLUMNS_GROUP_NAME + "/values";

    private static final EnumMap<TargetTableColumn, Integer> targetColumnToPoNIndex = new EnumMap<>(ImmutableMap.of(
            TargetTableColumn.CONTIG, 0, TargetTableColumn.START, 1, TargetTableColumn.END, 2));
    private static final int NUM_TARGET_COLUMNS = targetColumnToPoNIndex.size();

    private Lazy<List<Target>> targets;
    private Lazy<List<Target>> rawTargets;
    private Lazy<List<Target>> panelTargets;

    private Lazy<List<String>> targetNames;
    private Lazy<List<String>> rawTargetNames;
    private Lazy<List<String>> panelTargetNames;

    private Lazy<List<String>> sampleNames;
    private Lazy<List<String>> panelSampleNames;

    /*===============================================================================================================*
     * METHODS                                                                                                       *
     *===============================================================================================================*/

    /**
     * Create a new PoN interface to a HDF5 file.
     *
     * <p>DEV NOTE:  If you are adding attributes that are not RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * See the targetNames private attribute.  Otherwise, some operations will hang.</p>
     * @param file the underlying HDF5 file.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     */
    public HDF5PCACoveragePoN(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        targetNames = new Lazy<>(() -> readNames(file, TARGET_NAMES_PATH));
        rawTargetNames  = new Lazy<>(() -> readNames(file, RAW_TARGET_NAMES_PATH));
        panelTargetNames = new Lazy<>(() -> readNames(file, PANEL_TARGET_NAMES_PATH));
        targets  = new Lazy<>(() -> readTargets(file, TARGETS_PATH, TARGET_NAMES_PATH));
        rawTargets  = new Lazy<>(() -> readTargets(file, RAW_TARGETS_PATH, RAW_TARGET_NAMES_PATH));
        panelTargets = new Lazy<>(() -> readTargets(file, PANEL_TARGETS_PATH, PANEL_TARGET_NAMES_PATH));
        sampleNames = new Lazy<>(() -> readNames(file, SAMPLE_NAMES_PATH));
        panelSampleNames = new Lazy<>(() -> readNames(file, PANEL_SAMPLE_NAMES_PATH));
    }

    /**
     * Create a new PoN interface to a HDF5 file.  A version check is performed and a warning message logged if the
     * PoN version number is not up to date.
     *
     * @param file      the underlying HDF5 file.
     * @param logger    the logger to log the warning message with.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     */
    public HDF5PCACoveragePoN(final HDF5File file, final Logger logger) {
        this(file);
        if (getVersion() < CURRENT_PON_VERSION) {
            logger.warn("The version of the specified PoN (" + getVersion() + ") is older than the latest version " +
                    "(" + CURRENT_PON_VERSION + ").");
        }
    }

    @Override
    public double getVersion() {
        return file.readDouble(VERSION_PATH);
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

    @Override
    public double[] getTargetFactors() {
        final double[] values = file.readDoubleArray(TARGET_FACTORS_PATH);
        if (values.length != targetNames.get().size()) {
            throw new GATKException(String.format("Wrong number of elements in the target factors recovered " +
                    "from file '%s': number of target factors found in file (%d) != number of target names (%d)",
                    file.getFile(), values.length, targetNames.get().size()));
        }
        return values;
    }

    @Override
    public double[] getTargetVariances() {
        final double[] values = file.readDoubleArray(TARGET_VARIANCES_PATH);
        if (values.length != panelTargets.get().size()) {
            throw new GATKException(String.format("Wrong number of elements in the target variances recovered " +
                    "from file '%s': number of target variances found in file (%d) != number of panel-target names (%d)",
                    file.getFile(), values.length, panelTargets.get().size()));
        }
        return values;
    }

    @Override
    public RealMatrix getNormalizedCounts() {
        // Note the check uses sample names as number of rows and targets as number of columns.  This is due to the
        //  transposed storage.  The returned matrix is still targets (rows) x samples (columns).
        return readMatrixAndCheckDimensions(NORMALIZED_COUNTS_PATH, sampleNames.get().size(), targetNames.get().size())
                .transpose();
    }

    @Override
    public RealMatrix getLogNormalizedCounts() {
        // Note the check uses sample names as number of rows and targets as number of columns.  This is due to the
        //  transposed storage.  The returned matrix is still targets (rows) x samples (columns).
        return readMatrixAndCheckDimensions(LOG_NORMALIZED_COUNTS_PATH, getPanelSampleNames().size(), getPanelTargetNames().size())
                .transpose();
    }

    @Override
    public RealMatrix getLogNormalizedPInverseCounts() {
        return readMatrixAndCheckDimensions(LOG_NORMALIZED_PINV_PATH, getPanelSampleNames().size(),
                getPanelTargetNames().size());
    }

    @Override
    public RealMatrix getReducedPanelCounts() {
        // Note the check is using sample names as number of rows and targets as number of columns.  This is due to the
        //  transposed storage.  The returned matrix is still targets (rows) x pseudo-samples (columns).
        return readMatrixAndCheckDimensions(REDUCED_PANEL_COUNTS_PATH,
                c -> c <= getPanelSampleNames().size(),
                r -> r == panelTargetNames.get().size()).transpose();
    }

    @Override
    public RealMatrix getReducedPanelPInverseCounts() {
        return readMatrixAndCheckDimensions(REDUCED_PANEL_PINV_PATH,
                r -> r <= getPanelSampleNames().size(),
                c -> c == panelTargetNames.get().size());
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
        return panelTargetNames.get();
    }

    @Override
    public List<String> getSampleNames() {
        return sampleNames.get();
    }

    @Override
    public List<String> getPanelSampleNames() {
        return panelSampleNames.get();
    }

    /**
     * Write all of the coverage PoN fields to HDF5.
     */
    public static void write(final File outFile,
                             final HDF5File.OpenMode openMode,
                             final List<Target> rawTargets,
                             final ReadCountCollection normalizedCounts,
                             final ReadCountCollection logNormalizedCounts,
                             final double[] targetFactors,
                             final double[] targetVariances,
                             final ReductionResult reduction) {
        Utils.nonNull(outFile);
        Utils.nonNull(normalizedCounts);
        Utils.nonNull(logNormalizedCounts);
        Utils.nonNull(rawTargets);
        Utils.nonNull(targetFactors);
        Utils.nonNull(targetVariances);
        Utils.nonNull(reduction);
        try (final HDF5File file = new HDF5File(outFile, openMode)) {
            logger.info("Creating " + outFile.getAbsolutePath() + "...");
            final HDF5PCACoveragePoN pon = new HDF5PCACoveragePoN(file);

            logger.info("Setting version number (" + CURRENT_PON_VERSION + ")...");
            pon.setVersion(CURRENT_PON_VERSION);

            final List<Target> targets = normalizedCounts.targets();
            final List<Target> panelTargets = logNormalizedCounts.targets();
            logger.info("Setting targets ...");
            pon.setTargets(targets);
            logger.info("Setting raw targets ...");
            pon.setRawTargets(rawTargets);
            logger.info("Setting reduced panel targets ...");
            pon.setPanelTargets(panelTargets);

            logger.info("Setting target factors (" + targetFactors.length + ") ...");
            pon.setTargetFactors(targetFactors);
            logger.info("Setting target variances...");
            pon.setTargetVariances(targetVariances);

            logger.info("Setting normalized counts (" + normalizedCounts.counts().getRowDimension() +
                    " x " + normalizedCounts.counts().getColumnDimension() + ") (T)...");
            pon.setNormalizedCounts(normalizedCounts.counts());
            logger.info("Setting log-normalized counts (" + logNormalizedCounts.counts().getRowDimension() +
                    " x " + logNormalizedCounts.counts().getColumnDimension() + ") (T) ...");
            pon.setLogNormalizedCounts(logNormalizedCounts.counts());
            logger.info("Setting log-normalized pseudoinverse (" + reduction.getPseudoInverse().getRowDimension() +
                    " x " + reduction.getPseudoInverse().getColumnDimension() + ") ...");
            pon.setLogNormalizedPInverseCounts(reduction.getPseudoInverse());
            logger.info("Setting reduced panel counts (" + reduction.getReducedCounts().getRowDimension() +
                    " x " + reduction.getReducedCounts().getColumnDimension() + ") (T) ...");
            pon.setReducedPanelCounts(reduction.getReducedCounts());
            logger.info("Setting reduced panel pseudoinverse (" + reduction.getReducedPseudoInverse().getRowDimension() +
                    " x " + reduction.getReducedPseudoInverse().getColumnDimension() + ") ...");
            pon.setReducedPanelPInverseCounts(reduction.getReducedPseudoInverse());

            final List<String> targetNames = normalizedCounts.targets().stream().map(Target::getName).collect(Collectors.toList());
            final List<String> rawTargetNames = rawTargets.stream().map(Target::getName).collect(Collectors.toList());
            final List<String> panelTargetNames = logNormalizedCounts.targets().stream().map(Target::getName).collect(Collectors.toList());
            logger.info("Setting target names ...");
            pon.setTargetNames(targetNames);
            logger.info("Setting raw target names ...");
            pon.setRawTargetNames(rawTargetNames);
            logger.info("Setting reduced target names ...");
            pon.setPanelTargetNames(panelTargetNames);

            final List<String> sampleNames = normalizedCounts.columnNames();
            final List<String> panelSampleNames = logNormalizedCounts.columnNames();
            logger.info("Setting sample names ...");
            pon.setSampleNames(sampleNames);
            logger.info("Setting reduced sample names ...");
            pon.setPanelSampleNames(panelSampleNames);
        }
    }

    /**
     * Reads names from a path.
     * @param reader    the source HDF5 reader.
     * @param namesPath the path containing the names to be read.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<String> readNames(final HDF5File reader, final String namesPath) {
        final String[] values = reader.readStringArray(namesPath);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    /**
     * Reads targets from a path.
     * @param reader        the source HDF5 reader.
     * @param targetsPath   the path containing the targets to be read.
     * @param namesPath     the path containing the names of the targets to be read.
     * @return never {@code null}.
     * @throws GATKException if there was any problem reading the contents of the underlying HDF5 file.
     */
    private static List<Target> readTargets(final HDF5File reader, final String targetsPath, final String namesPath) {
        final String[][] values = reader.readStringMatrix(targetsPath, NUM_TARGET_COLUMNS_PATH);
        final List<String> targetNamesToRender = readNames(reader, namesPath);
        return renderPoNTargets(values, targetNamesToRender, reader);
    }

    private static List<Target> renderPoNTargets(final String[][] values, final List<String> targetNamesToRender, final HDF5File reader) {
        if (values.length != targetNamesToRender.size()) {
            throw new GATKException(String.format("Wrong number of elements in the targets recovered " +
                    "from file '%s': number of targets found in file (%d) != number of target names (%d)",
                    reader.getFile(), values.length, targetNamesToRender.size()));
        }

        final int numTargetCols = (int) reader.readDouble(NUM_TARGET_COLUMNS_PATH);
        final List<Target> result = new ArrayList<>(values.length);
        for (int i= 0; i < values.length; i++) {
            if (values[i].length != numTargetCols) {
                throw new GATKException(String.format("Wrong number of column elements in the targets recovered " +
                        "from file '%s': number of columns found in file (%d) != number of target columns (%d)",
                        reader.getFile(), values[i].length, numTargetCols));
            }
            result.add(new Target(targetNamesToRender.get(i), new SimpleInterval(values[i][0], Integer.parseInt(values[i][1]), Integer.parseInt(values[i][2]))));
        }
        return result;
    }

    //PRIVATE SETTERS (write values to HDF5 file and set internally held, lazily loaded fields to them)
    //These are private to prevent fields from being written individually, which could leave the file in a bad state

    private void setVersion(final double version) {
        file.makeDouble(VERSION_PATH, version);
    }

    private void setTargets(final List<Target> targets) {
        writeTargets(TARGETS_PATH, targets);
        this.targets  = new Lazy<>(() -> readTargets(file, TARGETS_PATH, TARGET_NAMES_PATH));
    }

    private void setRawTargets(final List<Target> targets) {
        writeTargets(RAW_TARGETS_PATH, targets);
        rawTargets = new Lazy<>(() -> readTargets(file, RAW_TARGETS_PATH, RAW_TARGET_NAMES_PATH));
    }

    private void setPanelTargets(final List<Target> targets) {
        writeTargets(PANEL_TARGETS_PATH, targets);
        panelTargets = new Lazy<>(() -> readTargets(file, PANEL_TARGETS_PATH, PANEL_TARGET_NAMES_PATH));
    }

    private void setTargetFactors(final double[] targetFactors) {
        writeTargetFactors(targetFactors);
    }

    private void setTargetVariances(final double[] targetVariances) {
        writeTargetVariances(targetVariances);
    }

    private void setNormalizedCounts(final RealMatrix counts) {
        file.makeDoubleMatrix(NORMALIZED_COUNTS_PATH, counts.transpose().getData());
    }

    private void setLogNormalizedCounts(final RealMatrix counts) {
        file.makeDoubleMatrix(LOG_NORMALIZED_COUNTS_PATH, counts.transpose().getData());
    }

    private void setLogNormalizedPInverseCounts(final RealMatrix counts) {
        file.makeDoubleMatrix(LOG_NORMALIZED_PINV_PATH, counts.getData());
    }

    private void setReducedPanelCounts(final RealMatrix counts) {
        file.makeDoubleMatrix(REDUCED_PANEL_COUNTS_PATH, counts.transpose().getData());
    }

    private void setReducedPanelPInverseCounts(final RealMatrix counts) {
        file.makeDoubleMatrix(REDUCED_PANEL_PINV_PATH, counts.getData());
    }

    private void setTargetNames(final List<String> names) {
        writeNames(TARGET_NAMES_PATH, names);
        targetNames = new Lazy<>(() -> readNames(file, TARGET_NAMES_PATH));
    }

    private void setRawTargetNames(final List<String> names) {
        writeNames(RAW_TARGET_NAMES_PATH, names);
        rawTargetNames = new Lazy<>(() -> readNames(file, RAW_TARGET_NAMES_PATH));
    }

    private void setPanelTargetNames(final List<String> names) {
        writeNames(PANEL_TARGET_NAMES_PATH, names);
        panelTargetNames = new Lazy<>(() -> readNames(file, PANEL_TARGET_NAMES_PATH));
    }

    private void setSampleNames(final List<String> names) {
        writeNames(SAMPLE_NAMES_PATH, names);
        sampleNames = new Lazy<>(() -> readNames(file, SAMPLE_NAMES_PATH));
    }

    private void setPanelSampleNames(final List<String> names) {
        writeNames(PANEL_SAMPLE_NAMES_PATH, names);
        panelSampleNames = new Lazy<>(() -> readNames(file, PANEL_SAMPLE_NAMES_PATH));
    }

    private void writeTargets(final String fullPath, List<Target> targets) {
        final String[][] targetValues = new String[targets.size()][NUM_TARGET_COLUMNS];
        for (int i = 0; i < targets.size(); i++) {
            final Target target = targets.get(i);
            targetValues[i][targetColumnToPoNIndex.get(TargetTableColumn.CONTIG)] = target.getContig();
            targetValues[i][targetColumnToPoNIndex.get(TargetTableColumn.START)] = String.valueOf(target.getStart());
            targetValues[i][targetColumnToPoNIndex.get(TargetTableColumn.END)] = String.valueOf(target.getEnd());
        }
        file.makeStringMatrix(fullPath, targetValues, NUM_TARGET_COLUMNS_PATH);
    }

    private void writeTargetFactors(final double[] targetFactors) {
        file.makeDoubleArray(TARGET_FACTORS_PATH, targetFactors);
    }

    private void writeTargetVariances(final double[] targetVariances) {
        file.makeDoubleArray(TARGET_VARIANCES_PATH, targetVariances);
    }

    private void writeNames(final String path, final List<String> names) {
        file.makeStringArray(path, names.toArray(new String[names.size()]));
    }

    //TODO: https://github.com/broadinstitute/gatk-protected/issues/637 move below methods to hdf5-java-bindings repo

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
}
