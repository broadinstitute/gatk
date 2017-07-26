package org.broadinstitute.hellbender.tools.pon.coverage.coveragehdf5;

import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * NOTE: There is code duplication in here, since some of the older codeis going to be removed in the future.
 */
public class HDF5Coverage {
    private static final String RAW_TARGETS_GROUP_NAME = "/raw_targets";
    private static final String SAMPLE_NAMES_GROUP_NAME = "/sample_names";

    /** Where the intervals live. */
    private static final String RAW_TARGETS_PATH = RAW_TARGETS_GROUP_NAME + "/block0_values";

    /** Where the target names live.*/
    private static final String RAW_TARGET_NAMES_PATH = RAW_TARGETS_GROUP_NAME + "/index";

    /** Where the sample names live.*/
    private static final String SAMPLE_NAMES_PATH = SAMPLE_NAMES_GROUP_NAME + "/index";

    /** Where the values associated with each target live.  This is usually coverage, but does not have to be. */
    private static final String TARGET_VALUES_GROUP_NAME = "/target_values";
    private static final String TARGET_VALUES_PATH = TARGET_VALUES_GROUP_NAME + "/values";

    private final HDF5File file;
    private Lazy<List<Target>> rawTargets;
    private Lazy<List<String>> rawTargetNames;
    private Lazy<List<String>> sampleNames;

    private static final Logger logger = LogManager.getLogger(HDF5Coverage.class);

    public HDF5Coverage(HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        rawTargetNames  = new Lazy<>(() -> readNames(file, RAW_TARGET_NAMES_PATH));
        sampleNames  = new Lazy<>(() -> readNames(file, SAMPLE_NAMES_PATH));
        rawTargets  = new Lazy<>(() -> readTargets(file, RAW_TARGETS_PATH, RAW_TARGET_NAMES_PATH));
    }

    public double[][] getTargetValues() {
        final double[][] values = file.readDoubleMatrix(TARGET_VALUES_PATH);
        if (values.length != rawTargetNames.get().size()) {
            throw new GATKException(String.format("Wrong number of elements in the target values recovered " +
                            "from file '%s': number of target values found in file (%d) != number of target names (%d)",
                    file.getFile(), values.length, rawTargetNames.get().size()));
        }
        return values;
    }

    public List<String> getRawTargetNames() {
        return rawTargetNames.get();
    }

    public List<Target> getRawTargets() {
        return rawTargets.get();
    }

    public List<String> getSampleNames() {
        return sampleNames.get();
    }

    private static List<String> readNames(final HDF5File reader, final String namesPath) {
        final String[] values = reader.readStringArray(namesPath);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    private static List<Target> readTargets(final HDF5File reader, final String targetsPath, final String namesPath) {
        final List<Locatable> targetIntervals = IntervalHelper.readIntervals(reader, targetsPath);
        final List<String> targetNamesToRender = readNames(reader, namesPath);

        return IntStream.range(0, targetIntervals.size()).mapToObj(i -> new Target(targetNamesToRender.get(i), targetIntervals.get(i))).collect(Collectors.toList());
    }

    private void writeTargets(final String fullPath, List<Target> targets) {
        IntervalHelper.writeIntervals(file, fullPath, targets);
    }

    private void writeTargetValues(final double[][] targetValues) {
        file.makeDoubleMatrix(TARGET_VALUES_PATH, targetValues);
    }

    private void writeNames(final String path, final List<String> names) {
        file.makeStringArray(path, names.toArray(new String[names.size()]));
    }

    /**
     *  Create (or modify) HDF5 file.
     *
     * @param outFile Path to the final HDF5 file.  Not {@code null}
     * @param openMode See {@link HDF5File.OpenMode}.  Must be {@link HDF5File.OpenMode} CREATE or READ_WRITE.  Not {@code null}
     * @param rawTargets the intervals and names for each target.  Not {@code null}
     * @param values a T x S matrix where T is the number of targets and S is the number of samples.  Not {@code null}
     * @param sampleNames the column names for each of the value columns.  Not {@code null}
     */
    public static void write(final File outFile,
                             final HDF5File.OpenMode openMode,
                             final List<Target> rawTargets, final double[][] values, final List<String> sampleNames) {
        if (!openMode.equals(HDF5File.OpenMode.CREATE) && !openMode.equals(HDF5File.OpenMode.READ_WRITE)) {
            throw new GATKException("Users should not try to fix this error.  Please contact the GATK forum.  The write method can only be used for HDF5 CREATE or READ_WRITE.");
        }
        Utils.nonNull(outFile);
        Utils.nonNull(openMode);
        Utils.nonNull(rawTargets);
        Utils.nonNull(values);
        Utils.nonNull(sampleNames);

        if (values.length != rawTargets.size()) {
            throw new GATKException("This is likely an issue to be solved by a GATK developer.  The shape of the values array (" + values.length + " x " + values[0].length + ") does not match the number of targets (" + rawTargets.size() + ").");
        }
        if (values[0].length != sampleNames.size()) {
            throw new GATKException("This is likely an issue to be solved by a GATK developer.  The shape of the values array (" + values.length + " x " + values[0].length + ") does not match the number of samples/columns (" + sampleNames.size() + ").");
        }

        try (final HDF5File file = new HDF5File(outFile, openMode)) {
            final HDF5Coverage hdf5CoverageFile = new HDF5Coverage(file);
            hdf5CoverageFile.writeTargets(RAW_TARGETS_PATH, rawTargets);
            hdf5CoverageFile.writeTargetValues(values);
            hdf5CoverageFile.writeNames(RAW_TARGET_NAMES_PATH, rawTargets.stream().map(Target::getName).collect(Collectors.toList()));
            hdf5CoverageFile.writeNames(SAMPLE_NAMES_PATH, sampleNames);
        }
    }

    /**
     *  See {@link HDF5Coverage::write}, but only creates new files.  Will fail if file exists.
     *
     *  This is the recommended method for creating HDF5 files.
     *
     * @param outFile
     * @param rawTargets
     * @param values
     * @param sampleNames
     */
    public static void create(final File outFile, final List<Target> rawTargets, final double[][] values, final List<String> sampleNames) {
        if (outFile.exists()) {
            throw new UserException.BadInput(outFile.getAbsolutePath() + " already exists.  This only allows creation of new files.");
        }
        write(outFile, HDF5File.OpenMode.CREATE, rawTargets, values, sampleNames);
    }

    private static final class IntervalHelper {
        //writing intervals as a string matrix is expensive,
        //so we instead store a map from integer indices to contig strings and
        //store (index, start, end) in a double matrix

        private static final String INTERVAL_CONTIG_NAMES_SUB_PATH = "/indexed_contig_names";
        private static final String INTERVAL_MATRIX_SUB_PATH = "/transposed_index_start_end";

        private enum IntervalField {
            CONTIG_INDEX(0),
            START (1),
            END (2);

            private final int index;

            IntervalField(final int index) {
                this.index = index;
            }
        }
        private static final int NUM_INTERVAL_FIELDS = IntervalField.values().length;

        private static List<Locatable> readIntervals(final HDF5File file,
                                                     final String path) {
            final String[] contigNames = file.readStringArray(path + INTERVAL_CONTIG_NAMES_SUB_PATH);
            final double[][] matrix = file.readDoubleMatrix(path + INTERVAL_MATRIX_SUB_PATH);
            final int numIntervals = matrix[0].length;
            return IntStream.range(0, numIntervals).boxed()
                    .map(i -> (new SimpleInterval(
                            contigNames[(int) matrix[IntervalField.CONTIG_INDEX.index][i]],
                            (int) matrix[IntervalField.START.index][i],
                            (int) matrix[IntervalField.END.index][i])))
                    .collect(Collectors.toList());
        }

        private static <T extends Locatable> void writeIntervals(final HDF5File file,
                                           final String path,
                                           final List<T> intervals) {
            final String[] contigNames = intervals.stream().map(Locatable::getContig).distinct().toArray(String[]::new);
            file.makeStringArray(path + INTERVAL_CONTIG_NAMES_SUB_PATH, contigNames);
            final Map<String, Double> contigNamesToIndexMap = IntStream.range(0, contigNames.length).boxed()
                    .collect(Collectors.toMap(i -> contigNames[i], i -> (double) i));
            final double[][] matrix = new double[NUM_INTERVAL_FIELDS][intervals.size()];
            for (int i = 0; i < intervals.size(); i++) {
                final Locatable interval = intervals.get(i);
                matrix[IntervalField.CONTIG_INDEX.index][i] = contigNamesToIndexMap.get(interval.getContig());
                matrix[IntervalField.START.index][i] = interval.getStart();
                matrix[IntervalField.END.index][i] = interval.getEnd();
            }
            file.makeDoubleMatrix(path + INTERVAL_MATRIX_SUB_PATH, matrix);
        }
    }
}
