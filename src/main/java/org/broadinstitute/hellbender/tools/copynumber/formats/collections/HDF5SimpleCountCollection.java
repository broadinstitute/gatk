package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Helper class for {@link SimpleCountCollection} used to read/write HDF5.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class HDF5SimpleCountCollection implements SampleMetadata {
    private static final String SAMPLE_NAME_PATH = "/sample_metadata/sample_name";
    private static final String INTERVALS_GROUP_NAME = "/intervals";
    private static final String COUNTS_PATH = "/counts/values";

    private final HDF5File file;
    private final Lazy<String> sampleName;
    private final Lazy<List<SimpleInterval>> intervals;
    private final Lazy<RealMatrix> counts;

    HDF5SimpleCountCollection(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        sampleName = new Lazy<>(() -> file.readStringArray(SAMPLE_NAME_PATH)[0]);
        intervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, INTERVALS_GROUP_NAME));
        counts = new Lazy<>(() -> new Array2DRowRealMatrix(file.readDoubleMatrix(COUNTS_PATH)));
    }

    public SampleMetadata getSampleMetadata() {
        return new SimpleSampleMetadata(getSampleName());
    }

    @Override
    public String getSampleName() {
        return sampleName.get();
    }

    List<SimpleInterval> getIntervals() {
        return intervals.get();
    }

    /**
     * @return single-row matrix containing the counts
     */
    RealMatrix getCounts() {
        return counts.get();
    }

    /**
     * @param intervals note that no particular sort order is assumed or checked for here,
     *                  but this package-protected method should only be called by {@link SimpleCountCollection#writeHDF5},
     *                  which enforces the order specified by {@link SampleLocatableCollection}
     */
    static void write(final File outFile,
                      final String sampleName,
                      final List<SimpleInterval> intervals,
                      final double[] counts) {
        Utils.nonNull(outFile);
        Utils.nonNull(sampleName);
        Utils.nonEmpty(intervals);
        Utils.nonNull(counts);

        Utils.validateArg(intervals.size() == counts.length, "Number of intervals and counts must match.");
        Utils.validateArg(intervals.stream().distinct().count() == intervals.size(), "Intervals must all be unique.");
        Utils.validateArg(Arrays.stream(counts).noneMatch(c -> c < 0), "Counts must all be non-negative integers.");

        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            final HDF5SimpleCountCollection hdf5CountCollection = new HDF5SimpleCountCollection(file);
            hdf5CountCollection.writeName(SAMPLE_NAME_PATH, sampleName);
            hdf5CountCollection.writeIntervals(intervals);
            hdf5CountCollection.writeCounts(counts);
        }
    }

    private <T extends SimpleInterval> void writeIntervals(final List<T> intervals) {
        HDF5Utils.writeIntervals(file, INTERVALS_GROUP_NAME, intervals);
    }

    private void writeName(final String path, final String name) {
        file.makeStringArray(path, name);
    }

    private void writeCounts(final double[] counts) {
        file.makeDoubleMatrix(COUNTS_PATH, new double[][]{counts});
    }
}
