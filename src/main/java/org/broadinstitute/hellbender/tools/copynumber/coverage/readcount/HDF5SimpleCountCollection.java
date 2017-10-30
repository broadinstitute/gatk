package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

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

final class HDF5SimpleCountCollection implements SampleMetadata {
    private static final String SAMPLE_NAME_PATH = "/sample_metadata/sample_name";
    private static final String INTERVALS_GROUP_NAME = "/intervals";
    private static final String COUNTS_PATH = "/counts/values";

    private final HDF5File file;
    private Lazy<List<SimpleInterval>> intervals;

    HDF5SimpleCountCollection(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        intervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, INTERVALS_GROUP_NAME));
    }

    public SampleMetadata getSampleMetadata() {
        return new SimpleSampleMetadata(getSampleName());
    }

    @Override
    public String getSampleName() {
        return file.readStringArray(SAMPLE_NAME_PATH)[0];
    }

    List<SimpleInterval> getIntervals() {
        return intervals.get();
    }

    /**
     * @return single-row matrix containing the counts
     */
    RealMatrix getCounts() {
        return new Array2DRowRealMatrix(file.readDoubleMatrix(COUNTS_PATH));
    }

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
