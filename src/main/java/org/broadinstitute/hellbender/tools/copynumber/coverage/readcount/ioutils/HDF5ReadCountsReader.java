package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ioutils;

import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class HDF5ReadCountsReader {
    private static final String SAMPLE_NAME_PATH = "/sample_name/value";
    private static final String INTERVALS_GROUP_NAME = "/intervals";
    private static final String READ_COUNTS_PATH = "/read_counts/values";

    private final HDF5File file;
    private Lazy<List<SimpleInterval>> intervals;

    public HDF5ReadCountsReader(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        intervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, INTERVALS_GROUP_NAME));
    }

    public String getSampleName() {
        return file.readStringArray(SAMPLE_NAME_PATH)[0];
    }

    public List<SimpleInterval> getIntervals() {
        return intervals.get();
    }

    /**
     * @return single-row matrix containing the read counts
     */
    public RealMatrix getReadCounts() {
        return new Array2DRowRealMatrix(file.readDoubleMatrix(READ_COUNTS_PATH));
    }

    public static void write(final File outFile,
                             final String sampleName,
                             final List<SimpleInterval> intervals,
                             final double[][] readCounts) {
        Utils.nonNull(outFile);
        Utils.nonNull(sampleName);
        Utils.nonEmpty(intervals);
        Utils.nonNull(readCounts);

        Utils.validateArg(readCounts.length == 1, "Read-count matrix must contain only a single row.");
        Utils.validateArg(intervals.size() == readCounts[0].length, "Number of intervals and read counts must match.");
        Utils.validateArg(Arrays.stream(readCounts[0]).noneMatch(c -> c < 0), "Read counts must all be non-negative integers.");
        Utils.validateArg(intervals.stream().distinct().count() == intervals.size(), "Intervals must all be unique.");

        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            final HDF5ReadCountsReader hdf5ReadCountCollection = new HDF5ReadCountsReader(file);
            hdf5ReadCountCollection.writeName(SAMPLE_NAME_PATH, sampleName);
            hdf5ReadCountCollection.writeIntervals(intervals);
            hdf5ReadCountCollection.writeReadCounts(readCounts);
        }
    }

    private void writeIntervals(final List<SimpleInterval> intervals) {
        HDF5Utils.writeIntervals(file, INTERVALS_GROUP_NAME, intervals);
    }

    private void writeName(final String path, final String name) {
        file.makeStringArray(path, name);
    }

    private void writeReadCounts(final double[][] readCounts) {
        file.makeDoubleMatrix(READ_COUNTS_PATH, readCounts);
    }
}
