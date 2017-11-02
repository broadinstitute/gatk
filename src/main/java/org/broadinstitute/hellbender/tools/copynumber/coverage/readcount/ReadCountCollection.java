package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;

/**
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public abstract class ReadCountCollection<RECORD extends ReadCountData> implements SampleMetadata {
    private final List<SimpleInterval> intervalList;
    private final Lazy<RealMatrix> totalReadCountMatrix;
    private final Lazy<Map<SimpleInterval, RECORD>> readCountDataMap;
    private final String sampleName;

    public ReadCountCollection(final File file, final List<SimpleInterval> intervalList) {
        Utils.nonNull(file);

        this.intervalList = intervalList;
        this.totalReadCountMatrix = new Lazy(initializeRealMatrix(file));
        this.readCountDataMap = new Lazy(initializeReadCountDataMap(file));
        this.sampleName = readSampleName(file);
    }

    /**
     * Read in the total read count matrix
     */
    abstract Supplier<RealMatrix> initializeRealMatrix(final File file);

    /**
     * Read in the read count data map
     */
    abstract Supplier<Map<SimpleInterval, RECORD>> initializeReadCountDataMap(final File file);

    /**
     * Read in sample name from read count file
     */
    abstract String readSampleName(final File file);

    /**
     * Get a record for a particular interval
     *
     * @param interval an interval
     * @return a record containing read count information
     */
    public RECORD getCountsForInterval(final SimpleInterval interval) {
        return readCountDataMap.get().get(interval);
    }

    /**
     * Get all intervals in this read count collection
     *
     * @return list of intervals
     */
    public List<SimpleInterval> getIntervals() {
        return intervalList;
    }

    /**
     * @return an aggregate matrix of read counts
     */
    public final RealMatrix getTotalReadCountMatrix() {
        return totalReadCountMatrix.get();
    }

    @Override
    public String getSampleName() {
        return sampleName;
    }
}
