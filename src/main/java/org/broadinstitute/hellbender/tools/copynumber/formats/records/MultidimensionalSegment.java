package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public class MultidimensionalSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPointsCopyRatio;
    private final int numPointsAlleleFraction;
    private final double meanLog2CopyRatio;

    public MultidimensionalSegment(final SimpleInterval interval,
                                   final int numPointsCopyRatio,
                                   final int numPointsAlleleFraction,
                                   final double meanLog2CopyRatio) {
        Utils.nonNull(interval);
        Utils.validateArg(numPointsCopyRatio > 0 || numPointsAlleleFraction > 0,
                String.format("Number of copy-ratio points or number of allele-fraction points must be positive: %s", interval));
        this.interval = interval;
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        this.meanLog2CopyRatio = meanLog2CopyRatio;
    }

    public MultidimensionalSegment(final SimpleInterval interval,
                                   final List<CopyRatio> denoisedLog2CopyRatios,
                                   final List<AllelicCount> allelicCounts) {
        Utils.nonNull(interval);
        Utils.nonNull(denoisedLog2CopyRatios);
        Utils.nonNull(allelicCounts);
        this.interval = interval;
        numPointsCopyRatio = denoisedLog2CopyRatios.size();
        numPointsAlleleFraction = allelicCounts.size();
        meanLog2CopyRatio = new CopyRatioSegment(interval, denoisedLog2CopyRatios).getMeanLog2CopyRatio();
    }

    public MultidimensionalSegment(final SimpleInterval interval,
                                   final Comparator<Locatable> comparator,
                                   final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector,
                                   final OverlapDetector<AllelicCount> allelicCountOverlapDetector) {
        this(
                interval,
                copyRatioMidpointOverlapDetector.getOverlaps(interval).stream()
                        .sorted(comparator)
                        .collect(Collectors.toList()),
                allelicCountOverlapDetector.getOverlaps(interval).stream()
                        .sorted(comparator)
                        .collect(Collectors.toList()));
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public int getNumPointsCopyRatio() {
        return numPointsCopyRatio;
    }

    public int getNumPointsAlleleFraction() {
        return numPointsAlleleFraction;
    }

    public double getMeanLog2CopyRatio() {
        return meanLog2CopyRatio;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final MultidimensionalSegment that = (MultidimensionalSegment) o;

        return numPointsCopyRatio == that.numPointsCopyRatio &&
                numPointsAlleleFraction == that.numPointsAlleleFraction &&
                Double.compare(that.meanLog2CopyRatio, meanLog2CopyRatio) == 0 &&
                interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPointsCopyRatio;
        result = 31 * result + numPointsAlleleFraction;
        temp = Double.doubleToLongBits(meanLog2CopyRatio);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "MultidimensionalSegment{" +
                "interval=" + interval +
                ", numPointsCopyRatio=" + numPointsCopyRatio +
                ", numPointsAlleleFraction=" + numPointsAlleleFraction +
                ", meanLog2CopyRatio=" + meanLog2CopyRatio +
                '}';
    }
}
