package org.broadinstitute.hellbender.tools.exome.copyratio;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.broadinstitute.hellbender.tools.exome.ReadCountRecord;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the copy-ratio model containing the set of coverages and the grouping of coverages into segments.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioData implements DataCollection {
    private final int numSegments;
    private final int numTargets;
    private final double coverageMin;
    private final double coverageMax;

    private final List<List<IndexedCoverage>> indexedCoveragesPerSegment = new ArrayList<>();

    public CopyRatioData(final SegmentedGenome segmentedGenome) {
        final TargetCollection<ReadCountRecord.SingleSampleRecord> targetCoverages = segmentedGenome.getGenome().getTargets();
        Utils.validateArg(targetCoverages.targetCount() > 0, "Cannot construct CopyRatioData with no target-coverage data.");
        //construct list of coverages (in order corresponding to that of segments in SegmentedGenome;
        //this may not be in genomic order, depending on how the segments are sorted in the segment file,
        //so we cannot simply take the list of coverages in the order from TargetCollection.targets()
        final List<SimpleInterval> segments = segmentedGenome.getSegments();
        numSegments = segments.size();
        final List<Double> coverages = segments.stream()
                .flatMap(s -> targetCoverages.targets(s).stream())
                .map(ReadCountRecord.SingleSampleRecord::getCount)
                .collect(Collectors.toList());
        numTargets = coverages.size();
        coverageMin = Collections.min(coverages);
        coverageMax = Collections.max(coverages);
        //partition coverages with target indices by segment
        int targetIndex = 0;
        for (final SimpleInterval segment : segments) {
            final List<ReadCountRecord.SingleSampleRecord> targetCoveragesInSegment = targetCoverages.targets(segment);
            final List<IndexedCoverage> indexedCoveragesInSegment = new ArrayList<>();
            for (final ReadCountRecord.SingleSampleRecord targetCoverage : targetCoveragesInSegment) {
                final double coverage = targetCoverage.getCount();
                indexedCoveragesInSegment.add(new IndexedCoverage(coverage, targetIndex++));
            }
            indexedCoveragesPerSegment.add(indexedCoveragesInSegment);
        }
    }

    public int getNumSegments() {
        return numSegments;
    }

    public int getNumTargets() {
        return numTargets;
    }

    public double getCoverageMin() {
        return coverageMin;
    }

    public double getCoverageMax() {
        return coverageMax;
    }

    public List<IndexedCoverage> getIndexedCoveragesInSegment(final int segment) {
        return Collections.unmodifiableList(indexedCoveragesPerSegment.get(segment));
    }

    //estimate global variance empirically by taking average of all per-segment variances
    public double estimateVariance() {
        return IntStream.range(0, numSegments)
                .mapToDouble(s -> new Variance().evaluate(Doubles.toArray(getIndexedCoveragesInSegment(s).stream().map(IndexedCoverage::getCoverage).collect(Collectors.toList()))))
                .average().getAsDouble();
    }

    //estimate segment means empirically by taking averages of coverages in each segment
    public CopyRatioState.SegmentMeans estimateSegmentMeans() {
        final List<Double> means = IntStream.range(0, numSegments).boxed()
                .map(s -> new Mean().evaluate(Doubles.toArray(getIndexedCoveragesInSegment(s).stream().map(IndexedCoverage::getCoverage).collect(Collectors.toList()))))
                .collect(Collectors.toList());
        return new CopyRatioState.SegmentMeans(means);
    }

    public final class IndexedCoverage {
        private final double coverage;
        private final int targetIndex;

        public IndexedCoverage(final double coverage, final int targetIndex) {
            this.coverage = coverage;
            this.targetIndex = targetIndex;
        }

        public double getCoverage() {
            return coverage;
        }

        public int getTargetIndex() {
            return targetIndex;
        }
    }
}
