package org.broadinstitute.hellbender.tools.copynumber.models;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SampleLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the copy-ratio model containing the copy-ratio data grouped by segment.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioSegmentedData implements DataCollection {
    private final CopyRatioCollection copyRatios;
    private final List<SimpleInterval> segments;
    private final double minLog2CopyRatioValue;
    private final double maxLog2CopyRatioValue;

    private final List<IndexedCopyRatio> indexedCopyRatios;
    private final List<IndexRange> indexRangesPerSegment;

    CopyRatioSegmentedData(final CopyRatioCollection copyRatios,
                           final List<SimpleInterval> segments) {
        this.copyRatios = Utils.nonNull(copyRatios);
        this.segments = Utils.nonEmpty(segments).stream().sorted(SampleLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());

        final List<Double> log2CopyRatioValues = copyRatios.getLog2CopyRatioValues();
        minLog2CopyRatioValue = log2CopyRatioValues.stream().min(Double::compareTo).orElse(Double.NaN);
        maxLog2CopyRatioValue = log2CopyRatioValues.stream().max(Double::compareTo).orElse(Double.NaN);

        indexedCopyRatios = new ArrayList<>(copyRatios.size());
        indexRangesPerSegment = new ArrayList<>(segments.size());

        //construct list of lists of copy ratios with an index in order corresponding to that of segments;
        //segment assignment is based on midpoint of copy-ratio interval
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = copyRatios.getMidpointOverlapDetector();
        int index = 0;
        for (int segmentIndex = 0; segmentIndex < segments.size(); segmentIndex++) {
            final SimpleInterval segment = segments.get(segmentIndex);
            final List<CopyRatio> copyRatiosInSegment = copyRatioMidpointOverlapDetector.getOverlaps(segment).stream()
                    .sorted(SampleLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                    .collect(Collectors.toList());
            final int segmentStartIndex = index;
            final int si = segmentIndex;
            IntStream.range(0, copyRatiosInSegment.size()).boxed()
                    .map(i -> new IndexedCopyRatio(copyRatiosInSegment.get(i), segmentStartIndex + i, si))
                    .forEach(indexedCopyRatios::add);
            indexRangesPerSegment.add(new IndexRange(segmentStartIndex, segmentStartIndex + copyRatiosInSegment.size()));
            index += copyRatiosInSegment.size();
        }
    }

    CopyRatioCollection getCopyRatios() {
        return copyRatios;
    }

    List<SimpleInterval> getSegments() {
        return Collections.unmodifiableList(segments);
    }

    int getNumSegments() {
        return segments.size();
    }

    int getNumPoints() {
        return copyRatios.size();
    }

    String getSampleName() {
        return copyRatios.getSampleName();
    }

    double getMinLog2CopyRatioValue() {
        return minLog2CopyRatioValue;
    }

    double getMaxLog2CopyRatioValue() {
        return maxLog2CopyRatioValue;
    }

    List<IndexedCopyRatio> getIndexedCopyRatios() {
        return Collections.unmodifiableList(indexedCopyRatios);
    }

    List<IndexedCopyRatio> getIndexedCopyRatiosInSegment(final int segmentIndex) {
        return Collections.unmodifiableList(indexedCopyRatios.subList(
                indexRangesPerSegment.get(segmentIndex).from, indexRangesPerSegment.get(segmentIndex).to));
    }

    //estimate global variance empirically by taking average of all per-segment variances
    double estimateVariance() {
        return IntStream.range(0, segments.size())
                .mapToDouble(s -> new Variance().evaluate(Doubles.toArray(
                        getIndexedCopyRatiosInSegment(s).stream()
                                .map(IndexedCopyRatio::getLog2CopyRatioValue)
                                .collect(Collectors.toList()))))
                .filter(v -> !Double.isNaN(v))
                .average().orElse(Double.NaN);
    }

    //estimate segment means empirically by taking averages of log2 copy ratios in each segment
    CopyRatioState.SegmentMeans estimateSegmentMeans() {
        final List<Double> means = IntStream.range(0, segments.size()).boxed()
                .map(s -> new Mean().evaluate(Doubles.toArray(
                        getIndexedCopyRatiosInSegment(s).stream()
                                .map(IndexedCopyRatio::getLog2CopyRatioValue)
                                .collect(Collectors.toList()))))
                .collect(Collectors.toList());
        return new CopyRatioState.SegmentMeans(means);
    }

    static final class IndexedCopyRatio extends CopyRatio {
        private final int index;
        private final int segmentIndex;

        private IndexedCopyRatio(final CopyRatio copyRatio,
                                 final int index,
                                 final int segmentIndex) {
            super(copyRatio.getInterval(), copyRatio.getLog2CopyRatioValue());
            this.index = index;
            this.segmentIndex = segmentIndex;
        }

        int getIndex() {
            return index;
        }

        int getSegmentIndex() {
            return segmentIndex;
        }
    }
}
