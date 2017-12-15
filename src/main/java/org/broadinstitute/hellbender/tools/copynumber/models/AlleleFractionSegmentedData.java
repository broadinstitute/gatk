package org.broadinstitute.hellbender.tools.copynumber.models;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SampleLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
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
 * {@link DataCollection} for the allele-fraction model containing the het alt and ref counts grouped by segment.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionSegmentedData implements DataCollection {
    private final AllelicCountCollection allelicCounts;
    private final List<SimpleInterval> segments;

    private final List<IndexedAllelicCount> indexedAllelicCounts;
    private final List<IndexRange> indexRangesPerSegment;

    AlleleFractionSegmentedData(final AllelicCountCollection allelicCounts,
                                final List<SimpleInterval> segments) {
        this.allelicCounts = Utils.nonNull(allelicCounts);
        this.segments = Utils.nonEmpty(segments).stream().sorted(SampleLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());

        indexedAllelicCounts = new ArrayList<>(allelicCounts.size());
        indexRangesPerSegment = new ArrayList<>(segments.size());

        final OverlapDetector<AllelicCount> allelicCountOverlapDetector = allelicCounts.getOverlapDetector();
        int startIndex = 0;
        for (int segmentIndex = 0; segmentIndex < segments.size(); segmentIndex++) {
            final SimpleInterval segment = segments.get(segmentIndex);
            final List<AllelicCount> allelicCountsInSegment = allelicCountOverlapDetector.getOverlaps(segment).stream()
                    .sorted(SampleLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                    .collect(Collectors.toList());
            final int segmentStartIndex = startIndex;
            final int si = segmentIndex;
            IntStream.range(0, allelicCountsInSegment.size()).boxed()
                    .map(i -> new IndexedAllelicCount(allelicCountsInSegment.get(i), segmentStartIndex + i, si))
                    .forEach(indexedAllelicCounts::add);
            indexRangesPerSegment.add(new IndexRange(segmentStartIndex, segmentStartIndex + allelicCountsInSegment.size()));
            startIndex += allelicCountsInSegment.size();
        }
    }

    AllelicCountCollection getAllelicCounts() {
        return allelicCounts;
    }

    List<SimpleInterval> getSegments() {
        return Collections.unmodifiableList(segments);
    }

    int getNumSegments() {
        return segments.size();
    }

    int getNumPoints() {
        return allelicCounts.size();
    }

    List<IndexedAllelicCount> getIndexedAllelicCounts() {
        return Collections.unmodifiableList(indexedAllelicCounts);
    }

    List<IndexedAllelicCount> getIndexedAllelicCountsInSegment(final int segmentIndex) {
        return Collections.unmodifiableList(indexedAllelicCounts.subList(
                indexRangesPerSegment.get(segmentIndex).from, indexRangesPerSegment.get(segmentIndex).to));
    }

    static final class IndexedAllelicCount extends AllelicCount {
        private final int index;
        private final int segmentIndex;

        private IndexedAllelicCount(final AllelicCount allelicCount,
                                    final int index,
                                    final int segmentIndex) {
            super(allelicCount.getInterval(), allelicCount.getRefReadCount(), allelicCount.getAltReadCount(), allelicCount.getRefNucleotide(), allelicCount.getAltNucleotide());
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
