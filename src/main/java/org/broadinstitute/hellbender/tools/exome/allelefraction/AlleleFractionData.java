package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.SegmentedModel;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Data for the Allele Fraction Model containing the set of het alt and ref counts
 * and the grouping of hets into segments.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionData implements DataCollection {
    private final List<AllelicCount> allelicCounts; // allelic counts indexed by het
    private final List<Integer> hetIndices; // the numbers 0, 1, 2 . . . N_hets
    private final List<Integer> startHetsPerSegment = new ArrayList<>();
    private final List<Integer> numHetsPerSegment = new ArrayList<>();

    public AlleleFractionData(final SegmentedModel segmentedModel) {
        allelicCounts = new ArrayList<>();
        final List<SimpleInterval> segmentIntervals = segmentedModel.getSegments();
        final TargetCollection<AllelicCount> alleleCounts = segmentedModel.getGenome().getSNPs();

        int startHet = 0;
        for (final SimpleInterval segment : segmentIntervals) {
            startHetsPerSegment.add(startHet);
            final List<AllelicCount> countsInSegment = alleleCounts.targets(segment);
            numHetsPerSegment.add(countsInSegment.size());
            startHet += countsInSegment.size();
            allelicCounts.addAll(countsInSegment);
        }

        hetIndices = IntStream.range(0, allelicCounts.size()).boxed().collect(Collectors.toList());
    }

    public List<AllelicCount> getAllelicCounts() { return allelicCounts; }

    public List<AllelicCount> countsInSegment(final int segment) {
        final int startInclusive = startHetsPerSegment.get(segment);
        final int endExclusive = startInclusive + numHetsPerSegment.get(segment);
        return Collections.unmodifiableList(allelicCounts.subList(startInclusive, endExclusive));
    }

    public List<Integer> hetsInSegment(final int segment) {
        final int startInclusive = startHetsPerSegment.get(segment);
        final int endExclusive = startInclusive + numHetsPerSegment.get(segment);
        return Collections.unmodifiableList(hetIndices.subList(startInclusive, endExclusive));
    }

    public int numHetsInSegment(final int segment) {
        return numHetsPerSegment.get(segment);
    }

    public int numSegments() { return startHetsPerSegment.size(); }

    public AllelicCount count(final int het) { return allelicCounts.get(het); }

    public int altCount(final int het) { return allelicCounts.get(het).getAltReadCount(); }

    public int refCount(final int het) { return allelicCounts.get(het).getRefReadCount(); }

    public int readCount(final int het) { return altCount(het) + refCount(het); }
}
