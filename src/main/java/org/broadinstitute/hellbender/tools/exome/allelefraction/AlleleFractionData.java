package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the allele-fraction model containing the set of het alt and ref counts
 * and the grouping of hets into segments.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionData implements DataCollection {
    private final int numSegments;
    private final AllelicPanelOfNormals allelicPoN;
    private final List<AllelicCount> allelicCounts; // allelic counts indexed by het
    private final List<Integer> hetIndices; // the numbers 0, 1, 2 . . . N_hets
    private final List<Integer> startHetsPerSegment = new ArrayList<>();
    private final List<Integer> numHetsPerSegment = new ArrayList<>();

    public AlleleFractionData(final SegmentedGenome segmentedGenome) {
        this(segmentedGenome, AllelicPanelOfNormals.EMPTY_PON);
    }

    public AlleleFractionData(final SegmentedGenome segmentedGenome, final AllelicPanelOfNormals allelicPoN) {
        numSegments = segmentedGenome.getSegments().size();
        this.allelicPoN = allelicPoN;
        allelicCounts = new ArrayList<>();
        final List<SimpleInterval> segmentIntervals = segmentedGenome.getSegments();
        final TargetCollection<AllelicCount> alleleCounts = segmentedGenome.getGenome().getSNPs();

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

    public AllelicPanelOfNormals getPoN() { return allelicPoN; }

    public List<AllelicCount> getAllelicCounts() { return Collections.unmodifiableList(allelicCounts); }

    public List<AllelicCount> getCountsInSegment(final int segment) {
        final int startInclusive = startHetsPerSegment.get(segment);
        final int endExclusive = startInclusive + numHetsPerSegment.get(segment);
        return Collections.unmodifiableList(allelicCounts.subList(startInclusive, endExclusive));
    }

    public List<Integer> getHetsInSegment(final int segment) {
        final int startInclusive = startHetsPerSegment.get(segment);
        final int endExclusive = startInclusive + numHetsPerSegment.get(segment);
        return Collections.unmodifiableList(hetIndices.subList(startInclusive, endExclusive));
    }

    public int getNumHetsInSegment(final int segment) {
        return numHetsPerSegment.get(segment);
    }

    public int getNumSegments() { return numSegments; }

    public AllelicCount getAllelicCount(final int het) { return allelicCounts.get(het); }

    public int getAltCount(final int het) { return allelicCounts.get(het).getAltReadCount(); }

    public int getRefCount(final int het) { return allelicCounts.get(het).getRefReadCount(); }

    public int getReadCount(final int het) { return getAltCount(het) + getRefCount(het); }
}
