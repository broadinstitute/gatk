package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Holds histograms of alt depth=1 sites for reference contexts.
 */
public class DepthOneHistograms {
    private final Map<String, Map<Pair<Nucleotide, ReadOrientation>, Histogram<Integer>>> map;
    private final int maxDepth;

    public DepthOneHistograms(final int maxDepth) {
        this.maxDepth = maxDepth;

        map = new HashMap<>(F1R2FilterConstants.NUM_KMERS);

        // Initialize, for each reference context, the (Alt Allele, Artifact Type) -> Histogram map
        F1R2FilterConstants.ALL_KMERS.forEach(context -> {
            map.put(context, new HashMap<>((Nucleotide.STANDARD_BASES.size() - 1) * ReadOrientation.SIZE));

            for (Nucleotide altAllele : Nucleotide.STANDARD_BASES) {
                // Skip e.g. AGT -> AGT because G is not an alt allele

                if (altAllele == F1R2FilterUtils.getMiddleBase(context)) {
                    continue;
                }

                for (ReadOrientation artifactType : ReadOrientation.values()) {
                    map.get(context).put(new ImmutablePair<>(altAllele, artifactType),
                            F1R2FilterUtils.createAltHistogram(context, altAllele, artifactType, maxDepth));
                }
            }
        });
    }

    public Histogram<Integer> get(final String referenceContext, final Nucleotide altAllele, final ReadOrientation orientation) {
        final Pair<Nucleotide, ReadOrientation> key = new ImmutablePair<>(altAllele, orientation);
        return map.get(referenceContext).get(key);
    }

    public void increment(final String referenceContext, final Nucleotide altAllele,
                          final ReadOrientation orientation, final int depth) {
        final int cappedDepth = Math.min(depth, maxDepth);
        map.get(referenceContext).get(new ImmutablePair<>(altAllele, orientation)).increment(cappedDepth);
    }

    public List<Histogram<Integer>> getHistograms() {
        return map.values().stream().flatMap(s -> s.values().stream()).collect(Collectors.toList());
    }
}
