package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by tsato on 3/14/18.
 */
public class ReadOrientationFilterConstants {
    public static final int REF_CONTEXT_PADDING_ON_EACH_SIDE = 1;

    static final int MIDDLE_INDEX = REF_CONTEXT_PADDING_ON_EACH_SIDE;

    static final int REFERENCE_CONTEXT_SIZE = 2 * REF_CONTEXT_PADDING_ON_EACH_SIDE + 1; // aka 3

    public static final List<Nucleotide> REGULAR_BASES = Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T);

    // the list of all possible kmers, where k = REFERENCE_CONTEXT_SIZE
    static final List<String> ALL_KMERS = SequenceUtil.generateAllKmers(REFERENCE_CONTEXT_SIZE).stream()
            .map(String::new).collect(Collectors.toList());

    // Each of these K-mers represents itself and its reverse complement
    static final List<String> CANONICAL_KMERS = ALL_KMERS.stream()
            .map(context -> new TreeSet<>(Arrays.asList(context, SequenceUtil.reverseComplement(context))))
            .distinct()
            .map(s -> s.first().compareTo(s.last()) < 0 ? s.first() : s.last())
            .collect(Collectors.toList());

    // If the posterior probability of neither F1R2 nor F2R1 is above this value, do not annotate the format fields with
    // information required for the filter.
    public static final double POSTERIOR_EMISSION_THRESHOLD = 0.3;

    // Delimiter of the alt histogram label
    static final String FIELD_SEPARATOR = "_";

    static final String binName = "depth";

    // We combine all sites of depths above this value in the last bin of the histogram
    static final int maxDepthForHistograms = 200;

    static final Integer[] bins = IntStream.rangeClosed(1, ReadOrientationFilterConstants.maxDepthForHistograms).boxed().toArray( Integer[]::new );

    static String contextToLabel(final String context, final Nucleotide altAllele, final ArtifactType type){
        return String.join(FIELD_SEPARATOR, context, altAllele.toString(), type.toString());
    }

    static Triple<String, Nucleotide, ArtifactType> labelToContext(final String label){
        final String[] parts = label.split(FIELD_SEPARATOR);
        Utils.validate(parts.length == 3, "Invalid label: " + label);
        return new ImmutableTriple<>(parts[0], Nucleotide.valueOf(parts[1]), ArtifactType.valueOf(parts[2]));
    }

    static final int numSubHistograms = (REGULAR_BASES.size() - 1) * (ArtifactType.values().length);

    static Histogram<Integer> initializeRefHistogram(final String refContext){
        final Histogram<Integer> h = new Histogram<>(binName, refContext);
        h.prefillBins(bins);
        return h;
    }

    static Histogram<Integer> initializeAltHistogram(final String refContext, final Nucleotide altAllele, final ArtifactType type){
        final Histogram<Integer> h = new Histogram<>(binName, contextToLabel(refContext, altAllele, type));
        h.prefillBins(bins);
        return h;
    }
}
