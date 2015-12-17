package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * Represents a segmented copy-number model of target-coverage and SNP-allele-count data.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentedModel {
    private final List<SimpleInterval> segments;
    private final Genome genome;

    /**
     * Constructs a segmented model from a segment file and a genome.
     * @param segmentFile   segment file (only "Sample", "Chromosome", "Start", "End" columns are used)
     * @param genome        genome containing target coverages and SNP allele counts to be segmented
     */
    public SegmentedModel(final File segmentFile, final Genome genome) {
        this(SegmentUtils.readIntervalsFromSegmentFile(segmentFile), genome);
    }

    /**
     * Constructs a segmented model from a list of segments and a genome.
     * Assumes the list of segments is well formed and non-shared.
     * @param segments  list of segments, cannot be {@code null}
     * @param genome    genome containing target coverages and SNP allele counts to be segmented, cannot be {@code null}
     */
    public SegmentedModel(final List<SimpleInterval> segments, final Genome genome) {
        Utils.nonNull(segments, "List of segments cannot be null.");
        Utils.nonNull(genome, "Genome cannot be null.");
        this.segments = Collections.unmodifiableList(segments);
        this.genome = genome;
    }

    /**
     * Returns an unmodifiable and immutable view of the list of segments.
     * @return an unmodifiable and immutable view of the list of segments, never {@code null}
     */
    public List<SimpleInterval> getSegments() {    return segments;    }

    /**
     * Returns the Genome held internally.
     * @return  the Genome held internally
     */
    public Genome getGenome() { return genome;  }

    /**
     * Returns a new SegmentedModel with the small segments merged, using the algorithm implemented in
     * {@link SegmentMergeUtils}.
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      a new SegmentedModel with the small segments merged
     */
    public SegmentedModel mergeSmallSegments(final int targetNumberThreshold) {
        final List<SimpleInterval> mergedSegments =
                SegmentMergeUtils.mergeSmallSegments(segments, genome, targetNumberThreshold);
        return new SegmentedModel(mergedSegments, genome);
    }

    /**
     * Writes the list of segments with number of targets and snps in each segment to file.
     * @param outFile       output file
     */
    void writeSegmentFileWithNumTargetsAndNumSNPs(final File outFile) {
        SegmentUtils.writeSegmentFileWithNumTargetsAndNumSNPs(outFile, segments, genome);
    }
}
