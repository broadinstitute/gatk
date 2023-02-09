package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;

public class IntervalOverlapCalculator {

    final SAMSequenceDictionary dictionary;
    final SVIntervalTree<Object> tree;

    public IntervalOverlapCalculator(final Collection<SVInterval> intervals, final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
        tree = new SVIntervalTree<>();
        intervals.stream().forEach(interval -> tree.put(interval, null));
    }

    public static IntervalOverlapCalculator create(final GATKPath path,
                                                   final SAMSequenceDictionary dictionary,
                                                   final IntervalSetRule intervalSetRule,
                                                   final IntervalMergingRule intervalMergingRule,
                                                   final int regionPadding) {
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final GenomeLocSortedSet includeSet = IntervalUtils.loadIntervals(Collections.singletonList(path.toString()), intervalSetRule, intervalMergingRule, regionPadding, parser);
        final List<SVInterval> intervals = convertGenomeLocsToSVIntervals(includeSet.toList());
        Utils.validate(!intervals.isEmpty(), "Intervals are empty: " + path);
        return new IntervalOverlapCalculator(intervals, dictionary);
    }

    private static List<SVInterval> convertGenomeLocsToSVIntervals( final List<GenomeLoc> genomeLocIntervals) {
        final List<SVInterval> convertedIntervals = new ArrayList<>(genomeLocIntervals.size());
        for ( final GenomeLoc genomeLoc : genomeLocIntervals ) {
            if ( genomeLoc.isUnmapped() ) {
                throw new UserException("Unmapped intervals cannot be converted to SVIntervals");
            }
            convertedIntervals.add(new SVInterval(genomeLoc.getContigIndex(), genomeLoc.getStart(), genomeLoc.getStop()));
        }
        return convertedIntervals;
    }

    public int getEndpointOverlapCount(final SVCallRecord record) {
        // Breakend overlap
        final SVInterval intervalA = new SVInterval(dictionary.getSequenceIndex(record.getContigA()), record.getPositionA(), record.getPositionA());
        final SVInterval intervalB = new SVInterval(dictionary.getSequenceIndex(record.getContigB()), record.getPositionB(), record.getPositionB());
        final int overlapsA = tree.overlappers(intervalA).hasNext() ? 1 : 0;
        final int overlapsB = tree.overlappers(intervalB).hasNext() ? 1 : 0;
        return overlapsA + overlapsB;
    }

    public Double getOverlapFraction(final SVCallRecord record) {
        // Total overlap as fraction of variant length
        final GATKSVVCFConstants.StructuralVariantAnnotationType type = record.getType();
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.BND
                || type == GATKSVVCFConstants.StructuralVariantAnnotationType.CTX) {
            return getEndpointOverlapCount(record) / 2.0;
        } else {
            if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INS || record.isDispersedDup()) {
                return getEndpointOverlapCount(record) > 0 ? 1. : 0.;
            }
        }
        final Integer length = record.getLength();
        if (length == null) {
            throw new IllegalArgumentException("Expected record length to be defined for " + record.getId());
        }
        return totalOverlap(record) / (double) length;
    }

    private long totalOverlap(final SVCallRecord record) {
        Utils.validate(record.isIntrachromosomal(), "Record must be intra-chromosomal");
        final SVInterval interval = new SVInterval(dictionary.getSequenceIndex(record.getContigA()), record.getPositionA(), record.getPositionB());
        final Iterator<SVIntervalTree.Entry<Object>> iter = tree.overlappers(interval);
        long overlap = 0;
        while (iter.hasNext()) {
            overlap += interval.overlapLen(iter.next().getInterval());
        }
        return overlap;
    }
}
