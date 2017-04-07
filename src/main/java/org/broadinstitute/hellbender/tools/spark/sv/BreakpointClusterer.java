package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.*;
import java.util.function.Function;

/**
 * A class that acts as a filter for breakpoint evidence.
 * It passes only that evidence that is part of a putative cluster.
 */
public final class BreakpointClusterer implements Function<BreakpointEvidence, Iterator<BreakpointEvidence>> {
    private final int minEvidenceCount;
    private final int staleEventDistance;
    private final SortedMap<BreakpointEvidence, Boolean> locMap = new TreeMap<>();
    private final List<Map.Entry<BreakpointEvidence, Boolean>> reportableEntries;
    private static final Iterator<BreakpointEvidence> noEvidence = Collections.emptyIterator();
    private int currentContig = -1;

    public BreakpointClusterer( final int minEvidenceCount, final int staleEventDistance ) {
        this.minEvidenceCount = minEvidenceCount;
        this.staleEventDistance = staleEventDistance;
        this.reportableEntries = new ArrayList<>(2 * minEvidenceCount);
    }

    @Override
    public Iterator<BreakpointEvidence> apply( final BreakpointEvidence evidence ) {
        if ( evidence.getContigIndex() != currentContig ) {
            currentContig = evidence.getContigIndex();
            locMap.clear();
        }

        locMap.put(evidence, true);

        final int locusStart = evidence.getEventStartPosition();
        final int locusEnd = evidence.getContigEnd();
        final int staleEnd = locusStart - staleEventDistance;
        int evidenceCount = 0;
        reportableEntries.clear();
        final Iterator<Map.Entry<BreakpointEvidence, Boolean>> itr = locMap.entrySet().iterator();
        while ( itr.hasNext() ) {
            final Map.Entry<BreakpointEvidence, Boolean> entry = itr.next();
            final BreakpointEvidence evidence2 = entry.getKey();
            final int contigEnd = evidence2.getContigEnd();
            if ( contigEnd <= staleEnd ) itr.remove();
            else if ( evidence2.getEventStartPosition() >= locusEnd ) break;
            else if ( contigEnd > locusStart ) {
                evidenceCount += 1;
                if ( entry.getValue() ) reportableEntries.add(entry);
            }
        }

        if ( evidenceCount >= minEvidenceCount ) {
            return reportableEntries.stream()
                    .map(entry -> {
                        entry.setValue(false);
                        return entry.getKey();
                    })
                    .iterator();
        }
        return noEvidence;
    }
}
