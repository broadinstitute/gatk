package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.*;
import java.util.function.Function;

/**
 * Class to fully sort a stream of nearly sorted BreakpointEvidences.
 */
public final class WindowSorter implements Function<BreakpointEvidence, Iterator<BreakpointEvidence>> {
    private final SortedSet<BreakpointEvidence> recordSet = new TreeSet<>();
    private final List<BreakpointEvidence> reportableEvidence = new ArrayList<>();
    private final int windowSize;
    private int currentContig = -1;

    public WindowSorter( final int windowSize ) {
        this.windowSize = windowSize;
    }

    @Override
    public Iterator<BreakpointEvidence> apply( final BreakpointEvidence evidence ) {
        reportableEvidence.clear();
        if ( evidence.getContigIndex() != currentContig ) {
            reportableEvidence.addAll(recordSet);
            recordSet.clear();
            currentContig = evidence.getContigIndex();
        } else {
            final int reportableEnd = evidence.getEventStartPosition() - windowSize;
            final Iterator<BreakpointEvidence> itr = recordSet.iterator();
            while ( itr.hasNext() ) {
                final BreakpointEvidence evidence2 = itr.next();
                if ( evidence2.getEventStartPosition() >= reportableEnd ) break;
                reportableEvidence.add(evidence2);
                itr.remove();
            }
        }
        recordSet.add(evidence);
        return reportableEvidence.iterator();
    }
}
