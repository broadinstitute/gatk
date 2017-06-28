package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;

import java.util.*;

/**
 * A class that acts as a filter for breakpoint evidence.
 * It passes only that evidence that is part of a putative cluster.
 */
public final class BreakpointClusterer implements Iterator<BreakpointEvidence> {
    private final int minEvidenceCount;
    private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;

    public BreakpointClusterer( final int minEvidenceCount, final Iterator<BreakpointEvidence> evidenceItr ) {
        this.minEvidenceCount = minEvidenceCount;
        this.evidenceTree = new SVIntervalTree<>();
        buildTree(evidenceItr);
    }

    @Override
    public boolean hasNext() {
        if ( listItr != null && listItr.hasNext() ) {
            return true;
        }
        listItr = null;
        while ( treeItr.hasNext() ) {
            final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = treeItr.next();
            if ( hasEnoughOverlappers(entry.getInterval()) ) {
                listItr = entry.getValue().iterator();
                return true;
            }
        }
        return false;
    }

    @Override
    public BreakpointEvidence next() {
        if ( !hasNext() ) {
            throw new NoSuchElementException("No next element.");
        }
        return listItr.next();
    }

    private void buildTree( final Iterator<BreakpointEvidence> evidenceItr ) {
        while ( evidenceItr.hasNext() ) {
            final BreakpointEvidence evidence = evidenceItr.next();
            final SVInterval location = evidence.getLocation();
            final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = evidenceTree.find(location);
            if ( entry != null ) {
                entry.getValue().add(evidence);
            } else {
                final List<BreakpointEvidence> valueList = new ArrayList<>(1);
                valueList.add(evidence);
                evidenceTree.put(location, valueList);
            }
        }
        treeItr = evidenceTree.iterator();
    }

    private boolean hasEnoughOverlappers( final SVInterval interval ) {
        final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> itr = evidenceTree.overlappers(interval);
        int count = 0;
        while ( itr.hasNext() ) {
            count += itr.next().getValue().size();
            if ( count >= minEvidenceCount ) {
                return true;
            }
        }
        return false;
    }
}
