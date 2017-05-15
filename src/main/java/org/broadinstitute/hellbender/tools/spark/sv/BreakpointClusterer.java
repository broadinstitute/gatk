package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * A class that acts as a filter for breakpoint evidence.
 * It passes only that evidence that is part of a putative cluster.
 */
public final class BreakpointClusterer implements Iterator<BreakpointEvidence> {
    private final int minEvidenceCount;
    private final int minCoherentEvidenceCount;
    private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
    private final ReadMetadata readMetadata;
    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;

    public BreakpointClusterer(final ReadMetadata readMetadata, final int minEvidenceCount, final int minCoherentEvidenceCount,  final Iterator<BreakpointEvidence> evidenceItr ) {
        this.minEvidenceCount = minEvidenceCount;
        this.minCoherentEvidenceCount = minCoherentEvidenceCount;
        this.evidenceTree = new SVIntervalTree<>();
        this.readMetadata = readMetadata;
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

    @VisibleForTesting
    protected boolean hasEnoughOverlappers( final SVInterval interval ) {
        final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> itr = evidenceTree.overlappers(interval);
        final SVIntervalTree<List<BreakpointEvidence>> targetIntervalTree = new SVIntervalTree<>();
        int evidenceIntervalOverlapCount = 0;
        while ( itr.hasNext() ) {
            final List<BreakpointEvidence> evidenceForInterval = itr.next().getValue();
            evidenceIntervalOverlapCount += evidenceForInterval.size();
            if ( evidenceIntervalOverlapCount >= minEvidenceCount ) {
                return true;
            }

            for (final BreakpointEvidence evidence : evidenceForInterval) {
                if (evidence.hasDistalTargets()) {
                    for (final SVInterval target : evidence.getDistalTargets(readMetadata)) {
                        if (targetIntervalTree.find(target) == null) {
                            targetIntervalTree.put(target, new ArrayList<>(1));
                        }
                        targetIntervalTree.find(target).getValue().add(evidence);
                    }
                }
            }
        }

        for (final SVIntervalTree.Entry<List<BreakpointEvidence>> targetTreeEntries : targetIntervalTree) {
            final SVInterval target = targetTreeEntries.getInterval();
            final int coherentEvidenceCount =
                    Utils.stream(targetIntervalTree.overlappers(target))
                            .mapToInt(overlapper -> overlapper.getValue().size()).sum();
            if (coherentEvidenceCount >= minCoherentEvidenceCount) {
                return true;
            }
        }

        return false;
    }
}
