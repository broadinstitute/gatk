package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * A class that acts as a filter for breakpoint evidence.
 * It passes only that evidence that is part of a putative cluster.
 */
public final class BreakpointDensityFilter implements Iterator<BreakpointEvidence> {
    private final ReadMetadata readMetadata;
    private final int minEvidenceWeight;
    private final int minCoherentEvidenceWeight;
    private final PartitionCrossingChecker partitionCrossingChecker;
    private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;

    public BreakpointDensityFilter( final Iterator<BreakpointEvidence> evidenceItr,
                                    final ReadMetadata readMetadata,
                                    final int minEvidenceWeight,
                                    final int minCoherentEvidenceWeight,
                                    final PartitionCrossingChecker partitionCrossingChecker ) {
        this.readMetadata = readMetadata;
        this.minEvidenceWeight = minEvidenceWeight;
        this.minCoherentEvidenceWeight = minCoherentEvidenceWeight;
        this.partitionCrossingChecker = partitionCrossingChecker;
        this.evidenceTree = buildTree(evidenceItr);
        this.treeItr = evidenceTree.iterator();
        this.listItr = null;
    }

    @Override
    public boolean hasNext() {
        if ( listItr != null && listItr.hasNext() ) {
            return true;
        }
        listItr = null;
        boolean result = false;
        while ( !result && treeItr.hasNext() ) {
            final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = treeItr.next();
            final SVInterval curInterval = entry.getInterval();
            if ( isValidated(entry.getValue()) || hasEnoughOverlappers(curInterval) ) {
                entry.getValue().forEach(ev -> ev.setValidated(true));
                result = true;
            } else if ( partitionCrossingChecker.onBoundary(curInterval) ) {
                result = true;
            }
            if ( result ) {
                listItr = entry.getValue().iterator();
            }
        }
        return result;
    }

    @Override
    public BreakpointEvidence next() {
        if ( !hasNext() ) {
            throw new NoSuchElementException("No next element.");
        }
        return listItr.next();
    }

    private static SVIntervalTree<List<BreakpointEvidence>> buildTree( final Iterator<BreakpointEvidence> evidenceItr ) {
        SVIntervalTree<List<BreakpointEvidence>> tree = new SVIntervalTree<>();
        while ( evidenceItr.hasNext() ) {
            final BreakpointEvidence evidence = evidenceItr.next();
            addToTree(tree, evidence.getLocation(), evidence);
        }
        return tree;
    }

    private boolean isValidated( final List<BreakpointEvidence> evList ) {
        for ( final BreakpointEvidence ev : evList ) {
            if ( ev.isValidated() ) return true;
        }
        return false;
    }

    @VisibleForTesting boolean hasEnoughOverlappers( final SVInterval interval ) {
        final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> itr = evidenceTree.overlappers(interval);
        final SVIntervalTree<List<BreakpointEvidence>> targetIntervalTree = new SVIntervalTree<>();
        int weight = 0;
        while ( itr.hasNext() ) {
            final List<BreakpointEvidence> evidenceForInterval = itr.next().getValue();
            weight += evidenceForInterval.stream().mapToInt(BreakpointEvidence::getWeight).sum();
            if ( weight >= minEvidenceWeight ) {
                return true;
            }

            for (final BreakpointEvidence evidence : evidenceForInterval) {
                if (evidence.hasDistalTargets()) {
                    for (final SVInterval target : evidence.getDistalTargets(readMetadata)) {
                        addToTree(targetIntervalTree, target, evidence);
                    }
                }
            }
        }

        for (final SVIntervalTree.Entry<List<BreakpointEvidence>> targetTreeEntries : targetIntervalTree) {
            final SVInterval target = targetTreeEntries.getInterval();
            final int coherentEvidenceWeight =
                    Utils.stream(targetIntervalTree.overlappers(target))
                            .mapToInt(overlapper ->
                                    overlapper.getValue().stream().mapToInt(BreakpointEvidence::getWeight).sum())
                            .sum();
            if (coherentEvidenceWeight >= minCoherentEvidenceWeight) {
                return true;
            }
        }
        return false;
    }

    private static void addToTree( final SVIntervalTree<List<BreakpointEvidence>> tree,
                                   final SVInterval interval,
                                   final BreakpointEvidence evidence ) {
        final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = tree.find(interval);
        if ( entry != null ) {
            entry.getValue().add(evidence);
        } else {
            final List<BreakpointEvidence> valueList = new ArrayList<>(1);
            valueList.add(evidence);
            tree.put(interval, valueList);
        }
    }
}
