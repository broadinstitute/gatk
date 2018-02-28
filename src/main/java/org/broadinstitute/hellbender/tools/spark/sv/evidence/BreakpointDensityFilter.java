package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * A class that acts as a filter for breakpoint evidence.
 * It passes only that evidence that is part of a putative cluster.
 */
public final class BreakpointDensityFilter implements Iterator<BreakpointEvidence> {
    private final ReadMetadata readMetadata;
    private final double minEvidenceWeight;
    private final double minCoherentEvidenceWeight;
    private final PartitionCrossingChecker partitionCrossingChecker;
    private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
    private final int minEvidenceMapq;
    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;

    public BreakpointDensityFilter( final Iterator<BreakpointEvidence> evidenceItr,
                                    final ReadMetadata readMetadata,
                                    final double minEvidenceWeightPerCoverage,
                                    final double minCoherentEvidenceWeightPerCoverage,
                                    final PartitionCrossingChecker partitionCrossingChecker,
                                    final int minEvidenceMapq) {
        this.readMetadata = readMetadata;
        this.minEvidenceWeight = minEvidenceWeightPerCoverage * readMetadata.getCoverage();
        this.minCoherentEvidenceWeight = minCoherentEvidenceWeightPerCoverage * readMetadata.getCoverage();
        this.partitionCrossingChecker = partitionCrossingChecker;
        this.evidenceTree = buildTree(evidenceItr);
        this.treeItr = evidenceTree.iterator();
        this.listItr = null;
        this.minEvidenceMapq = minEvidenceMapq;
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
        PairedStrandedIntervalTree<BreakpointEvidence> targetIntervalTree = new PairedStrandedIntervalTree<>();
        int weight = 0;
        while ( itr.hasNext() ) {
            final List<BreakpointEvidence> evidenceForInterval = itr.next().getValue();
            weight += evidenceForInterval.stream().mapToInt(BreakpointEvidence::getWeight).sum();
            if ( weight >= minEvidenceWeight) {
                return true;
            }

            for (final BreakpointEvidence evidence : evidenceForInterval) {
                if (evidence.hasDistalTargets(readMetadata, minEvidenceMapq)) {
                    final List<StrandedInterval> distalTargets = evidence.getDistalTargets(readMetadata, minEvidenceMapq);
                    for (int i = 0; i < distalTargets.size(); i++) {
                        targetIntervalTree.put(
                                new PairedStrandedIntervals(
                                        new StrandedInterval(evidence.getLocation(), evidence.isEvidenceUpstreamOfBreakpoint()),
                                        distalTargets.get(i)),
                                evidence
                                );
                    }
                }
            }
        }

        final Iterator<Tuple2<PairedStrandedIntervals, BreakpointEvidence>> targetLinkIterator = targetIntervalTree.iterator();
        while (targetLinkIterator.hasNext()) {
            Tuple2<PairedStrandedIntervals, BreakpointEvidence> next = targetLinkIterator.next();
            final int coherentEvidenceWeight = (int) Utils.stream(targetIntervalTree.overlappers(next._1())).count();
            if (coherentEvidenceWeight >= minCoherentEvidenceWeight) {
                return true;
            }

        }

        return false;
    }

    private static <T> void addToTree( final SVIntervalTree<List<T>> tree,
                                   final SVInterval interval,
                                   final T value ) {
        final SVIntervalTree.Entry<List<T>> entry = tree.find(interval);
        if ( entry != null ) {
            entry.getValue().add(value);
        } else {
            final List<T> valueList = new ArrayList<>(1);
            valueList.add(value);
            tree.put(interval, valueList);
        }
    }
}
