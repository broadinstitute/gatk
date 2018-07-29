package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;


public class EvidenceOverlapChecker {
    private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
    private final ReadMetadata readMetadata;
    private final int minEvidenceMapQ;

    EvidenceOverlapChecker(final Iterator<BreakpointEvidence> evidenceItr, final ReadMetadata readMetadata,
                           final int minEvidenceMapQ) {
        this.readMetadata = readMetadata;
        this.minEvidenceMapQ = minEvidenceMapQ;
        evidenceTree = new SVIntervalTree<>();
        Utils.stream(evidenceItr).forEach(
                evidence -> addToTree(evidenceTree, evidence.getLocation(), evidence)
        );
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

    EvidenceOverlapChecker.OverlapperIterator overlappers(final BreakpointEvidence evidence) {
        return new OverlapperIterator(evidence, evidenceTree);
    }

    // returns iterator to overlapping evidence, paired with a Boolean that is true if the evidence is also coherent
    OverlapAndCoherenceIterator overlappersWithCoherence(final BreakpointEvidence evidence) {
        return new OverlapAndCoherenceIterator(evidence, evidenceTree, readMetadata, minEvidenceMapQ);
    }

    Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> getTreeIterator() {
        return evidenceTree.iterator();
    }

    static class OverlapperIterator  implements Iterator<BreakpointEvidence> {
        private final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
        private Iterator<BreakpointEvidence> listItr;

        OverlapperIterator(final BreakpointEvidence evidence,
                           final SVIntervalTree<List<BreakpointEvidence>> evidenceTree) {
            treeItr = evidenceTree.overlappers(evidence.getLocation());
            listItr = null;
        }

        @Override
        public boolean hasNext() {
            if ( listItr != null && listItr.hasNext() ) {
                return true;
            }
            while(treeItr.hasNext()) {
                listItr = treeItr.next().getValue().iterator();
                if(listItr.hasNext()) {
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
    }

    /**
     * Implements iterator of BreakpointEvidence that overlaps a specified piece of evidence. Each call to next()
     * returns overlapping BreakpointEvidence paired with a Boolean that is true if the overlapper is also coherent
     * with the evidence.
     */
    static class OverlapAndCoherenceIterator implements Iterator<ImmutablePair<BreakpointEvidence, Boolean>> {
        private final ReadMetadata readMetadata;
        private final int minEvidenceMapQ;
        private final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
        private Iterator<BreakpointEvidence> listItr;
        private final boolean checkCoherence;
        private final Boolean isForwardStrand;
        private final List<StrandedInterval> distalTargets;

        OverlapAndCoherenceIterator(final BreakpointEvidence evidence,
                                    final SVIntervalTree<List<BreakpointEvidence>> evidenceTree,
                                    final ReadMetadata readMetadata,
                                    final int minEvidenceMapQ) {
            this.readMetadata = readMetadata;
            this.minEvidenceMapQ = minEvidenceMapQ;
            treeItr = evidenceTree.overlappers(evidence.getLocation());
            isForwardStrand = evidence.isEvidenceUpstreamOfBreakpoint();
            checkCoherence = (isForwardStrand != null) && evidence.hasDistalTargets(readMetadata, minEvidenceMapQ);
            distalTargets = checkCoherence ?
                    new ArrayList<>(evidence.getDistalTargets(readMetadata, minEvidenceMapQ)) : null;
        }

        @Override
        public boolean hasNext() {
            if ( listItr != null && listItr.hasNext() ) {
                return true;
            }
            while(treeItr.hasNext()) {
                listItr = treeItr.next().getValue().iterator();
                if(listItr.hasNext()) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public ImmutablePair<BreakpointEvidence, Boolean> next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException("No next element.");
            }
            final BreakpointEvidence overlapper = listItr.next();
            Boolean isCoherent = false;
            if(checkCoherence) {
                Boolean overlapperStrand = overlapper.isEvidenceUpstreamOfBreakpoint();
                if(overlapperStrand == isForwardStrand && overlapper.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
                    for(final StrandedInterval distalTarget : distalTargets) {
                        for(final StrandedInterval overlapperDistalTarget : overlapper.getDistalTargets(readMetadata, minEvidenceMapQ)) {
                            if(distalTarget.getStrand() == overlapperDistalTarget.getStrand()
                                    && distalTarget.getInterval().overlaps(overlapperDistalTarget.getInterval())) {
                                isCoherent = true;
                                break;
                            }
                        }
                    }
                }
            }
            return new ImmutablePair<>(overlapper, isCoherent);
        }
    }
}
