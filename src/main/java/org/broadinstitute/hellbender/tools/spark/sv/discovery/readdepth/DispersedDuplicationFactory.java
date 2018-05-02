package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;
import java.util.stream.Collectors;

/**
 * Calls large tandem duplication variants
 */
public class DispersedDuplicationFactory extends LargeSimpleSVFactory {

    public DispersedDuplicationFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                                       final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                                       final SVIntervalTree<VariantContext> structuralVariantCallTree,
                                       final SVIntervalTree<GATKRead> contigTree,
                                       final StructuralVariationDiscoveryArgumentCollection.StructuralVariantIntervalsForCNV arguments,
                                       final SAMSequenceDictionary dictionary) {
        super(intrachromosomalLinkTree, interchromosomalLinkTree, structuralVariantCallTree, contigTree, arguments, dictionary);
    }

    @Override
    protected LargeSimpleSV getNewSV(final int start,
                                     final int end,
                                     final int contigId,
                                     final String contig,
                                     final int readPairEvidence,
                                     final int splitReadEvidence,
                                     final int readPairCounterEvidence,
                                     final int splitReadCounterEvidence,
                                     final IntrachromosomalBreakpointPair breakpoints,
                                     final Collection<EvidenceTargetLink> supportingEvidence) {
        return new LargeSimpleSV(SimpleSVType.TYPES.DUP, LargeSimpleSV.EvidenceType.READ_PAIR, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints, supportingEvidence);
    }

    @Override
    protected boolean isCounterEvidenceOrientation(final EvidenceTargetLink link) {
        return !hasInnieOrientation(link) && !hasOutieOrientation(link);
    }

    @Override
    protected int countSupportingEvidenceReadPairs(final Collection<EvidenceTargetLink> links) {
        final int numFalseTrue = countReadPairs(links.stream()
                .filter(link -> !link.getPairedStrandedIntervals().getLeft().getStrand()
                        && link.getPairedStrandedIntervals().getRight().getStrand())
                .collect(Collectors.toList()));
        final int numTrueFalse = countReadPairs(links) - numFalseTrue;
        return Math.min(numFalseTrue, numTrueFalse);
    }

    @Override
    protected int countSupportingEvidenceSplitReads(final Collection<EvidenceTargetLink> links) {
        final int numFalseTrue = countSplitReads(links.stream()
                .filter(link -> !link.getPairedStrandedIntervals().getLeft().getStrand()
                        && link.getPairedStrandedIntervals().getRight().getStrand())
                .collect(Collectors.toList()));
        final int numTrueFalse = countSplitReads(links) - numFalseTrue;
        return Math.min(numFalseTrue, numTrueFalse);
    }

    /**
     * Checks if the link matches "outie" read pair orientation, i.e. -/+
     */
    @Override
    protected boolean hasSupportingEvidenceOrientation(final EvidenceTargetLink link) {
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        return intervals.getLeft().getStrand() != intervals.getRight().getStrand(); // +/- or -/+
    }
}
