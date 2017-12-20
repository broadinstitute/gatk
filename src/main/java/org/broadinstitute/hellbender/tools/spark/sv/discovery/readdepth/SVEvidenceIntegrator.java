package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.DiscoverVariantsFromReadDepthArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Integrates split read, read pair, structural variant call, and assembled breakpoint pair evidence to generate a sequence graph.
 */
public final class SVEvidenceIntegrator {

    private static final Logger logger = LogManager.getLogger(SVEvidenceIntegrator.class);
    private final SVIntervalTree<Object> highCoverageIntervalTree;
    private final SVIntervalTree<Object> blacklistTree;
    private final Collection<VariantContext> structuralVariantCalls;
    private final Collection<EvidenceTargetLink> evidenceTargetLinks;
    private final Collection<BreakpointPair> pairedBreakpoints;
    private final Collection<EventRecord> eventRecordCollection;
    private final SVMultiscaleIntervalTree<SVCopyNumberInterval> copyNumberIntervalTree;
    private final SAMSequenceDictionary dictionary;
    private final DiscoverVariantsFromReadDepthArgumentCollection arguments;

    public SVEvidenceIntegrator(final Collection<VariantContext> breakpoints,
                                final Collection<VariantContext> structuralVariantCalls,
                                final Collection<EvidenceTargetLink> evidenceTargetLinks,
                                final Collection<EventRecord> eventRecordCollection,
                                final SVMultiscaleIntervalTree<SVCopyNumberInterval> copyNumberIntervalTree,
                                final Collection<SVInterval> highCoverageIntervals,
                                final Collection<SVInterval> blacklist,
                                final SAMSequenceDictionary dictionary,
                                final DiscoverVariantsFromReadDepthArgumentCollection arguments) {
        Utils.nonNull(breakpoints, "Breakpoints collection cannot be null");
        Utils.nonNull(structuralVariantCalls, "Structural variant calls collection cannot be null");
        Utils.nonNull(evidenceTargetLinks, "Evidence target links collection cannot be null");
        Utils.nonNull(eventRecordCollection, "Event record collection cannot be null");
        Utils.nonNull(copyNumberIntervalTree, "Copy number posteriors tree cannot be null");
        Utils.nonNull(highCoverageIntervals, "High coverage intervals collection cannot be null");
        Utils.nonNull(blacklist, "Blacklist intervals collection cannot be null");
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        Utils.nonNull(arguments, "Parameter arguments collection cannot be null");

        this.dictionary = dictionary;
        this.arguments = arguments;
        this.structuralVariantCalls = structuralVariantCalls;
        this.evidenceTargetLinks = evidenceTargetLinks;
        this.eventRecordCollection = eventRecordCollection;

        pairedBreakpoints = getBreakpointPairs(breakpoints);
        highCoverageIntervalTree = SVIntervalUtils.buildIntervalTreeWithNullValues(highCoverageIntervals);
        blacklistTree = SVIntervalUtils.buildIntervalTreeWithNullValues(blacklist);
        this.copyNumberIntervalTree = copyNumberIntervalTree;
    }

    /**
     * Returns true if the two BNDs in vc1 and vc2 are a valid breakpoint pair, as indicated by their MATEID attributes
     */
    private static boolean isBreakpointPair(final VariantContext vc1, final VariantContext vc2) {
        return vc1.hasAttribute(GATKSVVCFConstants.BND_MATEID_STR) &&
                vc1.getAttributeAsString(GATKSVVCFConstants.BND_MATEID_STR, "").equals(vc2.getID()) &&
                vc2.getAttributeAsString(GATKSVVCFConstants.BND_MATEID_STR, "").equals(vc1.getID());
    }

    /**
     * Returns paired breakpoints from BND records
     */
    private Collection<BreakpointPair> getBreakpointPairs(final Collection<VariantContext> breakpoints) {
        final Map<String, VariantContext> unpairedVariants = new HashMap<>();
        final Collection<BreakpointPair> pairedBreakpoints = new ArrayList<>(breakpoints.size() / 2);
        final Iterator<VariantContext> breakpointIter = breakpoints.iterator();
        while (breakpointIter.hasNext()) {
            final VariantContext vc1 = breakpointIter.next();
            if (!vc1.hasAttribute(GATKSVVCFConstants.BND_MATEID_STR)) continue;
            final String mate = vc1.getAttributeAsString(GATKSVVCFConstants.BND_MATEID_STR, "");
            if (unpairedVariants.containsKey(mate)) {
                final VariantContext vc2 = unpairedVariants.remove(mate);
                if (isBreakpointPair(vc1, vc2)) {
                    pairedBreakpoints.add(new BreakpointPair(vc1, vc2, dictionary));
                } else {
                    throw new UserException.BadInput("Variant mate attributes did not match: " + vc1 + "\t" + vc2);
                }
            } else {
                unpairedVariants.put(vc1.getID(), vc1);
            }
        }
        if (!unpairedVariants.isEmpty()) {
            logger.warn("There were " + unpairedVariants.size() + " unpaired breakpoint variants with a " + GATKSVVCFConstants.BND_MATEID_STR + " attribute.");
        }
        return pairedBreakpoints;
    }

    /**
     * Gets list of breakpoint edges from structural variant calls and adds them to the provided tree
     */
    private Collection<CoordinateSVGraphEdge> getStructuralVariantCallEdges(final SVIntervalTree<CoordinateSVGraphEdge> edgeTree, final int treeOverlapPadding) {
        final Collection<CoordinateSVGraphEdge> edges = new ArrayList<>();
        for (final VariantContext variantContext : structuralVariantCalls) {
            if (variantContext.getStructuralVariantType() == StructuralVariantType.DEL
                    || variantContext.getStructuralVariantType() == StructuralVariantType.DUP
                    || variantContext.getStructuralVariantType() == StructuralVariantType.INV) {
                final int contig = dictionary.getSequenceIndex(variantContext.getContig());
                final int start = variantContext.getStart();
                final int end = variantContext.getEnd();
                final String strandInfo = variantContext.getAttributeAsString("STRANDS", "");
                if (strandInfo.isEmpty()) {
                    throw new UserException.BadInput("VCF record did not contain STRANDS field");
                }
                if (strandInfo.length() != 2) {
                    throw new UserException.BadInput("Expected VCF record STRANDS field to be 2 characters long but found: \"" + strandInfo + "\"");
                }
                final boolean strand1 = strandInfo.charAt(0) == '+' ? true : false;
                final boolean strand2 = strandInfo.charAt(1) == '+' ? true : false;
                final CoordinateSVGraphEdge edge = new CoordinateSVGraphEdge(contig, start, strand1, contig, end, strand2, false, new SVGraphEdgeEvidence(variantContext), dictionary);
                if (!hasEdgeOverlap(edgeTree, edge, treeOverlapPadding)) {
                    edges.add(edge);
                    edgeTree.put(edge.getIntervalA(), edge);
                    edgeTree.put(edge.getIntervalB(), edge);
                }
            }
        }
        return edges;
    }

    /**
     * Gets list of breakpoint edges from BND records and adds them to the provided tree. An edge is not created if a similar edge already exists in the tree.
     */
    private Collection<CoordinateSVGraphEdge> getBreakpointPairEdges(final SVIntervalTree<CoordinateSVGraphEdge> edgeTree, final int treeOverlapPadding) {
        final Collection<CoordinateSVGraphEdge> edges = new ArrayList<>();
        for (final BreakpointPair breakpointPair : pairedBreakpoints) {
            if (breakpointPair.getContigA() == breakpointPair.getContigB()) { //Only intrachromosomal events
                final CoordinateSVGraphEdge edge = new CoordinateSVGraphEdge(breakpointPair.getContigA(), breakpointPair.getPositionA(), breakpointPair.isStrandA(),
                        breakpointPair.getContigB(), breakpointPair.getPositionB(), breakpointPair.isStrandB(), false, new SVGraphEdgeEvidence(breakpointPair), dictionary);
                if (!hasEdgeOverlap(edgeTree, edge, treeOverlapPadding)) {
                    edges.add(edge);
                    edgeTree.put(edge.getIntervalA(), edge);
                    edgeTree.put(edge.getIntervalB(), edge);
                }
            }
        }
        return edges;
    }

    /**
     * Gets list of breakpoint edges from evidence target links and adds them to the provided tree. An edge is not created if a similar edge already exists in the tree.
     */
    private Collection<CoordinateSVGraphEdge> getEvidenceTargetLinkEdges(final SVIntervalTree<CoordinateSVGraphEdge> edgeTree, final int treeOverlapPadding) {
        final Collection<CoordinateSVGraphEdge> edges = new ArrayList<>();
        for (final EvidenceTargetLink link : evidenceTargetLinks) {
            if (link.getReadPairs() >= arguments.minLinkReadPairs) {
                final boolean leftStrand = link.getPairedStrandedIntervals().getLeft().getStrand();
                final boolean rightStrand = link.getPairedStrandedIntervals().getRight().getStrand();
                final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
                final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
                if (leftInterval.getContig() == rightInterval.getContig()) {    //Only intrachromosomal events

                    final int leftPosition = leftStrand ? leftInterval.getStart() : leftInterval.getEnd();
                    final int rightPosition = rightStrand ? rightInterval.getStart() : rightInterval.getEnd();
                    final CoordinateSVGraphEdge edge = new CoordinateSVGraphEdge(leftInterval.getContig(), leftPosition, leftStrand, rightInterval.getContig(), rightPosition, rightStrand, false, new SVGraphEdgeEvidence(link), dictionary);
                    if (!hasEdgeOverlap(edgeTree, edge, treeOverlapPadding)) {
                        edgeTree.put(edge.getIntervalA(), edge);
                        edgeTree.put(edge.getIntervalB(), edge);
                        edges.add(edge);
                    }
                    /*
                    for (int i = leftInterval.getStart(); i <= leftInterval.getEnd(); i += leftInterval.getLength()) {
                        for (int j = rightInterval.getStart(); j <= rightInterval.getEnd(); j += rightInterval.getLength()) {
                            final CoordinateSVGraphEdge edge = new CoordinateSVGraphEdge(leftInterval.getContig(), i, leftStrand, rightInterval.getContig(), j, rightStrand, false, new SVGraphEdgeEvidence(link), dictionary);
                            edgeTree.put(edge.getIntervalA(), edge);
                            edgeTree.put(edge.getIntervalB(), edge);
                            edges.add(edge);
                        }
                    }
                    */
                }
            }
        }
        return edges;
    }

    private Collection<CoordinateSVGraphEdge> getEventRecordEdges() {
        final Collection<CoordinateSVGraphEdge> edges = new ArrayList<>(eventRecordCollection.size());
        for (final EventRecord record : eventRecordCollection) {
            final int contig = dictionary.getSequenceIndex(record.getContig());
            final int start = record.getStart();
            final int end = record.getEnd();
            if (record.getType().equals(StructuralVariantType.DEL.name())) {
                edges.add(new CoordinateSVGraphEdge(contig, start, true, contig, end, false, false, new SVGraphEdgeEvidence(record), dictionary));
            } else if (record.getType().equals(StructuralVariantType.DUP.name())) {
                edges.add(new CoordinateSVGraphEdge(contig, start, false, contig, end, true, false, new SVGraphEdgeEvidence(record), dictionary));
            } else if (record.getType().equals(StructuralVariantType.INV.name())) {
                edges.add( new CoordinateSVGraphEdge(contig, start, true, contig, end, true, false, new SVGraphEdgeEvidence(record), dictionary));
                edges.add(new CoordinateSVGraphEdge(contig, start, false, contig, end, false, false, new SVGraphEdgeEvidence(record), dictionary));
            }
        }
        return edges;
    }

    /**
     * Filters edges that do not have CNV calls or whose ends overlap blacklist or high coverage intervals
     */
    private List<CoordinateSVGraphEdge> filterEdges(final Collection<CoordinateSVGraphEdge> edges, final int filterPadding) {
        return edges.stream()
                .filter(edge -> copyNumberIntervalTree.hasOverlapper(edge.getInterval()))
                .filter(edge -> !edgeEndpointsHaveTreeOverlap(edge, blacklistTree, filterPadding))
                .filter(edge -> !edgeEndpointsHaveTreeOverlap(edge, highCoverageIntervalTree, filterPadding))
                .collect(Collectors.toList());
    }

    /**
     * Returns true if either endpoint of the edge (give or take padding) overlaps with an interval in the provided tree.
     */
    private boolean edgeEndpointsHaveTreeOverlap(final CoordinateSVGraphEdge edge, final SVIntervalTree<Object> tree, final int endpointPadding) {
        return tree.hasOverlapper(SVIntervalUtils.getPaddedInterval(edge.getIntervalA(), endpointPadding, dictionary))
                || tree.hasOverlapper(SVIntervalUtils.getPaddedInterval(edge.getIntervalB(), endpointPadding, dictionary));
    }


    /**
     * Creates graph from structural variant calls, BNDs, and evidence target links
     */
    public SVGraph buildGraph() {
        final Collection<CoordinateSVGraphEdge> edges = new ArrayList<>();
        final SVIntervalTree<CoordinateSVGraphEdge> edgeTree = new SVIntervalTree<>();

        edges.addAll(getBreakpointPairEdges(edgeTree, arguments.insertSize));
        edges.addAll(getStructuralVariantCallEdges(edgeTree, arguments.insertSize));
        edges.addAll(getEvidenceTargetLinkEdges(edgeTree, arguments.insertSize));
        edges.addAll(getEventRecordEdges());

        final List<CoordinateSVGraphEdge> filteredEdges = filterEdges(edges, arguments.insertSize);
        return new SVGraph(filteredEdges, dictionary);
    }

    /**
     * Determines if a similar edges already exists in the provided tree, so as to not introduce duplicate evidence.
     */
    private boolean hasEdgeOverlap(final SVIntervalTree<CoordinateSVGraphEdge> tree, final CoordinateSVGraphEdge edge, final int padding) {
        final SVInterval intervalA = SVIntervalUtils.getPaddedInterval(edge.getIntervalA(), padding, dictionary);
        final SVInterval intervalB = SVIntervalUtils.getPaddedInterval(edge.getIntervalB(), padding, dictionary);
        final Stream<SVIntervalTree.Entry<CoordinateSVGraphEdge>> overlapStreamA = Utils.stream(tree.overlappers(intervalA));
        final Stream<SVIntervalTree.Entry<CoordinateSVGraphEdge>> overlapStreamB = Utils.stream(tree.overlappers(intervalB));
        return Stream.concat(overlapStreamA, overlapStreamB).map(entry -> entry.getValue())
                .anyMatch(overlapper -> overlapper.isStrandA() == edge.isStrandA() && overlapper.isStrandB() == edge.isStrandB()
                        && overlapper.getIntervalA().overlaps(intervalA) && overlapper.getIntervalB().overlaps(intervalB));
    }

}
