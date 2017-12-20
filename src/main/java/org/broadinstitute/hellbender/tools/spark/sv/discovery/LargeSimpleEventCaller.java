package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.*;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class LargeSimpleEventCaller {

    private static final Logger logger = LogManager.getLogger(LargeSimpleEventCaller.class);

    private final OverlapDetector<CopyRatio> readDepthOverlapDetector;
    private final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    private final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    private final SVIntervalTree<GATKRead> contigTree;
    private final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector;
    private final SimpleSVFactory tandemDuplicationFactory;
    private final SimpleSVFactory deletionFactory;
    private final SimpleSVFactory inversionFactory;
    private final ReadMetadata readMetadata;
    private final SAMSequenceDictionary dict;
    private final List<IntrachromosomalBreakpointPair> pairedBreakpoints;
    private final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments;

    public LargeSimpleEventCaller(final List<VariantContext> breakpoints,
                                  final List<GATKRead> assembledContigs,
                                  final List<EvidenceTargetLink> evidenceTargetLinks,
                                  final CopyRatioCollection readDepthData,
                                  final CalledCopyRatioSegmentCollection copyRatioSegments,
                                  final ReadMetadata readMetadata,
                                  final SAMSequenceDictionary dict,
                                  final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments) {
        this.readMetadata = readMetadata;
        this.dict = dict;
        this.arguments = arguments;

        final Map<String, VariantContext> unpairedVariants = new HashMap<>();
        pairedBreakpoints = new ArrayList<>(breakpoints.size() / 2);
        final Iterator<VariantContext> breakpointIter = breakpoints.iterator();
        while (breakpointIter.hasNext()) {
            final VariantContext vc1 = breakpointIter.next();
            if (!vc1.hasAttribute("MATEID")) continue;
            final String mate = vc1.getAttributeAsString("MATEID", "");
            if (unpairedVariants.containsKey(mate)) {
                final VariantContext vc2 = unpairedVariants.remove(mate);
                if (isBreakpointPair(vc1, vc2)) {
                    final int contig = readMetadata.getContigID(vc1.getContig());
                    final VariantContext first;
                    final VariantContext second;
                    if (vc1.getStart() < vc2.getStart()) {
                        first = vc1;
                        second = vc2;
                    } else {
                        first = vc2;
                        second = vc1;
                    }
                    final int start = first.getStart();
                    final int end = second.getStart();
                    final List<String> firstContigs = first.getAttributeAsStringList("CTG_NAMES", "");
                    final List<String> secondContigs = second.getAttributeAsStringList("CTG_NAMES", "");
                    pairedBreakpoints.add(new IntrachromosomalBreakpointPair(contig, start, end, firstContigs, secondContigs));
                } else {
                    throw new IllegalStateException("Variant mate attributes did not match: " + vc1 + "\t" + vc2);
                }
            } else {
                unpairedVariants.put(vc1.getID(), vc1);
            }
        }
        if (!unpairedVariants.isEmpty()) {
            logger.warn("There were " + unpairedVariants.size() + " unpaired breakpoint variants with a MATEID");
        }

        final List<EvidenceTargetLink> intrachromosomalEvidenceLinkTargets = evidenceTargetLinks.stream()
                .filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().getContig() == link.getPairedStrandedIntervals().getRight().getInterval().getContig()).collect(Collectors.toList());
        final List<EvidenceTargetLink> interchromosomalEvidenceLinkTargets = evidenceTargetLinks.stream()
                .filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().getContig() != link.getPairedStrandedIntervals().getRight().getInterval().getContig()).collect(Collectors.toList());

        intrachromosomalLinkTree = buildEvidenceIntervalTree(intrachromosomalEvidenceLinkTargets, 0, false);
        interchromosomalLinkTree = buildEvidenceIntervalTree(interchromosomalEvidenceLinkTargets, 0, true);
        contigTree = buildContigIntervalTree(assembledContigs, readMetadata);

        readDepthOverlapDetector = readDepthData.getOverlapDetector();
        copyRatioSegmentOverlapDetector = copyRatioSegments.getOverlapDetector();

        tandemDuplicationFactory = new TandemDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector,
                readDepthOverlapDetector, readMetadata, dict);
        deletionFactory = new DeletionFactory(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector,
                readDepthOverlapDetector, readMetadata, dict);
        inversionFactory = new CopyNeutralInversionFactory(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector,
                readDepthOverlapDetector, readMetadata, dict);

    }

    private static boolean isBreakpointPair(final VariantContext vc1, final VariantContext vc2) {
        return vc1.getAttributeAsString("MATEID", "").equals(vc2.getID()) &&
                vc2.getAttributeAsString("MATEID", "").equals(vc1.getID());
    }

    private static SVIntervalTree<Integer> buildIntervalTree(final GenomeLocSortedSet intervals) {
        final SVIntervalTree<Integer> tree = new SVIntervalTree<>();
        for (final GenomeLoc loc : intervals.toList()) {
            final int start = loc.getStart();
            final int end = loc.getEnd();
            final int contig = loc.getContigIndex();
            tree.put(new SVInterval(contig, start, end), contig);
        }
        return tree;
    }

    private static SVIntervalTree<GATKRead> buildContigIntervalTree(final List<GATKRead> contigs, final ReadMetadata readMetadata) {
        final SVIntervalTree<GATKRead> tree = new SVIntervalTree<>();
        for (final GATKRead read : contigs) {
            if (read.isUnmapped()) continue;
            final int start = read.getStart();
            final int end = read.getEnd();
            final int contig = readMetadata.getContigID(read.getContig());
            tree.put(new SVInterval(contig, start, end), read);
        }
        return tree;
    }

    public static void writeTandemDuplicationEventsAsBedFile(final String outputPath, final Collection<SimpleSV> events, final ReadMetadata readMetadata) {
        try (final OutputStream outputStream = BucketUtils.createFile(outputPath)) {
            for (final SimpleSV event : events) {
                outputStream.write((event.getString(readMetadata) + "\n").getBytes());
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing output BED file", e);
        }
    }

    static Optional<SimpleSV> getHighestScoringEventFromStream(final Stream<SimpleSV> events, final double counterEvidencePseudocount) {
        return events.max(Comparator.comparingDouble(event -> event.getScore(counterEvidencePseudocount)));
    }

    public SVIntervalTree<EvidenceTargetLink> getIntrachromosomalLinkTree() {
        return intrachromosomalLinkTree;
    }

    public SVIntervalTree<EvidenceTargetLink> getInterchromosomalLinkTree() {
        return interchromosomalLinkTree;
    }

    private SVIntervalTree<EvidenceTargetLink> buildEvidenceIntervalTree(final List<EvidenceTargetLink> links, final int padding, final boolean separateLeftRightIntervals) {
        final SVIntervalTree<EvidenceTargetLink> linkTree = new SVIntervalTree<>();
        for (final EvidenceTargetLink link : links) {
            final SVInterval linkLeftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval linkRightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            if (separateLeftRightIntervals) {
                final SVInterval eventIntervalLeft = new SVInterval(linkLeftInterval.getContig(), linkLeftInterval.getStart(), linkLeftInterval.getEnd());
                final SVInterval paddedEventIntervalLeft = IntervalUtils.getPaddedInterval(eventIntervalLeft, padding, dict);
                linkTree.put(paddedEventIntervalLeft, link);
                final SVInterval eventIntervalRight = new SVInterval(linkRightInterval.getContig(), linkRightInterval.getStart(), linkRightInterval.getEnd());
                final SVInterval paddedEventIntervalRight = IntervalUtils.getPaddedInterval(eventIntervalRight, padding, dict);
                linkTree.put(paddedEventIntervalRight, link);
            } else if (linkLeftInterval.getContig() == linkRightInterval.getContig()) {
                final SVInterval eventInterval = new SVInterval(linkLeftInterval.getContig(), linkLeftInterval.getStart(), linkRightInterval.getEnd());
                final SVInterval paddedEventInterval = IntervalUtils.getPaddedInterval(eventInterval, padding, dict);
                linkTree.put(paddedEventInterval, link);
            }
        }
        return linkTree;
    }

    private Stream<SimpleSV> getCandidateEventsOnInterval(final SVInterval leftInterval, final SVInterval rightInterval, final SVInterval callInterval, final SVIntervalTree<SimpleSV> disallowedIntervals, final IntrachromosomalBreakpointPair breakpoints, final boolean applyFilter) {
        return getEventsOnInterval(leftInterval, rightInterval, callInterval, breakpoints, applyFilter).stream()
                .filter(event -> !IntervalUtils.hasReciprocalOverlapInTree(event.getInterval(), disallowedIntervals, arguments.MAX_CANDIDATE_EVENT_RECIPROCAL_OVERLAP));
    }

    Optional<SimpleSV> getHighestScoringEventOnInterval(final SVInterval leftInterval, final SVInterval rightInterval, final SVInterval callInterval, final SVIntervalTree<SimpleSV> disallowedIntervals, final IntrachromosomalBreakpointPair breakpoints, final boolean applyFilter) {
        return getHighestScoringEventFromStream(getCandidateEventsOnInterval(leftInterval, rightInterval, callInterval, disallowedIntervals, breakpoints, applyFilter), arguments.COUNTEREVIDENCE_PSEUDOCOUNT);
    }

    public List<SimpleSV> getEvents() {
        //Search breakpoint pairs
        final SVIntervalTree<SimpleSV> calledEventTree = new SVIntervalTree<>();
        for (final IntrachromosomalBreakpointPair breakpoints : pairedBreakpoints) {
            final SVInterval leftBreakpointInterval = new SVInterval(breakpoints.getContig(), breakpoints.getInterval().getStart(), breakpoints.getInterval().getStart());
            final SVInterval rightBreakpointInterval = new SVInterval(breakpoints.getContig(), breakpoints.getInterval().getEnd(), breakpoints.getInterval().getEnd());
            final SVInterval leftInterval = IntervalUtils.getPaddedInterval(leftBreakpointInterval, arguments.BREAKPOINT_INTERVAL_PADDING, dict);
            final SVInterval rightInterval = IntervalUtils.getPaddedInterval(rightBreakpointInterval, arguments.BREAKPOINT_INTERVAL_PADDING, dict);
            final Optional<SimpleSV> event = getHighestScoringEventOnInterval(leftInterval, rightInterval, breakpoints.getInterval(), calledEventTree, breakpoints, true);
            if (event.isPresent()) {
                calledEventTree.put(event.get().getInterval(), event.get());
            }
        }

        //Search links
        final Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> linkIter = intrachromosomalLinkTree.iterator();
        while (linkIter.hasNext()) {
            final EvidenceTargetLink link = linkIter.next().getValue();
            final SVInterval leftInterval = IntervalUtils.getPaddedInterval(link.getPairedStrandedIntervals().getLeft().getInterval(), arguments.EVIDENCE_LINK_INTERVAL_PADDING, dict);
            final SVInterval rightInterval = IntervalUtils.getPaddedInterval(link.getPairedStrandedIntervals().getRight().getInterval(), arguments.EVIDENCE_LINK_INTERVAL_PADDING, dict);
            final int callStart = (leftInterval.getStart() + leftInterval.getEnd()) / 2;
            final int callEnd = (rightInterval.getStart() + rightInterval.getEnd()) / 2;
            final SVInterval callInterval = new SVInterval(leftInterval.getContig(), callStart, callEnd);
            final Optional<SimpleSV> event = getHighestScoringEventOnInterval(leftInterval, rightInterval, callInterval, calledEventTree, null, true);
            if (event.isPresent()) {
                calledEventTree.put(event.get().getInterval(), event.get());
            }
        }

        //Search called segments
        final Iterator<CalledCopyRatioSegment> segmentIter = copyRatioSegmentOverlapDetector.getAll().iterator();
        while (segmentIter.hasNext()) {
            final CalledCopyRatioSegment segment = segmentIter.next();
            final int contigId = readMetadata.getContigID(segment.getContig());
            final SVInterval leftSegmentInterval = new SVInterval(contigId, segment.getInterval().getStart(), segment.getInterval().getStart());
            final SVInterval rightSegmentInterval = new SVInterval(contigId, segment.getInterval().getEnd(), segment.getInterval().getEnd());
            final SVInterval leftInterval = IntervalUtils.getPaddedInterval(leftSegmentInterval, arguments.COPY_RATIO_SEGMENT_PADDING, dict);
            final SVInterval rightInterval = IntervalUtils.getPaddedInterval(rightSegmentInterval, arguments.COPY_RATIO_SEGMENT_PADDING, dict);
            final int callStart = (leftInterval.getStart() + leftInterval.getEnd()) / 2;
            final int callEnd = (rightInterval.getStart() + rightInterval.getEnd()) / 2;
            final SVInterval callInterval = new SVInterval(contigId, callStart, callEnd);
            final Optional<SimpleSV> event = getHighestScoringEventOnInterval(leftInterval, rightInterval, callInterval, calledEventTree, null, true);
            if (event.isPresent()) {
                calledEventTree.put(event.get().getInterval(), event.get());
            }
        }

        return Utils.stream(calledEventTree.iterator()).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList());
    }

    List<SimpleSV> getEventsOnInterval(final SVInterval leftInterval,
                                       final SVInterval rightInterval,
                                       final SVInterval callInterval,
                                       IntrachromosomalBreakpointPair breakpoints,
                                       boolean applyFilter) {

        if (leftInterval.getContig() != rightInterval.getContig()) return Collections.emptyList();
        if (callInterval.getLength() < arguments.MIN_SV_SIZE || callInterval.getLength() > arguments.MAX_SV_SIZE)
            return Collections.emptyList();

        final List<SimpleSV> events = new ArrayList<>(3);

        final SimpleSV tandemDuplication = tandemDuplicationFactory.create(leftInterval, rightInterval, callInterval, applyFilter, breakpoints);
        if (tandemDuplication != null) events.add(tandemDuplication);

        final SimpleSV deletion = deletionFactory.create(leftInterval, rightInterval, callInterval, applyFilter, breakpoints);
        if (deletion != null) events.add(deletion);

        final SimpleSV inversion = inversionFactory.create(leftInterval, rightInterval, callInterval, applyFilter, breakpoints);
        if (inversion != null) events.add(inversion);

        return events;
    }

}
