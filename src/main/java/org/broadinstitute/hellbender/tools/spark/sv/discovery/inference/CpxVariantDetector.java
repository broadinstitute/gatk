package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * This deals with the special case where a contig has multiple (> 2) alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromMultipleAlignments()}.
 * For cases where the contig's alignment shows signature that the assembly doesn't
 * paint the full picture we could decide to emit all BND records,
 * but that could be dealt with later.
 */
public final class CpxVariantDetector {
    private static final boolean DEBUG_OUTPUT = false; // for debugging, once the prototyping code is merged, this and its usage will be deleted


    public void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                   final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final String outputPath = svDiscoveryInputData.outputPath;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        final JavaRDD<AnnotatedContig> annotatedContigs =
                assemblyContigs.map(tig -> new AnnotatedContig(tig, referenceSequenceDictionaryBroadcast.getValue()));

        SVVCFWriter.writeVCF(annotatedContigs.map(tig -> tig.toVariantContext(referenceBroadcast.getValue())).collect(),
                outputPath, referenceSequenceDictionaryBroadcast.getValue(), toolLogger);

        if (DEBUG_OUTPUT) {
            try {
                // for easier view when debugging, will be taken out in the final commit.
                Files.write(Paths.get(Paths.get(outputPath).getParent().toAbsolutePath().toString() + "/cpxEvents.txt"),
                        () -> annotatedContigs
                                .sortBy(tig -> tig.tigWithInsMappings.getSourceContig().contigName, true, 1)
                                .map(annotatedContig -> (CharSequence) annotatedContig.toString())
                                .collect().iterator());
            } catch (final IOException ioe) {
                throw new UserException.CouldNotCreateOutputFile("Could not save filtering results to file", ioe);
            }
        }
    }

    /**
     * Fundamental class of complex variant interpretation and alt haplotype extraction.
     *
     * <p>
     *     The general strategy is to walk on the affected reference region
     *     as guided by the alignments along the 5-3 direction of assembly contig.
     *     The walking has two basic types: sliding (alignment block),
     *     and jumping (two chimeric alignment joint on the read).
     *     By walking along the contig, we extract the jumps, then
     *     segment the affected region with segments bounded by the jumps,
     *     and walk along the alignments one more time to figure out
     *     how the segments are arranged on the sample/alt-haplotype.
     * </p>
     */
    @DefaultSerializer(AnnotatedContig.Serializer.class)
    private static final class AnnotatedContig {
        private static final AssemblyContigWithFineTunedAlignments.Serializer contigSerializer =
                new AssemblyContigWithFineTunedAlignments.Serializer();
        private static final BasicInfo.Serializer basicInfoSerializer = new BasicInfo.Serializer();
        private static final Jump.Serializer jumpSerializer = new Jump.Serializer();
        private static final ReferenceSegmentsAndEventDescription.Serializer eventsSerializer =
                new ReferenceSegmentsAndEventDescription.Serializer();

        final AssemblyContigWithFineTunedAlignments tigWithInsMappings;
        final BasicInfo basicInfo;
        private List<Jump> jumps;
        private List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations;
        private ReferenceSegmentsAndEventDescription referenceSegmentsAndEventDescription;
        private byte[] altSeq;

        /**
         * These getters may return null if the corresponding fields are not properly initialized yet.
         */
        List<Jump> getJumps() {
            return jumps;
        }
        List<SimpleInterval> getEventPrimaryChromosomeSegmentingLocations() { return eventPrimaryChromosomeSegmentingLocations;}
        ReferenceSegmentsAndEventDescription getSegmentsAndDescription() { return referenceSegmentsAndEventDescription; }
        byte[] getAltSeq() {
            return altSeq;
        }

        AnnotatedContig(final AssemblyContigWithFineTunedAlignments tigWithInsMappings, final SAMSequenceDictionary refSequenceDictionary) {

            final AlignedContig sourceTig = tigWithInsMappings.getSourceContig();
            final List<AlignmentInterval> deOverlappedAlignmentConfiguration =
                    deOverlapAlignments(sourceTig.alignmentIntervals, refSequenceDictionary);
            final AlignedContig contig = new AlignedContig(sourceTig.contigName, sourceTig.contigSequence,
                    deOverlappedAlignmentConfiguration, sourceTig.hasEquallyGoodAlnConfigurations);
            this.tigWithInsMappings = new AssemblyContigWithFineTunedAlignments(contig, tigWithInsMappings.getInsertionMappings());

            this.basicInfo = new BasicInfo(contig);

            annotate(refSequenceDictionary);
        }

        private static List<AlignmentInterval> deOverlapAlignments(final List<AlignmentInterval> originalAlignments,
                                                                   final SAMSequenceDictionary refSequenceDictionary) {
            final List<AlignmentInterval> result = new ArrayList<>(originalAlignments.size());
            final Iterator<AlignmentInterval> iterator = originalAlignments.iterator();
            AlignmentInterval one = iterator.next();
            while (iterator.hasNext()) {
                final AlignmentInterval two = iterator.next();
                // TODO: 11/5/17 an edge case is possible where the best configuration contains two alignments,
                //       one of which contains a large gap, and since the gap split happens after the configuration scoring,
                //       (that gap split happens after scoring is due to how MQ and AS are used in the scoring step, gap-split alignment cannot use originating alignment's values, but it takes time to recompute)
                //       one of the alignment from the gap split may be contained in the other original alignment, leading to problems;
                //       here we skip the alignment that is BEFORE the child alignment from the gap-split,
                //       IFF that alignment contains the child alignment in terms of their spans on the read/contig
                //       if you are concerned about the first child alignment from the same gapped alignment being skipped,
                //       don't worry, that is impossible because child alignments of the same gapped alignment cannot overlap on the read.
                if (two.alnModType.equals(ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT)) {
                    final int overlapOnRead = AlignmentInterval.overlapOnContig(one, two);
                    if (overlapOnRead >= two.getSizeOnRead())
                        continue;
                }
                final List<AlignmentInterval> deoverlapped = ContigAlignmentsModifier.removeOverlap(one, two, refSequenceDictionary);
                result.add(deoverlapped.get(0));
                one = deoverlapped.get(1);
            }
            result.add(one);
            return result;
        }

        /**
         * The process of annotation is the processing of making sense of what happened.
         *
         * <p>
         *     Segment affected reference region by jumping locations, by
         *     <ul>
         *         <li>
         *             extract jumping locations on reference from de-overlapped alignments
         *         </li>
         *         <li>
         *             extract the segmenting locations on event primary chromosome from jumps
         *         </li>
         *     </ul>
         * </p>
         *
         * <p>
         *     Then make sense of event, i.e. provide interpretation.
         *     And finally extract corresponding alt haplotype sequence
         * </p>
         */
        private void annotate(final SAMSequenceDictionary refSequenceDictionary) {
            try {
                jumps = extractJumpsOnReference(tigWithInsMappings.getSourceContig().alignmentIntervals);

                eventPrimaryChromosomeSegmentingLocations =
                        extractSegmentingRefLocationsOnEventPrimaryChromosome(jumps, basicInfo, refSequenceDictionary);

                referenceSegmentsAndEventDescription =
                        segmentReferenceAndInterpret(basicInfo, tigWithInsMappings.getSourceContig().alignmentIntervals, jumps,
                                eventPrimaryChromosomeSegmentingLocations);

                altSeq = extractAltHaplotypeSeq(tigWithInsMappings, referenceSegmentsAndEventDescription.referenceSegments, basicInfo);
            } catch (final GATKException | IllegalArgumentException likelyNewEdgeCase) {
                throw new GATKException(toString(), likelyNewEdgeCase);
            }
        }

        public VariantContext toVariantContext(final ReferenceMultiSource reference) throws IOException {

            final SimpleInterval affectedRefRegion = getAffectedRefRegion();
            final CpxVariantType cpxVariant = new CpxVariantType(affectedRefRegion, typeSpecificExtraAttributes());
            final SimpleInterval pos = new SimpleInterval(affectedRefRegion.getContig(), affectedRefRegion.getStart(), affectedRefRegion.getStart());

            final VariantContextBuilder vcBuilder = new VariantContextBuilder()
                    .chr(affectedRefRegion.getContig()).start(affectedRefRegion.getStart()).stop(affectedRefRegion.getEnd())
                    .alleles(AnnotatedVariantProducer.produceAlleles(pos, reference, cpxVariant))
                    .id(cpxVariant.getInternalVariantId())
                    .attribute(SVTYPE, cpxVariant.toString())
                    .attribute(VCFConstants.END_KEY, affectedRefRegion.getEnd())
                    .attribute(SVLEN, cpxVariant.getSVLength())
                    .attribute(SEQ_ALT_HAPLOTYPE, new String(altSeq));

            cpxVariant.getTypeSpecificAttributes().forEach(vcBuilder::attribute);

            // evidence used for producing the novel adjacency

            return vcBuilder.make();
        }

        private SimpleInterval getAffectedRefRegion() {
            final SimpleInterval pos, end;
            if (referenceSegmentsAndEventDescription.referenceSegments.isEmpty()) {
                if( eventPrimaryChromosomeSegmentingLocations.size() != 2 )
                    throw new GATKException("Inference seems to be wrong: " +
                            "affected reference region is not segmented but segmenting locations are not the expected two\n" +
                    toString());
                pos = eventPrimaryChromosomeSegmentingLocations.get(0);
                end = eventPrimaryChromosomeSegmentingLocations.get(1);
            } else {
                pos = referenceSegmentsAndEventDescription.referenceSegments.get(0);
                end = referenceSegmentsAndEventDescription.referenceSegments.get(referenceSegmentsAndEventDescription.referenceSegments.size() - 1);
            }
            return new SimpleInterval(pos.getContig(), pos.getStart(), end.getEnd());
        }

        private Map<String, String> typeSpecificExtraAttributes() {
            final Map<String, String> attributes = new HashMap<>();
            attributes.put(CPX_EVENT_ALT_ARRANGEMENTS,
                            referenceSegmentsAndEventDescription.getAltArrangementString());
            if ( ! referenceSegmentsAndEventDescription.referenceSegments.isEmpty() ) {
                attributes.put(CPX_SV_REF_SEGMENTS,
                        String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR,
                                referenceSegmentsAndEventDescription.referenceSegments.stream()
                                        .map(SimpleInterval::toString).collect(Collectors.toList())));
            }
            return attributes;
        }

        private Map<String, String> dummy() {

            final Map<String, String> attributeMap = new HashMap<>();
            attributeMap.put(GATKSVVCFConstants.TOTAL_MAPPINGS,    "1");
            attributeMap.put(GATKSVVCFConstants.HQ_MAPPINGS,       tigWithInsMappings.getSourceContig().alignmentIntervals.stream().mapToInt(ai -> ai.mapQual).anyMatch( mq -> mq < CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD ) ? "0" : "1");
            attributeMap.put(GATKSVVCFConstants.MAPPING_QUALITIES, tigWithInsMappings.getSourceContig().alignmentIntervals.stream().map(ai -> String.valueOf(ai.mapQual)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
            attributeMap.put(GATKSVVCFConstants.ALIGN_LENGTHS,     tigWithInsMappings.getSourceContig().alignmentIntervals.stream().map(ai -> String.valueOf(ai.getSizeOnRead())).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
            attributeMap.put(GATKSVVCFConstants.MAX_ALIGN_LENGTH,  String.valueOf(tigWithInsMappings.getSourceContig().alignmentIntervals.stream().mapToInt(AlignmentInterval::getSizeOnRead).max().orElse(0)));
            attributeMap.put(GATKSVVCFConstants.CONTIG_NAMES,      tigWithInsMappings.getSourceContig().contigName);

            // TODO: 12/11/17 integrate these with those that survived the alignment filtering step?
            // known insertion mappings from filtered out alignments
            if ( !tigWithInsMappings.getInsertionMappings().isEmpty() ) {
                attributeMap.put(INSERTED_SEQUENCE_MAPPINGS,
                        String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, tigWithInsMappings.getInsertionMappings()));
            }

            return attributeMap;
        }

        @Override
        public String toString() {
            return "\nContig:\t" + tigWithInsMappings.toString() +
                    "\nBasicInfo:\t" + basicInfo.toString() +
                    "\nJumps:\t" + ( jumps == null ? "NULL" : jumps.toString()) +
                    "\nSeg.Boundaries:\t" + ( eventPrimaryChromosomeSegmentingLocations == null ? "NULL" : eventPrimaryChromosomeSegmentingLocations.toString()) +
                    "\n" + (referenceSegmentsAndEventDescription == null ? "NULL" : referenceSegmentsAndEventDescription.toString()) +
                    "\nAlt.Seg.:\t" + (altSeq == null ? "NULL" : new String(altSeq));
        }

        AnnotatedContig(final Kryo kryo, final Input input) {
            tigWithInsMappings = contigSerializer.read(kryo, input, AssemblyContigWithFineTunedAlignments.class);
            basicInfo = basicInfoSerializer.read(kryo, input, BasicInfo.class);
            final int numJumps = input.readInt();
            jumps = new ArrayList<>(numJumps);
            for (int i = 0; i < numJumps; ++i) {
                jumps.add(jumpSerializer.read(kryo, input, Jump.class));
            }
            final int numSegmentingLocs = input.readInt();
            eventPrimaryChromosomeSegmentingLocations = new ArrayList<>(numJumps);
            for (int i = 0; i < numSegmentingLocs; ++i) {
                eventPrimaryChromosomeSegmentingLocations.add(kryo.readObject(input, SimpleInterval.class));
            }
            referenceSegmentsAndEventDescription = eventsSerializer.read(kryo, input, ReferenceSegmentsAndEventDescription.class);
            altSeq = new byte[input.readInt()];
            input.read(altSeq);
        }

        void serialize(final Kryo kryo, final Output output) {
            contigSerializer.write(kryo, output, tigWithInsMappings);
            basicInfoSerializer.write(kryo, output, basicInfo);
            output.writeInt(jumps.size());
            for (final Jump jump : jumps)
                jumpSerializer.write(kryo, output, jump);
            output.writeInt(eventPrimaryChromosomeSegmentingLocations.size());
            for (final SimpleInterval segmentingLoc : eventPrimaryChromosomeSegmentingLocations)
                kryo.writeObject(output, segmentingLoc);
            eventsSerializer.write(kryo, output, referenceSegmentsAndEventDescription);
            output.writeInt(altSeq.length);
            output.write(altSeq);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AnnotatedContig> {
            @Override
            public void write(final Kryo kryo, final Output output, final AnnotatedContig alignedContig) {
                alignedContig.serialize(kryo, output);
            }

            @Override
            public AnnotatedContig read(final Kryo kryo, final Input input, final Class<AnnotatedContig> clazz) {
                return new AnnotatedContig(kryo, input);
            }
        }
    }

    // =================================================================================================================

    @DefaultSerializer(BasicInfo.Serializer.class)
    private static final class BasicInfo {

        final String eventPrimaryChromosome; // where the head and tail alignments are mapped to, mappings to other chromosomes will be considered insertion/MEI
        final boolean forwardStrandRep;     // if the signaling assembly contig is a forward strand representation
        final SimpleInterval alpha;         // length-1 interval for the starting ref location of the head/tail alignment if the signaling assembly contig is a '+'/'-' representation
        final SimpleInterval omega;         // length-1 interval for the ending   ref location of the tail/head alignment if the signaling assembly contig is a '+'/'-' representation

        BasicInfo(final AlignedContig contig) {
            final AlignmentInterval head = contig.getHeadAlignment();
            final AlignmentInterval tail = contig.getTailAlignment();
            if (head == null || tail == null)
                throw new GATKException("Head or tail alignment is null from contig:\n" + contig.toString());

            eventPrimaryChromosome = head.referenceSpan.getContig();
            forwardStrandRep = head.forwardStrand;
            if (forwardStrandRep) {
                alpha = new SimpleInterval(head.referenceSpan.getContig(), head.referenceSpan.getStart(), head.referenceSpan.getStart());
                omega = new SimpleInterval(tail.referenceSpan.getContig(), tail.referenceSpan.getEnd(), tail.referenceSpan.getEnd());
            } else {
                alpha = new SimpleInterval(tail.referenceSpan.getContig(), tail.referenceSpan.getStart(), tail.referenceSpan.getStart());
                omega = new SimpleInterval(head.referenceSpan.getContig(), head.referenceSpan.getEnd(), head.referenceSpan.getEnd());
            }
        }

        SimpleInterval getRefRegionBoundedByAlphaAndOmega() {
            return new SimpleInterval(eventPrimaryChromosome, alpha.getStart(), omega.getEnd());
        }

        @Override
        public String toString() {
            return "primary chr: " + eventPrimaryChromosome + "\tstrand rep:" + (forwardStrandRep ? '+' : '-') +
                    "\talpha: " + alpha.toString() +  "\tomega: " + omega.toString();
        }

        BasicInfo(final Kryo kryo, final Input input) {
            eventPrimaryChromosome = input.readString();
            forwardStrandRep = input.readBoolean();
            alpha = kryo.readObject(input, SimpleInterval.class);
            omega = kryo.readObject(input, SimpleInterval.class);
        }

        void serialize(final Kryo kryo, final Output output) {
            output.writeString(eventPrimaryChromosome);
            output.writeBoolean(forwardStrandRep);
            kryo.writeObject(output, alpha);
            kryo.writeObject(output, omega);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<BasicInfo> {
            @Override
            public void write(final Kryo kryo, final Output output, final BasicInfo alignedContig) {
                alignedContig.serialize(kryo, output);
            }

            @Override
            public BasicInfo read(final Kryo kryo, final Input input, final Class<BasicInfo> clazz) {
                return new BasicInfo(kryo, input);
            }
        }
    }

    // =================================================================================================================

    /**
     * A jump has a starting and landing ref location.
     *
     * <p>
     * A jump can be:
     * <ul>
     *     <li>gapped--meaning a part of read is uncovered by neighboring alignments;</li>
     *     <li>connected--meaning neighboring alignments leave no base on the read uncovered, but does not overlap on the read;</li>
     *     <li>retracting--meaning neighboring alignments overlap on the read, pointing to homology between their ref span</li>
     * </ul>
     * Among them, retracting jumps are the most difficult to deal with, mainly due to how to have a consistent
     * homology-yielding scheme.
     *
     * Here we DO NOT output jumps that are retracting: i.e. we enforce the homology-yielding convention as implemented in
     * {@link ContigAlignmentsModifier#yieldHomologousSequenceToAlignmentTwo(AlignmentInterval, AlignmentInterval, SAMSequenceDictionary)}.
     *
     * </p>
     */
    @DefaultSerializer(Jump.Serializer.class)
    private static final class Jump {

        final SimpleInterval start;
        final SimpleInterval landing;
        final StrandSwitch strandSwitch;
        final int gapSize; // jump is gapped when this size is > 0

        Jump(final AlignmentInterval one, final AlignmentInterval two) {
            Utils.validateArg(AlignmentInterval.overlapOnContig(one, two) <=0,
                    "assumption that input alignments DO NOT overlap is violated.");

            strandSwitch = ChimericAlignment.determineStrandSwitch(one, two);

            switch (strandSwitch){
                case NO_SWITCH:
                    if (one.forwardStrand) {
                        start = new SimpleInterval(one.referenceSpan.getContig(), one.referenceSpan.getEnd(), one.referenceSpan.getEnd());
                        landing = new SimpleInterval(two.referenceSpan.getContig(), two.referenceSpan.getStart(), two.referenceSpan.getStart());
                    } else {
                        start = new SimpleInterval(one.referenceSpan.getContig(), one.referenceSpan.getStart(), one.referenceSpan.getStart());
                        landing = new SimpleInterval(two.referenceSpan.getContig(), two.referenceSpan.getEnd(), two.referenceSpan.getEnd());
                    }
                    break;
                case FORWARD_TO_REVERSE:
                    start = new SimpleInterval(one.referenceSpan.getContig(), one.referenceSpan.getEnd(), one.referenceSpan.getEnd());
                    landing = new SimpleInterval(two.referenceSpan.getContig(), two.referenceSpan.getEnd(), two.referenceSpan.getEnd());
                    break;
                case REVERSE_TO_FORWARD:
                    start = new SimpleInterval(one.referenceSpan.getContig(), one.referenceSpan.getStart(), one.referenceSpan.getStart());
                    landing = new SimpleInterval(two.referenceSpan.getContig(), two.referenceSpan.getStart(), two.referenceSpan.getStart());
                    break;
                default: throw new NoSuchElementException("seeing a strand switch that doesn't make sense");
            }

            gapSize = Math.max(0, two.startInAssembledContig - one.endInAssembledContig - 1);
        }

        boolean isGapped() {
            return gapSize > 0;
        }

        @Override
        public String toString() {
            return "Jump start: " + start.toString() + "\tjump landing: " + landing.toString() +
                    "\t" + strandSwitch.name() + "\t" + (gapSize > 0 ? "G" : "C");
        }

        Jump(final Kryo kryo, final Input input) {
            start = kryo.readObject(input, SimpleInterval.class);
            landing = kryo.readObject(input, SimpleInterval.class);
            strandSwitch = StrandSwitch.values()[input.readInt()];
            gapSize = input.readInt();
        }

        void serialize(final Kryo kryo, final Output output) {
            kryo.writeObject(output, start);
            kryo.writeObject(output, landing);
            output.writeInt(strandSwitch.ordinal());
            output.writeInt(gapSize);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<Jump> {
            @Override
            public void write(final Kryo kryo, final Output output, final Jump alignedContig) {
                alignedContig.serialize(kryo, output);
            }

            @Override
            public Jump read(final Kryo kryo, final Input input, final Class<Jump> clazz) {
                return new Jump(kryo, input);
            }
        }
    }

    /**
     * Each pair of neighboring reference locations are meant to be used closed, i.e. [a, b].
     * The returned list of {@link Jump}'s are ordered along the alignments of the contig.
     */
    private static List<Jump> extractJumpsOnReference(final List<AlignmentInterval> alignmentConfiguration) {

        final List<Jump> unsortedJumps = new ArrayList<>(alignmentConfiguration.size() - 1);
        final Iterator<AlignmentInterval> iterator = alignmentConfiguration.iterator();
        AlignmentInterval one = iterator.next();
        while(iterator.hasNext()) {
            final AlignmentInterval two = iterator.next();
            unsortedJumps.add( new Jump(one, two) );
            one = two;
        }

        return unsortedJumps;
    }

    /**
     * Given {@code jumps} extracted from the chimeric alignments,
     * filter out those starting or landing locations that are
     * disjoint from region {@link BasicInfo#getRefRegionBoundedByAlphaAndOmega()}.
     * What's left, and returned, will be used as boundaries for constructing segments.
     */
    private static List<SimpleInterval> extractSegmentingRefLocationsOnEventPrimaryChromosome(
            final List<Jump> jumps,
            final BasicInfo basicInfo,
            final SAMSequenceDictionary refSequenceDictionary) {
        final SimpleInterval regionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();
        return jumps.stream()
                .flatMap(jump -> Stream.of(jump.start, jump.landing))
                .filter(loc -> !alignmentIsDisjointFromAlphaOmega(loc, regionBoundedByAlphaAndOmega))
                .sorted((one, two) -> IntervalUtils.compareLocatables(one, two, refSequenceDictionary)) // alignments are sorted
                .distinct()
                .collect(Collectors.toList());
    }

    // =================================================================================================================

    /**
     * This struct contains two key pieces of information that provides interpretation of the event:
     * <p>
     *     Ordered list of reference segments on the event primary chromosome that
     *     are bounded by segmenting locations extracted above from {@link Jump}'s.
     *
     *     The order in which these segments are stored in the list
     *     is the same as how they are tiled on the reference, and
     *     two neighboring segments always share a boundary base,
     *     e.g. {(chr1, 1, 100), (chr1, 100, 150), (chr1, 150, 180), ...}
     *     (this shared base is because the segmenting locations are extracted from jumps,
     *      but didn't keep track of the order of the originating jumps.)
     * </p>
     *
     * <p>
     *     Description of how the reference segments are arranged on the sample, including orientation changes.
     *     Each description must be one of the three:
     *     <ul>
     *         <li>
     *             a signed integer, indicating which reference segment is placed,
     *             and if negative, meaning the segment's orientation is flipped on the sample.
     *             The absolute value is 1 + index into the above segments.
     *         </li>
     *         <li>
     *             a string that conforms to the format by {@link SimpleInterval#toString()}
     *             for describing a string of sequence that could map to that particular location,
     *             which is disjoint from the region returned by {@link BasicInfo#getRefRegionBoundedByAlphaAndOmega()}.
     *         </li>
     *         <li>
     *             a string literal of the form {@link #UNMAPPED_INSERTION%d},
     *             used for indicating a part of the assembly contig is uncovered by any (high quality) mappings,
     *             with %d indicating how many bases unmapped
     *         </li>
     *     </ul>
     * </p>
     */
    @DefaultSerializer(ReferenceSegmentsAndEventDescription.Serializer.class)
    private static final class ReferenceSegmentsAndEventDescription {

        // used for indicating a part of the assembly contig is uncovered by any (high quality) mappings.
        public static final String UNMAPPED_INSERTION = "UINS";

        final List<SimpleInterval> referenceSegments;
        final List<String> eventDescriptions;

        ReferenceSegmentsAndEventDescription(final List<SimpleInterval> referenceSegments, final List<String> eventDescriptions) {
            this.referenceSegments = referenceSegments;
            this.eventDescriptions = eventDescriptions;
        }

        String getAltArrangementString() {
            return String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, eventDescriptions);
        }

        @Override
        public String toString() {
            return "Segments:\t" + referenceSegments.toString() +
                    "\nArrangements:\t" + eventDescriptions.toString();
        }

        ReferenceSegmentsAndEventDescription(final Kryo kryo, final Input input) {
            final int numSegments = input.readInt();
            referenceSegments = new ArrayList<>(numSegments);
            for (int i = 0; i < numSegments; ++i)
                referenceSegments.add(kryo.readObject(input, SimpleInterval.class));
            final int numDescriptions = input.readInt();
            eventDescriptions = new ArrayList<>(numDescriptions);
            for (int i = 0; i < numDescriptions; ++i)
                eventDescriptions.add(input.readString());
        }

        void serialize(final Kryo kryo, final Output output) {
            output.writeInt(referenceSegments.size());
            for (final SimpleInterval segment : referenceSegments)
                kryo.writeObject(output, segment);
            output.writeInt(eventDescriptions.size());
            for (final String description: eventDescriptions)
                output.writeString(description);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ReferenceSegmentsAndEventDescription> {
            @Override
            public void write(final Kryo kryo, final Output output, final ReferenceSegmentsAndEventDescription x) {
                x.serialize(kryo, output);
            }

            @Override
            public ReferenceSegmentsAndEventDescription read(final Kryo kryo, final Input input,
                                                             final Class<ReferenceSegmentsAndEventDescription> clazz) {
                return new ReferenceSegmentsAndEventDescription(kryo, input);
            }
        }
    }

    private static ReferenceSegmentsAndEventDescription segmentReferenceAndInterpret(final BasicInfo basicInfo,
                                                                                     final List<AlignmentInterval> contigAlignments,
                                                                                     final List<Jump> jumps,
                                                                                     final List<SimpleInterval> segmentingLocations) {

        if (segmentingLocations.size() == 1) {
            return onlyHeadAndTailShareOneBaseOnRef(basicInfo, contigAlignments, segmentingLocations.get(0));
        } else {
            return multipleAlignmentsInAlphaAndOmegaRegion(basicInfo, contigAlignments, jumps, segmentingLocations);
        }
    }

    /**
     * This is the case where a contig has all of the middle alignment(s) mapped to some disjoint places
     * (different chromosome, or same chromosome but not in the region returned by
     * {@link BasicInfo#getRefRegionBoundedByAlphaAndOmega()}, aka valid region),
     * AND the first and last alignment share a single base on their ref span
     * (hence we have only one jump location on the event primary chromosome).
     * Hence the {@code segmentingLocations} is a one single-base entry list.
     */
    private static ReferenceSegmentsAndEventDescription onlyHeadAndTailShareOneBaseOnRef(final BasicInfo basicInfo,
                                                                                         final List<AlignmentInterval> contigAlignments,
                                                                                         final SimpleInterval singleBase) {
        if ( singleBase.size() != 1)
            throw new UnhandledCaseSeen(
                    "run into unseen case where only one reference segmenting location is found but its size is not 1:\n"
                            + contigAlignments.toString());

        final SimpleInterval refRegionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();
        final boolean allMiddleAlignmentsDisjointFromAlphaOmega =
                contigAlignments
                        .subList(1, contigAlignments.size() - 1).stream()
                        .allMatch(ai -> alignmentIsDisjointFromAlphaOmega(ai.referenceSpan, refRegionBoundedByAlphaAndOmega));
        if ( ! allMiddleAlignmentsDisjointFromAlphaOmega )
            throw new UnhandledCaseSeen("run into unseen case where only one reference segmenting location is found" +
                    " but some middle alignments are overlapping alpha-omega region:\t" + refRegionBoundedByAlphaAndOmega + "\n"
                    + contigAlignments.toString());


        final List<SimpleInterval> segments = new ArrayList<>(Collections.singletonList(singleBase));

        final List<String> description =
                contigAlignments
                        .subList(1, contigAlignments.size() - 1).stream()
                        .map(ai -> {
                            if (basicInfo.forwardStrandRep)
                                return ai.forwardStrand ? ai.referenceSpan.toString() : "-"+ai.referenceSpan.toString();
                            else
                                return ai.forwardStrand ? "-"+ai.referenceSpan.toString() : ai.referenceSpan.toString();
                        }).collect(Collectors.toList());

        // the single base overlap from head and tail
        description.add(0, "1");
        description.add("1");
        return new ReferenceSegmentsAndEventDescription(segments, description);

    }

    /**
     * For dealing with case where at least one middle alignments of the assembly contig are overlapping with
     * the region returned by {@link BasicInfo#getRefRegionBoundedByAlphaAndOmega()}.
     * This is contrary to the case handled in {@link #onlyHeadAndTailShareOneBaseOnRef(BasicInfo, List, SimpleInterval)}.
     */
    private static ReferenceSegmentsAndEventDescription multipleAlignmentsInAlphaAndOmegaRegion(final BasicInfo basicInfo,
                                                                                                final List<AlignmentInterval> contigAlignments,
                                                                                                final List<Jump> jumps,
                                                                                                final List<SimpleInterval> segmentingLocations) {

        final List<SimpleInterval> segments = extractRefSegments(basicInfo, segmentingLocations);
        final List<String> descriptions = makeInterpretation(basicInfo, contigAlignments, jumps, segments);

        return new ReferenceSegmentsAndEventDescription(segments, descriptions);
    }

    private static List<SimpleInterval> extractRefSegments(final BasicInfo basicInfo,
                                                           final List<SimpleInterval> segmentingLocations) {

        final String eventPrimaryChromosome = basicInfo.eventPrimaryChromosome;

        final List<SimpleInterval> segments = new ArrayList<>(segmentingLocations.size() - 1);
        final Iterator<SimpleInterval> iterator = segmentingLocations.iterator();
        SimpleInterval leftBoundary = iterator.next();
        while (iterator.hasNext()) {
            final SimpleInterval rightBoundary = iterator.next();
            // there shouldn't be a segment constructed if two segmenting locations are adjacent to each other on the reference
            // this could happen when (in the simplest case), two alignments are separated by a mapped insertion (hence 3 total alignments),
            // and the two alignments' ref span are connected
            if (rightBoundary.getStart() - leftBoundary.getEnd() > 1) {
                segments.add(new SimpleInterval(eventPrimaryChromosome, leftBoundary.getStart(), rightBoundary.getStart()));
            }
            leftBoundary = rightBoundary;
        }
        return segments;
    }

    private static List<String> makeInterpretation(final BasicInfo basicInfo,
                                                   final List<AlignmentInterval> contigAlignments,
                                                   final List<Jump> jumps,
                                                   final List<SimpleInterval> segments) {

        // using overlap with alignments ordered along the '+' strand representation of
        // the signaling contig to make sense of how the reference segments are ordered,
        // including orientations--using signs of the integers, on the sample;
        final List<AlignmentInterval> alignmentIntervalList =
                basicInfo.forwardStrandRep ? contigAlignments
                                           : ImmutableList.copyOf(contigAlignments).reverse();
        final Iterator<Integer> jumpGapSizeIterator =
                basicInfo.forwardStrandRep ? jumps.stream().map(jump ->jump.gapSize).iterator()
                                           : ImmutableList.copyOf(jumps).reverse().stream().map(jump ->jump.gapSize).iterator();
        Integer currentJumpGapSize = jumpGapSizeIterator.next();
        boolean jumpIsLast = false;

        final SimpleInterval regionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();

        final List<String> descriptions = new ArrayList<>(2*segments.size()); //ini. cap. a guess
        final List<Tuple2<SimpleInterval, Integer>> insertionsMappedToDisjointRegionsAndInsertionLocations = new ArrayList<>();
        for (final AlignmentInterval alignment : alignmentIntervalList) {

            if (alignmentIsDisjointFromAlphaOmega(alignment.referenceSpan, regionBoundedByAlphaAndOmega)) {
                // disjoint alignment won't overlap any segments, so note down once where to insert, then move to next alignment
                final int indexABS = descriptions.size();
                insertionsMappedToDisjointRegionsAndInsertionLocations.add(new Tuple2<>(alignment.referenceSpan,
                        basicInfo.forwardStrandRep == alignment.forwardStrand ? indexABS : -1*indexABS));
                if ( currentJumpGapSize > 0 )
                    descriptions.add(ReferenceSegmentsAndEventDescription.UNMAPPED_INSERTION + "-" + currentJumpGapSize);
            } else {
                // depending on the representation and the current alignment's orientation, traverse segments in different order
                final int start, stop, step;
                if (basicInfo.forwardStrandRep == alignment.forwardStrand) {
                    start = 0;
                    stop = segments.size();
                    step = 1;
                } else {
                    start = segments.size()-1;
                    stop = -1;
                    step = -1;
                }
                // N*M overlaps
                for ( int i = start; i != stop; i += step ) {
                    final SimpleInterval currentSegment = segments.get(i);
                    // if current segment is contained in current alignment, note it down
                    if ( alignment.referenceSpan.contains(currentSegment) ) {
                        if (basicInfo.forwardStrandRep) // +1 below on i for 1-based description, no magic
                            descriptions.add( String.valueOf((alignment.forwardStrand ? 1 : -1) * (i+1)) );
                        else
                            descriptions.add( String.valueOf((alignment.forwardStrand ? -1 : 1) * (i+1)) );
                    }

                }
                // if the current alignment is associated with a gapped jump,
                // we need to signal that an unmapped insertion is present,
                // but only under two cases:
                //  1) no more segments to explore
                //  2) the next segment IS NOT contained in the current alignment's ref span
                if ( currentJumpGapSize > 0 && !jumpIsLast){
                    descriptions.add(ReferenceSegmentsAndEventDescription.UNMAPPED_INSERTION + "-" + currentJumpGapSize);
                }
            }

            if (jumpGapSizeIterator.hasNext()) // last alignment has no leaving jump so need to guard against that
                currentJumpGapSize = jumpGapSizeIterator.next();
            else
                jumpIsLast = true;
        }

        // post-hoc treatment of insertions that map to disjoint regions
        // go in reverse order because inserting into "descriptions", and insertion invalidates indices/iterators
        final ImmutableList<Tuple2<SimpleInterval, Integer>> reverse =
                ImmutableList.copyOf(insertionsMappedToDisjointRegionsAndInsertionLocations).reverse();
        for (final Tuple2<SimpleInterval, Integer> pair : reverse){
            final int index = pair._2;
            if (index > 0) {
                descriptions.add(pair._2, pair._1.toString());
            } else {
                descriptions.add(-1*pair._2, "-"+pair._1.toString());
            }
        }
        return descriptions;
    }

    private static boolean alignmentIsDisjointFromAlphaOmega(final SimpleInterval alignmentRefSpan,
                                                             final SimpleInterval regionBoundedByAlphaAndOmega) {
        return !alignmentRefSpan.overlaps(regionBoundedByAlphaAndOmega);
    }

    // =================================================================================================================

    /**
     * Extract alt haplotype sequence from the {@code tigWithInsMappings} to accompany the interpreted events.
     */
    private static byte[] extractAltHaplotypeSeq(final AssemblyContigWithFineTunedAlignments tigWithInsMappings,
                                                 final List<SimpleInterval> segments,
                                                 final BasicInfo basicInfo) {

        final AlignmentInterval head = tigWithInsMappings.getSourceContig().getHeadAlignment();
        final AlignmentInterval tail = tigWithInsMappings.getSourceContig().getTailAlignment();
        if (head == null || tail == null)
            throw new GATKException("Head or tail alignment is null from contig:\n" + tigWithInsMappings.getSourceContig().toString());

        if (segments.isEmpty()) { // case where middle alignments all map to disjoint locations
            final int start = head.endInAssembledContig;
            final int end = tail.startInAssembledContig;
            final byte[] altSeq = Arrays.copyOfRange(tigWithInsMappings.getSourceContig().contigSequence, start - 1, end);
            if ( ! basicInfo.forwardStrandRep ) {
                SequenceUtil.reverseComplement(altSeq);
            }
            return altSeq;
        }

        final SimpleInterval firstSegment, lastSegment;
        if ( basicInfo.forwardStrandRep ) {
            firstSegment = segments.get(0);
            lastSegment = segments.get(segments.size() - 1);
        } else {
            firstSegment = segments.get(segments.size() - 1);
            lastSegment = segments.get(0);
        }

        final int start, end;

        if ( !firstSegment.overlaps(head.referenceSpan) ) {
            // if first segment doesn't overlap with head alignment,
            // it must be the case that the base (and possibly following bases) immediately after (or before if reverse strand) the head alignment's ref span is deleted
            final boolean firstSegmentNeighborsHeadAlignment = basicInfo.forwardStrandRep ? (firstSegment.getStart() - head.referenceSpan.getEnd() == 1)
                                                                                          : (head.referenceSpan.getStart() - firstSegment.getEnd() == 1);
            if ( ! firstSegmentNeighborsHeadAlignment )
                throw new UnhandledCaseSeen("1st segment is not overlapping with head alignment but it is not immediately before/after the head alignment either\n"
                        + tigWithInsMappings.toString());
            start = head.endInAssembledContig;
        } else {
            final SimpleInterval intersect = firstSegment.intersect(head.referenceSpan);
            if (intersect.size() == 1) {
                start = head.endInAssembledContig;
            } else {
                start = head.readIntervalAlignedToRefSpan(intersect)._1;
            }
        }

        if ( !lastSegment.overlaps(tail.referenceSpan) ) {
            final boolean expectedCase = basicInfo.forwardStrandRep ? (tail.referenceSpan.getStart() - lastSegment.getEnd() == 1)
                                                                    : (lastSegment.getStart() - tail.referenceSpan.getEnd() == 1);
            if ( ! expectedCase )
                throw new UnhandledCaseSeen(tigWithInsMappings.toString());
            end = tail.startInAssembledContig;
        } else {
            final SimpleInterval intersect = lastSegment.intersect(tail.referenceSpan);
            if (intersect.size() == 1) {
                end = tail.startInAssembledContig;
            } else {
                end = tail.readIntervalAlignedToRefSpan(intersect)._2;
            }
        }

        // note from 1-based inclusive coordinate to C-style coordinate
        final byte[] altSeq = Arrays.copyOfRange(tigWithInsMappings.getSourceContig().contigSequence, start - 1, end);
        if ( ! basicInfo.forwardStrandRep ) {
            SequenceUtil.reverseComplement(altSeq);
        }

        return altSeq;
    }

    // =================================================================================================================


    private static final class UnhandledCaseSeen extends GATKException.ShouldNeverReachHereException {
        private static final long serialVersionUID = 0L;
        UnhandledCaseSeen( final String s ) {
            super(s);
        }

        UnhandledCaseSeen( final String s, final Throwable throwable ) {
            super(s, throwable);
        }

        UnhandledCaseSeen( final Throwable throwable) {this("Seeing unhandled case", throwable);}
    }
}
