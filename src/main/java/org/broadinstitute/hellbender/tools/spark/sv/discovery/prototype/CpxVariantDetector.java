package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.StrandSwitch;
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

/**
 * This deals with the special case where a contig's multiple (> 2) alignments has head and tail mapped to the same chr.
 * For the case where the head and tail mapped to different chromosome, we could decide to emit all BND records, but
 * that could be dealt with later.
 */
final class CpxVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {


    @Override
    public void inferSvAndWriteVCF(final String vcfOutputFileName, final String sampleId, final JavaRDD<AlignedContig> localAssemblyContigs,
                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                   final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                   final Logger toolLogger) {

        final JavaRDD<AnnotatedContig> annotatedContigs =
                localAssemblyContigs.map(tig -> new AnnotatedContig(tig, broadcastSequenceDictionary.getValue()));

        // extract corresponding alt haplotype sequence
        final JavaPairRDD<String, byte[]> contigNameAndAltSeq =
                annotatedContigs.mapToPair(annotatedContig -> new Tuple2<>(annotatedContig.contig.contigName, annotatedContig.extractAltHaplotypeSeq()));
        // output VCF

        // debug
        try {
            Files.write(Paths.get(Paths.get(vcfOutputFileName).getParent().toAbsolutePath().toString() + "/cpxEvents.txt"),
                    () -> annotatedContigs
                            .mapToPair(annotatedContig -> new Tuple2<>(annotatedContig.contig.contigName, annotatedContig))
                            .join(contigNameAndAltSeq)
                            .sortByKey()
                            .values()
                            .map(pair -> pair._1.toString() + "\nAlt. Hap.Seq.:\t" + new String(pair._2))
                            .map(s -> (CharSequence) s)
                            .collect().iterator());
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile("Could not save filtering results to file", ioe);
        }
    }

    @DefaultSerializer(AnnotatedContig.Serializer.class)
    private static final class AnnotatedContig {
        private static final AlignedContig.Serializer contigSerializer = new AlignedContig.Serializer();
        private static final BasicInfo.Serializer basicInfoSerializer = new BasicInfo.Serializer();
        private static final Jump.Serializer jumpSerializer = new Jump.Serializer();
        private static final ReferenceSegmentsAndEventDescription.Serializer eventsSerializer = new ReferenceSegmentsAndEventDescription.Serializer();
        final AlignedContig contig;
        final BasicInfo basicInfo;
        private List<Jump> jumps;
        private List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations;
        private ReferenceSegmentsAndEventDescription referenceSegmentsAndEventDescription;

        /**
         * These two getters may return null if the corresponding fields are not properly initialized yet.
         */
        List<Jump> getJumps() {
            return jumps;
        }
        List<SimpleInterval> getEventPrimaryChromosomeSegmentingLocations() { return eventPrimaryChromosomeSegmentingLocations;}
        ReferenceSegmentsAndEventDescription getSegmentsAndDescription() { return referenceSegmentsAndEventDescription; }

        /**
         * determine event primary chromosome and strand representation
         * add jumping locations on reference annotation then extract the segmenting locations on event primary chromosome
         * segment affected reference regions by jumping locations
         * then make sense of event, i.e. provide interpretation
         * finally extract corresponding alt haplotype sequence
         */
        AnnotatedContig(final AlignedContig contig, final SAMSequenceDictionary refSequenceDictionary) {
            this.contig = new AlignedContig(contig.contigName, contig.contigSequence,
                    deOverlapAlignments(contig.alignmentIntervals, refSequenceDictionary), contig.hasEquallyGoodAlnConfigurations);;
            this.basicInfo = new BasicInfo(this.contig);

            annotate(refSequenceDictionary);
        }

        byte[] extractAltHaplotypeSeq() {
            final Tuple2<Integer, Integer> alphaAndOmegaPosOnRead = basicInfo.getAlphaAndOmegaPosOnRead();
            final byte[] altSeq = Arrays.copyOfRange(contig.contigSequence, alphaAndOmegaPosOnRead._1 - 1, alphaAndOmegaPosOnRead._2);
            if ( ! basicInfo.forwardStrandRep )
                SequenceUtil.reverseComplement(altSeq);
            return altSeq;
        }


        private static List<AlignmentInterval> deOverlapAlignments(final List<AlignmentInterval> originalAlignments,
                                                                   final SAMSequenceDictionary refSequenceDictionary) {
            final List<AlignmentInterval> result = new ArrayList<>(originalAlignments.size());
            final Iterator<AlignmentInterval> iterator = originalAlignments.iterator();
            AlignmentInterval one = iterator.next();
            while (iterator.hasNext()) {
                final AlignmentInterval two = iterator.next();
                // TODO: 11/5/17 an edge case is possible where the best configuration contains two alignments, one of which contains a large gap, and since the gap split happens after the configuration scoring, one of the alignment from the gap split may be contained in the other original alignment, leading to problems; here we first skip it
                if (two.alnModType.equals(AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT)) {
                    final int overlapOnRead = AlignmentInterval.overlapOnContig(one, two);
                    if (overlapOnRead >= two.getSpanOnRead())
                        continue;
                }
                final List<AlignmentInterval> deoverlapped = ContigAlignmentsModifier.removeOverlap(one, two, refSequenceDictionary);
                result.add(deoverlapped.get(0));
                one = deoverlapped.get(1);
            }
            result.add(one);
            return result;
        }

        private void annotate(final SAMSequenceDictionary refSequenceDictionary) {
            try {
                jumps = extractJumpsOnReference(contig.alignmentIntervals);

                eventPrimaryChromosomeSegmentingLocations =
                        extractSegmentingRefLocationsOnEventPrimaryChromosome(jumps, basicInfo, refSequenceDictionary);

                referenceSegmentsAndEventDescription =
                        segmentReferenceAndInterpret(basicInfo, contig.alignmentIntervals, jumps,
                                eventPrimaryChromosomeSegmentingLocations);
            } catch (final GATKException ex) {
                throw new GATKException(toString(), ex);
            }
        }

        @Override
        public String toString() {
            return "\nContig:\t" + contig.toString() +
                    "\n" + basicInfo.toString() +
                    "\nJumps:\t" + ( jumps == null ? "NULL" : jumps.toString()) +
                    "\nSeg.Boundaries:\t" + ( eventPrimaryChromosomeSegmentingLocations == null ? "NULL" : eventPrimaryChromosomeSegmentingLocations.toString()) +
                    "\n" + (referenceSegmentsAndEventDescription == null ? "NULL" : referenceSegmentsAndEventDescription.toString());
        }

        AnnotatedContig(final Kryo kryo, final Input input) {
            contig = contigSerializer.read(kryo, input, AlignedContig.class);
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
        }

        void serialize(final Kryo kryo, final Output output) {
            contigSerializer.write(kryo, output, contig);
            basicInfoSerializer.write(kryo, output, basicInfo);
            output.writeInt(jumps.size());
            for (final Jump jump : jumps)
                jumpSerializer.write(kryo, output, jump);
            output.writeInt(eventPrimaryChromosomeSegmentingLocations.size());
            for (final SimpleInterval segmentingLoc : eventPrimaryChromosomeSegmentingLocations)
                kryo.writeObject(output, segmentingLoc);
            eventsSerializer.write(kryo, output, referenceSegmentsAndEventDescription);
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

        final String eventPrimaryChromosome; // where the head and tail alignments are mapped to, mappings to other chromosomes of a particular assembly contig will be considered insertion/MEI
        final boolean forwardStrandRep;     // if the signaling assembly contig is a forward strand representation
        final SimpleInterval alpha;         // the starting ref location of the head/tail alignment if the signaling assembly contig is a '+'/'-' representation
        final SimpleInterval omega;         // the ending   ref location of the tail/head alignment if the signaling assembly contig is a '+'/'-' representation
        final int alphaOnRead;              // position on read corresponding to alpha
        final int omegaOnRead;              // position on read corresponding to omega

        BasicInfo(final AlignedContig contig) {
            final AlignmentInterval head = contig.alignmentIntervals.get(0),
                                    tail = contig.alignmentIntervals.get(contig.alignmentIntervals.size()-1);

            eventPrimaryChromosome = head.referenceSpan.getContig();
            forwardStrandRep = head.forwardStrand;
            if (forwardStrandRep) {
                alpha = new SimpleInterval(head.referenceSpan.getContig(), head.referenceSpan.getStart(), head.referenceSpan.getStart());
                omega = new SimpleInterval(tail.referenceSpan.getContig(), tail.referenceSpan.getEnd(), tail.referenceSpan.getEnd());
            } else {
                alpha = new SimpleInterval(tail.referenceSpan.getContig(), tail.referenceSpan.getStart(), tail.referenceSpan.getStart());
                omega = new SimpleInterval(head.referenceSpan.getContig(), head.referenceSpan.getEnd(), head.referenceSpan.getEnd());
            }

            alphaOnRead = head.startInAssembledContig;
            omegaOnRead = tail.endInAssembledContig;
        }

        SimpleInterval getRefRegionBoundedByAlphaAndOmega() {
            return new SimpleInterval(eventPrimaryChromosome, alpha.getStart(), omega.getEnd());
        }

        Tuple2<Integer, Integer> getAlphaAndOmegaPosOnRead() {
            return new Tuple2<>(alphaOnRead, omegaOnRead);
        }

        @Override
        public String toString() {
            return "BasicInfo:\tprimary chr: " + eventPrimaryChromosome + "\tstrand rep:" + (forwardStrandRep ? '+' : '-') +
                    "\talpha: " + alpha.toString() +  "\tomega: " + omega.toString() +
                    "\talphaOnRead: " + alphaOnRead + "\tomegaOnRead: " + omegaOnRead;
        }

        BasicInfo(final Kryo kryo, final Input input) {
            eventPrimaryChromosome = input.readString();
            forwardStrandRep = input.readBoolean();
            alpha = kryo.readObject(input, SimpleInterval.class);
            omega = kryo.readObject(input, SimpleInterval.class);
            alphaOnRead = input.readInt();
            omegaOnRead = input.readInt();
        }

        void serialize(final Kryo kryo, final Output output) {
            output.writeString(eventPrimaryChromosome);
            output.writeBoolean(forwardStrandRep);
            kryo.writeObject(output, alpha);
            kryo.writeObject(output, omega);
            output.writeInt(alphaOnRead);
            output.writeInt(omegaOnRead);
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
     *     <li>gapped--meaning a part of read is uncovered by neighboring AI's;</li>
     *     <li>connected--meaning neighboring--but not overlapping on the read--AI's leave no base on the read uncovered;</li>
     *     <li>retracting--meaning neighboring AI's overlap on the read, pointing to homology between their ref span</li>
     * </ul>
     * Among them, retracting jumps are the most difficult to deal with, mainly due to how to have a consistent
     * homology-yielding scheme.
     *
     * Here we DO NOT output jumps that are retracting: i.e. we enforce the homology-yielding convention as implemented in
     * {@link ContigAlignmentsModifier#homologyYieldingStrategy(AlignmentInterval, AlignmentInterval, SAMSequenceDictionary)}.
     *
     * </p>
     */
    @DefaultSerializer(Jump.Serializer.class)
    private static final class Jump {

        final SimpleInterval start;
        final SimpleInterval landing;
        final StrandSwitch strandSwitch;
        final boolean gapped;

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

            gapped = two.startInAssembledContig > one.endInAssembledContig + 1;

        }

        @Override
        public String toString() {
            return "Jump start: " + start.toString() + "\tjump landing: " + landing.toString() +
                    "\t" + strandSwitch.name() + "\t" + (gapped ? "G" : "C");
        }

        Jump(final Kryo kryo, final Input input) {
            start = kryo.readObject(input, SimpleInterval.class);
            landing = kryo.readObject(input, SimpleInterval.class);
            strandSwitch = StrandSwitch.values()[input.readInt()];
            gapped = input.readBoolean();
        }

        void serialize(final Kryo kryo, final Output output) {
            kryo.writeObject(output, start);
            kryo.writeObject(output, landing);
            output.writeInt(strandSwitch.ordinal());
            output.writeBoolean(gapped);
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
     * The returned list of {@link Jump}'s are ordered along the alignments of the contig
     * if the contig is '+' strand representation, or reversed if the contig is '-' strand representation.
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

    private static List<SimpleInterval> extractSegmentingRefLocationsOnEventPrimaryChromosome(
            final List<Jump> jumps,
            final BasicInfo basicInfo,
            final SAMSequenceDictionary refSequenceDictionary) {
        final SimpleInterval regionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();
        return jumps.stream()
                .flatMap(jump -> Stream.of(jump.start, jump.landing))
                .filter(loc -> !alignmentIsDisjointFromAlphaOmega(loc, regionBoundedByAlphaAndOmega))
                .distinct()
                .sorted((one, two) -> IntervalUtils.compareLocatables(one, two, refSequenceDictionary))
                .collect(Collectors.toList());
    }

    // =================================================================================================================

    /**
     * This struct contains two pieces of information that provides interpretation of the event:
     * <p>
     *     Ordered list of reference segments on the event primary chromosome that are bounded
     *     by segmenting locations extracted above from {@link Jump}'s.
     *     The segments are always ordered in contiguously increasing order (e.g. 1, 2, 3, ...).
     *     Note that the resulting segments have a property that neighboring segments always share a boundary base.
     *     This is because the segmenting locations are extracted from jumps, but didn't keep track of the order of
     *     the originating jumps.
     * </p>
     *
     * <p>
     *     Description of how the reference segments are arranged on the sample, including orientation changes.
     *     Each description is either
     *     <ul>
     *         <li>
     *             a signed integer, indicating which reference segment is placed,
     *             and if negative, meaning the segment's orientation is flipped.
     *             The absolute value in the list is 1 + index into the above segments.
     *         </li>
     *         <li>
     *             a string that conforms to the format by {@link SimpleInterval#toString()}
     *             for describing a string of sequence that could map to that particular location,
     *             which is disjoint to the region bounded by {@link BasicInfo#alpha} and {@link BasicInfo#omega}.
     *         </li>
     *         <li>
     *             {@link #UNMAPPED_INSERTION},
     *             used for indicating a part of the assembly contig is uncovered by any mappings.
     *         </li>
     *     </ul>
     * </p>
     */
    @DefaultSerializer(ReferenceSegmentsAndEventDescription.Serializer.class)
    private static final class ReferenceSegmentsAndEventDescription {

        // used for indicating a part of the assembly contig is uncovered by any mappings.
        public static final String UNMAPPED_INSERTION = "UMIS";

        final List<SimpleInterval> referenceSegments;
        final List<String> eventDescriptions;

        ReferenceSegmentsAndEventDescription(final List<SimpleInterval> referenceSegments, final List<String> eventDescriptions) {
            this.referenceSegments = referenceSegments;
            this.eventDescriptions = eventDescriptions;
        }

        @Override
        public String toString() {
            return "Segments:\t" + referenceSegments.toString() +
                    "\nEvents:\t" + eventDescriptions.toString();
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

        final String eventPrimaryChromosome = basicInfo.eventPrimaryChromosome;

        if (segmentingLocations.size() == 1) {
            final boolean caseSeenBefore =
                    contigAlignments.stream()
                            .filter(ai -> ai.referenceSpan.getContig().equals(eventPrimaryChromosome)).count() == 2;
            if (caseSeenBefore) {
                return livingAbroad(basicInfo, contigAlignments, jumps, segmentingLocations);
            } else {
                throw new UnhandledCaseSeen("run into unseen case:\n");
            }
        } else {
            return homelandMover(basicInfo, contigAlignments, jumps, segmentingLocations);
        }
    }

    /**
     * This is the case where a contig has the middle alignment(s) mapped to some disjoint places
     * (different chromosome, or same chromosome but not in the region bounded by alpha and omega),
     * and the first and last alignment share a single base on their ref span
     * (hence we have only one jump location on the event primary chromosome)
     */
    private static ReferenceSegmentsAndEventDescription livingAbroad(final BasicInfo basicInfo,
                                                                     final List<AlignmentInterval> contigAlignments,
                                                                     final List<Jump> jumps,
                                                                     final List<SimpleInterval> segmentingLocations) {


        final SimpleInterval singleBase = segmentingLocations.get(0);
        final List<SimpleInterval> segments = new ArrayList<>(Arrays.asList(singleBase, singleBase));

        final List<String> description =
                contigAlignments
                        .subList(1, contigAlignments.size() - 1).stream()
                        .map(ai -> {
                            if (basicInfo.forwardStrandRep)
                                return ai.forwardStrand ? ai.referenceSpan.toString() : "-"+ai.referenceSpan.toString();
                            else
                                return ai.forwardStrand ? "-"+ai.referenceSpan.toString() : ai.referenceSpan.toString();
                        }).collect(Collectors.toList());

        description.add(0, "1");
        description.add("1");
        return new ReferenceSegmentsAndEventDescription(segments, description);

    }

    private static ReferenceSegmentsAndEventDescription homelandMover(final BasicInfo basicInfo,
                                                                      final List<AlignmentInterval> contigAlignments,
                                                                      final List<Jump> jumps,
                                                                      final List<SimpleInterval> segmentingLocations) {

        final List<SimpleInterval> segments = extractRefSegments(basicInfo, segmentingLocations);
        final List<String> descriptions = makeInterpretation(basicInfo, contigAlignments, jumps, segments);

        return new ReferenceSegmentsAndEventDescription(segments, descriptions);
    }

    // TODO: 11/6/17 fix 2-base long segment problem here
    private static List<SimpleInterval> extractRefSegments(final BasicInfo basicInfo,
                                                           final List<SimpleInterval> segmentingLocations) {

        final String eventPrimaryChromosome = basicInfo.eventPrimaryChromosome;

        final List<SimpleInterval> segments = new ArrayList<>(segmentingLocations.size() - 1);
        final Iterator<SimpleInterval> iterator = segmentingLocations.iterator();
        SimpleInterval leftBoundary = iterator.next();
        while (iterator.hasNext()) {
            final SimpleInterval rightBoundary = iterator.next();
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

        // using overlap with alignments ordered along the '+' strand representation of the signaling contig to
        // make sense of how the reference segments are ordered, including orientations--using signs of the integers,
        // on the sample;
        final List<AlignmentInterval> alignmentIntervalList =
                basicInfo.forwardStrandRep ? contigAlignments
                                           : ImmutableList.copyOf(contigAlignments).reverse();
        final Iterator<Boolean> jumpIterator =
                basicInfo.forwardStrandRep ? jumps.stream().map(jump -> jump.gapped).iterator()
                                           : ImmutableList.copyOf(jumps).reverse().stream().map(jump -> jump.gapped).iterator();
        Boolean currentJumpIsGapped = jumpIterator.next();
        boolean jumpIsLast = false;

        final SimpleInterval regionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();

        final List<String> descriptions = new ArrayList<>(2*segments.size()); //ini. cap. a guess
        final List<Tuple2<SimpleInterval, Integer>> insertionMappedToDisjointRegionAndWhereToInsert = new ArrayList<>();
        for (final AlignmentInterval alignment : alignmentIntervalList) {

            if (alignmentIsDisjointFromAlphaOmega(alignment.referenceSpan, regionBoundedByAlphaAndOmega)) {
                // disjoint alignment won't overlap any segments, so note down once where to insert, then move to next alignment
                final int indexABS = descriptions.size();
                insertionMappedToDisjointRegionAndWhereToInsert.add(new Tuple2<>(alignment.referenceSpan,
                        basicInfo.forwardStrandRep == alignment.forwardStrand ? indexABS : -1*indexABS));
                if ( currentJumpIsGapped )
                    descriptions.add(ReferenceSegmentsAndEventDescription.UNMAPPED_INSERTION);
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
                int i = start;
                for ( ; i != stop; i += step) {
                    final SimpleInterval currentSegment = segments.get(i);
                    // if current segment is contained in current alignment, note it down
                    if (alignmentContainsSegment(alignment.referenceSpan, currentSegment)) {
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
                if ( currentJumpIsGapped && !jumpIsLast){
                    descriptions.add(ReferenceSegmentsAndEventDescription.UNMAPPED_INSERTION);
                }
            }

            if (jumpIterator.hasNext()) // last alignment has no leaving jump so need to guard against that
                currentJumpIsGapped = jumpIterator.next();
            else
                jumpIsLast = true;
        }

        // post-hoc treatment of insertions that map to disjoint regions
        // go in reverse order because inserting into "descriptions", and insertion invalidates indices/iterators
        final ImmutableList<Tuple2<SimpleInterval, Integer>> reverse =
                ImmutableList.copyOf(insertionMappedToDisjointRegionAndWhereToInsert).reverse();
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

    // avoid case where segment and alignment ref span overlap only on one boundary base
    private static boolean alignmentContainsSegment(final SimpleInterval alignmentRefSpan, final SimpleInterval segment) {
        return alignmentRefSpan.overlaps(segment)
                &&
                alignmentRefSpan.intersect(segment).size() == segment.size();
    }

    private static boolean alignmentIsDisjointFromAlphaOmega(final SimpleInterval alignmentRefSpan,
                                                             final SimpleInterval regionBoundedByAlphaAndOmega) {
        return !alignmentRefSpan.overlaps(regionBoundedByAlphaAndOmega);
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
