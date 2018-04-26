package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * One of the two fundamental classes (the other is {@link CpxVariantCanonicalRepresentation})
 * for complex variant interpretation and alt haplotype extraction.
 *
 * <p>
 *     The general strategy is to walk on the affected reference region
 *     as guided by the alignments along the 5-3 direction of assembly contig.
 *     The walking has two basic types:
 *     <ul>
 *         <li>
 *             sliding (alignment block), and
 *         </li>
 *         <li>
 *             jumping (two split alignments joint on the read/contig) (see {@link Jump}).
 *         </li>
 *     </ul>
 *     By walking along the contig, we extract the jumps, then
 *     segment the affected region with segments bounded by the jumps,
 *     and walk along the alignments one more time to figure out
 *     how the segments are arranged on the sample/alt-haplotype.
 * </p>
 *
 * <p>
 *     This particular class is a prep step for {@link CpxVariantCanonicalRepresentation}, in that
 *     it extracts
 *     <ul>
 *         <li>
 *             basic information of the evidence contig (as defined in {@link BasicInfo}),
 *         </li>
 *         <li>
 *             how the various jumps are connected, and
 *         </li>
 *         <li>
 *            where the affected reference region are segmented
 *         </li>
 *     </ul>
 * </p>
 */
@DefaultSerializer(CpxVariantInducingAssemblyContig.Serializer.class)
final class CpxVariantInducingAssemblyContig {
    private static final AssemblyContigWithFineTunedAlignments.Serializer contigSerializer =
            new AssemblyContigWithFineTunedAlignments.Serializer();
    private static final BasicInfo.Serializer basicInfoSerializer =
            new BasicInfo.Serializer();
    private static final Jump.Serializer jumpSerializer =
            new Jump.Serializer();

    private final AssemblyContigWithFineTunedAlignments contigWithFineTunedAlignments;
    private final BasicInfo basicInfo;
    private final List<Jump> jumps;
    private final List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations;  // a list of 1bp "point" locations that will be used to segment the affected reference region

    @VisibleForTesting
    CpxVariantInducingAssemblyContig(final AssemblyContigWithFineTunedAlignments contigWithFineTunedAlignments,
                                     final BasicInfo basicInfo,
                                     final List<Jump> jumps,
                                     final List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations) {
        this.contigWithFineTunedAlignments = contigWithFineTunedAlignments;
        this.basicInfo = basicInfo;
        this.jumps = jumps;
        this.eventPrimaryChromosomeSegmentingLocations = eventPrimaryChromosomeSegmentingLocations;
    }

    /**
     * The process of initializing the fields is the processing of making sense of what happened:
     *
     * <p>
     *     Segment affected reference region by jumping locations, via
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
    CpxVariantInducingAssemblyContig(@Nonnull final AssemblyContigWithFineTunedAlignments contigWithFineTunedAlignments,
                                     @Nonnull final SAMSequenceDictionary refSequenceDictionary) {

        this.contigWithFineTunedAlignments = contigWithFineTunedAlignments;

        try {
            basicInfo = new BasicInfo(contigWithFineTunedAlignments.getSourceContig());

            jumps = extractJumpsOnReference(contigWithFineTunedAlignments.getAlignments());

            eventPrimaryChromosomeSegmentingLocations =
                    extractSegmentingRefLocationsOnEventPrimaryChromosome(jumps,
                            basicInfo,
                            refSequenceDictionary);
        } catch (final GATKException | IllegalArgumentException likelyNewEdgeCase) {
            throw new GATKException(toString(), likelyNewEdgeCase);
        }
    }

    /**
     * Each pair of neighboring reference locations are meant to be used closed, i.e. [a, b].
     * The returned list of {@link Jump}'s are ordered along the alignments of the contig.
     */
    @VisibleForTesting
    static List<Jump> extractJumpsOnReference(final List<AlignmentInterval> alignmentConfiguration) {

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
    @VisibleForTesting
    static List<SimpleInterval> extractSegmentingRefLocationsOnEventPrimaryChromosome(final List<Jump> jumps,
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

    static boolean alignmentIsDisjointFromAlphaOmega(final SimpleInterval alignmentRefSpan,
                                                     final SimpleInterval regionBoundedByAlphaAndOmega) {
        return !alignmentRefSpan.overlaps(regionBoundedByAlphaAndOmega);
    }

    // =============================================================================================================

    AssemblyContigWithFineTunedAlignments getPreprocessedTig() {
        return contigWithFineTunedAlignments;
    }
    BasicInfo getBasicInfo() {
        return basicInfo;
    }
    List<Jump> getJumps() {
        return jumps;
    }
    List<SimpleInterval> getEventPrimaryChromosomeSegmentingLocations() { return eventPrimaryChromosomeSegmentingLocations;}

    /**
     * @return two-base boundaries of contig alignments that are in the valid region, i.e. [alpha, omega]
     */
    @VisibleForTesting
    Set<SimpleInterval> getTwoBaseBoundaries() {
        // only look at alignments that are not disjoint from [alpha, omega]
        final SimpleInterval validRegion = new SimpleInterval(basicInfo.eventPrimaryChromosome, basicInfo.alpha.getStart(), basicInfo.omega.getEnd());

        final String chr = basicInfo.eventPrimaryChromosome;

        final Set<SimpleInterval> result = new HashSet<>(2 * contigWithFineTunedAlignments.getAlignments().size());
        for (final AlignmentInterval aln : contigWithFineTunedAlignments.getAlignments()) {
            final SimpleInterval referenceSpan = aln.referenceSpan;
            if (!validRegion.contains(referenceSpan))
                continue;

            final SimpleInterval left = new SimpleInterval(chr, referenceSpan.getStart(), referenceSpan.getStart()+1);
            final SimpleInterval right = new SimpleInterval(chr, referenceSpan.getEnd()-1, referenceSpan.getEnd());
            result.add(left);
            result.add(right);
        }

        return result;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("CpxVariantInducingAssemblyContig{");
        sb.append("contigWithFineTunedAlignments=").append(contigWithFineTunedAlignments);
        sb.append(", basicInfo=").append(basicInfo);
        sb.append(", jumps=").append(jumps);
        sb.append(", eventPrimaryChromosomeSegmentingLocations=").append(eventPrimaryChromosomeSegmentingLocations);
        sb.append('}');
        return sb.toString();
    }

    private CpxVariantInducingAssemblyContig(final Kryo kryo, final Input input) {
        contigWithFineTunedAlignments = contigSerializer.read(kryo, input, AssemblyContigWithFineTunedAlignments.class);

        basicInfo = basicInfoSerializer.read(kryo, input, BasicInfo.class);

        final int numJumps = input.readInt();
        jumps = new ArrayList<>(numJumps);
        for (int i = 0; i < numJumps; ++i) {
            jumps.add(jumpSerializer.read(kryo, input, Jump.class));
        }

        final int numSegmentingLocs = input.readInt();
        eventPrimaryChromosomeSegmentingLocations = new ArrayList<>(numJumps);
        for (int i = 0; i < numSegmentingLocs; ++i) {
            String chr = input.readString();
            int start = input.readInt();
            int end = input.readInt();
            eventPrimaryChromosomeSegmentingLocations.add( new SimpleInterval(chr, start, end));
        }
    }

    public void serialize(final Kryo kryo, final Output output) {
        contigSerializer.write(kryo, output, contigWithFineTunedAlignments);

        basicInfoSerializer.write(kryo, output, basicInfo);

        output.writeInt(jumps.size());
        for (final Jump jump : jumps)
            jumpSerializer.write(kryo, output, jump);

        output.writeInt(eventPrimaryChromosomeSegmentingLocations.size());
        for (final SimpleInterval segmentingLoc : eventPrimaryChromosomeSegmentingLocations) {
            output.writeString(segmentingLoc.getContig());
            output.writeInt(segmentingLoc.getStart());
            output.writeInt(segmentingLoc.getEnd());
        }
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<CpxVariantInducingAssemblyContig> {
        @Override
        public void write(final Kryo kryo, final Output output, final CpxVariantInducingAssemblyContig alignedContig) {
            alignedContig.serialize(kryo, output);
        }

        @Override
        public CpxVariantInducingAssemblyContig read(final Kryo kryo, final Input input, final Class<CpxVariantInducingAssemblyContig> clazz) {
            return new CpxVariantInducingAssemblyContig(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final CpxVariantInducingAssemblyContig that = (CpxVariantInducingAssemblyContig) o;

        if (!contigWithFineTunedAlignments.equals(that.contigWithFineTunedAlignments)) return false;
        if (!basicInfo.equals(that.basicInfo)) return false;
        if (!jumps.equals(that.jumps)) return false;
        return eventPrimaryChromosomeSegmentingLocations.equals(that.eventPrimaryChromosomeSegmentingLocations);
    }

    @Override
    public int hashCode() {
        int result = contigWithFineTunedAlignments.hashCode();
        result = 31 * result + basicInfo.hashCode();
        result = 31 * result + jumps.hashCode();
        result = 31 * result + eventPrimaryChromosomeSegmentingLocations.hashCode();
        return result;
    }

    // =============================================================================================================
    @DefaultSerializer(BasicInfo.Serializer.class)
    static final class BasicInfo {

        final String eventPrimaryChromosome; // where the head and tail alignments are mapped to, mappings to other chromosomes will be considered insertion/MEI
        final boolean forwardStrandRep;     // if the signaling assembly contig is a forward strand representation
        final SimpleInterval alpha;         // length-1 interval for the starting ref location of the head/tail alignment if the signaling assembly contig is a '+'/'-' representation
        final SimpleInterval omega;         // length-1 interval for the ending   ref location of the tail/head alignment if the signaling assembly contig is a '+'/'-' representation

        @VisibleForTesting
        BasicInfo(final String eventPrimaryChromosome, final boolean forwardStrandRep, final SimpleInterval alpha, final SimpleInterval omega) {
            this.eventPrimaryChromosome = eventPrimaryChromosome;
            this.forwardStrandRep = forwardStrandRep;
            this.alpha = alpha;
            this.omega = omega;
        }

        @VisibleForTesting
        BasicInfo(final AlignedContig contig) {
            if (contig.isUnmapped())
                throw new GATKException("Trying to extract information from unmapped read: " + contig.toString());
            final AlignmentInterval head = contig.getHeadAlignment();
            final AlignmentInterval tail = contig.getTailAlignment();

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
            String chr = input.readString();
            int start = input.readInt();
            int end = input.readInt();
            alpha = new SimpleInterval(chr, start, end);
            chr = input.readString();
            start = input.readInt();
            end = input.readInt();
            omega = new SimpleInterval(chr, start, end);
        }

        void serialize(final Kryo kryo, final Output output) {
            output.writeString(eventPrimaryChromosome);
            output.writeBoolean(forwardStrandRep);
            output.writeString(alpha.getContig());
            output.writeInt(alpha.getStart());
            output.writeInt(alpha.getEnd());
            output.writeString(omega.getContig());
            output.writeInt(omega.getStart());
            output.writeInt(omega.getEnd());
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

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final BasicInfo basicInfo = (BasicInfo) o;

            if (forwardStrandRep != basicInfo.forwardStrandRep) return false;
            if (!eventPrimaryChromosome.equals(basicInfo.eventPrimaryChromosome)) return false;
            if (!alpha.equals(basicInfo.alpha)) return false;
            return omega.equals(basicInfo.omega);
        }

        @Override
        public int hashCode() {
            int result = eventPrimaryChromosome.hashCode();
            result = 31 * result + (forwardStrandRep ? 1 : 0);
            result = 31 * result + alpha.hashCode();
            result = 31 * result + omega.hashCode();
            return result;
        }
    }

    /**
     * A simple helper struct for complex event interpretation.
     *
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
     * {@link CpxVariantInterpreter#yieldOverlapToAlignmentTwo(AlignmentInterval, AlignmentInterval, SAMSequenceDictionary)}.
     *
     * </p>
     */
    @DefaultSerializer(Jump.Serializer.class)
    static final class Jump {

        final SimpleInterval start;
        final SimpleInterval landing;
        final StrandSwitch strandSwitch;
        final int gapSize; // jump is gapped when this size is > 0

        @VisibleForTesting
        Jump(final SimpleInterval start, final SimpleInterval landing, final StrandSwitch strandSwitch, final int gapSize) {
            this.start = start;
            this.landing = landing;
            this.strandSwitch = strandSwitch;
            this.gapSize = gapSize;
        }

        @VisibleForTesting
        Jump(final AlignmentInterval one, final AlignmentInterval two) {
            Utils.validateArg(AlignmentInterval.overlapOnContig(one, two) <=0,
                    "assumption that input alignments DO NOT overlap is violated.");

            strandSwitch = SimpleChimera.determineStrandSwitch(one, two);

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

            gapSize = Math.max(0, two.startInAssembledContig - one.endInAssembledContig - 1); // -1 because we want bases in-between
        }

        boolean isGapped() {
            return gapSize > 0;
        }

        @Override
        public String toString() {
            return "Jump start: " + start.toString() + "\tjump landing: " + landing.toString() +
                    "\t" + strandSwitch.name() + "\t" + (gapSize > 0 ? "G" : "C");
        }

        private Jump(final Kryo kryo, final Input input) {

            String chr = input.readString();
            int beg = input.readInt();
            int end = input.readInt();
            start = new SimpleInterval(chr, beg, end);
            chr = input.readString();
            beg = input.readInt();
            end = input.readInt();
            landing = new SimpleInterval(chr, beg, end);

            strandSwitch = StrandSwitch.values()[input.readInt()];
            gapSize = input.readInt();
        }

        public void serialize(final Kryo kryo, final Output output) {
            output.writeString(start.getContig());
            output.writeInt(start.getStart());
            output.writeInt(start.getEnd());
            output.writeString(landing.getContig());
            output.writeInt(landing.getStart());
            output.writeInt(landing.getEnd());
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

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final Jump jump = (Jump) o;

            if (gapSize != jump.gapSize) return false;
            if (!start.equals(jump.start)) return false;
            if (!landing.equals(jump.landing)) return false;
            return strandSwitch == jump.strandSwitch;
        }

        @Override
        public int hashCode() {
            int result = start.hashCode();
            result = 31 * result + landing.hashCode();
            result = 31 * result + strandSwitch.ordinal();
            result = 31 * result + gapSize;
            return result;
        }
    }
}
