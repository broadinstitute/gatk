package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.TextMDCodec;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * Various types of read anomalies that provide evidence of genomic breakpoints.
 * There is a shallow hierarchy based on this class that classifies the anomaly as, for example, a split read,
 *   or a pair that has both reads on the same reference strand.
 * Each BreakpointEvidence object of type ReadEvidence object comes from examining a single read, and describing its
 *   funkiness, if any, by instantiating one of the subclasses that addresses that type of funkiness. Instances of
 *   BreakpointEvidence that are not of subtype ReadEvidence consolidate information from one or more pieces of
 *   ReadEvidence into a single object.
 */
@DefaultSerializer(BreakpointEvidence.Serializer.class)
public class BreakpointEvidence {
    private static final SVInterval.Serializer intervalSerializer = new SVInterval.Serializer();
    private final SVInterval location;
    private final int weight;
    private boolean validated; // this piece of evidence is consistent with enough other evidence to be taken seriously

    public BreakpointEvidence( final SVInterval location, final int weight, final boolean validated ) {
        this.location = location;
        this.weight = weight;
        this.validated = validated;
    }

    protected BreakpointEvidence( final Kryo kryo, final Input input ) {
        this.location = intervalSerializer.read(kryo, input, SVInterval.class);
        this.weight = input.readInt();
        this.validated = input.readBoolean();
    }

    public SVInterval getLocation() { return location; }
    public int getWeight() { return weight; }
    public boolean isValidated() { return validated; }
    public void setValidated( final boolean validated ) { this.validated = validated; }

    /**
     * If true: the evidence suggests a breakpoint at a reference location upstream of the interval's start coordinate
     * If false: the evidence suggests a breakpoint downstream of the interval's end coordinate
     */
    public Boolean isEvidenceUpstreamOfBreakpoint() { return null; }

    /**
     * Returns true if this piece of evidence specifies a possible distal target for the breakpoint.
     * @param readMetadata Read metadata for the library
     * @param minEvidenceMapQ The minimum mapping quality threshold (inclusive) for which a target (mate or SA mapping) should be created
     */
    public boolean hasDistalTargets(final ReadMetadata readMetadata, final int minEvidenceMapQ) {
        return false;
    }

    /**
     * Returns the distal interval implicated as a candidate adjacency to the breakpoint by this piece of evidence.
     * For example, in the case of a discordant read pair, this would be the region adjacent to the mate of the current
     * read. Returns null if the evidence does not specify or support a possible targeted region (for example, the case
     * of an read with an unmapped mate). Strands of the intervals indicate whether the distal target intervals are
     * upstream or downstream of their proposed breakpoints: true indicates that the breakpoint is upstream of the interval
     * start position; false indicates that the breakpoint is downstream of the interval end position.
     * @param minEvidenceMapq The minimum mapping quality (inclusive) of the target evidence; if the mapping quality is below
     *                        this value the evidence is treated as not having a distal target
     */
    public List<StrandedInterval> getDistalTargets(final ReadMetadata readMetadata, final int minEvidenceMapq) {
        return null;
    }

    @Override
    public String toString() {
        return location.toString() + "^" + weight;
    }

    /**
     * Returns string representation of BreakpointEvidence in tab-separated form:
     * Contig[begin:end]    weight  validated   EvidenceType    distal_targets
     * distal_targets is a (; separated) list of string representations (.toString()) for each distal target
     *     it is empty if there are no distal targets
     * Child classes may extend stringRep() by appending class-specific tab-separated info
     */
    public String stringRep(final ReadMetadata readMetadata, final int minEvidenceMapq) {
        final String dtString =  hasDistalTargets(readMetadata, minEvidenceMapq) ?
                getDistalTargets(readMetadata, minEvidenceMapq).stream()
                    .map(strandedInterval -> strandedInterval.getInterval().toString() + (strandedInterval.getStrand() ? "1" : "0"))
                    .collect(Collectors.joining(";"))
                : "";
        return location.toString() + "\t" + weight + "\t" + this.getClass().getSimpleName() + "\t" + dtString;
    }

    //* slicing equality -- just tests for equal fields */
    @VisibleForTesting boolean equalFields( final BreakpointEvidence that ) {
        return location.equals(that.location) && weight == that.weight && validated == that.validated;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BreakpointEvidence)) return false;
        BreakpointEvidence evidence = (BreakpointEvidence) o;
        return weight == evidence.weight &&
                validated == evidence.validated &&
                Objects.equals(location, evidence.location);
    }

    @Override
    public int hashCode() {

        return Objects.hash(location, weight, validated);
    }

    protected void serialize(final Kryo kryo, final Output output ) {
        intervalSerializer.write(kryo, output, location);
        output.writeInt(weight);
        output.writeBoolean(validated);
    }

    static SVInterval fixedWidthInterval( final int contigID,
                                                  final int contigOffset, final int offsetUncertainty ) {
        int width = 2 * offsetUncertainty;
        int start = contigOffset - offsetUncertainty;
        if ( start < 1 ) {
            width += start - 1;
            start = 1;
        }
        return new SVInterval(contigID, start, start + width);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<BreakpointEvidence> {
        @Override
        public void write( final Kryo kryo, final Output output, final BreakpointEvidence evidence ) {
            evidence.serialize(kryo, output);
        }

        @Override
        public BreakpointEvidence read( final Kryo kryo, final Input input, final Class<BreakpointEvidence> klass ) {
            return new BreakpointEvidence(kryo, input);
        }
    }

    @DefaultSerializer(ExternalEvidence.Serializer.class)
    public final static class ExternalEvidence extends BreakpointEvidence {
        public ExternalEvidence( final int contigID, final int start, final int end, final int weight ) {
            this(new SVInterval(contigID, start, end), weight);
        }

        public ExternalEvidence( final SVInterval interval, final int weight ) {
            super(interval, weight, false);
        }

        public ExternalEvidence( final Kryo kryo, final Input input ) {
            super(kryo, input);
        }

        @Override
        public String toString() { return super.toString() + "\tExternalEvidence"; }

        public final static class Serializer extends com.esotericsoftware.kryo.Serializer<ExternalEvidence> {
            @Override
            public void write( final Kryo kryo, final Output output, final ExternalEvidence externalEvidence ) {
                externalEvidence.serialize(kryo, output);
            }

            @Override
            public ExternalEvidence read( final Kryo kryo, final Input input, final Class<ExternalEvidence> klass ) {
                return new ExternalEvidence(kryo, input);
            }
        }
    }

    @DefaultSerializer(TemplateSizeAnomaly.Serializer.class)
    public final static class TemplateSizeAnomaly extends BreakpointEvidence {
        private final int readCount;

        public TemplateSizeAnomaly( final SVInterval interval, final int weight, final int readCount ) {
            super(interval, weight, false);
            this.readCount = readCount;
        }

        protected TemplateSizeAnomaly( final Kryo kryo, final Input input ) {
            super(kryo, input);
            readCount = input.readInt();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(readCount);
        }

        @Override
        public String toString() {
            return super.toString() + "\tTemplateSizeAnomaly\t" + readCount;
        }

        @Override
        public String stringRep(final ReadMetadata readMetadata, final int minEvidenceMapq) {
            return super.stringRep(readMetadata, minEvidenceMapq) + "\t" + readCount;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof TemplateSizeAnomaly)) return false;
            if (!super.equals(o)) return false;
            TemplateSizeAnomaly that = (TemplateSizeAnomaly) o;
            return readCount == that.readCount;
        }

        @Override
        public int hashCode() {

            return Objects.hash(super.hashCode(), readCount);
        }

        public final static class Serializer extends com.esotericsoftware.kryo.Serializer<TemplateSizeAnomaly> {
            @Override
            public void write( final Kryo kryo, final Output output, final TemplateSizeAnomaly templateSizeAnomaly ) {
                templateSizeAnomaly.serialize(kryo, output);
            }

            @Override
            public TemplateSizeAnomaly read( final Kryo kryo, final Input input, final Class<TemplateSizeAnomaly> klass ) {
                return new TemplateSizeAnomaly(kryo, input);
            }
        }

        public Integer getReadCount() {return readCount;}
    }

    @DefaultSerializer(ReadEvidence.Serializer.class)
    public static class ReadEvidence extends BreakpointEvidence {
        private static final int SINGLE_READ_WEIGHT = 1;
        private final String templateName; // QNAME of the read that was funky (i.e., the name of the fragment)
        private final TemplateFragmentOrdinal fragmentOrdinal; // which read we're talking about (first or last, for paired-end reads)
        private final boolean forwardStrand;
        private final String cigarString;
        private final int mappingQuality;
        private final int templateSize;
        private final String readGroup;

        /**
         * evidence offset and width is set to "the rest of the fragment" not covered by this read
         */
        protected ReadEvidence( final GATKRead read, final ReadMetadata metadata, final int readWeight ) {
            super(restOfFragmentInterval(read,metadata), readWeight, false);
            this.templateName = read.getName();
            if ( templateName == null ) throw new GATKException("Read has no name.");
            this.fragmentOrdinal = TemplateFragmentOrdinal.forRead(read);
            this.forwardStrand = ! read.isReverseStrand();
            this.cigarString = read.getCigar().toString();
            this.mappingQuality = read.getMappingQuality();
            this.templateSize = read.getFragmentLength();
            this.readGroup = read.getReadGroup();
        }

        /**
         * for use when the uncertainty in location has a fixed size
         */
        protected ReadEvidence( final GATKRead read, final ReadMetadata metadata, final int contigOffset,
                                final int offsetUncertainty, final boolean forwardStrand, final int readWeight ) {
            super(fixedWidthInterval(metadata.getContigID(read.getContig()),contigOffset,offsetUncertainty),
                  readWeight, false);
            this.templateName = read.getName();
            if ( templateName == null ) throw new GATKException("Read has no name.");
            this.fragmentOrdinal = TemplateFragmentOrdinal.forRead(read);
            this.forwardStrand = forwardStrand;
            this.cigarString = read.getCigar().toString();
            this.mappingQuality = read.getMappingQuality();
            this.templateSize = read.getFragmentLength();
            this.readGroup = read.getReadGroup();
        }

        /**
         * Directly construct ReadEvidence by supplying all fields. Used by testing
         */
        @VisibleForTesting ReadEvidence( final SVInterval interval, final int weight,
                                         final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                                         final boolean validated, final boolean forwardStrand,
                                         final String cigarString, final int mappingQuality,
                                         final int templateSize, final String readGroup) {
            super(interval, weight, validated);
            this.templateName = templateName;
            this.fragmentOrdinal = fragmentOrdinal;
            this.forwardStrand = forwardStrand;
            this.cigarString = cigarString;
            this.mappingQuality = mappingQuality;
            this.templateSize = templateSize;
            this.readGroup = readGroup;
        }

        /**
         * a technical constructor for use in Kryo (de-)serialization.
         * this creates an object by reading a Kryo-serialized stream.
         * it will be called by subclasses in their own constructors from Kryo streams (as super(kryo, input)).
         */
        protected ReadEvidence( final Kryo kryo, final Input input ) {
            super(kryo, input);
            this.templateName = input.readString();
            this.fragmentOrdinal = TemplateFragmentOrdinal.values()[input.readByte()];
            this.forwardStrand = input.readBoolean();
            this.cigarString = input.readString();
            this.mappingQuality = input.readInt();
            this.templateSize = input.readInt();
            this.readGroup = input.readString();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeString(templateName);
            output.writeByte(fragmentOrdinal.ordinal());
            output.writeBoolean(forwardStrand);
            output.writeString(cigarString);
            output.writeInt(mappingQuality);
            output.writeInt(templateSize);
            output.writeString(readGroup);
        }

        public boolean getForwardStrand() { return forwardStrand; }

        public String getTemplateName() {
            return templateName;
        }

        public TemplateFragmentOrdinal getFragmentOrdinal() {
            return fragmentOrdinal;
        }

        public String getCigarString() { return cigarString; }

        public int getMappingQuality() { return mappingQuality; }

        public int getTemplateSize() { return templateSize; }

        public String getReadGroup() { return readGroup; }

        @Override
        public Boolean isEvidenceUpstreamOfBreakpoint() {
            return forwardStrand;
        }

        @Override
        public String toString() {
            return super.toString() + "\t" + templateName + fragmentOrdinal;
        }

        /**
         * All ReadEvidence has basic BreakpointEvidence.stringRep() with appended info:
         * templateName/fragmentOrdinal  forwardStrand  cigarString    mappingQuality
         * Child classes may extend stringRep() by appending class-specific tab-separated info
         */
        @Override
        public String stringRep(final ReadMetadata readMetadata, final int minEvidenceMapq) {
            return super.stringRep(readMetadata, minEvidenceMapq) + "\t" + templateName + fragmentOrdinal
                    + "\t" + (forwardStrand ? "1" : "0") + "\t" + templateSize + "\t" + cigarString
                    + "\t" + mappingQuality;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof ReadEvidence)) return false;
            if (!super.equals(o)) return false;
            ReadEvidence that = (ReadEvidence) o;
            return forwardStrand == that.forwardStrand &&
                    mappingQuality == that.mappingQuality &&
                    templateSize == that.templateSize &&
                    Objects.equals(templateName, that.templateName) &&
                    fragmentOrdinal == that.fragmentOrdinal &&
                    Objects.equals(cigarString, that.cigarString) &&
                    Objects.equals(readGroup, that.readGroup);
        }

        @Override
        public int hashCode() {

            return Objects.hash(super.hashCode(), templateName, fragmentOrdinal, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ReadEvidence> {
            @Override
            public void write( final Kryo kryo, final Output output, final ReadEvidence evidence ) {
                evidence.serialize(kryo, output);
            }

            @Override
            public ReadEvidence read( final Kryo kryo, final Input input, final Class<ReadEvidence> klass ) {
                return new ReadEvidence(kryo, input);
            }
        }

        private static SVInterval restOfFragmentInterval( final GATKRead read, final ReadMetadata metadata ) {
            final int templateLen = metadata.getFragmentLengthStatistics(read.getReadGroup()).getMaxNonOutlierFragmentSize();
            int width;
            int start;
            if ( read.isReverseStrand() ) {
                // we can get a little more precise about the interval by checking to see if there are any leading mismatches
                // in the read's alignment and trimming them off.
                int leadingMismatches = getLeadingMismatches(read, true);
                final int readStart = read.getStart() + leadingMismatches;
                width = readStart - (read.getUnclippedEnd() + 1 - templateLen);
                start = readStart - width;
                if ( start < 1 ) {
                    width += start - 1;
                    start = 1;
                }
            } else {
                int trailingMismatches = getLeadingMismatches(read, false);
                final int readEnd = read.getEnd() + 1 - trailingMismatches;
                width = read.getUnclippedStart() + templateLen - readEnd;
                start = readEnd;
            }
            return new SVInterval(metadata.getContigID(read.getContig()), start, start + width);
        }

        @VisibleForTesting static int getLeadingMismatches(final GATKRead read, final boolean fromStart) {
            int leadingMismatches = 0;
            if (read.hasAttribute("MD")) {
                final String mdString = read.getAttributeAsString("MD");
                final List<TextMDCodec.MDElement> mdElements = TextMDCodec.parseMDString(mdString);
                int idx = fromStart ? 0 : (mdElements.size() - 1);
                while (fromStart ? (idx < mdElements.size()) : (idx >= 0)) {
                    TextMDCodec.MDElement mdElement = mdElements.get(idx);
                    if (mdElement instanceof TextMDCodec.MatchMDElement && mdElement.getLength() > 0) {
                        break;
                    } else {
                        leadingMismatches += mdElement.getLength();
                    }
                    idx = idx + (fromStart ? 1 : -1);
                }
            }
            return leadingMismatches;
        }

    }

    @DefaultSerializer(SplitRead.Serializer.class)
    public static final class SplitRead extends ReadEvidence {
        @VisibleForTesting static final int UNCERTAINTY = 3;
        private static final String SA_TAG_NAME = "SA";
        private static final int SPLIT_READ_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;
        private final String tagSA;
        private boolean primaryAlignmentClippedAtStart;
        private boolean primaryAlignmentForwardStrand;
        private final List<SAMapping> saMappings;

        public SplitRead( final GATKRead read, final ReadMetadata metadata, final boolean primaryAlignmentClippedAtAlignmentStart ) {
            // todo: if reads have multiple SA tags.. we should have two pieces of evidence with the right strands
            super(read, metadata, primaryAlignmentClippedAtAlignmentStart ? read.getStart() : read.getEnd(),
                  UNCERTAINTY, !primaryAlignmentClippedAtAlignmentStart, SPLIT_READ_WEIGHT);
            this.primaryAlignmentForwardStrand = !read.isReverseStrand();
            if ( getCigarString() == null || getCigarString().isEmpty() ) {
                throw new GATKException("Read has no cigar string.");
            }
            this.primaryAlignmentClippedAtStart = primaryAlignmentClippedAtAlignmentStart;
            if (read.hasAttribute(SA_TAG_NAME)) {
                tagSA = read.getAttributeAsString(SA_TAG_NAME);
            } else {
                tagSA = null;
            }
            saMappings = parseSATag(tagSA);
        }

        private SplitRead( final Kryo kryo, final Input input ) {
            super(kryo, input);
            tagSA = input.readString();
            primaryAlignmentClippedAtStart = input.readBoolean();
            saMappings = parseSATag(tagSA);
        }

        /**
         * Directly construct SplitRead by supplying all fields. Used by testing
         */
        @VisibleForTesting SplitRead(final SVInterval interval, final int weight,
                          final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                          final boolean validated, final boolean forwardStrand,
                          final String cigarString, final int mappingQuality,
                          final int templateSize, final String readGroup,
                          final boolean primaryAlignmentClippedAtStart, boolean primaryAlignmentForwardStrand,
                          final String tagSA) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup);
            this.primaryAlignmentClippedAtStart = primaryAlignmentClippedAtStart;
            this.primaryAlignmentForwardStrand = primaryAlignmentForwardStrand;
            this.tagSA = tagSA;
            this.saMappings = parseSATag(tagSA);
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeString(tagSA);
            output.writeBoolean(primaryAlignmentClippedAtStart);
        }

        @Override
        public String toString() {
            return super.toString()+"\tSplit\t"+getCigarString()+"\t"+(tagSA == null ? " SA: None" : (" SA: " + tagSA));
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof SplitRead)) return false;
            if (!super.equals(o)) return false;
            SplitRead splitRead = (SplitRead) o;
            return primaryAlignmentClippedAtStart == splitRead.primaryAlignmentClippedAtStart &&
                    primaryAlignmentForwardStrand == splitRead.primaryAlignmentForwardStrand &&
                    Objects.equals(tagSA, splitRead.tagSA) &&
                    Objects.equals(saMappings, splitRead.saMappings);
        }

        @Override
        public int hashCode() {

            return Objects.hash(super.hashCode(), tagSA, primaryAlignmentClippedAtStart, primaryAlignmentForwardStrand, saMappings);
        }

        private List<SAMapping> parseSATag(final String tagSA) {
            if (tagSA == null) {
                return null;
            } else {
                final String[] saStrings = tagSA.split(";");
                final List<SAMapping> supplementaryAlignments = new ArrayList<>(saStrings.length);
                for (final String saString : saStrings) {
                    final String[] saFields = saString.split(",", -1);
                    if (saFields.length != 6) {
                        throw new GATKException("Could not parse SATag: "+ saString);
                    }
                    final String contigId = saFields[0];
                    final int pos = Integer.parseInt(saFields[1]);
                    final boolean strand = "+".equals(saFields[2]);
                    final String cigarString = saFields[3];
                    final int mapQ = Integer.parseInt(saFields[4]);
                    final int mismatches = Integer.parseInt(saFields[5]);
                    SAMapping saMapping = new SAMapping(contigId, pos, strand, cigarString, mapQ, mismatches);

                    supplementaryAlignments.add(saMapping);
                }
                return supplementaryAlignments;

            }

        }

        @Override
        public boolean hasDistalTargets(final ReadMetadata readMetadata, final int minEvidenceMapQ) {
            // todo: right now we are limiting distal target calculation to split reads that have only one supplementary mapping
            return saMappings != null && saMappings.size() == 1 && hasHighQualitySupplementaryMappings(readMetadata, minEvidenceMapQ);
        }

        private boolean hasHighQualitySupplementaryMappings(final ReadMetadata readMetadata, final int minEvidenceMapq) {
            boolean hqMappingFound = false;
            if (saMappings != null) {
                for (final SAMapping mapping : saMappings) {
                    final int mapQ = mapping.getMapq();
                    SVInterval saInterval = saMappingToSVInterval(readMetadata, mapping,
                            calculateDistalTargetStrand(mapping, !primaryAlignmentClippedAtStart, primaryAlignmentForwardStrand));
                    if (isHighQualityMapping(readMetadata, mapQ, saInterval, minEvidenceMapq)) {
                        hqMappingFound = true;
                        break;
                    }
                }
            }
            return hqMappingFound;
        }

        private boolean isHighQualityMapping(final ReadMetadata readMetadata, final int mapQ, final SVInterval saInterval, final int minEvidenceMapq) {
            return mapQ >= minEvidenceMapq &&
                    (saInterval.getContig() == getLocation().getContig() ||
                    (!readMetadata.ignoreCrossContigID(saInterval.getContig()) && !readMetadata.ignoreCrossContigID(getLocation().getContig())))
                    && ! saInterval.overlaps(getLocation());
        }

        /**
         * In the case of a split read with an SA tag, the distal target interval is a StrandedInterval where the interval
         * location is the possible breakpoint interval given by the clipping location on the supplementary alignment,
         * and the strand indicates whether or not the rest of the supplementary alignment is upstream or downstream of
         * of the breakpoint. For example, for a split read spanning a deletion:
         *
         *    --C|xxxxxxxxxxx|C->
         *      +++         ^^^
         *
         *  C     = Clipping location
         *  --C   = Evidence
         *  C->   = Target evidence given by the SA tag
         *  +++   = Evidence getLocation
         *  |xxx| = Deletion
         *  ^^^   = Distal target interval
         *
         *  Since the clipping is at the left side of the distal target interval, the distal target strand will be (-)
         *  (because the putative breakpoint upstream of the end of the distal target interval), and the StrandedInterval
         *  implied by BreakpointEvidence.getLocation() and BreakpointEvidence.isEvidenceUpstreamOfBreakpoint()
         *  will have a strand of (+).
         */
        @Override
        public List<StrandedInterval> getDistalTargets(final ReadMetadata readMetadata, final int minEvidenceMapq) {
            if (hasDistalTargets(readMetadata, minEvidenceMapq)) {
                final List<StrandedInterval> supplementaryAlignments = new ArrayList<>(saMappings.size());
                for (final SAMapping saMapping : saMappings) {
                    final int mapQ = saMapping.getMapq();
                    final boolean strand = calculateDistalTargetStrand(saMapping, !primaryAlignmentClippedAtStart,
                                                                       primaryAlignmentForwardStrand);
                    final SVInterval saInterval = saMappingToSVInterval(readMetadata, saMapping, strand);
                    if (! isHighQualityMapping(readMetadata, mapQ, saInterval, minEvidenceMapq)) {
                        continue;
                    }
                    supplementaryAlignments.add(new StrandedInterval(saInterval, strand));
                }
                return supplementaryAlignments;
            } else {
                return null;
            }
        }

        @VisibleForTesting
        static boolean calculateDistalTargetStrand(final SAMapping saMapping, final boolean primaryEvidenceRightClipped,
                                                   final boolean primaryAlignmentForwardStrand) {
            final boolean primaryAlignmentRightClippedOnRead = (primaryAlignmentForwardStrand == primaryEvidenceRightClipped);

            // todo: assuming that we only have one SA tag and it goes with the clipped end of the SplitRead evidence
            if (primaryAlignmentRightClippedOnRead) {
                return !saMapping.isForwardStrand();
            } else {
                return saMapping.isForwardStrand();
            }

        }

        // todo: for now, taking the entire location of the supplementary alignment plus the uncertainty on each end
        // A better solution might be to find the location of the actual clip on the other end of the reference,
        // but that would be significantly more complex and possibly computationally expensive
        private SVInterval saMappingToSVInterval(final ReadMetadata readMetadata, final SAMapping saMapping,
                                                 final boolean saEvidenceDownstreamOfBreakpoint) {
            final int contigId = readMetadata.getContigID(saMapping.getContigName());
            final Cigar saCigar = TextCigarCodec.decode(saMapping.getCigar());
            final int pos = saEvidenceDownstreamOfBreakpoint ? saMapping.getStart() + saCigar.getPaddedReferenceLength() : saMapping.getStart();

            return new SVInterval( contigId,
                    Math.max(1, pos - UNCERTAINTY),
                    pos + UNCERTAINTY + 1);

        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SplitRead> {
            @Override
            public void write( final Kryo kryo, final Output output, final SplitRead splitRead ) {
                splitRead.serialize(kryo, output);
            }

            @Override
            public SplitRead read( final Kryo kryo, final Input input, final Class<SplitRead> klass ) {
                return new SplitRead(kryo, input);
            }
        }

        final static class SAMapping {
            private final String contigName;
            private final int start;
            private final boolean forwardStrand;
            private final String cigar;
            private final int mapq;
            private final int mismatches;

            public SAMapping(final String contigName, final int start, final boolean strand, final String cigar, final int mapq, final int mismatches) {
                this.contigName = contigName;
                this.start = start;
                this.forwardStrand = strand;
                this.cigar = cigar;
                this.mapq = mapq;
                this.mismatches = mismatches;
            }

            public String getContigName() {
                return contigName;
            }

            public int getStart() {
                return start;
            }

            public boolean isForwardStrand() {
                return forwardStrand;
            }

            public String getCigar() {
                return cigar;
            }

            public int getMapq() {
                return mapq;
            }

            public int getMismatches() {
                return mismatches;
            }
        }
    }

    @DefaultSerializer(LargeIndel.Serializer.class)
    public static final class LargeIndel extends ReadEvidence {
        private static final int UNCERTAINTY = 4;
        private static final int LARGE_INDEL_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;

        LargeIndel( final GATKRead read, final ReadMetadata metadata, final int contigOffset ) {
            super(read, metadata, contigOffset, UNCERTAINTY, true, LARGE_INDEL_WEIGHT);
            if ( getCigarString() == null || getCigarString().isEmpty() )
                throw new GATKException("Read has no cigar string.");
        }

        private LargeIndel( final Kryo kryo, final Input input ) {
            super(kryo, input);
        }

        /**
         * Directly construct LargeIndel by supplying all fields. Used by testing
         */
        @VisibleForTesting LargeIndel(final SVInterval interval, final int weight,
                             final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                             final boolean validated, final boolean forwardStrand,
                             final String cigarString, final int mappingQuality,
                             final int templateSize, final String readGroup) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup);
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
        }

        @Override
        public String toString() {return super.toString() + "\tIndel\t" + getCigarString();}

        @Override
        public boolean equals(Object o) {
            return (this == o) || ((o instanceof LargeIndel) && super.equals(o));
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LargeIndel> {
            @Override
            public void write( final Kryo kryo, final Output output, final LargeIndel largeIndel ) {
                largeIndel.serialize(kryo, output);
            }

            @Override
            public LargeIndel read( final Kryo kryo, final Input input, final Class<LargeIndel> klass ) {
                return new LargeIndel(kryo, input);
            }
        }
    }

    @DefaultSerializer(MateUnmapped.Serializer.class)
    public static final class MateUnmapped extends ReadEvidence {
        private static final int MATE_UNMAPPED_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;

        MateUnmapped( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata, MATE_UNMAPPED_WEIGHT);
        }

        private MateUnmapped( final Kryo kryo, final Input input ) { super(kryo, input); }

        /**
         * Directly construct MateUnmapped by supplying all fields. Used by testing
         */
        @VisibleForTesting MateUnmapped(final SVInterval interval, final int weight,
                             final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                             final boolean validated, final boolean forwardStrand,
                             final String cigarString, final int mappingQuality,
                             final int templateSize, final String readGroup) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup);
        }

        @Override
        public String toString() {
            return super.toString() + "\tUnmappedMate";
        }

        @Override
        public boolean equals(Object o) {
            return (this == o) || ((o instanceof MateUnmapped) && super.equals(o));
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<MateUnmapped> {
            @Override
            public void write( final Kryo kryo, final Output output, final MateUnmapped mateUnmapped ) {
                mateUnmapped.serialize(kryo, output);
            }

            @Override
            public MateUnmapped read( final Kryo kryo, final Input input, final Class<MateUnmapped> klass ) {
                return new MateUnmapped(kryo, input);
            }
        }
    }

    public static abstract class DiscordantReadPairEvidence extends ReadEvidence {
        protected final SVInterval target;
        protected final boolean targetForwardStrand;
        protected final int targetQuality;

        // even if we have access to and use the mate cigar, we still don't really know the exact breakpoint interval
        // specified by the mate since there could be unclipped mismatches at the ends of the alignment. This constant
        // tries to correct for that.
        public static final int MATE_ALIGNMENT_LENGTH_UNCERTAINTY = 2;

        public DiscordantReadPairEvidence(final GATKRead read, final ReadMetadata metadata, final int weight) {
            super(read, metadata, weight);
            target = getMateTargetInterval(read, metadata);
            targetForwardStrand = getMateForwardStrand(read);
            if (read.hasAttribute("MQ")) {
                targetQuality = read.getAttributeAsInteger("MQ");
            } else {
                targetQuality = Integer.MAX_VALUE;
            }
        }

        public DiscordantReadPairEvidence(final Kryo kryo, final Input input) {
            super(kryo, input);
            target = intervalSerializer.read(kryo, input, SVInterval.class);
            targetForwardStrand = input.readBoolean();
            targetQuality = input.readInt();
        }

        /**
         * Directly construct DiscordantReadPairEvidence by supplying all fields. Used by testing
         */
        private DiscordantReadPairEvidence(final SVInterval interval, final int weight,
                                  final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                                  final boolean validated, final boolean forwardStrand,
                                  final String cigarString, final int mappingQuality,
                                  final int templateSize, final String readGroup,
                                  final SVInterval target, final boolean targetForwardStrand, final int targetQuality) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup);
            this.target = target;
            this.targetForwardStrand = targetForwardStrand;
            this.targetQuality = targetQuality;
        }

        @Override
        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            intervalSerializer.write(kryo, output, target);
            output.writeBoolean(targetForwardStrand);
            output.writeInt(targetQuality);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof DiscordantReadPairEvidence)) return false;
            if (!super.equals(o)) return false;
            DiscordantReadPairEvidence that = (DiscordantReadPairEvidence) o;
            return targetForwardStrand == that.targetForwardStrand &&
                    targetQuality == that.targetQuality &&
                    Objects.equals(target, that.target);
        }

        @Override
        public int hashCode() {

            return Objects.hash(super.hashCode(), target, targetForwardStrand, targetQuality);
        }

        @Override
        public boolean hasDistalTargets(final ReadMetadata readMetadata, final int minEvidenceMapQ) {
            return getLocation().isUpstreamOf(target)
                    && isTargetHighQuality(readMetadata, minEvidenceMapQ);
        }

        private boolean isTargetHighQuality(final ReadMetadata readMetadata, final int minEvidenceMapq) {
            return targetQuality >= minEvidenceMapq
                    && ! target.overlaps(getLocation()) && ! readMetadata.ignoreCrossContigID(target.getContig());
        }

        /**
         * In the case of a discordant read pair, the distal target interval is a StrandedInterval where the interval
         * location is the possible breakpoint interval given by the inferred rest-of-fragment interval for the mate read
         * (ie the region that would contain the mate's mated read if its fragment size had been drawn from the non-outlier
         * fragment size distribution), and the strand can be easily computed as the strand of the reference the mate
         * mapped to. For example in the case of a deletion, if this piece of BreakpointEvidence represents the forward
         * strand mapping read in a long-fragment size pair:
         *
         *    ---->   |xxxxxxxxxxx|   <----
         *         +++++        ^^^^^^
         *  ----> = Evidence
         *  +++++ = Evidence getLocation
         *  <---- = Mate
         *  |xxx| = Deletion
         *  ^^^^^ = Distal target interval
         *
         *  The distal target strand will be (-), and the StrandedInterval implied by BreakpointEvidence.getLocation() and
         *  BreakpointEvidence.isEvidenceUpstreamOfBreakpoint() will have a strand of (+).
         */
        @Override
        public List<StrandedInterval> getDistalTargets(final ReadMetadata readMetadata, final int minEvidenceMapq) {
            if (hasDistalTargets(readMetadata, minEvidenceMapq)) {
                return Collections.singletonList(new StrandedInterval(target, targetForwardStrand));
            } else {
                return Collections.emptyList();
            }
        }

        /**
         * Finds the coordinates implicated by the read's mate as being part of the breakpoint, ie. the coordinates
         * to the 3' end of the mate, where the breakpoint might lie. Given that we don't have the actual mate read here,
         * we make two small assumptions: first, that the length of the mate is equal to the length of the read we are looking at.
         * Second, since we don't have the mate's CIGAR we can't actually compute the end coordinates of the mate alignment,
         * which might be pushed away from where we think it is by a large indel.
         */
        protected SVInterval getMateTargetInterval(final GATKRead read, final ReadMetadata metadata) {
            final int mateContigIndex = metadata.getContigID(read.getMateContig());
            final int mateStartPosition = read.getMateStart();
            final boolean mateReverseStrand = read.mateIsReverseStrand();

            final int maxAllowableFragmentSize = metadata.getFragmentLengthStatistics(read.getReadGroup()).getMaxNonOutlierFragmentSize();

            final int mateAlignmentLength;
            // if the read has an MC attribute we don't have to assume the aligned read length of the mate
            if (read.hasAttribute("MC")) {
                mateAlignmentLength = TextCigarCodec.decode(read.getAttributeAsString("MC")).getPaddedReferenceLength();
            } else {
                mateAlignmentLength = read.getLength();
            }
            return new SVInterval(mateContigIndex,
                    Math.max(0, mateReverseStrand ?
                                    mateStartPosition - maxAllowableFragmentSize + mateAlignmentLength :
                                    mateStartPosition + mateAlignmentLength - MATE_ALIGNMENT_LENGTH_UNCERTAINTY),
                    mateReverseStrand ?
                            mateStartPosition + MATE_ALIGNMENT_LENGTH_UNCERTAINTY :
                            mateStartPosition + maxAllowableFragmentSize);
        }

        protected boolean getMateForwardStrand(final GATKRead read) {
            return ! read.mateIsReverseStrand();
        }

    }

    @DefaultSerializer(InterContigPair.Serializer.class)
    public static final class InterContigPair extends DiscordantReadPairEvidence {
        private static final int INTER_CONTIG_PAIR_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;

        InterContigPair( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata, INTER_CONTIG_PAIR_WEIGHT);
        }

        private InterContigPair( final Kryo kryo, final Input input ) {
            super(kryo, input);
        }

        /**
         * Directly construct InterContigPair by supplying all fields. Used by testing
         */
        @VisibleForTesting InterContigPair(final SVInterval interval, final int weight,
                                final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                                final boolean validated, final boolean forwardStrand,
                                final String cigarString, final int mappingQuality,
                                final int templateSize, final String readGroup,
                                final SVInterval target, final boolean targetForwardStrand, final int targetQuality) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                  mappingQuality, templateSize, readGroup, target, targetForwardStrand, targetQuality);
        }

        @Override
        public String toString() {
            return super.toString() + "\tIntercontigPair\t" + target;
        }


        @Override
        public boolean equals(Object o) {
            return (this == o) || ((o instanceof InterContigPair) && super.equals(o));
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<InterContigPair> {
            @Override
            public void write( final Kryo kryo, final Output output, final InterContigPair interContigPair ) {
                interContigPair.serialize(kryo, output);
            }

            @Override
            public InterContigPair read( final Kryo kryo, final Input input, final Class<InterContigPair> klass ) {
                return new InterContigPair(kryo, input);
            }
        }
    }

    @DefaultSerializer(OutiesPair.Serializer.class)
    public static final class OutiesPair extends DiscordantReadPairEvidence {
        private static final int OUTIES_PAIR_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;

        OutiesPair( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata, OUTIES_PAIR_WEIGHT);
        }

        private OutiesPair( final Kryo kryo, final Input input ) {
            super(kryo, input);
        }

        /**
         * Directly construct OutiesPair by supplying all fields. Used by testing
         */
        @VisibleForTesting OutiesPair(final SVInterval interval, final int weight,
                           final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                           final boolean validated, final boolean forwardStrand,
                           final String cigarString, final int mappingQuality,
                           final int templateSize, final String readGroup,
                           final SVInterval target, final boolean targetForwardStrand, final int targetQuality) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup, target, targetForwardStrand, targetQuality);
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
        }

        @Override
        public String toString() {
            return super.toString() + "\tOutiesPair\t" + target;
        }

        @Override
        public boolean equals(Object o) {
            return (this == o) || ((o instanceof OutiesPair) && super.equals(o));
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<OutiesPair> {
            @Override
            public void write( final Kryo kryo, final Output output, final OutiesPair mateUnmapped ) {
                mateUnmapped.serialize(kryo, output);
            }

            @Override
            public OutiesPair read( final Kryo kryo, final Input input, final Class<OutiesPair> klass ) {
                return new OutiesPair(kryo, input);
            }
        }
    }

    @DefaultSerializer(SameStrandPair.Serializer.class)
    public static final class SameStrandPair extends DiscordantReadPairEvidence {
        private static final int SAME_STRAND_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;

        SameStrandPair( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata, SAME_STRAND_WEIGHT);
        }

        private SameStrandPair( final Kryo kryo, final Input input ) {
            super(kryo, input);
        }

        /**
         * Directly construct SameStrandPair by supplying all fields. Used by testing
         */
        @VisibleForTesting SameStrandPair(final SVInterval interval, final int weight,
                               final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                               final boolean validated, final boolean forwardStrand,
                               final String cigarString, final int mappingQuality,
                               final int templateSize, final String readGroup,
                               final SVInterval target, final boolean targetForwardStrand, final int targetQuality) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup, target, targetForwardStrand, targetQuality);
        }

        @Override
        public String toString() {
            return super.toString() + "\tSameStrandPair\t" + target;
        }

        @Override
        public boolean equals(Object o) {
            return (this == o) || ((o instanceof SameStrandPair) && super.equals(o));
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SameStrandPair> {
            @Override
            public void write( final Kryo kryo, final Output output, final SameStrandPair sameStrandPair ) {
                sameStrandPair.serialize(kryo, output);
            }

            @Override
            public SameStrandPair read( final Kryo kryo, final Input input, final Class<SameStrandPair> klass ) {
                return new SameStrandPair(kryo, input);
            }
        }
    }

    @DefaultSerializer(WeirdTemplateSize.Serializer.class)
    public static final class WeirdTemplateSize extends DiscordantReadPairEvidence {
        private final int mateStartPosition;
        private final boolean mateReverseStrand;
        private static final int WEIRD_TEMPLATE_SIZE_WEIGHT = ReadEvidence.SINGLE_READ_WEIGHT;

        WeirdTemplateSize( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata, WEIRD_TEMPLATE_SIZE_WEIGHT);
            this.mateStartPosition = read.getMateStart();
            this.mateReverseStrand = read.mateIsReverseStrand();
        }

        private WeirdTemplateSize( final Kryo kryo, final Input input ) {
            super(kryo, input);
            mateStartPosition = input.readInt();
            mateReverseStrand = input.readBoolean();
        }

        /**
         * Directly construct WeirdTemplateSize by supplying all fields. Used by testing
         */
        @VisibleForTesting WeirdTemplateSize(final SVInterval interval, final int weight,
                                  final String templateName, final TemplateFragmentOrdinal fragmentOrdinal,
                                  final boolean validated, final boolean forwardStrand,
                                  final String cigarString, final int mappingQuality,
                                  final int templateSize, final String readGroup,
                                  final SVInterval target, final boolean targetForwardStrand, final int targetQuality) {
            super(interval, weight, templateName, fragmentOrdinal, validated, forwardStrand, cigarString,
                    mappingQuality, templateSize, readGroup, target, targetForwardStrand, targetQuality);
            mateStartPosition = target.getStart();
            mateReverseStrand = !targetForwardStrand;
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(mateStartPosition);
            output.writeBoolean(mateReverseStrand);
        }

        @Override
        public String toString() {
            return super.toString() + "\tTemplateSize\t" + target + "\t" + getTemplateSize();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof WeirdTemplateSize)) return false;
            if (!super.equals(o)) return false;
            WeirdTemplateSize that = (WeirdTemplateSize) o;
            return mateStartPosition == that.mateStartPosition &&
                    mateReverseStrand == that.mateReverseStrand;
        }

        @Override
        public int hashCode() {

            return Objects.hash(super.hashCode(), mateStartPosition, mateReverseStrand);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<WeirdTemplateSize> {
            @Override
            public void write( final Kryo kryo, final Output output, final WeirdTemplateSize weirdTemplateSize ) {
                weirdTemplateSize.serialize(kryo, output);
            }

            @Override
            public WeirdTemplateSize read( final Kryo kryo, final Input input, final Class<WeirdTemplateSize> klass ) {
                return new WeirdTemplateSize(kryo, input);
            }
        }
    }
}
