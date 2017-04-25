package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Various types of read anomalies that provide evidence of genomic breakpoints.
 * There is a shallow hierarchy based on this class that classifies the anomaly as, for example, a split read,
 *   or a pair that has both reads on the same reference strand.
 * Each BreakpointEvidence object comes from examining a single read, and describing its funkiness, if any, by
 *   instantiating one of the subclasses that addresses that type of funkiness.
 */
public class BreakpointEvidence {
    private final SVInterval eventLocation;
    private final String templateName; // QNAME of the read that was funky (i.e., the name of the fragment)
    private final TemplateEnd templateEnd; // which read we're talking about (first or last, for paired-end reads)

    // Typically UNPAIRED (single read from the template), PAIRED_FIRST (first read from the template), or
    // PAIRED_LAST (second read from the template).  The SAM format, however, describes other weird possibilities,
    // and to avoid lying, we also allow the pairedness to be unknown, or for a read to be paired, but neither
    // first nor last (interior).
    public enum TemplateEnd {
        UNPAIRED(""), PAIRED_UNKNOWN("/?"), PAIRED_FIRST("/1"), PAIRED_SECOND("/2"), PAIRED_INTERIOR("/0");

        TemplateEnd( final String value ) { this.value = value; }

        @Override
        public String toString() { return value; }

        private final String value;
    }

    /**
     * evidence offset and width is set to "the rest of the fragment" not covered by this read
     */
    protected BreakpointEvidence( final GATKRead read, final ReadMetadata metadata ) {
        final int templateLen = metadata.getStatistics(read.getReadGroup()).getMedianFragmentSize();
        int width;
        int start;
        if ( read.isReverseStrand() ) {
            final int readStart = read.getStart();
            width = readStart - (read.getUnclippedEnd() + 1 - templateLen);
            start = readStart - width;
            if ( start < 1 ) {
                width += start - 1;
                start = 1;
            }
        } else {
            final int readEnd = read.getEnd() + 1;
            width = read.getUnclippedStart() + templateLen - readEnd;
            start = readEnd;
        }
        this.eventLocation = new SVInterval(metadata.getContigID(read.getContig()), start, start+width);
        this.templateName = read.getName();
        if ( templateName == null ) throw new GATKException("Read has no name.");
        this.templateEnd = findTemplateEnd(read);
    }

    /**
     * for use when the uncertainty in location has a fixed size
     */
    protected BreakpointEvidence( final GATKRead read, final ReadMetadata metadata,
                        final int contigOffset, final int offsetUncertainty ) {
        int width = 2*offsetUncertainty;
        int start = contigOffset - offsetUncertainty;
        if ( start < 1 ) {
            width += start - 1;
            start = 1;
        }
        this.eventLocation = new SVInterval(metadata.getContigID(read.getContig()), start, start+width);
        this.templateName = read.getName();
        if ( templateName == null ) throw new GATKException("Read has no name.");
        this.templateEnd = findTemplateEnd(read);
    }

    /**
     * a technical constructor for use in Kryo (de-)serialization.
     * this creates an object by reading a Kryo-serialized stream.
     * it will be called by subclasses in their own constructors from Kryo streams (as super(kryo, input)).
     */
    protected BreakpointEvidence( final Kryo kryo, final Input input ) {
        final int contig = input.readInt();
        final int start = input.readInt();
        final int end = input.readInt();
        this.eventLocation = new SVInterval(contig, start, end);
        this.templateName = input.readString();
        this.templateEnd = TemplateEnd.values()[input.readByte()];
    }

    /**
     * to make a sentinel (a bit of evidence that serves no function other than to mark the end of a stream).
     * used by the MapPartitioner to flush pending evidence in the FindBreakEvidenceSpark.WindowSorter.
     */
    public BreakpointEvidence( final int contigIndex ) {
        this.eventLocation = new SVInterval(contigIndex, 0, 0);
        this.templateName = "sentinel";
        this.templateEnd = TemplateEnd.PAIRED_UNKNOWN;
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(eventLocation.getContig());
        output.writeInt(eventLocation.getStart());
        output.writeInt(eventLocation.getEnd());
        output.writeString(templateName);
        output.writeByte(templateEnd.ordinal());
    }

    public SVInterval getLocation() { return eventLocation; }
    public String getTemplateName() { return templateName; }
    public TemplateEnd getTemplateEnd() { return templateEnd; }

    @Override
    public String toString() {
        return eventLocation+" "+templateName + templateEnd;
    }

    private static TemplateEnd findTemplateEnd( final GATKRead read ) {
        return !read.isPaired() ? TemplateEnd.UNPAIRED :
               !read.isFirstOfPair() && !read.isSecondOfPair() ? TemplateEnd.PAIRED_UNKNOWN :
                read.isFirstOfPair() && !read.isSecondOfPair() ? TemplateEnd.PAIRED_FIRST :
               !read.isFirstOfPair() &&  read.isSecondOfPair() ? TemplateEnd.PAIRED_SECOND :
                       TemplateEnd.PAIRED_INTERIOR;
    }

    @DefaultSerializer(SplitRead.Serializer.class)
    public static final class SplitRead extends BreakpointEvidence {
        private static final int UNCERTAINTY = 2;
        private static final String SA_TAG_NAME = "SA";
        private final String cigar;
        private final String tagSA;

        public SplitRead( final GATKRead read, final ReadMetadata metadata, final boolean atStart ) {
            super(read, metadata, atStart ? read.getStart() : read.getEnd(), UNCERTAINTY);
            cigar = read.getCigar().toString();
            if ( cigar == null ) throw new GATKException("Read has no cigar string.");
            tagSA = read.getAttributeAsString(SA_TAG_NAME);
        }

        private SplitRead( final Kryo kryo, final Input input ) {
            super(kryo, input);
            cigar = input.readString();
            tagSA = input.readString();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeString(cigar);
            output.writeString(tagSA);
        }

        @Override
        public String toString() {
            return super.toString() + " Split " + cigar + (tagSA == null ? " SA: None" : (" SA: " + tagSA));
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
    }

    @DefaultSerializer(LargeIndel.Serializer.class)
    public static final class LargeIndel extends BreakpointEvidence {
        private static final int UNCERTAINTY = 4;
        private final String cigar;

        LargeIndel( final GATKRead read, final ReadMetadata metadata, final int contigOffset ) {
            super(read, metadata, contigOffset, UNCERTAINTY);
            cigar = read.getCigar().toString();
            if ( cigar == null ) throw new GATKException("Read has no cigar string.");
        }

        private LargeIndel( final Kryo kryo, final Input input ) {
            super(kryo, input);
            cigar = input.readString();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeString(cigar);
        }

        @Override
        public String toString() {
            return super.toString() + " Indel " + cigar;
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
    public static final class MateUnmapped extends BreakpointEvidence {

        MateUnmapped( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata);
        }

        private MateUnmapped( final Kryo kryo, final Input input ) { super(kryo, input); }

        @Override
        public String toString() {
            return super.toString() + " UnmappedMate";
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

    @DefaultSerializer(InterContigPair.Serializer.class)
    public static final class InterContigPair extends BreakpointEvidence {
        private final int mateContigIndex;
        private final int mateStartPosition;

        InterContigPair( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata);
            this.mateContigIndex = metadata.getContigID(read.getMateContig());
            this.mateStartPosition = read.getMateStart();
        }

        private InterContigPair( final Kryo kryo, final Input input ) {
            super(kryo, input);
            mateContigIndex = input.readInt();
            mateStartPosition = input.readInt();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(mateContigIndex);
            output.writeInt(mateStartPosition);
        }

        @Override
        public String toString() {
            return super.toString() + " IntercontigPair " + mateContigIndex;
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
    public static final class OutiesPair extends BreakpointEvidence {
        OutiesPair( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata);
        }

        private OutiesPair( final Kryo kryo, final Input input ) { super(kryo, input); }

        @Override
        public String toString() {
            return super.toString() + " OutiesPair";
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
    public static final class SameStrandPair extends BreakpointEvidence {
        SameStrandPair( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata);
        }

        private SameStrandPair( final Kryo kryo, final Input input ) { super(kryo, input); }

        @Override
        public String toString() {
            return super.toString() + " SameStrandPair";
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
    public static final class WeirdTemplateSize extends BreakpointEvidence {
        private final int templateSize;

        WeirdTemplateSize( final GATKRead read, final ReadMetadata metadata ) {
            super(read, metadata);
            this.templateSize = read.getFragmentLength();
        }

        private WeirdTemplateSize( final Kryo kryo, final Input input ) {
            super(kryo, input);
            templateSize = input.readInt();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(templateSize);
        }

        @Override
        public String toString() {
            return super.toString() + " TemplateSize " + templateSize;
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
