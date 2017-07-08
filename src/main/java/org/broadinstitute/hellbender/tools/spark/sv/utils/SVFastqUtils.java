package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.TemplateFragmentOrdinal;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Memory-economical utilities for producing a FASTQ file.
 */
public class SVFastqUtils {

    public static final char HEADER_PREFIX_CHR = '@';
    public static final char FRAGMENT_NUMBER_SEPARATOR_CHR = '/';
    public static final char HEADER_FIELD_SEPARATOR_CHR = '\t';
    public static final String HEADER_FIELD_SEPARATOR_REGEXP = "\\t"; // make sure its consistent with {@link #HEADER_FIELD_SEPARATOR_CHR}
    public static final char HEADER_FIELD_EQUAL_CHR = '=';
    private static final char HEADER_FIELD_ARRAY_SEPARATOR_CHR = ';';
    private static final char LINE_SEPARATOR_CHR = '+';
    private static final String MAPPING_FIELD_NAME = "mapping";
    private static final String UNMAPPED_STR = "unmapped";

    private static final String FRAGMENT_NUMBER_SEPARATOR_STR = "" + FRAGMENT_NUMBER_SEPARATOR_CHR;
    private static final String HEADER_FIELD_SEPARATOR_STR = "" + HEADER_FIELD_SEPARATOR_CHR;
    private static final String MAPPING_FIELD_EQUAL_TO = MAPPING_FIELD_NAME + HEADER_FIELD_EQUAL_CHR;
    private static final String UNMAPPED_DESCRIPTION_STR = MAPPING_FIELD_EQUAL_TO + UNMAPPED_STR;
    private static final String HEADER_FIELD_ARRAY_SEPARATOR_STR = "" + HEADER_FIELD_ARRAY_SEPARATOR_CHR;
    private static final Pattern MAPPING_DESCRIPTION_PATTERN = Pattern.compile("^" + MAPPING_FIELD_EQUAL_TO +
            "(.+):(\\d+);(" + Strand.PATTERN.pattern() + ");(.+)$");
    private static final Pattern FASTQ_READ_HEADER_PATTERN = Pattern.compile("^" + HEADER_PREFIX_CHR +
            "([^" + HEADER_FIELD_SEPARATOR_REGEXP +  "]+)(" + HEADER_FIELD_SEPARATOR_REGEXP + "(.*))?$");

    // TODO htsjdk has it own Strand annotation enum. These classes could be merged if that Strand would me updated
    // TODO so that one can get the enum constant char encoding; currently one can only do the transformation the other way.
    public enum Strand {
        POSITIVE('+'),
        NEGATIVE('-');

        public static final Pattern PATTERN = Pattern.compile("\\+|\\-");

        private final char charEncoding;

        Strand(final char ce) {
            charEncoding = ce;
        }

        static Strand decode(final String ce) {
            Utils.nonNull(ce);
            if (ce.length() == 1)
                return decode(ce.charAt(0));
            else
                throw new NoSuchElementException(String.format("there is no strand designation for encoding %s valid encodings are: %s.",
                        ce, Stream.of(values()).map(Strand::encodeAsString).collect(Collectors.joining(", "))));
        }

        static Strand decode(final char ce) {
            if (ce == POSITIVE.charEncoding)
                return POSITIVE;
            else if (ce == NEGATIVE.charEncoding)
                return NEGATIVE;
            else
                throw new NoSuchElementException("there is no strand designation for encoding " + ce + " valid encodings are: " +
                        Stream.of(values()).map(s -> "" + s.charEncoding).collect(Collectors.joining(", ")) + ".");
        }

        @Override
        public String toString() { return encodeAsString(); };

        String encodeAsString() { return "" + charEncoding; }
    }

    public static final class Mapping implements Locatable {
        
        private final String contig;

        private final int start;

        private final boolean forwardStrand;

        private final Cigar cigar;

        public Mapping(final GATKRead read) {
            Utils.nonNull(read);
            if (read.isUnmapped()) {
                contig = null;
                start = -1;
                forwardStrand = true;
                cigar = null;
            } else {
                contig = read.getContig();
                start = read.getStart();
                forwardStrand = !read.isReverseStrand();
                cigar = read.getCigar();
            }
        }

        public Mapping(final String mappingDescription) {
            Utils.nonNull(mappingDescription);
            if (mappingDescription.equals(UNMAPPED_DESCRIPTION_STR)) {
                contig = null;
                start = -1;
                forwardStrand = true;
                cigar = null;
            } else {
                final Matcher matcher = MAPPING_DESCRIPTION_PATTERN.matcher(mappingDescription);
                if (!matcher.find()) {
                    throw new IllegalArgumentException("invalid mapping description: '" + mappingDescription + "'");
                } else {
                    contig = matcher.group(1);
                    start = Integer.parseInt(matcher.group(2));
                    cigar = TextCigarCodec.decode(matcher.group(4));
                    forwardStrand = Strand.decode(matcher.group(3)) == Strand.POSITIVE;
                }
            }
        }

        public boolean isMapped() {
            return contig != null;
        }

        @Override
        public String getContig() {
            return contig;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return start;
        }

        public Cigar getCigar() {
            return cigar;
        }

        public boolean isForwardStrand() {
            return forwardStrand;
        }

        public static String toString(final GATKRead read) {
            Utils.nonNull(read);
            if (read.isUnmapped()) {
                return UNMAPPED_DESCRIPTION_STR;
            } else {
                return MAPPING_FIELD_EQUAL_TO + String.join(HEADER_FIELD_ARRAY_SEPARATOR_STR, read.getContig() + SimpleInterval.CONTIG_SEPARATOR + read.getStart(), (read.isReverseStrand() ? Strand.NEGATIVE : Strand.POSITIVE).encodeAsString(), read.getCigar().toString());
            }
        }

        public String toString() {
            if (isMapped()) {
                return MAPPING_FIELD_EQUAL_TO + String.join(HEADER_FIELD_ARRAY_SEPARATOR_STR, getContig() + SimpleInterval.CONTIG_SEPARATOR + getStart(), (forwardStrand ? Strand.POSITIVE : Strand.NEGATIVE).encodeAsString(), cigar.toString());
            } else {
                return UNMAPPED_DESCRIPTION_STR;
            }
        }
    }

    @DefaultSerializer(FastqRead.Serializer.class)
    public static final class FastqRead implements FermiLiteAssembler.BasesAndQuals {
        private final String header;
        private final byte[] bases;
        private final byte[] quals;

        public FastqRead(final GATKRead read) {
            this(read, true);
        }

        public FastqRead( final GATKRead read, final boolean includeMappingLocation ) {
            this.header = composeHeaderLine(read, includeMappingLocation);
            this.bases = read.getBases();
            this.quals = read.getBaseQualities();
            if (!read.isUnmapped() && read.isReverseStrand()) {
                SequenceUtil.reverseComplement(this.bases);
                SequenceUtil.reverseQualities(this.quals);
            }
        }

        private static String composeHeaderLine(final GATKRead read, final boolean includeMappingLocation) {
            Utils.nonNull(read);
            return HEADER_PREFIX_CHR + read.getName() + TemplateFragmentOrdinal.forRead(read)
                    + (includeMappingLocation ?  (HEADER_FIELD_SEPARATOR_STR + Mapping.toString(read)) : "");
        }

        private FastqRead(final String header, final byte[] bases, final byte[] quals) {
            this.header = header;
            this.bases = bases;
            this.quals = quals;
        }

        private FastqRead( final Kryo kryo, final Input input ) {
            header = input.readString();
            final int nBases = input.readInt();
            bases = new byte[nBases];
            input.readBytes(bases);
            quals = new byte[nBases];
            input.readBytes(quals);
        }

        /**
         * Returns the header line of this Fastq read starting with '@' followed by the read id and description separated by tab characters.
         * @return never {@code null}.
         */
        public String getHeader() {
            return header;
        }

        public String getId() {
            final String[] headerParts = header.split(HEADER_FIELD_SEPARATOR_STR);
            return headerParts[0].substring(1); // skip the '@'.
        }
        public String getName() { final String id = getId();
            final int fragmentNumberSeparatorIndex = id.lastIndexOf(FRAGMENT_NUMBER_SEPARATOR_CHR);
            return fragmentNumberSeparatorIndex >= 0 ? id.substring(0, fragmentNumberSeparatorIndex) : id;
        }
        public String getDescription() {
            final int tabIndex = header.indexOf(HEADER_FIELD_SEPARATOR_CHR);
            return tabIndex >= 0 ? header.substring(tabIndex + 1) : "";
        }

        @Override public byte[] getBases() { return bases; }
        @Override public byte[] getQuals() { return quals; }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeAscii(header);
            output.writeInt(bases.length);
            output.writeBytes(bases);
            output.writeBytes(quals);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<FastqRead> {
            @Override
            public void write( final Kryo kryo, final Output output, final FastqRead read ) {
                read.serialize(kryo, output);
            }

            @Override
            public FastqRead read( final Kryo kryo, final Input input, final Class<FastqRead> type ) {
                return new FastqRead(kryo, input);
            }
        }
    }

    public static List<FastqRead> readFastqFile( final String fileName ) {
        final int INITIAL_CAPACITY = 10000; // absolute guess, just something not too crazy small
        final List<FastqRead> reads = new ArrayList<>(INITIAL_CAPACITY);
        try ( final BufferedReader reader = new BufferedReader(new InputStreamReader(BucketUtils.openFile(fileName))) ) {
            String header;
            int lineNo = 0;
            while ( (header = reader.readLine()) != null ) {
                lineNo += 1;
                if (!FASTQ_READ_HEADER_PATTERN.matcher(header).find()) {
                    throw new GATKException("In FASTQ file "+fileName+" sequence identifier line does not start with @ on line "+lineNo);
                }
                final String callLine = reader.readLine();
                lineNo += 1;
                if ( callLine == null ) {
                    throw new GATKException("In FASTQ file "+fileName+" file truncated: missing calls.");
                }
                final String sepLine = reader.readLine();
                lineNo += 1;
                if ( sepLine == null ) {
                    throw new GATKException("In FASTQ file "+fileName+" file truncated: missing + line.");
                }
                if ( sepLine.length() < 1 || sepLine.charAt(0) != LINE_SEPARATOR_CHR ) {
                    throw new GATKException("In FASTQ file " + fileName + " separator line does not start with + on line " + lineNo);
                }
                final String qualLine = reader.readLine();
                lineNo += 1;
                if ( qualLine == null ) {
                    throw new UserException.BadInput("In FASTQ file "+fileName+" file truncated: missing quals.");
                }
                if ( callLine.length() != qualLine.length() ) {
                    throw new UserException.BadInput("In FASTQ file "+fileName+" there are "+qualLine.length()+
                            " quality scores on line "+lineNo+" but there are "+callLine.length()+" base calls.");
                }
                final byte[] quals = qualLine.getBytes();
                final byte[] calls = callLine.getBytes();
                reads.add(new FastqRead(header, calls, quals));
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Can't read "+fileName, ioe);
        }
        return reads;
    }

    /** Convert a read's name into a FASTQ record sequence ID */
    public static String readToFastqSeqId( final GATKRead read, final boolean includeMappingLocation ) {
        final String nameSuffix = TemplateFragmentOrdinal.forRead(read).nameSuffix();
        final String description;
        if ( includeMappingLocation ) {
            final Mapping mapping = new Mapping(read);
            description = HEADER_FIELD_SEPARATOR_CHR + mapping.toString();
        } else {
            description = "";
        }
        return read.getName() + nameSuffix + HEADER_FIELD_SEPARATOR_CHR + description;
    }

    /** Write a list of FASTQ records into a file. */
    public static void writeFastqFile(
            final String fileName,
            final Iterator<FastqRead> fastqReadItr ) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName)) ) {
            writeFastqStream(writer, fastqReadItr);
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    public static void writeFastqStream( final OutputStream writer, final Iterator<FastqRead> fastqReadItr )
        throws IOException {
        int index = 0;
        while ( fastqReadItr.hasNext() ) {
            final FastqRead read = fastqReadItr.next();
            final String header = read.getHeader();
            if (header.contains(" ")) {
                throw new IllegalStateException("Blank found: " + header);
            }
            if ( header == null ) writer.write(Integer.toString(++index).getBytes());
            else writer.write(header.getBytes());
            writer.write('\n');
            writer.write(read.getBases());
            writer.write('\n');
            writer.write(SVFastqUtils.LINE_SEPARATOR_CHR);
            writer.write('\n');
            final byte[] quals = read.getQuals();
            final int nQuals = quals.length;
            final byte[] fastqQuals = new byte[nQuals];
            for ( int idx = 0; idx != nQuals; ++idx ) {
                fastqQuals[idx] = (byte)SAMUtils.phredToFastq(quals[idx]);
            }
            writer.write(fastqQuals);
            writer.write('\n');
        }
    }
}
