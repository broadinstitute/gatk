package org.broadinstitute.hellbender.tools.spark.sv;

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
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Memory-economical utilities for producing a FASTQ file.
 */
public class SVFastqUtils {

    public static final class SATag implements Locatable {
        private final String contig;
        private final int start;
        private final boolean forwardStrand;
        private final Cigar cigar;
        private final int mapQ;
        private final int nMismatches;

        public SATag( final GATKRead read ) {
            this.contig = read.getContig();
            this.start = read.getStart();
            this.forwardStrand = !read.isReverseStrand();
            this.cigar = read.getCigar();
            this.mapQ = read.getMappingQuality();
            final Integer nm = read.getAttributeAsInteger("NM");
            this.nMismatches = nm == null ? -1 : nm;
        }

        public SATag( final String contig, final int start, final boolean forwardStrand, final Cigar cigar,
                      final int mapQ, final int nMismatches ) {
            this.contig = contig;
            this.start = start;
            this.forwardStrand = forwardStrand;
            this.cigar = cigar;
            this.mapQ = mapQ;
            this.nMismatches = nMismatches;
        }

        @Override public String getContig() { return contig; }
        @Override public int getStart() { return start; }
        @Override public int getEnd() { return start + cigar.getReferenceLength() - 1; }
        public boolean isForwardStrand() { return forwardStrand; }
        public Cigar getCigar() { return cigar; }
        public int getMapQ() { return mapQ; }
        public int getNMismatches() { return nMismatches; }

        @Override public String toString() {
            return contig+","+start+(forwardStrand ? ",+," : ",-,")+cigar+","+mapQ+","+nMismatches+";";
        }

        public static SATag fromString( final String saString ) {
            final String[] tokens = saString.split(",");
            try {
                return new SATag(tokens[0], Integer.parseInt(tokens[1]), tokens[2].charAt(0)=='+',
                        TextCigarCodec.decode(tokens[3]), Integer.parseInt(tokens[4]),
                        Integer.parseInt(tokens[5].substring(0, tokens[5].length() - 1)));
            } catch ( final Exception e ) {
                throw new GATKException("Not in SA tag format: "+saString, e);
            }
        }
    }

    @DefaultSerializer(FastqRead.Serializer.class)
    public static final class FastqRead implements FermiLiteAssembler.BasesAndQuals {
        private final String name;
        private final byte[] bases;
        private final byte[] quals;

        public FastqRead( final GATKRead read, final boolean includeMappingLocation ) {
            name = readToFastqSeqId(read, includeMappingLocation);
            if ( read.isReverseStrand() ) {
                bases = BaseUtils.simpleReverseComplement(read.getBases());
                quals = read.getBaseQualities();
                SequenceUtil.reverseQualities(quals);
            } else {
                bases = read.getBases();
                quals = read.getBaseQualities();
            }
        }

        public FastqRead( final String name, final byte[] bases, final byte[] quals ) {
            this.name = name;
            this.bases = bases;
            this.quals = quals;
        }

        private FastqRead( final Kryo kryo, final Input input ) {
            name = input.readString();
            final int nBases = input.readInt();
            bases = new byte[nBases];
            input.readBytes(bases);
            quals = new byte[nBases];
            input.readBytes(quals);
        }

        public String getName() { return name; }
        @Override public byte[] getBases() { return bases; }
        @Override public byte[] getQuals() { return quals; }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeAscii(name);
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
            String seqIdLine;
            int lineNo = 0;
            while ( (seqIdLine = reader.readLine()) != null ) {
                lineNo += 1;
                if ( seqIdLine.length() < 1 || seqIdLine.charAt(0) != '@' ) {
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
                if ( sepLine.length() < 1 || sepLine.charAt(0) != '+' ) {
                    throw new GATKException("In FASTQ file " + fileName + " separator line does not start with + on line " + lineNo);
                }
                final String qualLine = reader.readLine();
                lineNo += 1;
                if ( qualLine == null ) {
                    throw new GATKException("In FASTQ file "+fileName+" file truncated: missing quals.");
                }
                if ( callLine.length() != qualLine.length() ) {
                    throw new GATKException("In FASTQ file "+fileName+" there are "+qualLine.length()+
                            " quality scores on line "+lineNo+" but there are "+callLine.length()+" base calls.");
                }
                final byte[] quals = qualLine.getBytes();
                SAMUtils.fastqToPhred(quals);
                reads.add(new FastqRead(seqIdLine.substring(1), callLine.getBytes(), quals));
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Can't read "+fileName, ioe);
        }
        return reads;
    }

    /** Convert a read into a FASTQ record sequence ID */
    public static String readToFastqSeqId( final GATKRead read, final boolean includeMappingLocation ) {
        final String nameSuffix = read.isPaired() ? (read.isFirstOfPair() ? "/1" : "/2") : "";
        String mapLoc = "";
        if ( includeMappingLocation && !read.isUnmapped() ) {
            mapLoc = " mapping=" + new SATag(read);
        }
        return read.getName() + nameSuffix + mapLoc;
    }

    /** Convert a FASTQ record sequence ID into an SA tag, if present. */
    public static SATag seqIdToSATag( final String seqId ) {
        final int tagPos = seqId.indexOf("mapping=");
        if ( tagPos == -1 ) return null;
        return SATag.fromString(seqId.substring(tagPos + 8));
    }

    /** Extract template name from FASTQ sequence ID */
    public static String seqIdToTemplate( final String seqId ) {
        int tokLen = seqId.indexOf(' ');
        if ( tokLen == -1 ) tokLen = seqId.length();
        if ( tokLen > 2 && seqId.charAt(tokLen-2)=='/' ) tokLen -= 2;
        return seqId.substring(0, tokLen);
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
            writer.write('@');
            final String name = read.getName();
            if ( name == null ) writer.write(Integer.toString(++index).getBytes());
            else writer.write(name.getBytes());
            writer.write('\n');
            writer.write(read.getBases());
            writer.write('\n');
            writer.write('+');
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
