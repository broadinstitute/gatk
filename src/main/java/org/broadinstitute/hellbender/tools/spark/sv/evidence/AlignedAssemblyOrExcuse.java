package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextWriter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly.Connection;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly.Contig;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * An assembly with its contigs aligned to reference, or a reason that there isn't an assembly.
 */
@DefaultSerializer(AlignedAssemblyOrExcuse.Serializer.class)
public final class AlignedAssemblyOrExcuse {
    private final int assemblyId;
    private final String errorMessage;
    private final FermiLiteAssembly assembly;
    private final List<List<BwaMemAlignment>> contigAlignments;

    public AlignedAssemblyOrExcuse( final int assemblyId, final String errorMessage ) {
        this.assemblyId = assemblyId;
        this.errorMessage = errorMessage;
        this.assembly = null;
        this.contigAlignments = null;
    }

    public AlignedAssemblyOrExcuse( final int assemblyId, final FermiLiteAssembly assembly,
                                    final List<List<BwaMemAlignment>> contigAlignments ) {
        Utils.validate(assembly.getNContigs()==contigAlignments.size(),
                "Number of contigs in assembly doesn't match length of list of alignments.");
        Utils.validateArg(assembly.getContigs().stream().noneMatch(contig -> contig.getConnections()==null),
                "Some assembly has contigs that have null connections");
        this.assemblyId = assemblyId;
        this.errorMessage = null;
        this.assembly = assembly;
        this.contigAlignments = contigAlignments;
    }

    private AlignedAssemblyOrExcuse( final Kryo kryo, final Input input ) {
        this.assemblyId = input.readInt();
        this.errorMessage = input.readString();
        if ( errorMessage != null ) {
            this.assembly = null;
            this.contigAlignments = null;
        } else {
            final int nContigs = input.readInt();
            final List<Contig> contigs = new ArrayList<>(nContigs);
            for ( int idx = 0; idx != nContigs; ++idx ) {
                contigs.add(readContig(input));
            }
            for ( int idx = 0; idx != nContigs; ++idx ) {
                final int nConnections = input.readInt();
                if ( nConnections == 0 ) continue;
                final List<Connection> connections = new ArrayList<>(nConnections);
                for ( int connIdx = 0; connIdx != nConnections; ++connIdx ) {
                    connections.add(readConnection(input, contigs));
                }
                contigs.get(idx).setConnections(connections);
            }
            this.assembly = new FermiLiteAssembly(contigs);
            final List<List<BwaMemAlignment>> contigAlignments = new ArrayList<>(nContigs);
            for ( int idx = 0; idx != nContigs; ++idx ) {
                final int nAlignments = input.readInt();
                final List<BwaMemAlignment> alignments = new ArrayList<>(nAlignments);
                for ( int alnIdx = 0; alnIdx != nAlignments; ++alnIdx ) {
                    alignments.add(readAlignment(input));
                }
                contigAlignments.add(alignments);
            }
            this.contigAlignments = contigAlignments;
        }
    }

    public int getAssemblyId() {
        return assemblyId;
    }

    /**
     * Either this is null, or the assembly and list of alignments is null.
     */
    public String getErrorMessage() {
        return errorMessage;
    }

    public boolean isNotFailure() {
        return errorMessage==null;
    }

    public FermiLiteAssembly getAssembly() {
        return assembly;
    }

    /**
     * List is equal in length to the number of contigs in the assembly.
     */
    public List<List<BwaMemAlignment>> getContigAlignments() {
        return contigAlignments;
    }

    /**
     * write a SAM file containing records for each aligned contig
     */
    static void writeSAMFile( final String samFile,
                                     final SAMFileHeader header,
                                     final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList ) {
        try ( final OutputStream os = BucketUtils.createFile(samFile) ) {
            final SAMTextWriter writer = new SAMTextWriter(os);
            writer.setSortOrder(SAMFileHeader.SortOrder.queryname, true);
            writer.setHeader(header);

            final List<String> refNames = getRefNames(header);
            alignedAssemblyOrExcuseList.stream()
                    .filter(AlignedAssemblyOrExcuse::isNotFailure)
                    .flatMap(aa -> aa.toSAMStreamForAlignmentsOfThisAssembly(header,refNames))
                    .forEach(writer::addAlignment);

            writer.finish();
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write SAM file of aligned contigs.", ioe);
        }
    }

    private Stream<SAMRecord> toSAMStreamForAlignmentsOfThisAssembly(final SAMFileHeader header, final List<String> refNames) {
        Utils.validate(isNotFailure(), "Can't stream SAM records from a failed assembly.");
        return IntStream.range(0, contigAlignments.size()).boxed()
                .flatMap(contigIdx ->
                        toSAMStreamForOneContig(header,refNames, assemblyId,contigIdx,
                                assembly.getContig(contigIdx).getSequence(),contigAlignments.get(contigIdx)));
    }

    public static List<String> getRefNames(final SAMFileHeader header) {
        return header.getSequenceDictionary().getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList());
    }

    public static Stream<SAMRecord> toSAMStreamForOneContig(final SAMFileHeader header, final List<String> refNames,
                                                            final int assemblyId, final int contigIdx,
                                                            final byte[] contigSequence, final List<BwaMemAlignment> alignments) {
        if ( alignments.isEmpty() ) return Stream.empty();

        final String readName = formatContigName(assemblyId, contigIdx);
        final Map<BwaMemAlignment,String> saTagMap = BwaMemAlignmentUtils.createSATags(alignments,refNames);

        return alignments.stream()
                .map(alignment -> {
                    final SAMRecord samRecord =
                            BwaMemAlignmentUtils.applyAlignment(readName, contigSequence, null, null, alignment,
                                    refNames, header, false, false);
                    final String saTag = saTagMap.get(alignment);
                    if ( saTag != null ) samRecord.setAttribute("SA", saTag);
                    return samRecord;
                });
    }

    public static String formatContigName(final int assemblyId, final int contigIdx) {
        return formatAssemblyID(assemblyId) + ":" + formatContigID(contigIdx);
    }

    static String formatAssemblyID(final int assemblyId) {
        return String.format("asm%06d", assemblyId);
    }

    static String formatContigID(final int contigIdx) {
        return String.format("tig%05d", contigIdx);
    }

    /**
     * write a file describing each interval
     */
    public static void writeIntervalFile( final String intervalFile,
                                          final SAMFileHeader header,
                                          final List<SVInterval> intervals,
                                          final List<AlignedAssemblyOrExcuse> intervalDispositions ) {
        final Map<Integer, AlignedAssemblyOrExcuse> resultsMap = new HashMap<>();
        intervalDispositions.forEach(alignedAssemblyOrExcuse ->
                resultsMap.put(alignedAssemblyOrExcuse.getAssemblyId(), alignedAssemblyOrExcuse));

        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(intervalFile))) ) {
            final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
            final int nIntervals = intervals.size();
            for ( int intervalId = 0; intervalId != nIntervals; ++intervalId ) {
                final SVInterval interval = intervals.get(intervalId);
                final String seqName = contigs.get(interval.getContig()).getSequenceName();
                final AlignedAssemblyOrExcuse alignedAssemblyOrExcuse = resultsMap.get(intervalId);
                final String disposition;
                if ( alignedAssemblyOrExcuse == null ) {
                    disposition = "unknown";
                } else if ( alignedAssemblyOrExcuse.getErrorMessage() != null ) {
                    disposition = alignedAssemblyOrExcuse.getErrorMessage();
                } else {
                    disposition = "produced " + alignedAssemblyOrExcuse.getAssembly().getNContigs() + " contigs";
                }
                writer.write(intervalId + "\t" +
                        seqName + ":" + interval.getStart() + "-" + interval.getEnd() + "\t" +
                        disposition + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write intervals file " + intervalFile, ioe);
        }
    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(assemblyId);
        output.writeString(errorMessage);
        if ( errorMessage == null ) {
            final int nContigs = assembly.getNContigs();
            final Map<Contig, Integer> contigMap = new HashMap<>();
            output.writeInt(nContigs);
            for ( int idx = 0; idx != nContigs; ++idx ) {
                final Contig contig = assembly.getContig(idx);
                writeContig(contig, output);
                contigMap.put(contig, idx);
            }
            for ( final Contig contig : assembly.getContigs() ) {
                final List<Connection> connections = contig.getConnections();
                output.writeInt(connections.size());
                for ( final Connection connection : connections ) {
                    writeConnection(connection, contigMap, output);
                }
            }
            for ( final List<BwaMemAlignment> alignments : contigAlignments ) {
                output.writeInt(alignments.size());
                for ( final BwaMemAlignment alignment : alignments ) {
                    writeAlignment(alignment, output);
                }
            }
        }
    }

    private static void writeContig( final Contig contig, final Output output ) {
        output.writeInt(contig.getSequence().length);
        output.writeBytes(contig.getSequence());
        output.writeBytes(contig.getPerBaseCoverage());
        output.writeInt(contig.getNSupportingReads());
    }

    private static Contig readContig( final Input input ) {
        final int sequenceLen = input.readInt();
        final byte[] sequence = new byte[sequenceLen];
        input.readBytes(sequence);
        final byte[] perBaseCoverage = new byte[sequenceLen];
        input.readBytes(perBaseCoverage);
        final int nSupportingReads = input.readInt();
        return new Contig(sequence, perBaseCoverage, nSupportingReads);
    }

    private static void writeConnection( final Connection connection,
                                         final Map<Contig, Integer> contigMap,
                                         final Output output ) {
        output.writeInt(contigMap.get(connection.getTarget()));
        output.writeInt(connection.getOverlapLen());
        output.writeBoolean(connection.isRC());
        output.writeBoolean(connection.isTargetRC());
    }

    private static Connection readConnection( final Input input, final List<Contig> contigs ) {
        final Contig target = contigs.get(input.readInt());
        final int overlapLen = input.readInt();
        final boolean isRC = input.readBoolean();
        final boolean isTargetRC = input.readBoolean();
        return new Connection(target, overlapLen, isRC, isTargetRC);
    }

    private static void writeAlignment(final BwaMemAlignment alignment, final Output output) {
        output.writeInt(alignment.getSamFlag());
        output.writeInt(alignment.getRefId());
        output.writeInt(alignment.getRefStart());
        output.writeInt(alignment.getRefEnd());
        output.writeInt(alignment.getSeqStart());
        output.writeInt(alignment.getSeqEnd());
        output.writeInt(alignment.getMapQual());
        output.writeInt(alignment.getNMismatches());
        output.writeInt(alignment.getAlignerScore());
        output.writeInt(alignment.getSuboptimalScore());
        output.writeString(alignment.getCigar());
        output.writeString(alignment.getMDTag());
        output.writeString(alignment.getXATag());
        output.writeInt(alignment.getMateRefId());
        output.writeInt(alignment.getMateRefStart());
        output.writeInt(alignment.getTemplateLen());
    }

    private static BwaMemAlignment readAlignment(final Input input) {
        final int samFlag = input.readInt();
        final int refId = input.readInt();
        final int refStart = input.readInt();
        final int refEnd = input.readInt();
        final int seqStart = input.readInt();
        final int seqEnd = input.readInt();
        final int mapQual = input.readInt();
        final int nMismatches = input.readInt();
        final int alignerScore = input.readInt();
        final int suboptimalScore = input.readInt();
        final String cigar = input.readString();
        final String mdTag = input.readString();
        final String xaTag = input.readString();
        final int mateRefId = input.readInt();
        final int mateRefStart = input.readInt();
        final int templateLen = input.readInt();
        return new BwaMemAlignment(samFlag, refId, refStart, refEnd, seqStart, seqEnd,
                mapQual, nMismatches, alignerScore, suboptimalScore,
                cigar, mdTag, xaTag,
                mateRefId, mateRefStart, templateLen);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedAssemblyOrExcuse> {
        @Override
        public void write( final Kryo kryo, final Output output, final AlignedAssemblyOrExcuse alignedAssemblyOrExcuse ) {
            alignedAssemblyOrExcuse.serialize(kryo, output);
        }

        @Override
        public AlignedAssemblyOrExcuse read( final Kryo kryo, final Input input, final Class<AlignedAssemblyOrExcuse> klass ) {
            return new AlignedAssemblyOrExcuse(kryo, input);
        }
    }
}
