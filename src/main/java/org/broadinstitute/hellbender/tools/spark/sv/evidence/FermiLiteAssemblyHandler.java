package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMultiMap;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly.Contig;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly.Connection;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;

/** This LocalAssemblyHandler aligns assembly contigs with BWA, along with some optional writing of intermediate results. */
public final class FermiLiteAssemblyHandler implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
    private static final long serialVersionUID = 1L;
    private final String alignerIndexFile;
    private final int maxFastqSize;
    private final String fastqDir;
    private final boolean writeGFAs;
    private final boolean popVariantBubbles;
    private final boolean removeShadowedContigs;
    private final boolean expandAssemblyGraph;

    public FermiLiteAssemblyHandler( final String alignerIndexFile, final int maxFastqSize,
                                     final String fastqDir, final boolean writeGFAs,
                                     final boolean popVariantBubbles, final boolean removeShadowedContigs,
                                     final boolean expandAssemblyGraph ) {
        this.alignerIndexFile = alignerIndexFile;
        this.maxFastqSize = maxFastqSize;
        this.fastqDir = fastqDir;
        this.writeGFAs = writeGFAs;
        this.popVariantBubbles = popVariantBubbles;
        this.removeShadowedContigs = removeShadowedContigs;
        this.expandAssemblyGraph = expandAssemblyGraph;
    }

    @Override
    public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<SVFastqUtils.FastqRead>> intervalAndReads ) {
        final int intervalID = intervalAndReads._1();
        final String assemblyName = AlignedAssemblyOrExcuse.formatAssemblyID(intervalID);
        final List<SVFastqUtils.FastqRead> readsList = intervalAndReads._2();

        // bail if the assembly will be too large
        final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
        if ( fastqSize > maxFastqSize ) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- too big (" + fastqSize + " bytes).");
        }

        // record the reads in the assembly as a FASTQ, if requested
        if ( fastqDir != null ) {
            final String fastqName = String.format("%s/%s.fastq", fastqDir, assemblyName);
            final ArrayList<SVFastqUtils.FastqRead> sortedReads = new ArrayList<>(readsList);
            sortedReads.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            SVFastqUtils.writeFastqFile(fastqName, sortedReads.iterator());
        }

        // assemble the reads
        final FermiLiteAssembler assembler = new FermiLiteAssembler();
        if ( popVariantBubbles ) {
            assembler.setCleaningFlag(0x60);
        }
        final long timeStart = System.currentTimeMillis();
        final FermiLiteAssembly initialAssembly = assembler.createAssembly(readsList);
        final int secondsInAssembly = (int)((System.currentTimeMillis() - timeStart + 500)/1000);
        if ( initialAssembly.getNContigs() == 0 ) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- no contigs produced by assembler.");
        }

        // patch up the assembly to improve contiguity
        final FermiLiteAssembly assembly = reviseAssembly(initialAssembly, removeShadowedContigs, expandAssemblyGraph);

        // record the assembly as a GFA, if requested
        if ( fastqDir != null && writeGFAs ) {
            final String gfaName =  String.format("%s/%s.gfa", fastqDir, assemblyName);
            try ( final OutputStream os = BucketUtils.createFile(gfaName) ) {
                assembly.writeGFA(os);
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write "+gfaName, ioe);
            }
        }

        // align the assembled contigs to the genomic reference
        final AlignedAssemblyOrExcuse result;
        try ( final BwaMemAligner aligner = new BwaMemAligner(BwaMemIndexCache.getInstance(alignerIndexFile)) ) {
            aligner.setIntraCtgOptions();
            final List<byte[]> sequences =
                    assembly.getContigs().stream()
                            .map(Contig::getSequence)
                            .collect(SVUtils.arrayListCollector(assembly.getNContigs()));
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(sequences);
            result = new AlignedAssemblyOrExcuse(intervalID, assembly, secondsInAssembly, alignments);
        }

        return result;
    }

    @VisibleForTesting
    static FermiLiteAssembly reviseAssembly( final FermiLiteAssembly initialAssembly,
                                             final boolean removeShadowedContigs,
                                             final boolean expandAssemblyGraph ) {
        final FermiLiteAssembly unshadowedAssembly =
                removeShadowedContigs ? removeShadowedContigs(initialAssembly) : initialAssembly;
        return expandAssemblyGraph ? expandAssemblyGraph(removeUnbranchedConnections(unshadowedAssembly)) : unshadowedAssembly;
    }

    @VisibleForTesting
    static FermiLiteAssembly removeShadowedContigs( final FermiLiteAssembly assembly ) {
        final int kmerSize = 31;
        final double maxMismatchRate = .05;

        // make a map of assembled kmers
        final int capacity = assembly.getContigs().stream().mapToInt(tig -> tig.getSequence().length - kmerSize + 1).sum();
        final HopscotchMultiMap<SVKmerShort, ContigLocation, KmerLocation> kmerMap = new HopscotchMultiMap<>(capacity);
        assembly.getContigs().forEach(tig -> {
            int contigOffset = 0;
            final Iterator<SVKmer> contigItr = new SVKmerizer(tig.getSequence(), kmerSize, new SVKmerShort());
            while ( contigItr.hasNext() ) {
                final SVKmerShort kmer = (SVKmerShort)contigItr.next();
                final SVKmerShort canonicalKmer = kmer.canonical(kmerSize);
                final ContigLocation location = new ContigLocation(tig, contigOffset++, kmer.equals(canonicalKmer));
                kmerMap.add(new KmerLocation(canonicalKmer, location));
            }
        });

        // remove contigs that vary by a small number of SNPs from a sub-sequence of another contig
        final Set<Contig> contigsToRemove = new HashSet<>();
        assembly.getContigs().forEach(tig -> {
            final Set<ContigLocation> testedLocations = new HashSet<>();
            final byte[] tigBases = tig.getSequence();
            final int maxMismatches = (int)(tigBases.length * maxMismatchRate);
            int tigOffset = 0;
            final SVKmerizer kmerItr = new SVKmerizer(tig.getSequence(), kmerSize, new SVKmerShort(kmerSize));
            while ( kmerItr.hasNext() ) {
                final SVKmerShort contigKmer = (SVKmerShort)kmerItr.next();
                final SVKmerShort canonicalContigKmer = contigKmer.canonical(kmerSize);
                final boolean canonical = contigKmer.equals(canonicalContigKmer);
                final Iterator<KmerLocation> locItr = kmerMap.findEach(canonicalContigKmer);
                while ( locItr.hasNext() ) {
                    final ContigLocation contigLocation = locItr.next().getLocation();
                    final Contig tig2 = contigLocation.getContig();
                    if ( tig == tig2 || contigsToRemove.contains(tig2) ) continue;
                    final byte[] tig2Bases = tig2.getSequence();
                    final boolean isRC = canonical != contigLocation.isCanonical();
                    final int tig2Offset =
                            isRC ? tig2Bases.length - contigLocation.getOffset() - kmerSize : contigLocation.getOffset();
                    if ( tigOffset > tig2Offset ||
                            tigBases.length - tigOffset > tig2Bases.length - tig2Offset ) continue;
                    final int tig2Start = tig2Offset - tigOffset;
                    if ( testedLocations.add(new ContigLocation(tig2, tig2Start, isRC)) ) {
                        int nMismatches = 0;
                        if ( !isRC ) {
                            for ( int idx = 0; idx != tigBases.length; ++idx ) {
                                if ( tigBases[idx] != tig2Bases[tig2Start+idx] ) {
                                    if ( (nMismatches += 1) > maxMismatches ) break;
                                }
                            }
                        } else {
                            final int tig2RCOffset = tig2Bases.length - tig2Start - 1;
                            for ( int idx = 0; idx != tigBases.length; ++idx ) {
                                if ( tigBases[idx] != BaseUtils.simpleComplement(tig2Bases[tig2RCOffset-idx]) ) {
                                    if ( (nMismatches += 1) > maxMismatches ) break;
                                }
                            }
                        }
                        if ( nMismatches <= maxMismatches ) {
                            contigsToRemove.add(tig);
                            break;
                        }
                    }
                }
                tigOffset += 1;
            }
        });

        final List<Contig> contigList = new ArrayList<>(assembly.getContigs().size()-contigsToRemove.size());
        assembly.getContigs().stream()
                .filter(tig -> !contigsToRemove.contains(tig))
                .forEach(contigList::add);

        final Set<Contig> staleConnectionContigs = new HashSet<>(SVUtils.hashMapCapacity(contigsToRemove.size()));
        contigsToRemove.forEach(tig ->
                tig.getConnections().forEach(conn ->
                        staleConnectionContigs.add(conn.getTarget())));

        staleConnectionContigs.forEach(tig -> {
            final List<Connection> connections = new ArrayList<>(tig.getConnections().size() - 1);
            tig.getConnections().stream()
                    .filter(conn -> !contigsToRemove.contains(conn.getTarget()))
                    .forEach(connections::add);
            tig.setConnections(connections);
        });

        return new FermiLiteAssembly(contigList);
    }

    // join contigs that connect without any other branches
    @VisibleForTesting
    static FermiLiteAssembly removeUnbranchedConnections( final FermiLiteAssembly assembly ) {
        final int nContigs = assembly.getNContigs();
        final List<Contig> contigList = new ArrayList<>(nContigs);
        final Set<Contig> examined = new HashSet<>(SVUtils.hashMapCapacity(nContigs));

        // find contigs with a single predecessor contig that has a single successor (or vice versa) and join them
        assembly.getContigs().forEach(tig -> {
            if ( !examined.add(tig) ) return;
            Connection conn;
            Connection conn2;
            while ( (conn = getSolePredecessor(tig)) != null &&
                    !examined.contains(conn.getTarget()) &&
                    (conn2 = getSingletonConnection(conn.getTarget(), !conn.isTargetRC())) != null ) {
                examined.add(conn.getTarget());
                tig = joinContigsWithConnections(tig, conn, conn2);
            }
            while ( (conn = getSoleSuccessor(tig)) != null &&
                    !examined.contains(conn.getTarget()) &&
                    (conn2 = getSingletonConnection(conn.getTarget(), !conn.isTargetRC())) != null ) {
                examined.add(conn.getTarget());
                tig = joinContigsWithConnections(tig, conn, conn2);
            }
            contigList.add(tig);
        });
        return new FermiLiteAssembly(contigList);
    }

    private static Connection getSolePredecessor( final Contig contig ) {
        return getSingletonConnection(contig, true);
    }

    private static Connection getSoleSuccessor( final Contig contig ) {
        return getSingletonConnection(contig, false);
    }

    private static Connection getSingletonConnection( final Contig contig, final boolean isRC ) {
        Connection singleton = null;
        for ( Connection conn : contig.getConnections() ) {
            if ( conn.isRC() == isRC ) {
                if ( singleton == null ) singleton = conn;
                else return null;
            }
        }
        return singleton;
    }

    private static Contig joinContigsWithConnections( final Contig firstContig,
                                                      final Connection connection,
                                                      final Connection rcConnection ) {
        final Contig joinedContig = joinContigs(firstContig, Collections.singletonList(connection));
        final Contig lastContig = connection.getTarget();
        final int capacity = firstContig.getConnections().size() + lastContig.getConnections().size() - 2;
        final List<Connection> connections = new ArrayList<>(capacity);
        for ( final Connection conn : firstContig.getConnections() ) {
            if ( conn != connection ) {
                final Connection newConnection =
                        new Connection(conn.getTarget(), conn.getOverlapLen(), true, conn.isTargetRC());
                replaceConnection(conn.getTarget(), rcConnection(firstContig, conn),
                        rcConnection(joinedContig, newConnection));
                connections.add(newConnection);
            }
        }
        for ( final Connection conn : lastContig.getConnections() ) {
            if ( conn != rcConnection ) {
                final Connection newConnection =
                        new Connection(conn.getTarget(), conn.getOverlapLen(), false, conn.isTargetRC());
                replaceConnection(conn.getTarget(), rcConnection(lastContig, conn),
                        rcConnection(joinedContig, newConnection));
                connections.add(newConnection);
            }
        }
        joinedContig.setConnections(connections);
        return joinedContig;
    }

    private static Connection rcConnection( final Contig contig, final Connection connection ) {
        return new Connection(contig, connection.getOverlapLen(), !connection.isTargetRC(), !connection.isRC());
    }

    private static void replaceConnection( final Contig contig,
                                           final Connection oldConnection,
                                           final Connection newConnection ) {
        final List<Connection> oldConnections = contig.getConnections();
        final List<Connection> newConnections = new ArrayList<>(oldConnections.size());
        for ( final Connection conn : oldConnections ) {
            final Connection toAdd;
            if ( conn.getTarget() == oldConnection.getTarget() &&
                    conn.isRC() == oldConnection.isRC() &&
                    conn.isTargetRC() == oldConnection.isTargetRC() ) {
                toAdd = newConnection;
            } else {
                toAdd = conn;
            }
            newConnections.add(toAdd);
        }
        contig.setConnections(newConnections);
    }

    private static Contig joinContigs( final Contig firstContig, final List<Connection> path ) {
        final int nSupportingReads =
                path.stream()
                        .mapToInt(conn -> conn.getTarget().getNSupportingReads())
                        .reduce(firstContig.getNSupportingReads(), Integer::sum);
        final int newContigLen =
                path.stream()
                        .mapToInt(conn -> conn.getTarget().getSequence().length - conn.getOverlapLen())
                        .reduce(firstContig.getSequence().length, Integer::sum);
        final byte[] sequence = new byte[newContigLen];
        int dstIndex = firstContig.getSequence().length;
        System.arraycopy(firstContig.getSequence(), 0, sequence, 0, dstIndex);
        if ( path.get(0).isRC() )
            SequenceUtil.reverseComplement(sequence, 0, dstIndex);
        for ( final Connection conn : path ) {
            final byte[] contigSequence = conn.getTarget().getSequence();
            final int len = contigSequence.length - conn.getOverlapLen();
            if ( !conn.isTargetRC() ) {
                System.arraycopy(contigSequence, conn.getOverlapLen(), sequence, dstIndex, len);
                final int mismatchCount =
                        getMismatchCount(Arrays.copyOfRange(sequence,dstIndex-conn.getOverlapLen(),dstIndex),
                                         Arrays.copyOfRange(contigSequence,0,conn.getOverlapLen()));
                if ( mismatchCount > 0 ) {
                    throw new GATKException("ends don't match! " + mismatchCount + " mismatches.");
                }
            }
            else {
                System.arraycopy(contigSequence, 0, sequence, dstIndex, len);
                SequenceUtil.reverseComplement(sequence, dstIndex, len);
                final int mismatchCount =
                        getMismatchCount(Arrays.copyOfRange(sequence,dstIndex-conn.getOverlapLen(),dstIndex),
                                         BaseUtils.simpleReverseComplement(
                                            Arrays.copyOfRange(contigSequence,contigSequence.length-conn.getOverlapLen(),contigSequence.length)));
                if ( mismatchCount > 0 ) {
                    throw new GATKException("rc ends don't match! " + mismatchCount + " mismatches.");
                }
            }
            dstIndex += len;
        }
        return new Contig(sequence, null, nSupportingReads);
    }

    private static int getMismatchCount( final byte[] seq1, final byte[] seq2 ) {
        if ( seq1.length != seq2.length ) return -1;
        int count = 0;
        for ( int idx = 0; idx != seq1.length; ++idx ) {
            if ( seq1[idx] != seq2[idx] ) ++count;
        }
        return count;
    }

    // trace all paths, breaking only at cycles and at contigs that require phasing
    @VisibleForTesting
    static FermiLiteAssembly expandAssemblyGraph( final FermiLiteAssembly assembly ) {
        final int nContigs = assembly.getNContigs();
        final List<Contig> contigList = new ArrayList<>(nContigs);
        final Set<ContigStrand> visited = new HashSet<>();
        final Set<Contig> examined = new HashSet<>(SVUtils.hashMapCapacity(nContigs));

        // trace paths from sources and sinks
        assembly.getContigs().forEach(tig -> {
            if ( examined.contains(tig) ) return;
            if ( tig.getConnections().isEmpty() ) {
                contigList.add(tig);
                examined.add(tig);
            } else {
                final int nPredecessors = countPredecessors(tig);
                final int nSuccessors = tig.getConnections().size() - nPredecessors;
                if ( nPredecessors == 0 ) {
                    tracePaths(tig, false, contigList, examined, visited);
                } else if ( nSuccessors == 0 ) {
                    tracePaths(tig, true, contigList, examined, visited);
                }
            }
        });

        // anything not examined must be part of a smooth cycle -- just start anywhere to pick it up
        assembly.getContigs().forEach(tig -> {
            if ( !examined.contains(tig) ) {
                tracePaths(tig, false, contigList, examined, visited);
            }
        });

        return new FermiLiteAssembly(contigList);
    }

    private static int countPredecessors( final Contig contig ) {
        return contig.getConnections().stream().mapToInt(conn -> conn.isRC() ? 1 : 0).sum();
    }

    // called for each connected source and each connected sink that hasn't already been examined
    private static void tracePaths( final Contig contig,
                                    final boolean isRC,
                                    final List<Contig> contigList,
                                    final Set<Contig> examined,
                                    final Set<ContigStrand> visited ) {
        examined.add(contig);
        final LinkedList<Connection> path = new LinkedList<>();
        final ContigStrand contigStrand = new ContigStrand(contig, isRC);
        visited.add(contigStrand);
        for ( Connection connection : contig.getConnections() ) {
            if ( connection.isRC() == isRC && connection.getOverlapLen() >= 0 ) {
                extendPath(contig, connection, path, contigList, examined, visited);
            }
        }
        visited.remove(contigStrand);
    }

    private static void extendPath( final Contig firstContig,
                                    final Connection connection,
                                    final LinkedList<Connection> path,
                                    final List<Contig> contigList,
                                    final Set<Contig> examined,
                                    final Set<ContigStrand> visited ) {
        path.addLast(connection);
        final ContigStrand contigStrand = new ContigStrand(connection.getTarget(), connection.isTargetRC());
        boolean isCycle = !visited.add(contigStrand);
        boolean atEndOfPath = true;
        if ( !isCycle ) {
            final Contig target = connection.getTarget();
            examined.add(target);
            final int nPredecessors = countPredecessors(target);
            final int nSuccessors = target.getConnections().size() - nPredecessors;
            boolean needsPhasing = nPredecessors > 1 && nSuccessors > 1;
            if ( needsPhasing ) {
                tracePaths(target, connection.isTargetRC(), contigList, examined, visited);
            } else {
                for ( Connection conn : target.getConnections() ) {
                    if ( conn.isRC() == connection.isTargetRC() && conn.getOverlapLen() >= 0 ) {
                        extendPath(firstContig, conn, path, contigList, examined, visited);
                        atEndOfPath = false;
                    }
                }
            }
        }
        if ( atEndOfPath ) {
            contigList.add(joinContigs(firstContig, path));
        }
        if ( !isCycle ) visited.remove(contigStrand);
        path.removeLast();
    }

    private static final class ContigLocation {
        private final Contig contig;
        private final int offset;
        private final boolean canonical;
        public ContigLocation(final Contig contig, final int offset, final boolean canonical ) {
            this.contig = contig;
            this.offset = offset;
            this.canonical = canonical;
        }
        public Contig getContig() { return contig; }
        public int getOffset() { return offset; }
        public boolean isCanonical() { return canonical; }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof ContigLocation && equals((ContigLocation)obj);
        }

        public boolean equals( final ContigLocation that ) {
            if ( this == that ) return true;
            return contig == that.contig && offset == that.offset && canonical == that.canonical;
        }

        @Override public int hashCode() {
            return 47 * (47 * (contig.hashCode() + 47 * offset) + (canonical ? 31 : 5));
        }
    }

    private static final class KmerLocation implements Map.Entry<SVKmerShort, ContigLocation> {
        private final SVKmerShort kmer;
        private final ContigLocation location;

        public KmerLocation( final SVKmerShort kmer, final ContigLocation location ) {
            this.kmer = kmer;
            this.location = location;
        }

        public SVKmerShort getKmer() { return kmer; }
        public ContigLocation getLocation() { return location; }
        @Override public SVKmerShort getKey() { return kmer; }
        @Override public ContigLocation getValue() { return location; }
        @Override public ContigLocation setValue(final ContigLocation value ) {
            throw new UnsupportedOperationException("KmerLocation is immutable");
        }
    }

    private static final class ContigStrand {
        private final Contig contig;
        private final boolean isRC;

        public ContigStrand( final Contig contig, final boolean isRC ) {
            this.contig = contig;
            this.isRC = isRC;
        }

        public Contig getContig() { return contig; }
        public boolean isRC() { return isRC; }

        public ContigStrand rc() { return new ContigStrand(contig, !isRC); }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof ContigStrand && equals((ContigStrand) obj);
        }

        public boolean equals( final ContigStrand that ) {
            if ( this == that ) return true;
            return contig == that.contig && isRC == that.isRC;
        }

        @Override public int hashCode() {
            return isRC ? -contig.hashCode() : contig.hashCode();
        }
    }
}
