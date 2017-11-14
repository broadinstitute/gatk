package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer.Base;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils.FastqRead;
import scala.Tuple2;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) junk", summary = "complete crap",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class KmerAdjacencyBuilder extends CommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "input fastq",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME)
    private String fastqFile;

    @Argument(doc = "graph output",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String graphOutput;

    @Argument(doc = "kmer size", fullName = "kSize", optional = true)
    private int kSizeArg = 63;

    @Argument(doc = "minimum quality", fullName = "minQ", optional = true)
    private int minQArg = 10;

    @Argument(doc = "minimum kmer count", fullName = "minKCount", optional = true)
    private int minKCountArg = 4;

    @Override protected Object doWork() {
        final List<FastqRead> reads = SVFastqUtils.readFastqFile(fastqFile);
        final int kSize = kSizeArg;
        final int minQ = minQArg;
        final int minKCount = minKCountArg;

        // kmerize each read, counting the observations of each kmer.
        // trim each read at the first base having a quality less than minQ
        final int nKmers = reads.stream().mapToInt(read -> Math.min(0, read.getBases().length - kSize + 1)).sum();
        final Map<SVKmerLong, Integer> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        for ( final FastqRead read : reads ) {
            SVKmerizer.canonicalStream(trimmedRead(read, minQ), kSize, new SVKmerLong(kSize))
                    .forEach(kmer -> kmerCounts.merge((SVKmerLong)kmer, 1, Integer::sum));
        }

        // dump a histogram of kmer counts
        //final Map<Integer, Integer> kmerCountHistogram = new TreeMap<>();
        //kmerCounts.values().forEach( count -> kmerCountHistogram.merge(count, 1, Integer::sum));
        //kmerCountHistogram.forEach( (count, countCount) -> System.out.println(count + "\t" + countCount));

        // ignore kmers that appear less than minKCount times
        kmerCounts.entrySet().removeIf(entry -> entry.getValue() < minKCount);

        // build 61-mer adjacency map from 63-mers
        final Map<SVKmerLong, int[]> kmerAdjacencyMap = new HashMap<>(SVUtils.hashMapCapacity(kmerCounts.size()));
        kmerCounts.forEach((kmer, value) -> {
            final SVKmerLong kkk = kmer.removeFirstAndLastBase(kSize);
            final int[] counts = kmerAdjacencyMap.computeIfAbsent(kkk, key -> new int[8]);
            counts[kmer.firstBase(kSize).ordinal()] += value;
            counts[kmer.lastBase().ordinal() + 4] += value;
        });

        // the adjacency map won't include source and sink kmers
        // add them to the map now
        final int kSize2 = kSize - 2;
        final List<Tuple2<SVKmerLong, int[]>> sourcesAndSinks = new ArrayList<>();
        kmerAdjacencyMap.forEach((kmer, counts) -> {
            for ( final Base base : Base.values() ) {
                final int predCount = counts[base.ordinal()];
                if ( predCount > 0 ) {
                    final SVKmerLong predecessor = kmer.predecessor(base, kSize2);
                    final SVKmerLong canonicalPredecessor = predecessor.canonical(kSize2);
                    final int idx;
                    if ( predecessor.equals(canonicalPredecessor) ) { // if predecessor is in canonical form
                        idx = kmer.lastBase().ordinal() + 4;
                    }
                    else {
                        idx = kmer.lastBase().complement().ordinal();
                    }
                    final int[] oldCounts = kmerAdjacencyMap.get(canonicalPredecessor);
                    if ( oldCounts == null ) {
                        final int[] newCounts = new int[8];
                        newCounts[idx] = predCount;
                        sourcesAndSinks.add(new Tuple2<>(canonicalPredecessor, newCounts));
                    } else if ( oldCounts[idx] == 0 ) {
                        oldCounts[idx] = predCount;
                    }
                }
                final int succCount = counts[base.ordinal() + 4];
                if ( succCount > 0 ) {
                    final SVKmerLong successor = kmer.successor(base, kSize2);
                    final SVKmerLong canonicalSuccessor = successor.canonical(kSize2);
                    final int idx;
                    if ( successor.equals(canonicalSuccessor) ) { // if successor is in canonical form
                        idx = kmer.firstBase(kSize2).ordinal();
                    } else {
                        idx = kmer.firstBase(kSize2).complement().ordinal() + 4;
                    }
                    final int[] oldCounts = kmerAdjacencyMap.get(canonicalSuccessor);
                    if ( oldCounts == null ) {
                        final int[] newCounts = new int[8];
                        newCounts[idx] = succCount;
                        sourcesAndSinks.add(new Tuple2<>(canonicalSuccessor, newCounts));
                    } else if ( oldCounts[idx] == 0 ) {
                        oldCounts[idx] = succCount;
                    }
                }
            }
        });
        sourcesAndSinks.forEach(tuple2 ->
                kmerAdjacencyMap.merge(tuple2._1(), tuple2._2(), (oldCounts, newCounts) -> {
                    for ( int idx = 0; idx != oldCounts.length; ++idx ) {
                        oldCounts[idx] += newCounts[idx];
                    }
                    return oldCounts;
                }));

        // dump the adjacency map
        //for ( final Map.Entry<SVKmerLong, int[]> entry : kmerAdjacencyMap.entrySet() ) {
        //    final StringBuilder sb = new StringBuilder(entry.getKey().toString(kSize2));
        //    Arrays.stream(entry.getValue()).forEach(iii -> sb.append('\t').append(iii));
        //    System.out.println(sb);
        //}

        // build contigs
        final List<Contig> contigs = new ArrayList<>();
        final Map<SVKmerLong, ContigEnd> contigEnds = new HashMap<>();
        final Map<SVKmerLong, ContigLocus> locusMap = new HashMap<>(SVUtils.hashMapCapacity(kmerAdjacencyMap.size()));
        kmerAdjacencyMap.forEach( (kmer, counts) -> {
            if ( !contigEnds.containsKey(kmer) ) {
                final int contigId = contigs.size();
                Contig contig = null;
                final SVKmerLong predecessor = getSolePredecessor(kmer, kSize2, counts);
                if ( predecessor == null || getSuccessorCount(predecessor, kSize2, kmerAdjacencyMap) > 1 ) {
                    contig = buildContig(kmer, kSize2, counts, kmerAdjacencyMap, contigId, locusMap);
                } else {
                    final SVKmerLong successor = getSoleSuccessor(kmer, kSize2, counts);
                    if ( successor == null || getPredecessorCount(successor, kSize2, kmerAdjacencyMap) > 1 ) {
                        contig = buildContig(kmer.reverseComplement(kSize2), kSize2, counts, kmerAdjacencyMap, contigId, locusMap);
                    }
                }
                if ( contig != null ) {
                    contigs.add(contig);
                    contigEnds.put(contig.getFirst().getKmer(), contig.getFirst());
                    contigEnds.put(contig.getLast().getKmer(), contig.getLast());
                }
            }
        });

        // dump contigs
        dumpGFA( contigs, contigEnds, graphOutput + ".gfa" );
        dumpFASTA( contigs, graphOutput + ".fasta" );
        dumpDOT( contigs, contigEnds, graphOutput + ".dot" );
/*
        reads.forEach(read -> {
            System.out.println(read.getHeader());
            SVKmerizer.stream(trimmedRead(read,minQ), kSize2, new SVKmerLong()).forEach(kmer -> {
                final ContigLocus locus = locusMap.get(((SVKmerLong)kmer).canonical(kSize2));
                if ( locus == null ) {
                    System.out.println("  no mapping");
                } else if ( ((SVKmerLong)kmer).isCanonical(kSize2) == locus.isCanonical() ) {
                    System.out.println("  " + locus.getContigId() + ":" + locus.getOffset());
                } else {
                    System.out.println("  " + locus.getContigId() + ":" + ~locus.getOffset());
                }
            });
        });
*/
        return null;
    }

    private byte[] trimmedRead( final FastqRead read, final int minQ ) {
        final byte[] quals = read.getQuals();
        final int trimLen =
                IntStream.range(0, quals.length).filter(idx -> quals[idx] < minQ).findFirst().orElse(quals.length);
        return Arrays.copyOf(read.getBases(), trimLen);
    }

    private void dumpDOT( final List<Contig> contigs, final Map<SVKmerLong, ContigEnd> contigEnds, final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            writer.write("digraph {\n");
            for ( final Contig contig : contigs ) {
                final int id = contig.getFirst().getContigId();
                final double width = contig.getSequence().length() / 100.;
                writer.write("tig" + id + " [width=" + width + "]\n");
                writer.write("tig" + id + "RC [width=" + width + "]\n");
            }
            for ( final Contig contig : contigs ) {
                final int id = contig.getFirst().getContigId();
                for ( final Connection connection : contig.getFirst().getConnections() ) {
                    final ContigEnd target = contigEnds.get(connection.getKmer());
                    if ( target.getContigId() >= id ) {
                        final String label = " [label=" + connection.getWeight() + "]\n";
                        writer.write("tig" + id + "RC -> tig" + target.getContigId() + (target.isFirst() ? "" : "RC") + label);
                        writer.write("tig" + target.getContigId() + (target.isFirst() ? "RC" : "") + " -> tig" + id + label);
                    }
                }
                for ( final Connection connection : contig.getLast().getConnections() ) {
                    final ContigEnd target = contigEnds.get(connection.getKmer());
                    if ( target.getContigId() >= id ) {
                        final String label = " [label=" + connection.getWeight() + "]\n";
                        writer.write("tig" + id + " -> tig" + target.getContigId() + (target.isFirst() ? "" : "RC") + label);
                        writer.write("tig" + target.getContigId() + (target.isFirst() ? "RC" : "") + " -> tig" + id + "RC" + label);
                    }
                }
            }
            writer.write("}\n");
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly DOT file.", ioe);
        }
    }

    private void dumpFASTA( final List<Contig> contigs, final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            for ( final Contig contig : contigs ) {
                writer.write(">" + contig.getFirst().getContigId() + "\n");
                writer.write(contig.getSequence() + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly FASTA file.", ioe);
        }
    }

    private void dumpGFA( final List<Contig> contigs, final Map<SVKmerLong, ContigEnd> contigEnds, final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            for ( final Contig contig : contigs ) {
                final int contigId = contig.getFirst().getContigId();
                final int seqLen = contig.getSequence().length();
                writer.write("S\ttig" + contigId + "\t" + contig.getSequence() + "\tLN:i:" + seqLen + "\n");
                for ( final Connection connection : contig.getFirst().getConnections() ) {
                    final ContigEnd target = contigEnds.get(connection.getKmer());
                    if ( target.getContigId() >= contigId ) {
                        final String dir = target.isFirst() ? "+" : "-";
                        writer.write("L\ttig" + contigId + "\t-\ttig" + target.getContigId() + "\t" + dir + "\n");
                    }
                }
                for ( final Connection connection : contig.getLast().getConnections() ) {
                    final ContigEnd target = contigEnds.get(connection.getKmer());
                    if ( target.getContigId() >= contigId ) {
                        final String dir = target.isFirst() ? "+" : "-";
                        writer.write("L\ttig" + contigId + "\t+\ttig" + target.getContigId() + "\t" + dir + "\n");
                    }
                }
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly FASTA file.", ioe);
        }
    }

    private static SVKmerLong getSolePredecessor( final SVKmerLong kmer, final int kSize, final int[] counts ) {
        Base solePred = null;
        final int offset = kmer.isCanonical(kSize) ? 0 : 4;
        for ( final Base base : Base.values() ) {
            if ( counts[base.ordinal() + offset] > 0 ) {
                if ( solePred != null ) return null; // found a 2nd predecessor
                solePred = base;
            }
        }
        if ( solePred == null ) return null;
        if ( !kmer.isCanonical(kSize) ) solePred = solePred.complement();
        return kmer.predecessor(solePred, kSize);
    }

    private static SVKmerLong getSoleSuccessor( final SVKmerLong kmer, final int kSize, final int[] counts ) {
        Base soleSucc = null;
        final int offset = kmer.isCanonical(kSize) ? 4 : 0;
        for ( final Base base : Base.values() ) {
            if ( counts[base.ordinal() + offset] > 0 ) {
                if ( soleSucc != null ) return null; // found a 2nd successor
                soleSucc = base;
            }
        }
        if ( soleSucc == null ) return null;
        if ( !kmer.isCanonical(kSize) ) soleSucc = soleSucc.complement();
        return kmer.successor(soleSucc, kSize);
    }

    private static int getPredecessorCount( final SVKmerLong kmer, final int kSize,
                                            final Map<SVKmerLong, int[]> kmerAdjacencyMap ) {
        if ( !kmer.isCanonical(kSize) ) {
            return getSuccessorCount(kmer.reverseComplement(kSize), kSize, kmerAdjacencyMap);
        }
        return getPredecessorCount(kmerAdjacencyMap.get(kmer));
    }

    private static int getPredecessorCount( final int[] counts ) {
        return (counts[0] > 0 ? 1 : 0) + (counts[1] > 0 ? 1 : 0) + (counts[2] > 0 ? 1 : 0) + (counts[3] > 0 ? 1 : 0);
    }

    private static int getSuccessorCount( final SVKmerLong kmer, final int kSize,
                                          final Map<SVKmerLong, int[]> kmerAdjacencyMap ) {
        if ( !kmer.isCanonical(kSize) ) {
            return getPredecessorCount(kmer.reverseComplement(kSize), kSize, kmerAdjacencyMap);
        }
        return getSuccessorCount(kmerAdjacencyMap.get(kmer));
    }

    private static int getSuccessorCount( final int[] counts ) {
        return (counts[4] > 0 ? 1 : 0) + (counts[5] > 0 ? 1 : 0) + (counts[6] > 0 ? 1 : 0) + (counts[7] > 0 ? 1 : 0);
    }

    private static Contig buildContig( final SVKmerLong kmerArg, final int kSize, final int[] countsArg,
                                       final Map<SVKmerLong, int[]> kmerAdjacencyMap, final int contigId,
                                       final Map<SVKmerLong, ContigLocus> locusMap ) {
        final ContigEnd firstEnd = new ContigEnd(contigId, true, kmerArg, kSize, countsArg);
        final StringBuilder sb = new StringBuilder(kmerArg.toString(kSize));
        SVKmerLong kmer = kmerArg;
        SVKmerLong lastKmer = kmerArg;
        int[] counts = countsArg;
        int[] lastCounts = countsArg;
        int offset = 0;
        locusMap.put(kmer.canonical(kSize), new ContigLocus(contigId, offset++, kmer.isCanonical(kSize)));
        while ( (kmer = getSoleSuccessor(kmer, kSize, counts)) != null ) {
            counts = kmerAdjacencyMap.get(kmer.canonical(kSize));
            if ( (kmer.isCanonical(kSize) ? getPredecessorCount(counts) : getSuccessorCount(counts)) != 1 ) break;
            sb.append(kmer.lastBase().toString());
            locusMap.put(kmer.canonical(kSize), new ContigLocus(contigId, offset++, kmer.isCanonical(kSize)));
            lastKmer = kmer;
            lastCounts = counts;
        }
        final ContigEnd lastEnd = new ContigEnd(contigId, false, lastKmer, kSize, lastCounts);
        return new Contig(sb.toString(), firstEnd, lastEnd);
    }

    private static final class ContigLocus {
        private final int contigId;
        private final int offset;
        private final boolean canonical;

        public ContigLocus( final int contigId, final int offset, final boolean canonical ) {
            this.contigId = contigId;
            this.offset = offset;
            this.canonical = canonical;
        }

        public int getContigId() { return contigId; }
        public int getOffset() { return offset; }
        public boolean isCanonical() { return canonical; }
    }

    private static final class Connection {
        private final SVKmerLong kmer;
        private final int weight;

        public Connection( final SVKmerLong kmer, final int weight ) {
            this.kmer = kmer;
            this.weight = weight;
        }

        public SVKmerLong getKmer() { return kmer; }
        public int getWeight() { return weight; }
    }

    private static final class ContigEnd {
        private final SVKmerLong kmer;
        private final int contigId;
        private final boolean isFirst;
        private final boolean isCanonical;
        private final List<Connection> connections;

        public ContigEnd( final int contigId, final boolean isFirst,
                          final SVKmerLong kmer, final int kSize, final int[] counts ) {
            this.kmer = kmer.canonical(kSize);
            this.contigId = contigId;
            this.isFirst = isFirst;
            this.isCanonical = kmer.isCanonical(kSize);
            this.connections =
                    new ArrayList<>(isFirst==isCanonical?getPredecessorCount(counts):getSuccessorCount(counts));
            addConnections( kSize, counts );
        }

        public SVKmerLong getKmer() { return kmer; }
        public int getContigId() { return contigId; }
        public boolean isFirst() { return isFirst; }
        public boolean isCanonical() { return isCanonical; }
        public List<Connection> getConnections() { return Collections.unmodifiableList(connections); }

        private void addConnections( final int kSize, final int[] counts ) {
            for ( final Base base : Base.values() ) {
                final int weight;
                SVKmerLong connectingKmer = null;
                if ( isFirst ) {
                    if ( isCanonical ) {
                        if ( (weight = counts[base.ordinal()]) > 0 ) {
                            connectingKmer = kmer.predecessor(base, kSize).canonical(kSize);
                        }
                    } else {
                        if ( (weight = counts[4 + base.complement().ordinal()]) > 0 ) {
                            connectingKmer = kmer.successor(base.complement(), kSize).canonical(kSize);
                        }
                    }
                } else {
                    if ( isCanonical ) {
                        if ( (weight = counts[4 + base.ordinal()]) > 0 ) {
                            connectingKmer = kmer.successor(base, kSize).canonical(kSize);
                        }
                    } else {
                        if ( (weight = counts[base.complement().ordinal()]) > 0 ) {
                            connectingKmer = kmer.predecessor(base.complement(), kSize).canonical(kSize);
                        }
                    }
                }
                if ( weight > 0 ) {
                    connections.add(new Connection(connectingKmer, weight));
                }
            }
        }
    }

    private static final class Contig {
        private final String sequence;
        private final ContigEnd first;
        private final ContigEnd last;

        public Contig( final String sequence, final ContigEnd first, final ContigEnd last ) {
            this.sequence = sequence;
            this.first = first;
            this.last = last;
            if ( !first.isFirst() || last.isFirst() || first.getContigId() != last.getContigId() ) {
                throw new GATKException("something got mixed up when building the contig");
            }
        }

        public String getSequence() { return sequence; }
        public ContigEnd getFirst() { return first; }
        public ContigEnd getLast() { return last; }
    }
}
