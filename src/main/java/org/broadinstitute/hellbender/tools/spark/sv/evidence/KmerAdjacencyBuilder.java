package org.broadinstitute.hellbender.tools.spark.sv.evidence;

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
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

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
    private int minQArg = 7;

    @Argument(doc = "minimum kmer count", fullName = "minKCount", optional = true)
    private int minKCountArg = 3;

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
            SVKmerizer.canonicalStream(maskedSequence(read, minQ), kSize, new SVKmerLong(kSize))
                    .forEach(kmer -> kmerCounts.merge((SVKmerLong)kmer, 1, Integer::sum));
        }

        // dump a histogram of kmer counts
        //final Map<Integer, Integer> kmerCountHistogram = new TreeMap<>();
        //kmerCounts.values().forEach( count -> kmerCountHistogram.merge(count, 1, Integer::sum));
        //kmerCountHistogram.forEach( (count, countCount) -> System.out.println(count + "\t" + countCount));

        // ignore kmers that appear less than minKCount times
        kmerCounts.entrySet().removeIf(entry -> entry.getValue() < minKCount);

        // build 61-mer adjacency map from 63-mers (assuming default kmer size)
        final int kSize2 = kSize - 2;
        final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet =
                new HopscotchSet<>(3*SVUtils.hashMapCapacity(kmerCounts.size()));
        kmerCounts.forEach((kmer, value) -> {
            final SVKmerAdjacencies newAdj =
                    new SVKmerAdjacencies(kmer.removeFirstAndLastBase(kSize), kmer.firstBase(kSize), kmer.lastBase());
            final SVKmerAdjacencies oldAdj = kmerAdjacenciesSet.find(newAdj);
            if ( oldAdj == null ) kmerAdjacenciesSet.add(newAdj);
            else oldAdj.mergeAdjacencies(newAdj);

            final SVKmerAdjacencies newPrevAdj =
                    (SVKmerAdjacencies)new SVKmerAdjacencies(newAdj.predecessor(kmer.firstBase(kSize), kSize2),
                                                             null,
                                                             newAdj.lastBase()).canonical(kSize2);
            final SVKmerAdjacencies oldPrevAdj = kmerAdjacenciesSet.find(newPrevAdj);
            if ( oldPrevAdj == null ) kmerAdjacenciesSet.add(newPrevAdj);
            else oldPrevAdj.mergeAdjacencies(newPrevAdj);

            final SVKmerAdjacencies newNextAdj =
                    (SVKmerAdjacencies)new SVKmerAdjacencies(newAdj.successor(kmer.lastBase(), kSize2),
                                                             kmer.lastBase(), null).canonical(kSize2);
            final SVKmerAdjacencies oldNextAdj = kmerAdjacenciesSet.find(newNextAdj);
            if ( oldNextAdj == null ) kmerAdjacenciesSet.add(newNextAdj);
            else oldNextAdj.mergeAdjacencies(newNextAdj);
        });

        for ( final SVKmerAdjacencies adj : kmerAdjacenciesSet ) {
            Utils.validateArg(adj.isCanonical(kSize2), "non-canonical kmer");
            for ( final Base base : Base.values() ) {
                if ( adj.hasPredecessor(base) ) {
                    final SVKmerLong kmer = adj.predecessor(base, kSize2).canonical(kSize2);
                    if ( kmerAdjacenciesSet.find(kmer) == null ) {
                        if ( kmerAdjacenciesSet.find(adj.predecessor(base.complement(), kSize2).canonical(kSize2)) != null ) {
                            System.out.println("can't find predecessor, but found complement " + kmer.toString(kSize2));
                        } else {
                            System.out.println("can't find predecessor " + kmer.toString(kSize2));
                        }
                    }
                }
                if ( adj.hasSuccessor(base) ) {
                    final SVKmerLong kmer = adj.successor(base, kSize2).canonical(kSize2);
                    if ( kmerAdjacenciesSet.find(kmer) == null ) {
                        if ( kmerAdjacenciesSet.find(adj.successor(base.complement(), kSize2).canonical(kSize2)) != null ) {
                            System.out.println("can't find successor, but found complement " + kmer.toString(kSize2));
                        } else {
                            System.out.println("can't find successor " + kmer.toString(kSize2));
                        }
                    }
                }
            }
        }

        // build contigs
        final List<Contig> contigs = new ArrayList<>();
        final Map<SVKmerLong, Integer> contigEnds = new HashMap<>();
        kmerAdjacenciesSet.forEach( adj -> {
            if ( !contigEnds.containsKey(adj) && !contigEnds.containsKey(adj.reverseComplement(kSize2)) ) {
                Contig contig = null;
                if ( isContigStart(adj, kmerAdjacenciesSet, kSize2) ) {
                    contig = buildContig(adj, kSize2, kmerAdjacenciesSet);
                } else if ( isContigEnd(adj, kmerAdjacenciesSet, kSize2) ) {
                    contig = buildContig(adj.reverseComplement(kSize2), kSize2, kmerAdjacenciesSet);
                }
                if ( contig != null ) {
                    final int contigId = contigs.size();
                    contigs.add(contig);
                    contigEnds.put(contig.getFirst(), contigId);
                    contigEnds.put(contig.getLast().reverseComplement(kSize2), ~contigId);
                }
            }
        });

        // dump contigs
        dumpGFA( contigs, contigEnds, kSize2, graphOutput + ".gfa" );
        dumpFASTA( contigs, graphOutput + ".fasta" );
        dumpDOT( contigs, contigEnds, kSize2, graphOutput + ".dot" );
/*
        reads.forEach(read -> {
            System.out.println(read.getHeader());
            SVKmerizer.stream(maskedSequence(read,minQ), kSize2, new SVKmerLong()).forEach(kmer -> {
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

    private static byte[] maskedSequence(final FastqRead read, final int minQ) {
        final byte[] quals = read.getQuals();
        final byte[] calls = Arrays.copyOf(read.getBases(), quals.length);
        for ( int idx = 0; idx != quals.length; ++idx ) {
            if ( quals[idx] < minQ ) calls[idx] = 'N';
        }
        return calls;
    }

    private static void dumpDOT( final List<Contig> contigs,
                                 final Map<SVKmerLong, Integer> contigEnds,
                                 final int kSize,
                                 final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            writer.write("digraph {\n");
            final int nContigs = contigs.size();
            for ( int contigId = 0; contigId != nContigs; ++contigId ) {
                final Contig contig = contigs.get(contigId);
                final double width = contig.getSequence().length() / 100.;
                writer.write("tig" + contigId + " [width=" + width + "]\n");
                writer.write("tig" + contigId + "RC [width=" + width + "]\n");
            }
            for ( int contigId = 0; contigId != nContigs; ++contigId ) {
                final Contig contig = contigs.get(contigId);
                for ( final int strandId : contig.getFirst().getPredecessorContigs(contigEnds, kSize) ) {
                    final int targetId = strandId < 0 ? ~strandId : strandId;
                    final boolean targetRC = strandId < 0;
                    if ( targetId >= contigId ) {
                        writer.write("tig" + contigId + "RC -> tig" + targetId + (targetRC ? "RC\n" : "\n"));
                        writer.write("tig" + targetId + (targetRC ? "" : "RC") + " -> tig" + contigId + "\n");
                    }
                }
                for ( final int strandId : contig.getLast().getSuccessorContigs(contigEnds, kSize) ) {
                    final int targetId = strandId < 0 ? ~strandId : strandId;
                    final boolean targetRC = strandId < 0;
                    if ( targetId >= contigId ) {
                        writer.write("tig" + strandId + " -> tig" + targetId + (targetRC ? "RC\n" : "\n"));
                        writer.write("tig" + targetId + (targetRC ? "" : "RC") + " -> tig" + contigId + "RC\n");
                    }
                }
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly DOT file.", ioe);
        }
    }

    private static void dumpFASTA( final List<Contig> contigs, final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            final int nContigs = contigs.size();
            for ( int contigId = 0; contigId != nContigs; ++contigId ) {
                writer.write(">" + contigId + "\n");
                writer.write(contigs.get(contigId).getSequence() + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly FASTA file.", ioe);
        }
    }

    private static void dumpGFA( final List<Contig> contigs,
                                 final Map<SVKmerLong, Integer> contigEnds,
                                 final int kSize,
                                 final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            final int nContigs = contigs.size();
            for ( int contigId = 0; contigId != nContigs; ++contigId ) {
                final Contig contig = contigs.get(contigId);
                final int seqLen = contig.getSequence().length();
                writer.write("S\ttig" + contigId + "\t" + contig.getSequence() + "\tLN:i:" + seqLen + "\n");
                for ( final int strandId : contig.getFirst().getPredecessorContigs(contigEnds, kSize) ) {
                    final int targetId = strandId < 0 ? ~strandId : strandId;
                    if ( targetId >= contigId ) {
                        final String dir = strandId < 0 ? "-" : "+";
                        writer.write("L\ttig" + contigId + "\t-\ttig" + targetId + "\t" + dir + "\n");
                    }
                }
                for ( final int strandId : contig.getLast().getSuccessorContigs(contigEnds, kSize) ) {
                    final int targetId = strandId < 0 ? ~strandId : strandId;
                    if ( targetId >= contigId ) {
                        final String dir = strandId < 0 ? "-" : "+";
                        writer.write("L\ttig" + contigId + "\t+\ttig" + targetId + "\t" + dir + "\n");
                    }
                }
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly FASTA file.", ioe);
        }
    }

    private static SVKmerAdjacencies strandSensitiveLookup( final SVKmerLong kmer,
                                                            final int kSize,
                                                            final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet ) {
        final SVKmerLong canonicalKmer = kmer.canonical(kSize);
        final SVKmerAdjacencies kmerAdjacencies =
                kmerAdjacenciesSet.find(canonicalKmer);
        if ( kmerAdjacencies == null ) {
            throw new GATKException("can't find expected kmer in adjacencies set");
        }
        return kmer.equals(canonicalKmer) ? kmerAdjacencies : kmerAdjacencies.reverseComplement(kSize);
    }

    private static boolean isContigStart( final SVKmerAdjacencies kmerAdjacencies,
                                          final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet,
                                          final int kSize ) {
        final SVKmerLong predecessorKmer = kmerAdjacencies.getSolePredecessor(kSize);
        if ( predecessorKmer == null ) return true;
        final SVKmerAdjacencies predecessorAdjacencies =
                strandSensitiveLookup(predecessorKmer, kSize, kmerAdjacenciesSet);
        return predecessorAdjacencies.successorCount() > 1;
    }

    private static boolean isContigEnd( final SVKmerAdjacencies kmerAdjacencies,
                                        final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet,
                                        final int kSize ) {
        final SVKmerLong successorKmer = kmerAdjacencies.getSoleSuccessor(kSize);
        if ( successorKmer == null ) return true;
        final SVKmerAdjacencies successorAdjacencies =
                strandSensitiveLookup(successorKmer, kSize, kmerAdjacenciesSet);
        return successorAdjacencies.predecessorCount() > 1;
    }

    private static Contig buildContig( final SVKmerAdjacencies kmerAdjacencies, final int kSize,
                                       final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet ) {
        final StringBuilder contigSequence = new StringBuilder(kmerAdjacencies.toString(kSize));
        SVKmerAdjacencies currentAdjacencies = kmerAdjacencies;
        SVKmerLong successorKmer;
        while ( (successorKmer = currentAdjacencies.getSoleSuccessor(kSize)) != null ) {
            final SVKmerAdjacencies successorAdjacencies =
                    strandSensitiveLookup(successorKmer, kSize, kmerAdjacenciesSet);
            if ( successorAdjacencies.predecessorCount() > 1 ) break;
            contigSequence.append(successorAdjacencies.lastBase().name());
            currentAdjacencies = successorAdjacencies;
        }
        return new Contig(contigSequence.toString(), kmerAdjacencies, currentAdjacencies);
    }

    private static final class Contig {
        private final String sequence;
        private final SVKmerAdjacencies first;
        private final SVKmerAdjacencies last;

        public Contig( final String sequence, final SVKmerAdjacencies first, final SVKmerAdjacencies last ) {
            this.sequence = sequence;
            this.first = first;
            this.last = last;
        }

        public String getSequence() { return sequence; }
        public SVKmerAdjacencies getFirst() { return first; }
        public SVKmerAdjacencies getLast() { return last; }
    }
}
