package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer.ASCIICharSequence;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils.FastqRead;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;

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
    private static final int GAP_CONTIG_ID = Integer.MIN_VALUE;

    @Argument(doc = "input fastq",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME)
    private String fastqFile;

    @Argument(doc = "graph output",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String graphOutput;

    @Argument(doc = "kmer size", fullName = "kSize", optional = true)
    private int kSize = 39;

    @Argument(doc = "minimum quality", fullName = "minQ", optional = true)
    private int minQ = 7;

    @Argument(doc = "minimum reliable kmer count", fullName = "minKCount", optional = true)
    private int minKCount = 3;

    @Override protected Object doWork() {
        final List<FastqRead> reads = SVFastqUtils.readFastqFile(fastqFile);
        if ( reads.size()%2 != 0 ) {
            throw new UserException("FASTQ input must be interleaved pairs, but there are an odd number of reads.");
        }

        final int kSize2 = kSize - 2;
        final Set<SVKmerLong> goodKmers = countKmers(reads, kSize, minQ, minKCount);
        final Map<SVKmerLong, Integer> contigEnds = new HashMap<>();
        final List<Contig> contigs = buildContigs(buildAdjacenciesSet(goodKmers, kSize), kSize2, contigEnds);
        final List<List<SVInterval>> readPaths = pathReadPairs(reads, contigs, kSize2);
        final int nCapturedGaps = countCapturedGaps(readPaths, kSize2);
        final Set<SVKmerLong> gapKmers = fillGaps(readPaths, contigs, reads, kSize);
        final Set<SVKmerLong> allKmers =
                new HashSet<>(SVUtils.hashMapCapacity(goodKmers.size() + gapKmers.size()));
        allKmers.addAll(goodKmers);
        allKmers.addAll(gapKmers);

        final Map<SVKmerLong, Integer> patchedEnds = new HashMap<>();
        final List<Contig> patchedContigs = buildContigs(buildAdjacenciesSet(allKmers, kSize), kSize2, patchedEnds);
        final List<List<SVInterval>> patchedReadPaths = pathReadPairs(reads, patchedContigs, kSize2);
        final int nPatchedGaps = countCapturedGaps(patchedReadPaths, kSize2);

        System.out.println("Closed " + (nCapturedGaps-nPatchedGaps) + " of " + nCapturedGaps + " captured gaps.");
        printPaths(patchedReadPaths, patchedContigs, kSize2);

        // dump contigs
        dumpGFA( patchedContigs, patchedEnds, kSize - 2, graphOutput + ".gfa" );
        dumpFASTA( patchedContigs, graphOutput + ".fasta" );
        dumpDOT( patchedContigs, patchedEnds, kSize - 2, graphOutput + ".dot" );

        return null;
    }

    private static int countCapturedGaps( final List<List<SVInterval>> readPaths, final int kSize2 ) {
        int count = 0;
        for ( final List<SVInterval> path : readPaths ) {
            final int nIntervals = path.size() - 1;
            for ( int idx = 1; idx < nIntervals; ++idx ) {
                final SVInterval interval = path.get(idx);
                if ( interval.getContig() == GAP_CONTIG_ID && interval.getLength() < kSize2 ) {
                    count += 1;
                }
            }
        }
        return count;
    }

    private static Set<SVKmerLong> fillGaps( final List<List<SVInterval>> readPaths,
                                             final List<Contig> contigs,
                                             final List<FastqRead> reads,
                                             final int kSize ) {
//        System.out.println("Captured gaps:");
        final int kSize2 = kSize - 2;
        final int nKmers =
                readPaths.stream().mapToInt(path ->
                        path.stream().mapToInt(interval ->
                                interval.getContig() == GAP_CONTIG_ID ? interval.getLength() : 0).sum()).sum();
        final Map<SVKmerLong, Integer> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        final int nReads = readPaths.size();
        for ( int readIdx = 0; readIdx != nReads; ++readIdx ) {
            final List<SVInterval> path = readPaths.get(readIdx);
            final int nIntervals = path.size() - 1; // skip any final gap -- we only want captured gaps
            int readOffset = -1; // we want to get the base before the gap into the extended kmer
            for ( int idx = 0; idx < nIntervals; ++idx ) {
                final SVInterval interval = path.get(idx);
                if ( idx > 0 && interval.getContig() == GAP_CONTIG_ID && interval.getLength() < kSize2 ) {
                    final FastqRead read = reads.get(readIdx);
                    final int endOffset = readOffset + kSize - 1 + interval.getLength();
                    final CharSequence gapSequence =
                            new ASCIICharSequence(read.getBases()).subSequence(readOffset, endOffset);
                    SVKmerizer.canonicalStream(gapSequence, kSize, new SVKmerLong(kSize))
                            .forEach(kmer -> kmerCounts.merge(kmer, 1, Integer::sum));
/*
                    final StringBuilder sb = new StringBuilder();
                    sb.append(new ASCIICharSequence(read.getBases()).subSequence(readOffset+kSize2, endOffset-1));
                    sb.append(" ");
                    appendInterval(sb, path.get(idx-1), contigs, kSize2);
                    appendInterval(sb, path.get(idx+1), contigs, kSize2);
                    sb.append(' ').append(readIdx);
                    System.out.println(sb);
*/
                }
                readOffset += interval.getLength();
            }
        }

        kmerCounts.entrySet().removeIf(entry -> entry.getValue() <= 1);

        return kmerCounts.keySet();
    }

    private static void printPaths( final List<List<SVInterval>> readPaths,
                                    final List<Contig> contigs,
                                    final int kSize2 ) {
        final int nReads = readPaths.size();
        for ( int readId = 0; readId != nReads; readId += 2 ) {
            final StringBuilder sb = new StringBuilder();
            sb.append("Pair ").append(readId / 2).append('\t');
            describePath(sb, readPaths.get(readId), contigs, kSize2);
            sb.append(" | ");
            describePath(sb, readPaths.get(readId+1), contigs, kSize2);
            System.out.println(sb);
        }
    }

    private static void describePath( final StringBuilder sb, final List<SVInterval> path,
                                      final List<Contig> contigs, final int kSize2 ) {
        for ( final SVInterval interval : path ) {
            if ( interval.getContig() == GAP_CONTIG_ID ) {
                sb.append(interval.getLength()).append("X ");
            } else {
                appendInterval(sb, interval, contigs, kSize2);
            }
        }
    }

    private static void appendInterval( final StringBuilder sb, final SVInterval interval,
                                        final List<Contig> contigs, final int kSize2 ) {
        final int contigId = interval.getContig();
        if ( contigId < 0 ) sb.append('~').append(~contigId);
        else sb.append(contigId);
        sb.append(":").append(interval.getStart()).append('-').append(interval.getEnd()+kSize2-2).append('/');
        sb.append(contigs.get(contigId < 0 ? ~contigId : contigId).getSequence().length()-1).append(' ');
    }

    private static List<List<SVInterval>> pathReadPairs( final List<FastqRead> reads,
                                                         final List<Contig> contigs,
                                                         final int kSize2 ) {
        final int nKmers = contigs.stream().mapToInt(tig -> tig.getSequence().length()-kSize2+1).sum();
        final Map<SVKmerLong, SVLocation> contigKmerMap = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        final int nContigs = contigs.size();
        for ( int contigId = 0; contigId != nContigs; ++contigId ) {
            int contigOffset = 0;
            final Iterator<SVKmerLong> contigItr =
                    new SVKmerizer<>(contigs.get(contigId).getSequence(), kSize2, new SVKmerLong());
            while ( contigItr.hasNext() ) {
                final SVKmerLong kmer = contigItr.next();
                if ( kmer.isCanonical(kSize2) ) {
                    contigKmerMap.put(kmer, new SVLocation(contigId, contigOffset));
                } else {
                    final int contigLen = contigs.get(contigId).getSequence().length();
                    final int rcPosition = contigLen - contigOffset - kSize2;
                    contigKmerMap.put(kmer.reverseComplement(kSize2), new SVLocation(~contigId, rcPosition));
                }
                contigOffset += 1;
            }
        }

        final int nReads = reads.size();
        final List<List<SVInterval>> readPaths = new ArrayList<>(nReads);
        for ( int readId = 0; readId+1 < nReads; readId += 2 ) {
            final FastqRead read1 = reads.get(readId);
            readPaths.add(processRead(read1.getBases(), read1.getQuals(), kSize2, contigKmerMap, contigs));
            final FastqRead read2 = reads.get(readId + 1);
            final byte[] bases2 = Arrays.copyOf(read2.getBases(), read2.getBases().length);
            SequenceUtil.reverseComplement(bases2);
            final byte[] quals2 = Arrays.copyOf(read2.getQuals(), read2.getQuals().length);
            SequenceUtil.reverseQualities(quals2);
            readPaths.add(processRead(bases2, quals2, kSize2, contigKmerMap, contigs));
        }
        return readPaths;
    }

    private static List<SVInterval> processRead( final byte[] sequence,
                                                 final byte[] quals,
                                                 final int kSize2,
                                                 final Map<SVKmerLong, SVLocation> contigKmerMap,
                                                 final List<Contig> contigs ) {
        final List<SVInterval> readPath = new ArrayList<>();
        SVInterval currentSpan = null;
        int missCount = 0;
        final Iterator<SVKmerLong> readItr = new SVKmerizer<>(sequence, kSize2, new SVKmerLong());
        while ( readItr.hasNext() ) {
            final SVKmerLong kmer = readItr.next();
            final SVKmerLong kmerCanonical = kmer.canonical(kSize2);
            SVLocation location = contigKmerMap.get(kmerCanonical);
            if ( location != null && !kmer.equals(kmerCanonical) ) {
                int contigId = location.getContig();
                if ( contigId < 0 ) contigId = ~contigId;
                final int contigLen = contigs.get(contigId).getSequence().length();
                location = new SVLocation(~location.getContig(), contigLen - location.getPosition() - kSize2);
            }
            if ( location == null ) {
                if ( currentSpan != null ) {
                    readPath.add(currentSpan);
                    currentSpan = null;
                }
                missCount += 1;
            } else if ( currentSpan == null ) {
                if ( missCount > 0 ) {
                    readPath.add(new SVInterval(GAP_CONTIG_ID, 0, missCount));
                    missCount = 0;
                }
                final int position = location.getPosition();
                currentSpan = new SVInterval(location.getContig(), position, position + 1);
            } else if ( location.getContig() == currentSpan.getContig() &&
                        location.getPosition() == currentSpan.getEnd() ) {
                currentSpan =
                        new SVInterval(currentSpan.getContig(), currentSpan.getStart(), currentSpan.getEnd() + 1);
            } else {
                readPath.add(currentSpan);
                final int position = location.getPosition();
                currentSpan = new SVInterval(location.getContig(), position, position + 1);
            }
        }
        if ( missCount > 0 ) {
            readPath.add(new SVInterval(-1, 0, missCount));
        } else if ( currentSpan != null ) {
            readPath.add(currentSpan);
        }
        return readPath;
    }

    private static List<Contig> buildContigs( final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet,
                                              final int kSize2,
                                              final Map<SVKmerLong, Integer> contigEnds ) {
        // build contigs
        final List<Contig> contigs = new ArrayList<>();
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

        return contigs;
    }

    private static HopscotchSet<SVKmerAdjacencies> buildAdjacenciesSet( final Set<SVKmerLong> goodKmers,
                                                                        final int kSize ) {
        // build 61-mer adjacency map from 63-mers (assuming default kmer size)
        final int kSize2 = kSize - 2;
        final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet =
                new HopscotchSet<>(3*SVUtils.hashMapCapacity(goodKmers.size()));
        goodKmers.forEach(kmer -> {
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
                            newAdj.firstBase(kSize2), null).canonical(kSize2);
            final SVKmerAdjacencies oldNextAdj = kmerAdjacenciesSet.find(newNextAdj);
            if ( oldNextAdj == null ) kmerAdjacenciesSet.add(newNextAdj);
            else oldNextAdj.mergeAdjacencies(newNextAdj);
        });

        return kmerAdjacenciesSet;
    }

    private static Set<SVKmerLong> countKmers( final List<FastqRead> reads,
                                               final int kSize,
                                               final int minQ,
                                               final int minKCount ) {
        // kmerize each read, counting the observations of each kmer.
        // ignore kmers that contain a call with a quality less than minQ
        final int nKmers = reads.stream().mapToInt(read -> Math.min(0, read.getBases().length - kSize + 1)).sum();
        final Map<SVKmerLong, Integer> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        for ( final FastqRead read : reads ) {
            SVKmerizer.canonicalStream(maskedSequence(read, minQ), kSize, new SVKmerLong(kSize))
                    .forEach(kmer -> kmerCounts.merge(kmer, 1, Integer::sum));
        }

        // ignore kmers that appear less than minKCount times
        kmerCounts.entrySet().removeIf(entry -> entry.getValue() < minKCount);

        return kmerCounts.keySet();
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
                        if ( targetId != contigId ) {
                            writer.write("tig" + targetId + (targetRC ? "" : "RC") + " -> tig" + contigId + "\n");
                        }
                    }
                }
                for ( final int strandId : contig.getLast().getSuccessorContigs(contigEnds, kSize) ) {
                    final int targetId = strandId < 0 ? ~strandId : strandId;
                    final boolean targetRC = strandId < 0;
                    if ( targetId >= contigId ) {
                        writer.write("tig" + contigId + " -> tig" + targetId + (targetRC ? "RC\n" : "\n"));
                        if ( targetId != contigId ) {
                            writer.write("tig" + targetId + (targetRC ? "" : "RC") + " -> tig" + contigId + "RC\n");
                        }
                    }
                }
            }
            writer.write("}\n");
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
                    if ( targetId > contigId ) { // intentionally different than analogous comparison a few lines below
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
        final SVKmerAdjacencies kmerAdjacencies = kmerAdjacenciesSet.find(canonicalKmer);
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
