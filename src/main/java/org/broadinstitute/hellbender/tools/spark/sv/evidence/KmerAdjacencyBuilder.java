package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer.ASCIICharSequence;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils.FastqRead;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.*;

@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) junk", summary = "complete crap",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class KmerAdjacencyBuilder extends CommandLineProgram {
    private static final long serialVersionUID = 1L;
    private static final int GAP_CONTIG_ID = Integer.MIN_VALUE;
    private static final int MAX_PATCHED_QUALITY_SUM = 30;

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

    @Argument(doc = "kmer size", fullName = "bigKSize", optional = true)
    private int bigKSize = 111;

    @Argument(doc = "minimum quality", fullName = "minQ", optional = true)
    private int minQ = 7;

    @Argument(doc = "minimum reliable kmer count", fullName = "minKCount", optional = true)
    private int minKCount = 3;

    @Argument(doc = "minimum reliable kmer count", fullName = "minBigKCount", optional = true)
    private int minBigKCount = 4;

    @Override protected Object doWork() {
        final List<FastqRead> reads = squashNs(SVFastqUtils.readFastqFile(fastqFile));

        final int kSize2 = kSize - 2;
        final Set<SVKmerLong> goodKmers = countKmers(reads, kSize, minQ, minKCount, new SVKmerLong(kSize));
        final Map<SVKmer, Integer> contigEnds = new HashMap<>();
        final List<Contig> contigs = buildContigs(buildAdjacenciesSet(goodKmers, kSize), kSize2, contigEnds);
        final List<List<SVInterval>> readPaths = pathReads(reads, contigs, kSize2);
        final int nCapturedGaps = countCapturedGaps(readPaths, kSize2);
        goodKmers.addAll(fillGaps(readPaths, contigs, reads, kSize));

        final Map<SVKmer, Integer> patchedEnds = new HashMap<>();
        final List<Contig> patchedContigs = buildContigs(buildAdjacenciesSet(goodKmers, kSize), kSize2, patchedEnds);
        final List<List<SVInterval>> patchedReadPaths = pathReads(reads, patchedContigs, kSize2);
        final int nPatchedGaps = countCapturedGaps(patchedReadPaths, kSize2);

        System.out.println("Closed " + (nCapturedGaps-nPatchedGaps) + " of " + nCapturedGaps + " captured gaps.");

        final Set<SVKmerHuge> bigKmers = countBigKmers(reads, patchedReadPaths, bigKSize, minBigKCount, new SVKmerHuge(bigKSize));
        final Map<SVKmer, Integer> bigEnds = new HashMap<>();
        final List<Contig> bigContigs =
                buildContigs(buildAdjacenciesSet(bigKmers, bigKSize), bigKSize - 2, bigEnds);

        FermiLiteAssembly initialAssembly = convertAssembly(bigContigs, bigEnds, bigKSize - 2);
        final FermiLiteAssembly assembly =
                FermiLiteAssemblyHandler.reviseAssembly(initialAssembly, true, true);

        try ( final Writer writer = new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(graphOutput + ".gfa"))) ) {
            assembly.writeGFA(writer);
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write " + graphOutput + ".gfa", ioe);
        }

        // dump contigs
//        dumpGFA( bigContigs, bigEnds, bigKSize - 2, graphOutput + ".gfa" );
        dumpFASTA( bigContigs, graphOutput + ".fasta" );
        dumpDOT( bigContigs, bigEnds, bigKSize - 2, graphOutput + ".dot" );

        return null;
    }

    private static Set<SVKmerHuge> countBigKmers( final List<FastqRead> reads,
                                                  final List<List<SVInterval>> paths,
                                                  final int kSize,
                                                  final int minKCount,
                                                  final SVKmerHuge exemplar ) {
        final int nKmers = reads.stream().mapToInt(read -> Math.min(0, read.getBases().length - kSize + 1)).sum();
        final Map<SVKmerHuge, Integer> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));

        final int nReads = reads.size();
        for ( int readId = 0; readId != nReads; ++readId ) {
            final List<SVInterval> path = paths.get(readId);
            final int nIntervals = path.size();
            final SVInterval interval0 = path.get(0);
            if ( nIntervals == 1 && interval0.getContig() == GAP_CONTIG_ID ) continue;

            final SVInterval intervalN = path.get(nIntervals - 1);
            final byte[] calls = reads.get(readId).getBases();
            final int begin = interval0.getContig() == GAP_CONTIG_ID ? interval0.getLength() : 0;
            final int end = calls.length -
                    ((nIntervals > 1 && intervalN.getContig() == GAP_CONTIG_ID) ? intervalN.getLength() : 0);
            SVKmerizer.canonicalStream(new SVKmerizer.ASCIICharSubSequence(calls, begin, end), kSize, exemplar)
                    .forEach(kmer -> kmerCounts.merge(kmer, 1, Integer::sum));
        }

        final IntHistogram countHistogram = new IntHistogram(100);
        for ( final Integer count : kmerCounts.values() ) {
            countHistogram.addObservation(count);
        }
        System.out.println("Kmer count histogram:");
        for ( int idx = 1; idx != 101; ++idx ) {
            System.out.println(idx + "\t" + countHistogram.getNObservations(idx));
        }

        // ignore kmers that appear less than minKCount times
        kmerCounts.entrySet().removeIf(entry -> entry.getValue() < minKCount);

        Set<SVKmerHuge> result = new HashSet<>(SVUtils.hashMapCapacity(kmerCounts.size()));
        result.addAll(kmerCounts.keySet());
        return result;
    }

    public static FermiLiteAssembly convertAssembly( final List<Contig> contigs,
                                                     final Map<SVKmer, Integer> contigEnds,
                                                     final int kSize ) {
        final List<FermiLiteAssembly.Contig> fContigs = new ArrayList<>(contigs.size());
        for ( final Contig contig : contigs ) {
            final byte[] perBaseCoverage = new byte[contig.getSequence().length()];
            fContigs.add(new FermiLiteAssembly.Contig(contig.getSequence().getBytes(), perBaseCoverage, 0));
        }
        for ( int contigId = 0; contigId != contigs.size(); ++contigId ) {
            final Contig contig = contigs.get(contigId);
            List<Integer> predIds = contig.getFirst().getPredecessorContigs(contigEnds, kSize);
            List<Integer> succIds = contig.getLast().getSuccessorContigs(contigEnds, kSize);
            final List<FermiLiteAssembly.Connection> connections =
                    new ArrayList<>(predIds.size()+succIds.size());
            for ( final Integer id : predIds ) {
                final int targetId = id < 0 ? ~id : id;
                connections.add(new FermiLiteAssembly.Connection(fContigs.get(targetId),kSize-1, true, id >= 0));
            }
            for ( final Integer id : succIds ) {
                final int targetId = id < 0 ? ~id : id;
                connections.add(new FermiLiteAssembly.Connection(fContigs.get(targetId),kSize-1, false, id < 0));
            }
            fContigs.get(contigId).setConnections(connections);
        }
        return new FermiLiteAssembly(fContigs);
    }

    private static List<FastqRead> squashNs( final List<FastqRead> reads ) {
        for ( final FastqRead read : reads ) {
            final byte[] calls = read.getBases();
            final byte[] quals = read.getQuals();
            final int nCalls = calls.length;
            if ( calls[0] == 'N' || calls[0] == 'n' ) {
                calls[0] = calls[1];
                quals[0] = 1;
            }
            for ( int idx = 1; idx != nCalls; ++idx ) {
                if ( calls[idx] == 'N' || calls[idx] == 'n' ) {
                    calls[idx] = calls[idx - 1];
                    quals[idx] = 1;
                }
            }
        }
        return reads;
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
        final int kSize2 = kSize - 2;
        final int nKmers =
                readPaths.stream().mapToInt(path ->
                        path.stream().mapToInt(interval ->
                                interval.getContig() == GAP_CONTIG_ID ? interval.getLength()+2*kSize-2 : 0).sum()).sum();
        final Map<SVKmerLong, Integer> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        final int nReads = readPaths.size();
        for ( int readIdx = 0; readIdx != nReads; ++readIdx ) {
            final List<SVInterval> path = readPaths.get(readIdx);
            final int nIntervals = path.size() - 1; // skip any final gap -- we only want captured gaps
            int readOffset = 0;
            for ( int idx = 0; idx < nIntervals; ++idx ) {
                final SVInterval interval = path.get(idx);
                if ( idx > 0 && interval.getContig() == GAP_CONTIG_ID ) {
                    if ( interval.getLength() == kSize2 ) {
                        final SVInterval prevInterval = path.get(idx - 1);
                        final CharSequence prevSequence = getContigSequence(prevInterval.getContig(), contigs);
                        final int prevPosition = prevInterval.getEnd() + kSize2 - 1;
                        if ( prevPosition < prevSequence.length() ) {
                            final SVInterval nextInterval = path.get(idx + 1);
                            final int nextPosition = nextInterval.getStart() - 1;
                            if ( nextPosition >= 0 ) {
                                final CharSequence nextSequence = getContigSequence(nextInterval.getContig(), contigs);
                                if ( nextSequence.charAt(nextPosition) == prevSequence.charAt(prevPosition) ) {
                                    final FastqRead read = reads.get(readIdx);
                                    final byte[] calls = read.getBases();
                                    final byte[] quals = read.getQuals();
                                    final int readPosition = readOffset + kSize2 - 1;
                                    calls[readPosition] = (byte)nextSequence.charAt(nextPosition);
                                    quals[readPosition] = 1;
                                }
                            }
                        }
                    } else if ( interval.getLength() < kSize2 ) {
                        final FastqRead read = reads.get(readIdx);
                        final byte[] calls = read.getBases();
                        final int gapStart = Math.max(0, readOffset - 2);
                        final int gapEnd = Math.min(calls.length, readOffset + interval.getLength() + kSize - 1);
                        final CharSequence gapSequence =
                                new ASCIICharSequence(calls).subSequence(gapStart, gapEnd);
                        SVKmerizer.canonicalStream(gapSequence, kSize, new SVKmerLong(kSize))
                                .forEach(kmer -> kmerCounts.merge(kmer, 1, Integer::sum));
                    }
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
        int contigId = interval.getContig();
        if ( contigId < 0 ) {
            sb.append('~');
            contigId = ~contigId;
        }
        sb.append(contigId);
        final int end = interval.getEnd() - 1;
        final int last = contigs.get(contigId).getSequence().length() - kSize2;
        if ( interval.getStart() == 0 && end == last ) {
            sb.append(' ');
        } else {
            sb.append(":").append(interval.getStart()).append('-').append(end).append('/');
            sb.append(last).append(' ');
        }
    }

    private static List<List<SVInterval>> pathReads( final List<FastqRead> reads,
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

        final int[] contigLengths = new int[nContigs];
        for ( int idx = 0; idx != nContigs; ++idx ) {
            contigLengths[idx] = contigs.get(idx).getSequence().length();
        }

        final int nReads = reads.size();
        final List<List<SVInterval>> readPaths = new ArrayList<>(nReads);
        for ( final FastqRead read : reads ) {
            final byte[] sequence = read.getBases();
            final List<SVInterval> readPath = pathRead(sequence, kSize2, contigKmerMap, contigLengths);
            while ( readPath.size() > 1 && readPath.get(readPath.size() - 1).getContig() == GAP_CONTIG_ID ) {
                if ( !patchEndGap(sequence, readPath, kSize2, contigKmerMap, contigLengths) ) {
                    break;
                }
            }
            readPaths.add(readPath);
        }
        return readPaths;
    }

    private static boolean patchEndGap( final byte[] sequenceArg,
                                        final List<SVInterval> path,
                                        final int kSize2,
                                        final Map<SVKmerLong, SVLocation> contigKmerMap,
                                        final int[] contigLengths ) {
        final int subSequenceStart = sequenceArg.length - path.get(path.size() - 1).getLength() - kSize2 + 1;
        final byte[] sequence = Arrays.copyOfRange(sequenceArg, subSequenceStart, sequenceArg.length);
        final int patchOffset = kSize2 - 1;
        final byte call = (byte)Character.toUpperCase((char)(sequence[patchOffset] & 0xff));
        List<SVInterval> bestPathPatch = null;
        int bestMatchLen = 0;
        byte bestCall = call;
        for ( final byte newCall : new byte[] {'A', 'C', 'G', 'T'} ) {
            if ( call != newCall ) {
                sequence[kSize2 - 1] = newCall;
                List<SVInterval> pathPatch = pathRead(sequence, kSize2, contigKmerMap, contigLengths);
                final int matchLen = pathPatch.stream()
                        .mapToInt(interval -> interval.getContig() == GAP_CONTIG_ID ? 0 : interval.getLength())
                        .sum();
                if ( matchLen > bestMatchLen ) {
                    bestMatchLen = matchLen;
                    bestPathPatch = pathPatch;
                    bestCall = newCall;
                }
            }
        }
        sequence[patchOffset] = bestCall;
        if ( bestPathPatch != null ) {
            path.remove(path.size() - 1);
            final SVInterval last = path.get(path.size() - 1);
            final SVInterval first = bestPathPatch.get(0);
            if ( last.gapLen(first) == 0 ) {
                path.remove(path.size() - 1);
                bestPathPatch.set(0, last.join(first));
            }
            path.addAll(bestPathPatch);
            System.arraycopy(sequence, 0, sequenceArg, subSequenceStart, sequence.length);
        }
        return bestPathPatch != null;
    }
/*
    private static List<SVInterval> ecEnds( final byte[] sequenceArg,
                                            final byte[] quals,
                                            final int kSize2,
                                            final Map<SVKmerLong, SVLocation> contigKmerMap,
                                            final List<Contig> contigs,
                                            final int[] contigLengths ) {
        final byte[] sequence = Arrays.copyOf(sequenceArg, sequenceArg.length);
        final List<SVInterval> path = pathRead(sequence, kSize2, contigKmerMap, contigLengths);
        List<SVInterval> patchedPath = new ArrayList<>(path.size());
        patchedPath.addAll(path);
        int qualSum = 0;
        while ( patchedPath.size() > 1 && patchedPath.get(0).getContig() == GAP_CONTIG_ID && patchedPath.get(1).getStart() > 0 ) {
            final int gapLen = patchedPath.get(0).getLength();
            final int patchOffset = gapLen - 1;
            final int patchQual = quals[patchOffset] * 10 / 10;
            if ( qualSum + patchQual > MAX_PATCHED_QUALITY_SUM ) {
                break;
            }
            final SVInterval path1 = patchedPath.get(1);
            final int contigOffset = path1.getStart() - 1;
            sequence[patchOffset] =
                    (byte)getContigSequence(path1.getContig(), contigs).charAt(contigOffset);
            final CharSequence headSequence =
                    new SVKmerizer.ASCIICharSubSequence(sequence, 0, gapLen + kSize2 - 1);
            final List<SVInterval> gapPath = pathRead(headSequence, kSize2, contigKmerMap, contigLengths);
            patchedPath.remove(0);
            patchedPath.set(0, new SVInterval(path1.getContig(), path1.getStart() - 1, path1.getEnd()));
            gapPath.remove(gapPath.size() - 1);
            patchedPath.addAll(0, gapPath);
            qualSum += patchQual;
        }
        if ( (path.get(0).getContig() == GAP_CONTIG_ID ? path.get(0).getLength() : 0) -
                (patchedPath.get(0).getContig() == GAP_CONTIG_ID ? patchedPath.get(0).getLength() : 0) > 1 ) {
            System.arraycopy(sequence, 0, sequenceArg, 0, sequence.length);
            return patchedPath;
        }
        return path;
    }
*/
    public static final class StringRC implements CharSequence {
        final String sequence;
        final int offset;

        public StringRC( final String sequence ) {
            this.sequence = sequence;
            this.offset = sequence.length() - 1;
        }

        @Override public int length() { return sequence.length(); }

        @Override public char charAt( final int index ) {
            switch ( sequence.charAt(offset - index) ) {
                case 'a': case 'A': return 'T';
                case 'c': case 'C': return 'G';
                case 'g': case 'G': return 'C';
                case 't': case 'T': return 'A';
                default: throw new IllegalStateException("sequence contains bogus base call");
            }
        }

        @Override public CharSequence subSequence( final int start, final int end ) {
            return new StringRC(sequence.substring(offset - end, offset - start));
        }

        @Override public String toString() { return new StringBuilder(this).toString(); }
    }

    private static CharSequence getContigSequence( final int contigId, final List<Contig> contigs ) {
        if ( contigId < 0 ) {
            return new StringRC(contigs.get(~contigId).sequence);
        }
        return contigs.get(contigId).sequence;
    }

    private static List<SVInterval> pathRead( final byte[] sequence,
                                              final int kSize2,
                                              final Map<SVKmerLong, SVLocation> contigKmerMap,
                                              final int[] contigLengths ) {
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
                final int contigLen = contigLengths[contigId < 0 ? ~contigId : contigId];
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
                final int currentContigId = currentSpan.getContig();
                final int contigLen = contigLengths[currentContigId < 0 ? ~currentContigId : currentContigId];
                if ( currentSpan.getEnd() + kSize2 - 1 != contigLen || position > 0 ) {
                    readPath.add(new SVInterval(GAP_CONTIG_ID, 0, 0));
                }
                currentSpan = new SVInterval(location.getContig(), position, position + 1);
            }
        }
        if ( missCount > 0 ) {
            readPath.add(new SVInterval(GAP_CONTIG_ID, 0, missCount));
        } else if ( currentSpan != null ) {
            readPath.add(currentSpan);
        }

        return readPath;
    }

    private static List<Contig> buildContigs( final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet,
                                              final int kSize2,
                                              final Map<SVKmer, Integer> contigEnds ) {
        // build contigs
        final List<Contig> contigs = new ArrayList<>();
        kmerAdjacenciesSet.forEach( adj -> {
            if ( !contigEnds.containsKey(adj.getKmer()) && !contigEnds.containsKey(adj.getKmer().reverseComplement(kSize2)) ) {
                Contig contig = null;
                if ( isContigStart(adj, kmerAdjacenciesSet, kSize2) ) {
                    contig = buildContig(adj, kSize2, kmerAdjacenciesSet);
                } else if ( isContigEnd(adj, kmerAdjacenciesSet, kSize2) ) {
                    contig = buildContig(adj.reverseComplement(kSize2), kSize2, kmerAdjacenciesSet);
                }
                if ( contig != null ) {
                    final int contigId = contigs.size();
                    contigs.add(contig);
                    contigEnds.put(contig.getFirst().getKmer(), contigId);
                    contigEnds.put(contig.getLast().getKmer().reverseComplement(kSize2), ~contigId);
                }
            }
        });

        System.out.println("nKmers in adjacency set = " + kmerAdjacenciesSet.size());
        System.out.println("nKmers in contigs = " + contigs.stream().mapToInt(tig -> tig.getSequence().length() - kSize2 + 1).sum());

        return contigs;
    }

    private static <T extends SVKmer> HopscotchSet<SVKmerAdjacencies> buildAdjacenciesSet( final Set<T> goodKmers,
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
                    new SVKmerAdjacencies(newAdj.getKmer().predecessor(kmer.firstBase(kSize), kSize2),
                                          null,
                                          newAdj.getKmer().lastBase()).canonical(kSize2);
            final SVKmerAdjacencies oldPrevAdj = kmerAdjacenciesSet.find(newPrevAdj);
            if ( oldPrevAdj == null ) kmerAdjacenciesSet.add(newPrevAdj);
            else oldPrevAdj.mergeAdjacencies(newPrevAdj);

            final SVKmerAdjacencies newNextAdj =
                    new SVKmerAdjacencies(newAdj.getKmer().successor(kmer.lastBase(), kSize2),
                                          newAdj.getKmer().firstBase(kSize2),
                                          null).canonical(kSize2);
            final SVKmerAdjacencies oldNextAdj = kmerAdjacenciesSet.find(newNextAdj);
            if ( oldNextAdj == null ) kmerAdjacenciesSet.add(newNextAdj);
            else oldNextAdj.mergeAdjacencies(newNextAdj);
        });

        return kmerAdjacenciesSet;
    }

    private static Set<SVKmerLong> countKmers( final List<FastqRead> reads,
                                                         final int kSize,
                                                         final int minQ,
                                                         final int minKCount,
                                                         final SVKmerLong exemplar ) {
        // kmerize each read, counting the observations of each kmer.
        // ignore kmers that contain a call with a quality less than minQ
        final int nKmers = reads.stream().mapToInt(read -> Math.min(0, read.getBases().length - kSize + 1)).sum();
        final Map<SVKmerLong, Integer> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        for ( final FastqRead read : reads ) {
            SVKmerizer.canonicalStream(maskedSequence(read, minQ), kSize, exemplar)
                    .forEach(kmer -> kmerCounts.merge(kmer, 1, Integer::sum));
        }
/*
        final IntHistogram countHistogram = new IntHistogram(100);
        for ( final Integer count : kmerCounts.values() ) {
            countHistogram.addObservation(count);
        }
        System.out.println("Kmer count histogram:");
        for ( int idx = 1; idx != 101; ++idx ) {
            System.out.println(idx + "\t" + countHistogram.getNObservations(idx));
        }
*/
        // ignore kmers that appear less than minKCount times
        kmerCounts.entrySet().removeIf(entry -> entry.getValue() < minKCount);

        Set<SVKmerLong> result = new HashSet<>(SVUtils.hashMapCapacity(kmerCounts.size()));
        result.addAll(kmerCounts.keySet());
        return result;
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
                                 final Map<SVKmer, Integer> contigEnds,
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
                                 final Map<SVKmer, Integer> contigEnds,
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

    private static SVKmerAdjacencies strandSensitiveLookup( final SVKmer kmer,
                                                            final int kSize,
                                                            final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet ) {
        final SVKmer canonicalKmer = kmer.canonical(kSize);
        final SVKmerAdjacencies kmerAdjacencies = kmerAdjacenciesSet.find(new SVKmerAdjacencies(canonicalKmer));
        if ( kmerAdjacencies == null ) {
            throw new GATKException("can't find expected kmer in adjacencies set");
        }
        return kmer.equals(canonicalKmer) ? kmerAdjacencies : kmerAdjacencies.reverseComplement(kSize);
    }

    private static boolean isContigStart( final SVKmerAdjacencies kmerAdjacencies,
                                          final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet,
                                          final int kSize ) {
        final SVKmer predecessorKmer = kmerAdjacencies.getSolePredecessor(kSize);
        if ( predecessorKmer == null ) return true;
        final SVKmerAdjacencies predecessorAdjacencies =
                strandSensitiveLookup(predecessorKmer, kSize, kmerAdjacenciesSet);
        return predecessorAdjacencies.successorCount() > 1;
    }

    private static boolean isContigEnd( final SVKmerAdjacencies kmerAdjacencies,
                                        final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet,
                                        final int kSize ) {
        final SVKmer successorKmer = kmerAdjacencies.getSoleSuccessor(kSize);
        if ( successorKmer == null ) return true;
        final SVKmerAdjacencies successorAdjacencies =
                strandSensitiveLookup(successorKmer, kSize, kmerAdjacenciesSet);
        return successorAdjacencies.predecessorCount() > 1;
    }

    private static Contig buildContig( final SVKmerAdjacencies kmerAdjacencies, final int kSize,
                                       final HopscotchSet<SVKmerAdjacencies> kmerAdjacenciesSet ) {
        final StringBuilder contigSequence = new StringBuilder(kmerAdjacencies.getKmer().toString(kSize));
        SVKmerAdjacencies currentAdjacencies = kmerAdjacencies;
        SVKmer successorKmer;
        while ( (successorKmer = currentAdjacencies.getSoleSuccessor(kSize)) != null ) {
            final SVKmerAdjacencies successorAdjacencies =
                    strandSensitiveLookup(successorKmer, kSize, kmerAdjacenciesSet);
            if ( successorAdjacencies.predecessorCount() > 1 ) break;
            contigSequence.append(successorAdjacencies.getKmer().lastBase().name());
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
