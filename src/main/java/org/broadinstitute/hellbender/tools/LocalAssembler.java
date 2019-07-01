package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.PairWalker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.HopscotchSet;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.zip.GZIPOutputStream;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "experiment",
        oneLineSummary = "experiment",
        usageExample = "gatk LocalAssembler",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class LocalAssembler extends PairWalker {
    public static final byte QMIN = 25;
    public static final int MIN_THIN_OBS = 4;
    public static final int MIN_GAPFILL_COUNT = 3;
    public static final int TOO_MANY_TRAVERSALS = 100000;
    public static final int TOO_MANY_SCAFFOLDS = 50000;
    public static final int MIN_SV_SIZE = 50;

    @Argument(fullName=StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write outputs to this file name prefix", optional = true)
    public static String output;

    @Argument(fullName="assembly-name", doc="name of assembly used as a prefix for traversal names")
    public static String assemblyName;

    private final List<GATKRead> reads = new ArrayList<>();

    @Override public boolean requiresIntervals() { return true; }

    @Override public void apply( final GATKRead read, final GATKRead mate ) {
        trimOverruns(read, mate);
        reads.add(read);
        reads.add(mate);
    }

    @Override public void applyUnpaired( final GATKRead read ) {
        reads.add(read);
    }

    @Override public Object onTraversalSuccess() {
        super.onTraversalSuccess(); // flush any incomplete pairs

        final int regionSize = getTraversalIntervals().stream().mapToInt(SimpleInterval::size).sum();
        final KmerSet<KmerAdjacency> kmerAdjacencySet = new KmerSet<>(10 * regionSize);
        kmerizeReads(kmerAdjacencySet);
        List<ContigImpl> contigs = buildContigs(kmerAdjacencySet);
        connectContigs(contigs);

        removeThinContigs(contigs, kmerAdjacencySet);
        weldPipes(contigs);
        markComponents(contigs);

        if ( fillGaps(kmerAdjacencySet) ) {
            contigs = buildContigs(kmerAdjacencySet);
            connectContigs(contigs);
            removeThinContigs(contigs, kmerAdjacencySet);
            weldPipes(contigs);
            markComponents(contigs);
        }

        markCycles(contigs);

        final String outputFilePrefix = output != null ? output : assemblyName;
        final List<Path> readPaths = pathReads(kmerAdjacencySet);
        final Map<Contig,List<TransitPairCount>> contigTransitsMap =
                collectTransitPairCounts(contigs, readPaths);
        final String traversalsFilename = outputFilePrefix + ".traversals.fa.gz";
        try {
            final List<Traversal> allTraversals =
                    new ArrayList<>(traverseAllPaths(contigs, readPaths, contigTransitsMap));
            writeTraversals(allTraversals, traversalsFilename);
            try {
                final String scaffoldsFileName = outputFilePrefix + ".scaffolds.fa.gz";
                writeTraversals(createScaffolds(allTraversals), scaffoldsFileName);
            } catch ( final AssemblyTooComplexException x ) {
                logger.warn("Assembly too complex for scaffolding.");
            }
        } catch ( final AssemblyTooComplexException x ) {
            logger.warn("Assembly too complex.  Writing contigs as traversals in " +
                    traversalsFilename + ".");
            final Collection<Traversal> contigTraversals = new ArrayList<>(contigs.size());
            for ( final Contig contig : contigs ) {
                contigTraversals.add(new Traversal(Collections.singletonList(contig)));
            }
            writeTraversals(contigTraversals, traversalsFilename);
        }

        contigs.sort(Comparator.comparingInt(ContigImpl::getId));
        writeDOT(contigs, outputFilePrefix + ".assembly.dot");
        writeContigs(contigs, outputFilePrefix + ".contigs.txt.gz");
        writePaths(readPaths, outputFilePrefix + ".paths.txt.gz");
        writeReads(reads, outputFilePrefix + ".reads.fastq.gz");
        return null;
    }

    /** trim read pairs of base calls that have gone past the end of a short fragment */
    private void trimOverruns( final GATKRead read, final GATKRead mate ) {
        // if both mapped and they're on different strands
        if ( !read.isUnmapped() && !mate.isUnmapped() &&
                read.isReverseStrand() != mate.isReverseStrand() ) {
            // and both start within 1 base on the ref
            if ( Math.abs(read.getStart() - read.getMateStart()) <= 1 ) {
                // and both end within 1 base
                final int readRefLen = read.getCigar().getReferenceLength();
                final int mateRefLen = mate.getCigar().getReferenceLength();
                if ( Math.abs(readRefLen - mateRefLen) <= 1 ) {
                    if ( mate.isReverseStrand() ) {
                        trimClips(read, mate);
                    } else {
                        trimClips(mate, read);
                    }
                }
            }
        }
    }

    private void trimClips( final GATKRead fwd, final GATKRead rev ) {
        final List<CigarElement> fwdElements = fwd.getCigarElements();
        final List<CigarElement> revElements = rev.getCigarElements();
        final int lastElementIdx = fwdElements.size() - 1;
        final CigarElement fwdLastElement = fwdElements.get(lastElementIdx);
        final CigarElement revFirstElement = revElements.get(0);
        if ( fwdLastElement.getOperator() == CigarOperator.S &&
                revFirstElement.getOperator() == CigarOperator.S ) {
            final byte[] fwdBases = fwd.getBasesNoCopy();
            final int lastElementLen = fwdLastElement.getLength();
            fwd.setBases(Arrays.copyOfRange(fwdBases, 0, fwdBases.length - lastElementLen));
            final byte[] fwdQuals = fwd.getBaseQualitiesNoCopy();
            if ( fwdQuals.length > 0 ) {
                final int qualsLen = fwdQuals.length - lastElementLen;
                fwd.setBaseQualities(Arrays.copyOfRange(fwdQuals, 0, qualsLen));
            }
            final List<CigarElement> newFwdElements = new ArrayList<>(fwdElements);
            newFwdElements.set(lastElementIdx, new CigarElement(lastElementLen, CigarOperator.H));
            fwd.setCigar(new Cigar(newFwdElements));

            final byte[] revBases = rev.getBasesNoCopy();
            final int firstElementLen = revFirstElement.getLength();
            rev.setBases(Arrays.copyOfRange(revBases, firstElementLen, revBases.length));
            final byte[] revQuals = rev.getBaseQualitiesNoCopy();
            if ( revQuals.length > 0 ) {
                rev.setBaseQualities(Arrays.copyOfRange(revQuals, firstElementLen, revQuals.length));
            }
            final List<CigarElement> newRevElements = new ArrayList<>(revElements);
            newRevElements.set(0, new CigarElement(firstElementLen, CigarOperator.H));
            rev.setCigar(new Cigar(newRevElements));
        }
    }

    private void kmerizeReads( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        for ( final GATKRead read : reads ) {
            final byte[] calls = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            KmerAdjacency.kmerize(calls, quals, QMIN, kmerAdjacencySet);
        }
    }

    private static List<ContigImpl> buildContigs( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final List<ContigImpl> contigs = new ArrayList<>();
        ContigImpl.nContigs = 1;
        for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
            if ( kmerAdjacency.getContig() == null ) {
                ContigImpl contig = null;
                final KmerAdjacency predecessor = kmerAdjacency.getSolePredecessor();
                if ( predecessor == null ||
                        predecessor.getSuccessorCount() > 1 ||
                        predecessor == kmerAdjacency.rc() ) {
                    contig = new ContigImpl(kmerAdjacency);
                } else {
                    final KmerAdjacency successor = kmerAdjacency.getSoleSuccessor();
                    if ( successor == null ||
                            successor.getPredecessorCount() > 1 ||
                            successor == kmerAdjacency.rc() ) {
                        contig = new ContigImpl(kmerAdjacency.rc());
                    }
                }
                if ( contig != null ) {
                    contigs.add(contig);
                }
            }
        }

        // smooth circles?
        for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
            if ( kmerAdjacency.getContig() == null ) {
                contigs.add(new ContigImpl(kmerAdjacency));
            }
        }

        return contigs;
    }

    private static void connectContigs( final List<ContigImpl> contigs ) {
        final int nContigs = contigs.size();
        final KmerSet<ContigEndKmer> contigEnds = new KmerSet<>(2*nContigs);
        for ( int contigId = 0; contigId != nContigs; ++contigId ) {
            final ContigImpl contig = contigs.get(contigId);
            final KmerAdjacency fwdKmer = contig.getFirstKmer();
            final KmerAdjacency revKmer = contig.getLastKmer().rc();
            if ( fwdKmer == revKmer ) {
                contigEnds.add(new ContigEndKmer(fwdKmer.getKVal(), contig, ContigOrientation.BOTH));
            } else {
                contigEnds.add(new ContigEndKmer(fwdKmer.getKVal(), contig, ContigOrientation.FWD));
                contigEnds.add(new ContigEndKmer(revKmer.getKVal(), contig, ContigOrientation.REV));
            }
        }

        for ( int contigId = 0; contigId != nContigs; ++contigId ) {
            final Contig contig = contigs.get(contigId);

            final KmerAdjacency start = contig.getFirstKmer();
            final int predecessorCount = start.getPredecessorCount();
            if ( predecessorCount > 0 ) {
                final List<Contig> predecessors = contig.getPredecessors();
                final int mask = start.getPredecessorMask();
                for ( int call = 0; call != 4; ++call ) {
                    if ( (mask & (1 << call)) != 0 ) {
                        final long kVal =
                                KmerAdjacency.reverseComplement(start.getPredecessorVal(call));
                        final ContigEndKmer contigEndKmer = contigEnds.find(new Kmer(kVal));
                        switch ( contigEndKmer.getContigOrientation() ) {
                            case FWD:
                                predecessors.add(contigEndKmer.getContig().rc());
                                break;
                            case REV:
                                predecessors.add(contigEndKmer.getContig());
                                break;
                            case BOTH:
                                predecessors.add(contigEndKmer.getContig());
                                predecessors.add(contigEndKmer.getContig().rc());
                                break;
                        }
                    }
                }
            }

            final KmerAdjacency end = contig.getLastKmer();
            final int successorCount = end.getSuccessorCount();
            if ( successorCount > 0 ) {
                final List<Contig> successors = contig.getSuccessors();
                final int mask = end.getSuccessorMask();
                for ( int call = 0; call != 4; ++call ) {
                    if ( (mask & (1 << call)) != 0 ) {
                        final long kVal = end.getSuccessorVal(call);
                        final ContigEndKmer contigEndKmer = contigEnds.find(new Kmer(kVal));
                        switch ( contigEndKmer.getContigOrientation() ) {
                            case FWD:
                                successors.add(contigEndKmer.getContig());
                                break;
                            case REV:
                                successors.add(contigEndKmer.getContig().rc());
                                break;
                            case BOTH:
                                successors.add(contigEndKmer.getContig());
                                successors.add(contigEndKmer.getContig().rc());
                                break;
                        }
                    }
                }
            }
        }
    }

    private static void removeThinContigs( final List<ContigImpl> contigs,
                                           final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        contigs.sort(Comparator.comparingInt(ContigImpl::getMaxObservations));
        boolean contigRemoved;
        do {
            final int nContigs = contigs.size();
            final Map<Contig, CutData> cutDataMap = new HashMap<>(nContigs * 3);

            for ( final ContigImpl contig : contigs ) {
                if ( cutDataMap.containsKey(contig) ) {
                    continue;
                }

                cutDataMap.put(contig, new CutData());
                int children = 0;
                for ( final Contig nextContig : contig.getSuccessors() ) {
                    if ( !cutDataMap.containsKey(nextContig) ) {
                        findCuts(nextContig, contig, cutDataMap);
                        children += 1;
                    }
                }
                for ( final Contig nextContig : contig.getPredecessors() ) {
                    if ( !cutDataMap.containsKey(nextContig) ) {
                        findCuts(nextContig, contig, cutDataMap);
                        children += 1;
                    }
                }
                if ( children >= 2 ) {
                    contig.setCut(true);
                }
            }

            contigRemoved = false;
            final Iterator<ContigImpl> itr = contigs.iterator();
            while ( itr.hasNext() ) {
                final Contig contig = itr.next();
                if ( contig.getMaxObservations() < MIN_THIN_OBS && !contig.isCut() ) {
                    unlinkContig(contig, kmerAdjacencySet);
                    itr.remove();
                    contigRemoved = true;
                    break;
                }
            }
        } while ( contigRemoved );
        contigs.sort(Comparator.comparingInt(ContigImpl::getId));
    }

    private static CutData findCuts( final Contig contig,
                                     final Contig parent,
                                     final Map<Contig, CutData> cutDataMap ) {
        final CutData cutData = new CutData();
        cutDataMap.put(contig, cutData);
        for ( final Contig nextContig : contig.getSuccessors() ) {
            if ( nextContig == parent ) {
                continue;
            }
            CutData nextCutData = cutDataMap.get(nextContig);
            if ( nextCutData != null ) {
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.visitNum);
            } else {
                nextCutData = findCuts(nextContig, contig, cutDataMap);
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.minVisitNum);
                if ( nextCutData.minVisitNum >= cutData.visitNum ) {
                    contig.setCut(true);
                }
            }
        }
        for ( final Contig nextContig : contig.getPredecessors() ) {
            if ( nextContig == parent ) {
                continue;
            }
            CutData nextCutData = cutDataMap.get(nextContig);
            if ( nextCutData != null ) {
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.visitNum);
            } else {
                nextCutData = findCuts(nextContig, contig, cutDataMap);
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.minVisitNum);
                if ( nextCutData.minVisitNum >= cutData.visitNum ) {
                    contig.setCut(true);
                }
            }
        }
        return cutData;
    }

    private static void unlinkContig( final Contig contig,
                                      final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final KmerAdjacency firstKmer = contig.getFirstKmer();
        final int firstKmerFinalCall = firstKmer.getFinalCall();
        for ( final Contig predecessor : contig.getPredecessors() ) {
            if ( predecessor != contig && predecessor != contig.rc() ) {
                predecessor.getLastKmer().removeSuccessor(firstKmerFinalCall, kmerAdjacencySet);
                if ( !predecessor.getSuccessors().remove(contig) ) {
                    throw new GATKException("failed to find predecessor link");
                }
            }
        }

        final KmerAdjacency lastKmer = contig.getLastKmer();
        final int lastKmerInitialCall = lastKmer.getInitialCall();
        for ( final Contig successor : contig.getSuccessors() ) {
            if ( successor != contig && successor != contig.rc() ) {
                successor.getFirstKmer().removePredecessor(lastKmerInitialCall, kmerAdjacencySet);
                if ( !successor.getPredecessors().remove(contig) ) {
                    throw new GATKException("failed to find successor link");
                }
            }
        }

        KmerAdjacency nextKmer = firstKmer;
        KmerAdjacency kmer;
        do {
            kmer = nextKmer;
            nextKmer = kmer.getSoleSuccessor();
            kmerAdjacencySet.remove(kmer.canonical());
        } while ( kmer != lastKmer );
    }

    private static void updateKmerContig( final KmerAdjacency firstKmer,
                                          final KmerAdjacency lastKmer,
                                          final Contig contig ) {
        int offset = 0;
        for ( KmerAdjacency kmer = firstKmer; kmer != lastKmer; kmer = kmer.getSoleSuccessor() ) {
            if ( kmer == null ) {
                throw new GATKException("contig does not have a flat pipeline of kmers");
            }
            kmer.setContig(null, 0); // erase the old info first or we'll throw
            kmer.setContig(contig, offset++);
        }
        lastKmer.setContig(null, 0);
        lastKmer.setContig(contig, offset);
        if ( offset + Kmer.KSIZE != contig.size() ) {
            throw new GATKException("kmer chain length does not equal contig size");
        }
    }

    private static void weldPipes( final List<ContigImpl> contigs ) {
        for ( int contigId = 0; contigId != contigs.size(); ++contigId ) {
            final ContigImpl contig = contigs.get(contigId);
            if ( contig.getSuccessors().size() == 1 ) {
                final Contig successor = contig.getSuccessors().get(0);
                if ( successor != contig && successor != contig.rc() &&
                        successor.getPredecessors().size() == 1 ) {
                    contigs.set(contigId, join(contig, successor));
                    if ( !contigs.remove(successor.canonical()) ) {
                        throw new GATKException("successor linkage is messed up");
                    }
                    contigId -= 1; // reconsider the new contig -- there might be more joining possible
                    continue;
                }
            }
            if ( contig.getPredecessors().size() == 1 ) {
                final Contig predecessor = contig.getPredecessors().get(0);
                if ( predecessor != contig && predecessor != contig.rc() &&
                        predecessor.getSuccessors().size() == 1 ) {
                    contigs.set(contigId, join(predecessor, contig));
                    if ( !contigs.remove(predecessor.canonical()) ) {
                        throw new GATKException("predecessor linkage is messed up");
                    }
                    contigId -= 1; // reconsider
                }
            }
        }
    }

    private static ContigImpl join( final Contig predecessor, final Contig successor ) {
        if ( !checkOverlap(predecessor.getSequence(), successor.getSequence()) ) {
                throw new GATKException("sequences can't be joined");
        }
        final ContigImpl joinedContig = new ContigImpl(predecessor, successor);
        updateKmerContig(joinedContig.getFirstKmer(), joinedContig.getLastKmer(), joinedContig);
        return joinedContig;
    }

    private static boolean checkOverlap( final CharSequence seq1, final CharSequence seq2 ) {
        final int seq1Len = seq1.length();
        final CharSequence seq1SubSeq = seq1.subSequence(seq1Len - Kmer.KSIZE + 1, seq1Len);
        final CharSequence seq2SubSeq = seq2.subSequence(0, Kmer.KSIZE - 1);
        return seq1SubSeq.equals(seq2SubSeq);
    }

    private static void markComponents( final List<ContigImpl> contigs ) {
        for ( final ContigImpl contig : contigs ) {
            contig.setComponentId(0);
        }

        int componentId = 0;
        for ( final ContigImpl contig : contigs ) {
            if ( contig.getComponentId() == 0 ) {
                contig.setComponentId(++componentId);
                markSuccessorComponents(contig);
                markSuccessorComponents(contig.rc());
            }
        }
    }

    private static void markSuccessorComponents( final Contig contig ) {
        final int componentId = contig.getComponentId();
        for ( final Contig successor : contig.getSuccessors() ) {
            if ( successor.getComponentId() == 0 ) {
                successor.canonical().setComponentId(componentId);
                markSuccessorComponents(successor);
                markSuccessorComponents(successor.rc());
            }
        }
    }

    private static void markCycles( final List<ContigImpl> contigs ) {
        for ( final Contig contig : contigs ) {
            contig.setCyclic(false);
        }

        final int nContigs = contigs.size();
        final Deque<Contig> deque = new ArrayDeque<>(nContigs);
        final Map<Contig, CutData> cutDataMap = new HashMap<>(nContigs * 3);
        for ( final Contig contig : contigs ) {
            if ( !cutDataMap.containsKey(contig) ) {
                markCyclesRecursion(contig, deque, cutDataMap);
            }
        }
    }

    private static CutData markCyclesRecursion( final Contig contig,
                                                final Deque<Contig> deque,
                                                final Map<Contig, CutData> cutDataMap ) {
        final CutData cutData = new CutData();
        cutDataMap.put(contig, cutData);
        deque.addFirst(contig);

        for ( final Contig successor : contig.getSuccessors() ) {
            final CutData successorCutData = cutDataMap.get(successor);
            if ( successorCutData == null ) {
                final int recursionVisitNum =
                        markCyclesRecursion(successor, deque, cutDataMap).minVisitNum;
                cutData.minVisitNum = Math.min(cutData.minVisitNum, recursionVisitNum);
            } else {
                cutData.minVisitNum = Math.min(cutData.minVisitNum, successorCutData.visitNum);
            }
        }

        if ( cutData.visitNum == cutData.minVisitNum ) {
            Contig tig = deque.removeFirst();
            if ( tig == contig ) {
                cutDataMap.get(tig).visitNum = Integer.MAX_VALUE;

                // single-vertex component -- cyclic only if self-referential
                if ( tig.getSuccessors().contains(tig) ) {
                    tig.setCyclic(true);
                }
            } else {
                while ( true ) {
                    // kill cross-links
                    cutDataMap.get(tig).visitNum = Integer.MAX_VALUE;
                    tig.setCyclic(true);
                    if ( tig == contig ) break;
                    tig = deque.removeFirst();
                }
            }
        }
        return cutData;
    }

    private boolean fillGaps( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final Map<String, Integer> gapFillCounts = new HashMap<>();
        for ( final GATKRead read : reads ) {
            final Path path = new Path(read.getBasesNoCopy(), kmerAdjacencySet);
            final List<PathPart> parts = path.getParts();
            final int lastIdx = parts.size() - 1;
            for ( int idx = 1; idx < lastIdx; ++idx ) {
                final PathPart pathPart = parts.get(idx);
                if ( pathPart.isGap() ) {
                    final char prevCall = parts.get(idx - 1).getLastCall();
                    final char nextCall = parts.get(idx + 1).getFirstCall();
                    String gapFill = prevCall + pathPart.getSequence().toString() + nextCall;
                    final SequenceRC gapFillRC = new SequenceRC(gapFill);
                    if ( gapFillRC.compareTo(gapFill) < 0 ) {
                        gapFill = gapFillRC.toString();
                    }
                    gapFillCounts.merge(gapFill, 1, Integer::sum);
                }
            }
        }

        boolean newKmers = false;
        for ( final Map.Entry<String, Integer> entry : gapFillCounts.entrySet() ) {
            final int nObservations = entry.getValue();
            if ( nObservations >= MIN_GAPFILL_COUNT ) {
                KmerAdjacency.kmerize(entry.getKey(), nObservations, kmerAdjacencySet);
                newKmers = true;
            }
        }

        if ( newKmers ) {
            for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
                kmerAdjacency.setContig(null, 0);
            }
        }
        return newKmers;
    }

    private List<Path> pathReads( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final List<Path> readPaths = new ArrayList<>(reads.size());
        for ( final GATKRead read : reads ) {
            readPaths.add(new Path(read.getBasesNoCopy(), kmerAdjacencySet));
        }
        return readPaths;
    }

    private static Map<Contig,List<TransitPairCount>> collectTransitPairCounts(
            final List<ContigImpl> contigs,
            final List<Path> readPaths ) {
        final Map<Contig,List<TransitPairCount>> contigTransitsMap =
                new HashMap<>(3 * contigs.size());
        for ( final Path path : readPaths ) {
            final List<PathPart> parts = path.getParts();
            final int lastPart = parts.size() - 1;
            for ( int partIdx = 1; partIdx < lastPart; ++partIdx ) {
                final Contig prevContig = parts.get(partIdx - 1).getContig();
                if ( prevContig == null ) continue;
                final Contig curContig = parts.get(partIdx).getContig();
                if ( curContig == null ) {
                    partIdx += 1;
                    continue;
                }
                final Contig nextContig = parts.get(partIdx + 1).getContig();
                if ( nextContig == null ) {
                    partIdx += 2;
                    continue;
                }
                final TransitPairCount tpc = new TransitPairCount(prevContig, nextContig);
                final List<TransitPairCount> tpcList =
                        contigTransitsMap.computeIfAbsent(curContig, tig -> new ArrayList<>(4));
                final int idx = tpcList.indexOf(tpc);
                if ( idx != -1 ) {
                    tpcList.get(idx).observe();
                } else {
                    tpcList.add(tpc);
                    contigTransitsMap.computeIfAbsent(curContig.rc(), tig -> new ArrayList<>(4))
                            .add(tpc.getRC());
                }
            }
        }
        return contigTransitsMap;
    }

    private static Set<Traversal> traverseAllPaths(
            final List<ContigImpl> contigs,
            final List<Path> readPaths,
            final Map<Contig, List<TransitPairCount>> contigTransitsMap ) {
        final Set<Traversal> traversalSet = new HashSet<>();
        final List<Contig> contigsList = new ArrayList<>();
        for ( final Contig contig : contigs ) {
            // untransited contigs are sources, sinks, or large contigs that can't be crossed by a read
            // build traversals from these
            if ( !contigTransitsMap.containsKey(contig) ) {
                boolean done = false;
                for ( final Contig successor : contig.getSuccessors() ) {
                    traverse(successor, contig,
                            contigsList, readPaths, contigTransitsMap, traversalSet);
                    done = true;
                }
                for ( final Contig predecessor : contig.getPredecessors() ) {
                    traverse(predecessor.rc(), contig.rc(),
                            contigsList, readPaths, contigTransitsMap, traversalSet);
                    done = true;
                }
                if ( !done ) { // if there were no predecessors or successors, it stands alone
                    addTraversal(new Traversal(Collections.singletonList(contig)), traversalSet);
                }
            }
        }

        // look for transits that haven't been traced
        for ( final Map.Entry<Contig, List<TransitPairCount>> entry :
                contigTransitsMap.entrySet() ) {
            for ( final TransitPairCount tpc : entry.getValue() ) {
                if ( tpc.getCount() > 0 ) {
                    tpc.resetCount();
                    final Contig contig = entry.getKey();
                    final Set<Traversal> fwdTraversalSet = new HashSet<>();
                    traverse(tpc.getNextContig(), contig,
                            contigsList, readPaths, contigTransitsMap, fwdTraversalSet);
                    final Set<Traversal> revTraversalSet = new HashSet<>();
                    traverse(tpc.getPrevContig().rc(), contig.rc(),
                            contigsList, readPaths, contigTransitsMap, revTraversalSet);
                    for ( final Traversal revTraversal : revTraversalSet ) {
                        final Traversal revTraversalRC = revTraversal.rc();
                        for ( final Traversal fwdTraversal : fwdTraversalSet ) {
                            final Traversal combo = Traversal.combine(revTraversalRC, fwdTraversal);
                            addTraversal(combo, traversalSet);
                        }
                    }
                }
            }
        }

        return traversalSet;
    }

    private static void traverse( final Contig contig,
                                  final Contig predecessor,
                                  final List<Contig> contigsList,
                                  final List<Path> readPaths,
                                  final Map<Contig, List<TransitPairCount>> contigTransitsMap,
                                  final Set<Traversal> traversalSet ) {
        contigsList.add(predecessor);
        if ( contig.isCyclic() ) {
            final int cycleIdx = contigsList.indexOf(contig);
            if ( cycleIdx != -1 ) {
                traverseCycle(cycleIdx, contig, contigsList, readPaths, contigTransitsMap, traversalSet);
                contigsList.remove(contigsList.size() - 1);
                return;
            }
        }
        final List<TransitPairCount> transits = contigTransitsMap.get(contig);
        boolean done = false;
        if ( transits != null ) {
            for ( final TransitPairCount tpc : transits ) {
                if ( tpc.getPrevContig() == predecessor ) {
                    final Contig successor = tpc.getNextContig();
                    if ( predecessor == contig.rc() ) {
                        final int nContigs = contigsList.size();
                        if ( nContigs > 1 ) {
                            if ( successor.rc() == contigsList.get(nContigs - 2) ) {
                                continue;
                            }
                        }
                    }
                    tpc.resetCount();
                    traverse(successor, contig, contigsList, readPaths, contigTransitsMap, traversalSet);
                    done = true;
                }
            }
        }
        if ( !done ) {
            contigsList.add(contig);
            addTraversal(new Traversal(contigsList), traversalSet);
            contigsList.remove(contigsList.size() - 1);
        }
        contigsList.remove(contigsList.size() - 1);
    }

    private static void traverseCycle( final int cycleIdx,
                                       final Contig contig,
                                       final List<Contig> contigsList,
                                       final List<Path> readPaths,
                                       final Map<Contig, List<TransitPairCount>> contigTransitsMap,
                                       final Set<Traversal> traversalSet ) {
        contigsList.add(contig);
        final int startIdx = cycleIdx > 0 ? cycleIdx - 1 : 0;
        final List<List<Contig>> longestPaths =
                findLongestPaths(contigsList.subList(startIdx, contigsList.size()), readPaths);
        if ( longestPaths.isEmpty() ) {
            addTraversal(new Traversal(contigsList, true), traversalSet);
        } else {
            for ( final List<Contig> path : longestPaths ) {
                if ( path.isEmpty() ) {
                    addTraversal(new Traversal(contigsList, true), traversalSet);
                    continue;
                }
                final List<Contig> extendedContigsList =
                        new ArrayList<>(contigsList.size() + path.size());
                extendedContigsList.addAll(contigsList);
                if ( path.get(path.size() - 1).isCyclic() ) {
                    extendedContigsList.addAll(path);
                    addTraversal(new Traversal(extendedContigsList, true), traversalSet);
                } else {
                    for ( final Contig curContig : path ) {
                        if ( curContig.isCyclic() ) {
                            extendedContigsList.add(curContig);
                        } else {
                            final Contig prevContig =
                                    extendedContigsList.remove(extendedContigsList.size() - 1);
                            traverse(curContig, prevContig, extendedContigsList, readPaths,
                                    contigTransitsMap, traversalSet);
                        }
                    }
                }
                clearTransitPairs(contigTransitsMap, extendedContigsList);
            }
        }
        contigsList.remove(contigsList.size() - 1);
    }

    private static void clearTransitPairs(
            final Map<Contig, List<TransitPairCount>> contigTransitsMap,
            final List<Contig> contigsList ) {
        final int lastIdx = contigsList.size() - 1;
        for ( int idx = 1; idx < lastIdx; ++idx ) {
            final List<TransitPairCount> pairCounts = contigTransitsMap.get(contigsList.get(idx));
            if ( pairCounts != null ) {
                final Contig predecessor = contigsList.get(idx - 1);
                final Contig successor = contigsList.get(idx + 1);
                for ( final TransitPairCount tpc : pairCounts ) {
                    if ( tpc.getPrevContig() == predecessor && tpc.getNextContig() == successor ) {
                        tpc.resetCount();
                        break;
                    }
                }
            }
        }
    }

    private static void addTraversal( final Traversal traversal,
                                      final Set<Traversal> traversalSet ) {
        if ( !traversalSet.contains(traversal.rc()) ) {
            traversalSet.add(traversal);
            if ( traversalSet.size() >= TOO_MANY_TRAVERSALS ) {
                throw new AssemblyTooComplexException();
            }
        }
    }

    private static List<List<Contig>> findLongestPaths( final List<Contig> toMatch,
                                                        final List<Path> readPaths ) {
        final List<List<Contig>> results = new ArrayList<>();
        for ( final Path path : readPaths ) {
            testPath(path, toMatch, results);
            testPath(path.rc(), toMatch, results);
        }
        return results;
    }

    private static void testPath( final Path path,
                                  final List<Contig> toMatch,
                                  final List<List<Contig>> results ) {
        final List<PathPart> pathParts = path.getParts();
        final int nPathParts = pathParts.size();
        final List<Contig> pathContigs = new ArrayList<>(nPathParts);
        pathParts.forEach(pp -> pathContigs.add(pp.getContig()));
        final int matchIdx = Collections.indexOfSubList(pathContigs, toMatch);
        if ( matchIdx != -1 ) {
            final int suffixIdx = matchIdx + toMatch.size();
            if ( suffixIdx < nPathParts ) {
                resolveResult(grabParts(pathContigs, suffixIdx), results);
            }
        }
    }

    private static List<Contig> grabParts( final List<Contig> pathContigs, final int suffixIdx ) {
        final int nPathContigs = pathContigs.size();
        Contig prev = pathContigs.get(suffixIdx - 1);
        final List<Contig> result = new ArrayList<>(nPathContigs - suffixIdx);
        for ( int idx = suffixIdx; idx != nPathContigs; ++idx ) {
            final Contig tig = pathContigs.get(idx);
            if ( tig == null || !prev.getSuccessors().contains(tig) ) break;
            result.add(tig);
            prev = tig;
        }
        return result;
    }

    private static void resolveResult( final List<Contig> result,
                                       final List<List<Contig>> results ) {
        final int nResults = results.size();
        for ( int idx = 0; idx != nResults; ++idx ) {
            final List<Contig> test = results.get(idx);
            if ( isPrefix(result, test) ) return;
            if ( isPrefix(test, result) ) {
                results.set(idx, result);
                return;
            }
        }
        results.add(result);
    }

    private static boolean isPrefix( final List<Contig> list1, final List<Contig> list2 ) {
        final int list1Size = list1.size();
        final int list2Size = list2.size();
        if ( list1Size > list2Size ) return false;
        for ( int idx = 0; idx != list1Size; ++idx ) {
            if ( list1.get(idx) != list2.get(idx) ) return false;
        }
        return true;
    }

    private static Collection<Traversal> createScaffolds( final List<Traversal> allTraversals ) {
        removeTriviallyDifferentTraversals(allTraversals);

        final int nTraversals = allTraversals.size();
        final Map<Contig, List<Integer>> traversalsByFirstContig = new HashMap<>(3 * nTraversals);
        for ( int idx = 0; idx != nTraversals; ++idx ) {
            final Traversal traversal = allTraversals.get(idx);
            traversalsByFirstContig.compute(traversal.getFirstContig(),
                    ( k, v ) -> v == null ? new ArrayList<>(3) : v).add(idx);
            final Traversal rcTraversal = traversal.rc();
            traversalsByFirstContig.compute(rcTraversal.getFirstContig(),
                    ( k, v ) -> v == null ? new ArrayList<>(3) : v).add(~idx);
        }

        final List<Traversal> scaffolds = new ArrayList<>(nTraversals);
        final boolean[] touched = new boolean[nTraversals];
        for ( int idx = 0; idx != nTraversals; ++idx ) {
            if ( !touched[idx] ) {
                expandTraversal(idx, touched, traversalsByFirstContig, allTraversals, scaffolds);
            }
        }
        return scaffolds;
    }

    private static void expandTraversal( final int traversalIdx,
                                         final boolean[] touched,
                                         final Map<Contig, List<Integer>> traversalsByFirstContig,
                                         final List<Traversal> allTraversals,
                                         final List<Traversal> scaffolds ) {
        final Traversal traversal = allTraversals.get(traversalIdx);
        touched[traversalIdx] = true;
        final List<Traversal> downExtensions = new ArrayList<>();
        final Set<Contig> startingContigSet = new HashSet<>();
        walkTraversals(traversal, touched, startingContigSet, traversalsByFirstContig,
                        allTraversals, downExtensions);
        final List<Traversal> upExtensions = new ArrayList<>();
        walkTraversals(traversal.rc(), touched, startingContigSet, traversalsByFirstContig,
                        allTraversals, upExtensions);
        for ( final Traversal down : downExtensions ) {
            for ( final Traversal up : upExtensions ) {
                if ( scaffolds.size() >= TOO_MANY_SCAFFOLDS ) {
                    throw new AssemblyTooComplexException();
                }
                scaffolds.add(Traversal.combineOverlappers(up.rc(), down, traversal.getContigs().size()));
            }
        }
    }

    private static void walkTraversals( final Traversal traversal,
                                        final boolean[] touched,
                                        final Set<Contig> startingContigSet,
                                        final Map<Contig, List<Integer>> traversalsByFirstContig,
                                        final List<Traversal> allTraversals,
                                        final List<Traversal> extensions ) {
        final Contig firstContig = traversal.getFirstContig();
        final List<Integer> indexList;
        if ( startingContigSet.contains(firstContig) ||
                traversal.isInextensible() ||
                (indexList = traversalsByFirstContig.get(traversal.getLastContig())) == null ) {
            extensions.add(traversal);
            return;
        }
        startingContigSet.add(firstContig);
        for ( int idx : indexList ) {
            final Traversal extension;
            if ( idx >= 0 ) {
                extension = allTraversals.get(idx);
                touched[idx] = true;
            } else {
                final int rcIdx = ~idx;
                extension = allTraversals.get(rcIdx).rc();
                touched[rcIdx] = true;
            }
            walkTraversals(Traversal.combine(traversal, extension), touched, startingContigSet,
                            traversalsByFirstContig, allTraversals, extensions );
        }
        startingContigSet.remove(firstContig);
    }

    private static void removeTriviallyDifferentTraversals(
                                            final Collection<Traversal> allTraversals ) {
        if ( allTraversals.isEmpty() ) {
            return;
        }
        final TreeSet<Traversal> sortedTraversals = new TreeSet<>(new TraversalEndpointComparator());
        for ( final Traversal traversal : allTraversals ) {
            sortedTraversals.add(traversal);
            sortedTraversals.add(traversal.rc());
        }
        final Iterator<Traversal> traversalIterator = sortedTraversals.iterator();
        Traversal prevTraversal = traversalIterator.next();
        while ( traversalIterator.hasNext() ) {
            final Traversal curTraversal = traversalIterator.next();
            if ( isTriviallyDifferent(prevTraversal, curTraversal) ) {
                traversalIterator.remove();
            } else {
                prevTraversal = curTraversal;
            }
        }
        sortedTraversals.removeIf(Traversal::isRC);
        allTraversals.clear();
        allTraversals.addAll(sortedTraversals);
    }

    private static boolean isTriviallyDifferent( final Traversal traversal1,
                                                 final Traversal traversal2 ) {
        final Contig firstContig1 = traversal1.getFirstContig();
        final Contig lastContig1 = traversal1.getLastContig();
        final Contig firstContig2 = traversal2.getFirstContig();
        final Contig lastContig2 = traversal2.getLastContig();
        if ( firstContig1 != firstContig2 || lastContig1 != lastContig2 ) {
            return false;
        }
        final int interiorSize1 = traversal1.getSequenceLength() - firstContig1.size() - lastContig1.size();
        final int interiorSize2 = traversal2.getSequenceLength() - firstContig2.size() - lastContig2.size();

        // if the path lengths are so different that one could harbor an SV, they're not trivially different
        if ( Math.abs(interiorSize1 - interiorSize2) >= MIN_SV_SIZE ) {
            return false;
        }

        // if the paths are small enough that there can't be an SV's worth of differences, they're trivially different
        final int maxInteriorSize = Math.max(interiorSize1, interiorSize2);
        if ( maxInteriorSize < MIN_SV_SIZE ) {
            return true;
        }

        // dang, maybe there's enough material in common that there can't be an SV's worth of differences
        // run a longest common subsequence algorithm to figure out the length of the common material
        // DP matrix holds length of common material
        final List<Contig> contigs1 = traversal1.getContigs();
        final int rowLen = contigs1.size() - 1;
        final int[][] rowPair = new int[2][];
        rowPair[0] = new int[rowLen];
        rowPair[1] = new int[rowLen];
        int pairIdx = 0;
        final List<Contig> contigs2 = traversal2.getContigs();
        final int nRows = contigs2.size() - 1;
        for ( int idx2 = 1; idx2 != nRows; ++idx2 ) {
            final int[] curRow = rowPair[pairIdx];
            final int[] prevRow = rowPair[pairIdx ^ 1];
            pairIdx ^= 1;

            final int id2 = contigs2.get(idx2).getId();
            for ( int idx1 = 1; idx1 != rowLen; ++idx1 ) {
                final Contig tig1 = contigs1.get(idx1);
                if ( tig1.getId() == id2 ) {
                    // if the previous cells also contain a match we've already removed the K-1 bases upstream
                    final boolean extendMatch =
                            contigs1.get(idx1 -1).getId() == contigs2.get(idx2 - 1).getId();
                    curRow[idx1] = prevRow[idx1 - 1] + (extendMatch ? tig1.getNKmers() : tig1.size());
                } else {
                    curRow[idx1] = Math.max(curRow[idx1 - 1], prevRow[idx1]);
                }
            }
        }
        final int commonLen = rowPair[pairIdx ^ 1][rowLen - 1];
        return (maxInteriorSize - commonLen) < MIN_SV_SIZE;
    }

    private static class TraversalEndpointComparator implements Comparator<Traversal> {
        @Override
        public int compare( final Traversal traversal1, final Traversal traversal2 ) {
            int cmp = Integer.compare(traversal1.contigs.get(0).getId(),
                                      traversal2.contigs.get(0).getId());
            if ( cmp != 0 ) {
                return cmp;
            }
            final int last1 = traversal1.contigs.size() - 1;
            final int last2 = traversal2.contigs.size() - 1;
            cmp = Integer.compare(traversal1.contigs.get(last1).getId(),
                                  traversal2.contigs.get(last2).getId());
            if ( cmp != 0 ) {
                return cmp;
            }
            // among those starting and ending at the same place, sort least observed last
            return -Integer.compare(traversal1.getMinMaxObservations(), traversal2.getMinMaxObservations());
        }
    }

    private static void writeDOT( final List<ContigImpl> contigs, final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            writer.write("digraph {\n");
            for ( final Contig contig : contigs ) {
                final double width = contig.getSequence().length() / 100.;
                writer.write(contig + " [width=" + width + "]\n");
                writer.write( contig.rc() + " [width=" + width + "]\n");
            }
            for ( final Contig contig : contigs ) {
                for ( final Contig predecessor : contig.getPredecessors() ) {
                    final String predecessorName = predecessor.rc().toString();
                    writer.write(contig.rc() + " -> " + predecessorName + "\n");
                }
                for ( final Contig successor : contig.getSuccessors() ) {
                    final String successorName = successor.toString();
                    writer.write(contig + " -> " + successorName + "\n");
                }
            }
            writer.write("}\n");
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly DOT file.", ioe);
        }
    }

    private static BufferedWriter makeGZFile( final String fileName ) throws IOException {
        final GZIPOutputStream gzOS = new GZIPOutputStream(BucketUtils.createFile(fileName));
        return new BufferedWriter(new OutputStreamWriter(gzOS));
    }

    private static void writeContigs( final List<ContigImpl> contigs, final String fileName ) {
        try ( final BufferedWriter writer = makeGZFile(fileName) ) {
            for ( final Contig contig : contigs ) {
                final List<Contig> predecessors = contig.getPredecessors();
                final String predecessorDescription;
                if ( predecessors.size() == 0 ) {
                    predecessorDescription = "\tnone";
                } else {
                    final StringBuilder sb = new StringBuilder();
                    char prefix = '\t';
                    for ( final Contig predecessor : predecessors ) {
                        sb.append(prefix);
                        prefix = ',';
                        sb.append(predecessor);
                    }
                    predecessorDescription = sb.toString();
                }

                final List<Contig> successors = contig.getSuccessors();
                final String successorDescription;
                if ( successors.size() == 0 ) {
                    successorDescription = "\tnone";
                } else {
                    final StringBuilder sb = new StringBuilder();
                    char prefix = '\t';
                    for ( final Contig successor : successors ) {
                        sb.append(prefix);
                        prefix = ',';
                        sb.append(successor);
                    }
                    successorDescription = sb.toString();
                }

                final String contigName = contig.toString();
                final String component =
                        (contig.isCyclic() ? "(C)\t" : "\t") + contig.getComponentId();
                writer.write(contigName + component + predecessorDescription +
                        successorDescription + "\t" +
                        contig.getMaxObservations() + "\t" +
                        contig.getFirstKmer().getNObservations() + "\t" +
                        contig.getLastKmer().getNObservations() + "\t" +
                        contig.size() + "\t" +
                        contig.getSequence() + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write contigs file.", ioe);
        }
    }

    private static void writePaths( final List<Path> readPaths, final String fileName ) {
        try ( final BufferedWriter writer = makeGZFile(fileName) ) {
            final int nReads = readPaths.size();
            for ( int readId = 0; readId != nReads; ++readId ) {
                final Path path = readPaths.get(readId);
                final String pathDesc = path.toString();
                writer.write((readId + 1) + ": " + pathDesc + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write paths file.", ioe);
        }
    }

    private static void writeReads( final List<GATKRead> reads, final String fileName ) {
        try ( final BufferedWriter writer = makeGZFile(fileName) ) {
            for ( final GATKRead read : reads ) {
                writer.write("@" + read.getName());
                writer.write('\n');
                writer.write(new String(read.getBasesNoCopy()));
                writer.write("\n+\n");
                final byte[] quals = read.getBaseQualitiesNoCopy();
                final int nQuals = quals.length;
                final byte[] fastqQuals = new byte[nQuals];
                for ( int idx = 0; idx != nQuals; ++idx ) {
                    fastqQuals[idx] = (byte)SAMUtils.phredToFastq(quals[idx]);
                }
                writer.write(new String(fastqQuals));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly sam file.", ioe);
        }
    }

    private static void writeTraversals( final Collection<Traversal> traversals,
                                         final String fileName ) {
        try ( final BufferedWriter writer = makeGZFile(fileName) ) {
            int traversalNo = 0;
            for ( final Traversal traversal : traversals ) {
                writer.write(">");
                if ( assemblyName != null ) {
                    writer.write(assemblyName);
                    writer.write("_");
                }
                writer.write("t");
                writer.write(Integer.toString(++traversalNo));
                writer.write(" ");
                writer.write(traversal.getName());
                writer.newLine();
                writer.write(traversal.getSequence());
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly sam file.", ioe);
        }
    }

    public static class Kmer {
        public static final int KSIZE = 31;
        public static final long KMASK = (1L << 2*KSIZE) - 1L;
        private final long kVal;

        public Kmer( final long kVal ) { this.kVal = kVal; }

        public long getKVal() { return kVal; }
        public boolean isCanonical() { return isCanonical(kVal); }
        public int getInitialCall() { return (int)(kVal >> (KSIZE*2 - 2)) & 3; }
        public int getFinalCall() { return (int)kVal & 3; }

        public long getPredecessorVal( final int call ) {
            return (kVal >> 2) | ((long)call << (2 * (KSIZE - 1)));
        }
        public long getSuccessorVal( final int call ) { return ((kVal << 2) & KMASK) | call; }

        public static boolean isCanonical( final long val ) {
            return (val & (1L << KSIZE)) == 0L;
        }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof Kmer && kVal == ((Kmer)obj).kVal;
        }

        @Override public int hashCode() {
            return (int)(kVal ^ kVal >>> 32);
        }
    }

    public static final class KmerSet<KMER extends Kmer> extends HopscotchSet<KMER> {
        public KmerSet( final int capacity ) { super(capacity); }

        @Override
        protected int hashToIndex( final Object kmer ) {
            return (int)(((241 * ((Kmer)kmer).getKVal()) & Long.MAX_VALUE) % capacity());
        }
    }

    public static abstract class KmerAdjacency extends Kmer {
        public KmerAdjacency( final long kVal ) { super(kVal); }

        public abstract KmerAdjacency getSolePredecessor();
        public abstract int getPredecessorMask();
        public abstract int getPredecessorCount();
        public abstract void removePredecessor( final int callToRemove,
                                                final KmerSet<KmerAdjacency> kmerAdjacencySet );

        public abstract KmerAdjacency getSoleSuccessor();
        public abstract int getSuccessorMask();
        public abstract int getSuccessorCount();
        public abstract void removeSuccessor( final int callToRemove,
                                              final KmerSet<KmerAdjacency> kmerAdjacencySet );

        public abstract Contig getContig();
        public abstract int getContigOffset();
        public abstract void setContig( final Contig contig, final int contigOffset );

        public abstract int getNObservations();
        public abstract KmerAdjacency rc();
        public abstract KmerAdjacencyImpl canonical();

        public void observe( final KmerAdjacency predecessor, final KmerAdjacency successor ) {
            observe(predecessor, successor, 1);
        }

        public abstract void observe( final KmerAdjacency predecessor,
                                      final KmerAdjacency successor,
                                      final int count );

        @Override public String toString() {
            final StringBuilder sb = new StringBuilder(KSIZE);
            long currentVal = getKVal();
            for ( int idx = 0; idx != KSIZE; ++idx ) {
                sb.append("ACGT".charAt((int)currentVal & 3));
                currentVal >>= 2;
            }
            sb.reverse();
            return sb.toString();
        }

        public static void kmerize( final byte[] calls,
                                    final byte[] quals,
                                    final byte qMin,
                                    final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            int currentCount = 0;
            long currentKVal = 0;
            KmerAdjacency prevAdjacency = null;
            KmerAdjacency currentAdjacency = null;
            for ( int idx = 0; idx < calls.length; ++idx ) {
                if ( quals[idx] <  qMin ) {
                    if ( currentAdjacency != null ) {
                        currentAdjacency.observe(prevAdjacency, null);
                    }
                    currentCount = 0;
                    currentAdjacency = prevAdjacency = null;
                    continue;
                }
                currentKVal <<= 2;
                switch ( calls[idx] ) {
                    case 'A': case 'a': break;
                    case 'C': case 'c': currentKVal += 1; break;
                    case 'G': case 'g': currentKVal += 2; break;
                    case 'T': case 't': currentKVal += 3; break;
                    default:
                        if ( currentAdjacency != null ) {
                            currentAdjacency.observe(prevAdjacency, null);
                        }
                        currentCount = 0;
                        currentAdjacency = prevAdjacency = null;
                        continue;
                }
                if ( ++currentCount >= KSIZE ) {
                    final KmerAdjacency nextAdjacency = findOrAdd(currentKVal, kmerAdjacencySet);
                    if ( currentAdjacency != null ) {
                        currentAdjacency.observe(prevAdjacency, nextAdjacency);
                    }
                    prevAdjacency = currentAdjacency;
                    currentAdjacency = nextAdjacency;
                }
            }
            if ( currentAdjacency != null ) {
                currentAdjacency.observe(prevAdjacency, null);
            }
        }

        public static void kmerize( final String sequence,
                                    final int nObservations,
                                    final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            int currentCount = 0;
            long currentKVal = 0;
            int nObs = 0;
            KmerAdjacency prevAdjacency = null;
            KmerAdjacency currentAdjacency = null;
            final int nCalls = sequence.length();
            for ( int idx = 0; idx != nCalls; ++idx ) {
                currentKVal <<= 2;
                switch ( sequence.charAt(idx) ) {
                    case 'A': case 'a': break;
                    case 'C': case 'c': currentKVal += 1; break;
                    case 'G': case 'g': currentKVal += 2; break;
                    case 'T': case 't': currentKVal += 3; break;
                }
                if ( ++currentCount >= KSIZE ) {
                    final KmerAdjacency nextAdjacency = findOrAdd(currentKVal, kmerAdjacencySet);
                    if ( currentAdjacency != null ) {
                        currentAdjacency.observe(prevAdjacency, nextAdjacency, nObs);
                        nObs = nObservations;
                    }
                    prevAdjacency = currentAdjacency;
                    currentAdjacency = nextAdjacency;
                }
            }
            if ( currentAdjacency != null ) {
                currentAdjacency.observe(prevAdjacency, null, 0);
            }
        }

        // Lookup table for reverse-complementing each possible byte value.
        // Each pair of bits represents a base, so you have to reverse bits pairwise and then invert all bits.
        // This is most quickly and easily done with a lookup table.
        private static final long[] BYTEWISE_REVERSE_COMPLEMENT;
        static {
            BYTEWISE_REVERSE_COMPLEMENT = new long[256];
            for ( int bIn = 0; bIn != 256; ++bIn ) {
                BYTEWISE_REVERSE_COMPLEMENT[bIn] =
                        ~(((bIn & 3) << 6) | (((bIn >> 2) & 3) << 4) |
                                (((bIn >> 4) & 3) << 2) | ((bIn >> 6) & 3)) & 0xffL;
            }
        }

        public static long reverseComplement( long val ) {
            // process val one byte at a time
            long result = BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF]; // handle the low-order byte
            int nBytes = 8;
            while ( --nBytes != 0 ) { // pre-decrementing:  we'll go through the loop 7 times
                // rotate down by a byte
                val >>= 8;
                // rotate up by a byte and OR in the reverse complement of the next byte
                result = (result << 8) | BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF];
            }
            return result >>> (Long.SIZE - 2*KSIZE);
        }

        public static KmerAdjacency find( final long kVal,
                                          final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            if ( isCanonical(kVal) ) return kmerAdjacencySet.find(new Kmer(kVal & KMASK));
            final KmerAdjacency result = kmerAdjacencySet.find(new Kmer(reverseComplement(kVal)));
            return result == null ? null : result.rc();
        }

        public static KmerAdjacency findOrAdd( final long kVal,
                                               final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            if ( isCanonical(kVal) ) {
                return kmerAdjacencySet.findOrAdd(new Kmer(kVal & KMASK), kmer ->
                        new KmerAdjacencyImpl(((Kmer)kmer).getKVal()));
            }
            return kmerAdjacencySet.findOrAdd(new Kmer(reverseComplement(kVal)), kmer ->
                    new KmerAdjacencyImpl(((Kmer)kmer).getKVal())).rc();
        }
    }

    public static final class KmerAdjacencyRC extends KmerAdjacency {
        private final KmerAdjacencyImpl rc;
        private static final int[] NIBREV =
        // 0000,  0001,  0010,  0011,  0100,  0101,  0110,  0111,  1000,  1001,  1010,  1011,  1100,  1101,  1110,  1111
        {0b0000,0b1000,0b0100,0b1100,0b0010,0b1010,0b0110,0b1110,0b0001,0b1001,0b0101,0b1101,0b0011,0b1011,0b0111,0b1111};

        public KmerAdjacencyRC( final KmerAdjacencyImpl rc ) {
            super(reverseComplement(rc.getKVal()));
            this.rc = rc;
        }

        @Override public KmerAdjacency getSolePredecessor() {
            final KmerAdjacency successor = rc.getSoleSuccessor();
            return successor == null ? null : successor.rc();
        }
        @Override public int getPredecessorMask() { return NIBREV[rc.getSuccessorMask()]; }
        @Override public int getPredecessorCount() { return rc.getSuccessorCount(); }
        @Override
        public void removePredecessor( final int callToRemove,
                                       final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            rc.removeSuccessor(3 - callToRemove, kmerAdjacencySet);
        }

        @Override public KmerAdjacency getSoleSuccessor() {
            final KmerAdjacency predecessor = rc.getSolePredecessor();
            return predecessor == null ? null : predecessor.rc();
        }
        @Override public int getSuccessorMask() { return NIBREV[rc.getPredecessorMask()]; }
        @Override public int getSuccessorCount() { return rc.getPredecessorCount(); }
        @Override
        public void removeSuccessor( final int callToRemove,
                                     final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            rc.removePredecessor(3 - callToRemove, kmerAdjacencySet);
        }

        @Override public Contig getContig() {
            final Contig contig = rc.getContig();
            return contig == null ? null : contig.rc();
        }
        @Override public int getContigOffset() {
            final Contig contig = rc.getContig();
            return contig == null ? 0 : contig.size() - rc.getContigOffset() - KSIZE;
        }
        @Override public void setContig( final Contig contig, final int contigOffset ) {
            if ( contig == null ) rc.setContig(null, 0);
            else rc.setContig(contig.rc(), contig.size() - contigOffset - KSIZE);
        }

        @Override public int getNObservations() { return rc.getNObservations(); }
        @Override public KmerAdjacency rc() { return rc; }
        @Override public KmerAdjacencyImpl canonical() { return rc; }

        @Override public void observe( final KmerAdjacency predecessor,
                                       final KmerAdjacency successor,
                                       final int count ) {
            rc.observe(successor == null ? null : successor.rc(),
                    predecessor == null ? null : predecessor.rc(),
                    count);
        }
    }

    public static final class KmerAdjacencyImpl extends KmerAdjacency {
        private KmerAdjacency solePredecessor; // set to null if there are no predecessors, or multiple predecessors
        private KmerAdjacency soleSuccessor; // set to null if there are no successors, or multiple successors
        private int predecessorMask; // bit mask of observed kmers preceding this one
        private int successorMask; // bit mask observed kmers following this one
        private Contig contig; // the contig that contains this Kmer
        private int contigOffset; // the offset within the contig where this kmer is found
        private int nObservations; // the reads in which this kmer was observed
        private final KmerAdjacencyRC rc; // the reverse-complement of this kmer
        private static final int[] COUNT_FOR_MASK =
                //side sum for binary values from 0 -> 15
                //0000  0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
                {    0,    1,   1,   2,   1,   2,   2,   3,   1,   2,   2,   3,   2,   3,   3,   4 };

        public KmerAdjacencyImpl( final long kVal ) {
            super(kVal);
            this.rc = new KmerAdjacencyRC(this);
        }

        @Override public KmerAdjacency getSolePredecessor() { return solePredecessor; } // may return null
        @Override public int getPredecessorMask() { return predecessorMask; }
        @Override public int getPredecessorCount() { return COUNT_FOR_MASK[predecessorMask]; }
        @Override
        public void removePredecessor( final int callToRemove,
                                       final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            predecessorMask &= ~(1 << callToRemove);
            solePredecessor = null;
            if ( getPredecessorCount() == 1 ) {
                for ( int call = 0; call != 4; ++call ) {
                    if ( ((1 << call) & predecessorMask) != 0 ) {
                        solePredecessor = find(getPredecessorVal(call), kmerAdjacencySet);
                        break;
                    }
                }
            }
        }

        @Override public KmerAdjacency getSoleSuccessor() { return soleSuccessor; } // may return null
        @Override public int getSuccessorMask() { return successorMask; }
        @Override public int getSuccessorCount() { return COUNT_FOR_MASK[successorMask]; }
        @Override
        public void removeSuccessor( final int callToRemove,
                                     final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            successorMask &= ~(1 << callToRemove);
            soleSuccessor = null;
            if ( getSuccessorCount() == 1 ) {
                for ( int call = 0; call != 4; ++call ) {
                    if ( ((1 << call) & successorMask) != 0 ) {
                        soleSuccessor = find(getSuccessorVal(call), kmerAdjacencySet);
                        break;
                    }
                }
            }
        }

        @Override public Contig getContig() { return contig; }
        @Override public int getContigOffset() { return contigOffset; }
        @Override public void setContig( final Contig contig, final int contigOffset ) {
            if ( contig != null && this.contig != null ) {
                throw new GATKException("Internal error: overwriting kmer contig.");
            }
            this.contig = contig;
            this.contigOffset = contigOffset;
        }

        @Override public int getNObservations() { return nObservations; }
        @Override public KmerAdjacency rc() { return rc; }
        @Override public KmerAdjacencyImpl canonical() { return this; }

        @Override public void observe( final KmerAdjacency predecessor,
                                       final KmerAdjacency successor,
                                       final int count ) {
            if ( predecessor != null ) {
                if ( predecessor.getSuccessorVal(getFinalCall()) != getKVal() ) {
                    throw new GATKException("illegal predecessor");
                }
                final int initialCall = predecessor.getInitialCall();
                final int newPredecessorMask = 1 << initialCall;
                if ( (newPredecessorMask & predecessorMask) == 0 ) {
                    if ( predecessorMask == 0 ) {
                        solePredecessor = predecessor;
                        predecessorMask = newPredecessorMask;
                    } else {
                        solePredecessor = null;
                        predecessorMask |= newPredecessorMask;
                    }
                }
            }
            if ( successor != null ) {
                if ( successor.getPredecessorVal(getInitialCall()) != getKVal() ) {
                    throw new GATKException("illegal successor");
                }
                final int finalCall = successor.getFinalCall();
                final int newSuccessorMask = 1 << finalCall;
                if ( (newSuccessorMask & successorMask) == 0 ) {
                    if ( successorMask == 0 ) {
                        soleSuccessor = successor;
                        successorMask = newSuccessorMask;
                    } else {
                        soleSuccessor = null;
                        successorMask |= newSuccessorMask;
                    }
                }
            }
            nObservations += count;
        }
    }

    public enum ContigOrientation {
        FWD, // k-mer appears at the 5' end of the contig
        REV, // k-mer appears at the 5' end of the reverse-complemented contig
        BOTH // k-mer occurs on 5' end of the contig and its RC (can happen when the contig is a palindrome)
    }

    public static final class ContigEndKmer extends Kmer {
        private final Contig contig;
        private final ContigOrientation contigOrientation;

        public ContigEndKmer( final long kVal,
                              final Contig contig,
                              final ContigOrientation contigEnd ) {
            super(kVal);
            this.contig = contig;
            this.contigOrientation = contigEnd;
        }

        public Contig getContig() { return contig; }
        public ContigOrientation getContigOrientation() { return contigOrientation; }
    }

    public interface Contig {
        int getId();
        CharSequence getSequence();
        int getMaxObservations();
        KmerAdjacency getFirstKmer();
        KmerAdjacency getLastKmer();
        List<Contig> getPredecessors();
        List<Contig> getSuccessors();
        int getComponentId();
        int size();
        default int getNKmers() { return size() - Kmer.KSIZE + 1; }
        Contig rc();
        boolean isCyclic();
        void setCyclic( final boolean cyclic );
        boolean isCut();
        void setCut( final boolean cut );
        boolean isCanonical();
        ContigImpl canonical();
    }

    public static final class ContigImpl implements Contig {
        private static int nContigs;
        private final int id;
        private final CharSequence sequence;
        private final int maxObservations;
        private final KmerAdjacency firstKmer;
        private final KmerAdjacency lastKmer;
        private final List<Contig> predecessors;
        private final List<Contig> successors;
        private int componentId;
        private boolean cyclic;
        private boolean cut;
        private final Contig rc;

        public ContigImpl( final KmerAdjacency firstKmerAdjacency ) {
            this.id = nContigs++;
            final StringBuilder sb = new StringBuilder(firstKmerAdjacency.toString());
            int maxObservations = firstKmerAdjacency.getNObservations();
            KmerAdjacency lastKmerAdjacency = firstKmerAdjacency;
            for ( KmerAdjacency kmerAdjacency = firstKmerAdjacency.getSoleSuccessor();
                  kmerAdjacency != null;
                  kmerAdjacency = kmerAdjacency.getSoleSuccessor() ) {
                // if we've gone around a circle, or if we're branching backwards, or if we hit a palindrome u-turn
                if ( firstKmerAdjacency == kmerAdjacency ||
                        kmerAdjacency.getPredecessorCount() != 1 ||
                        kmerAdjacency == lastKmerAdjacency.rc() ) {
                    break;
                }
                sb.append("ACGT".charAt(kmerAdjacency.getFinalCall()));
                maxObservations = Math.max(maxObservations, kmerAdjacency.getNObservations());
                lastKmerAdjacency = kmerAdjacency;
            }
            this.sequence = sb.toString();
            this.maxObservations = maxObservations;
            this.firstKmer = firstKmerAdjacency;
            this.lastKmer = lastKmerAdjacency;
            this.predecessors = new ArrayList<>(firstKmer.getPredecessorCount());
            this.successors = new ArrayList<>(lastKmer.getSuccessorCount());
            this.rc = new ContigRCImpl(this);

            int offset = 0;
            for ( KmerAdjacency kmerAdjacency = firstKmerAdjacency;
                  kmerAdjacency != lastKmerAdjacency;
                  kmerAdjacency = kmerAdjacency.getSoleSuccessor() ) {
                kmerAdjacency.setContig(this, offset++);
            }
            lastKmerAdjacency.setContig(this, offset);
        }

        // create a new contig by joining two contigs
        public ContigImpl( final Contig predecessor, final Contig successor ) {
            if ( predecessor == successor || predecessor == successor.rc() ) {
                throw new GATKException("can't self-join");
            }
            this.id = nContigs++;
            final StringBuilder sb = new StringBuilder(predecessor.getSequence());
            final CharSequence successorSequence = successor.getSequence();
            sb.append(successorSequence.subSequence(Kmer.KSIZE - 1, successorSequence.length()));
            this.sequence = sb.toString();
            this.maxObservations =
                    Math.max(predecessor.getMaxObservations(), successor.getMaxObservations());
            this.firstKmer = predecessor.getFirstKmer();
            this.lastKmer = successor.getLastKmer();
            this.predecessors = new ArrayList<>(predecessor.getPredecessors().size());
            this.successors = new ArrayList<>(successor.getSuccessors().size());
            this.rc = new ContigRCImpl(this);

            // fix predecessor linkages to point to new contig
            for ( final Contig predPredecessor : predecessor.getPredecessors() ) {
                if ( predPredecessor == successor ) {
                    predecessors.add(this);
                } else if ( predPredecessor == predecessor.rc() ) {
                    predecessors.add(rc);
                } else {
                    predecessors.add(predPredecessor);
                    final List<Contig> successors = predPredecessor.getSuccessors();
                    successors.set(successors.indexOf(predecessor), this);
                }
            }

            // fix successor linkages to point to new contig
            for ( final Contig succSuccessor : successor.getSuccessors() ) {
                if ( succSuccessor == predecessor ) {
                    successors.add(this);
                } else if ( succSuccessor == successor.rc() ) {
                    successors.add(rc);
                } else {
                    successors.add(succSuccessor);
                    final List<Contig> predecessors = succSuccessor.getPredecessors();
                    predecessors.set(predecessors.indexOf(successor), this);
                }
            }
        }

        @Override public int getId() { return id; }
        @Override public CharSequence getSequence() { return sequence; }
        @Override public int getMaxObservations() { return maxObservations; }
        @Override public KmerAdjacency getFirstKmer() { return firstKmer; }
        @Override public KmerAdjacency getLastKmer() { return lastKmer; }
        @Override public List<Contig> getPredecessors() { return predecessors; }
        @Override public List<Contig> getSuccessors() { return successors; }
        @Override public int getComponentId() { return componentId; }
        public void setComponentId( final int id ) { this.componentId = id; }
        @Override public int size() { return sequence.length(); }
        @Override public Contig rc() { return rc; }
        @Override public boolean isCyclic() { return cyclic; }
        @Override public void setCyclic( final boolean cyclic ) { this.cyclic = cyclic; }
        @Override public boolean isCut() { return cut; }
        @Override public void setCut( final boolean cut ) { this.cut = cut; }
        @Override public boolean isCanonical() { return true; }
        @Override public ContigImpl canonical() { return this; }
        @Override public String toString() { return "c" + id; }
    }

    public static final class ContigRCImpl implements Contig {
        private final CharSequence sequence;
        private final List<Contig> predecessors;
        private final List<Contig> successors;
        private final ContigImpl rc;

        public ContigRCImpl( final ContigImpl contig ) {
            this.sequence = new SequenceRC(contig.getSequence());
            this.predecessors = new ListRC(contig.getSuccessors());
            this.successors = new ListRC(contig.getPredecessors());
            this.rc = contig;
        }

        @Override public int getId() { return ~rc.getId(); }
        @Override public CharSequence getSequence() { return sequence; }
        @Override public int getMaxObservations() { return rc.getMaxObservations(); }
        @Override public KmerAdjacency getFirstKmer() { return rc.getLastKmer().rc(); }
        @Override public KmerAdjacency getLastKmer() { return rc.getFirstKmer().rc(); }
        @Override public List<Contig> getPredecessors() { return predecessors; }
        @Override public List<Contig> getSuccessors() { return successors; }
        @Override public int getComponentId() { return rc.getComponentId(); }
        @Override public int size() { return sequence.length(); }
        @Override public Contig rc() { return rc; }
        @Override public boolean isCyclic() { return rc.isCyclic(); }
        @Override public void setCyclic( final boolean cyclic ) { rc.setCyclic(cyclic); }
        @Override public boolean isCut() { return rc.isCut(); }
        @Override public void setCut( final boolean cut ) { rc.setCut(cut); }
        @Override public boolean isCanonical() { return false; }
        @Override public ContigImpl canonical() { return rc; }
        @Override public String toString() { return rc.toString() + "RC"; }
    }

    public static final class SequenceRC implements CharSequence, Comparable<CharSequence> {
        private final int lenLess1;
        private final CharSequence sequence;

        public SequenceRC( final CharSequence sequence ) {
            this.lenLess1 = sequence.length() - 1;
            this.sequence = sequence;
        }

        @Override public int length() { return sequence.length(); }
        @Override public char charAt( final int index ) {
            final char result;
            switch ( Character.toUpperCase(sequence.charAt(lenLess1 - index)) ) {
                case 'A': result = 'T'; break;
                case 'C': result = 'G'; break;
                case 'G': result = 'C'; break;
                case 'T': result = 'A'; break;
                default: result = 'N'; break;
            }
            return result;
        }
        @Override public CharSequence subSequence( final int start, final int end ) {
            return new StringBuilder(end - start).append(this, start, end).toString();
        }
        @Override public String toString() { return new StringBuilder(this).toString(); }

        @Override public int compareTo( final CharSequence charSequence ) {
            final int len1 = length();
            final int len2 = charSequence.length();
            final int cmpLen = Math.min(len1, len2);
            for ( int idx = 0; idx != cmpLen; ++idx ) {
                final char char1 = charAt(idx);
                final char char2 = Character.toUpperCase(charSequence.charAt(idx));
                if ( char1 > char2 ) return 1;
                if ( char1 < char2 ) return -1;
            }
            return Integer.compare(len1, len2);
        }
    }

    public static final class ListRC extends AbstractList<Contig> {
        private final List<Contig> contigList;

        public ListRC( final List<Contig> contigList ) {
            this.contigList = contigList;
        }

        @Override public Contig get( final int index ) {
            return contigList.get(reflectIndex(index)).rc();
        }
        @Override public int size() { return contigList.size(); }
        @Override public Contig set( final int index, final Contig contig ) {
            return contigList.set(reflectIndex(index), contig.rc()).rc();
        }
        @Override public void add( final int index, final Contig contig ) {
            contigList.add(reflectIndex(index), contig.rc());
        }
        @Override public Contig remove( final int index ) {
            return contigList.remove(reflectIndex(index)).rc();
        }

        private int reflectIndex( final int index ) { return size() - 1 - index; }
    }

    public interface PathPart {
        Contig getContig();
        CharSequence getSequence();
        void extend( final char call );
        int getStart();
        int getStop();
        boolean isGap();
        int getLength();
        PathPart rc();
        char getFirstCall();
        char getLastCall();
        default boolean startsAtBeginning() { return getStart() == 0; }
        default boolean stopsAtEnd() { return getStop() + Kmer.KSIZE - 1 == getContig().size(); }
    }

    public static final class PathPartGap implements PathPart {
        private final StringBuilder sequence = new StringBuilder();

        public PathPartGap( final KmerAdjacency kmer ) { sequence.append(kmer.toString()); }
        private PathPartGap( final CharSequence sequence ) { this.sequence.append(sequence); }

        @Override public Contig getContig() { return null; }
        @Override public CharSequence getSequence() { return sequence.toString(); }
        @Override public void extend( final char call ) { sequence.append(call); }
        @Override public int getStart() { return 0; }
        @Override public int getStop() { return sequence.length(); }
        @Override public boolean isGap() { return true; }
        @Override public int getLength() { return sequence.length() - Kmer.KSIZE + 1; }
        @Override public PathPart rc() { return new PathPartGap(new SequenceRC(sequence)); }
        @Override public char getFirstCall() { return sequence.charAt(Kmer.KSIZE - 1); }
        @Override public char getLastCall() {
            return sequence.charAt(sequence.length() - Kmer.KSIZE + 1);
        }
    }

    public static final class PathPartContig implements PathPart {
        private final Contig contig;
        private final int start;
        private int stop;

        public PathPartContig( final Contig contig, final int start ) {
            this(contig, start, start+1);
        }
        public PathPartContig( final Contig contig, final int start, final int stop ) {
            this.contig = contig;
            this.start = start;
            this.stop = stop;
        }

        @Override public Contig getContig() { return contig; }
        @Override public String getSequence() { return null; }
        @Override public void extend( final char call ) { stop += 1; }
        @Override public int getStart() { return start; }
        @Override public int getStop() { return stop; }
        @Override public boolean isGap() { return false; }
        @Override public int getLength() { return stop - start; }
        @Override public PathPart rc() {
            final int revBase = contig.size() - Kmer.KSIZE + 1;
            return new PathPartContig(contig.rc(), revBase - stop, revBase - start);
        }
        @Override public char getFirstCall() {
            return getContig().getSequence().charAt(start + Kmer.KSIZE - 1);
        }
        @Override public char getLastCall() { return getContig().getSequence().charAt(stop - 1); }
    }

    public static final class Error {
        private final ContigImpl contig;
        private final int offset;
        private final char call;
        private final byte quality;

        public Error( final Contig contig, final int offset, final char call, final byte quality ) {
            this.contig = contig.canonical();
            this.offset = this.contig == contig ? offset : contig.size() - offset - 1;
            this.call = call;
            this.quality = quality;
        }

        public Contig getContig() { return contig; }
        public int getOffset() { return offset; }
        public char getCall() { return call; }
        public byte getQuality() { return quality; }
    }

    public static final class Path {
        private final List<PathPart> parts;

        public Path( final byte[] calls,
                     final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            parts = new ArrayList<>();
            long kVal = 0;
            int count = 0;
            PathPart currentPathPart = null;
            for ( int idx = 0; idx != calls.length; ++idx ) {
                final char call = (char)calls[idx];
                kVal <<= 2;
                switch ( call ) {
                    case 'C': case 'c': kVal += 1; break;
                    case 'G': case 'g': kVal += 2; break;
                    case 'T': case 't': kVal += 3; break;
                }
                if ( ++count >= Kmer.KSIZE ) {
                    final KmerAdjacency kmer = KmerAdjacencyImpl.find(kVal, kmerAdjacencySet);
                    // if we fail to look up the kmer
                    if ( kmer == null ) {
                        if ( currentPathPart == null ) {
                            // if there's no current path part, just create the 1st one as a PathPartGap
                            currentPathPart = new PathPartGap(new KmerAdjacencyImpl(kVal));
                            parts.add(currentPathPart);
                        } else if ( currentPathPart.isGap() ) {
                            // if the current path part is a PathPartGap, just extend it
                            currentPathPart.extend(call);
                        } else {
                            // new PathPartGap
                            currentPathPart = new PathPartGap(new KmerAdjacencyImpl(kVal));
                            parts.add(currentPathPart);
                        }
                    } else {
                        // we've found our kmer
                        final Contig contig = kmer.getContig();
                        if ( currentPathPart == null ) {
                            // we've looked up a kmer, but don't have a current path part -- create one
                            currentPathPart = new PathPartContig(contig, kmer.getContigOffset());
                            parts.add(currentPathPart);
                        } else if ( contig == currentPathPart.getContig() ) {
                            // our lookup is on the current path part's contig -- extend it
                            final int kmerOffset = kmer.getContigOffset();
                            final int curStop = currentPathPart.getStop();
                            if ( kmerOffset == curStop ) {
                                currentPathPart.extend(call);
                            } else if ( kmerOffset == 0 && contig.getNKmers() == curStop ) {
                                // cycle onto same contig
                                currentPathPart = new PathPartContig(contig, 0);
                                parts.add(currentPathPart);
                            } else {
                                // weird:  kmer is non-contiguous.  start a new path part after a zero-length gap
                                parts.add(zeroLengthGap(currentPathPart));
                                currentPathPart = new PathPartContig(contig, kmerOffset);
                                parts.add(currentPathPart);
                            }
                        } else {
                            final int kmerContigOffset = kmer.getContigOffset();
                            if ( currentPathPart.isGap() ) {
                                // squash captured gaps caused by single-base sequencing errors
                                final int gapLen = currentPathPart.getLength();
                                if ( gapLen == Kmer.KSIZE ) {
                                    final int prevPartIdx = parts.size() - 2;
                                    if ( prevPartIdx >= 0 ) {
                                        final PathPart prevPart = parts.get(prevPartIdx);
                                        final Contig prevPartContig = prevPart.getContig();
                                        final int prevPartStart = prevPart.getStart();
                                        final int prevPartStop = prevPart.getStop();
                                        final int prevPartMaxStop =
                                                prevPartContig.size() - Kmer.KSIZE + 1;
                                        final int newStop = kmerContigOffset + 1;
                                        if ( prevPartContig == kmer.getContig() ) {
                                            if ( kmerContigOffset - prevPartStop == gapLen ) {
                                                currentPathPart =
                                                        new PathPartContig(prevPartContig, prevPartStart, newStop);
                                                parts.set(prevPartIdx, currentPathPart);
                                                parts.remove(prevPartIdx + 1);
                                                continue;
                                            }
                                        } else if ( prevPartMaxStop - prevPartStop + kmerContigOffset == gapLen ) {
                                            parts.set(prevPartIdx,
                                                    new PathPartContig(prevPartContig, prevPartStart, prevPartMaxStop));
                                            currentPathPart = new PathPartContig(kmer.getContig(), 0, newStop);
                                            parts.set(prevPartIdx + 1, currentPathPart);
                                            continue;
                                        }
                                    }
                                }
                            } else if ( !currentPathPart.stopsAtEnd() || kmerContigOffset != 0 ) {
                                // not an end-to-end join across contigs -- record a zero-length gap
                                parts.add(zeroLengthGap(currentPathPart));
                            }
                            // we're jumping to a new contig.  start a new path part
                            currentPathPart = new PathPartContig(contig, kmerContigOffset);
                            parts.add(currentPathPart);
                        }
                    }
                }
            }
        }

        private static PathPart zeroLengthGap( final PathPart currentPathPart ) {
            final int currentStop = currentPathPart.getStop();
            final CharSequence currentSequence = currentPathPart.getContig().getSequence();
            final CharSequence almostAKmer =
                    currentSequence.subSequence(currentStop, currentStop + Kmer.KSIZE - 1);
            return new PathPartGap(almostAKmer);
        }

        // RCing constructor
        private Path( final Path that ) {
            this.parts = new ArrayList<>();
            final List<PathPart> thoseParts = that.parts;
            for ( int idx = thoseParts.size() - 1; idx >= 0; --idx ) {
                parts.add(thoseParts.get(idx).rc());
            }
        }

        public List<PathPart> getParts() { return parts; }
        public Path rc() { return new Path(this); }

        @Override public String toString() {
            if ( parts.size() == 0 ) return "";
            final StringBuilder sb = new StringBuilder();
            String prefix = "";
            final PathPart firstPart = parts.get(0);
            final PathPart lastPart = parts.get(parts.size() - 1);
            for ( final PathPart pp : parts ) {
                sb.append(prefix);
                prefix = ", ";
                if ( pp.isGap() ) {
                    sb.append("NoKmer(").append(pp.getLength()).append(")");
                } else {
                    final Contig contig = pp.getContig();
                    sb.append(contig);
                    final int maxStop = contig.size() - Kmer.KSIZE + 1;
                    if ( (pp != firstPart && pp.getStart() != 0) ||
                         (pp != lastPart && pp.getStop() != maxStop) ) {
                        sb.append('(').append(pp.getStart()).append('-')
                                .append(pp.getStop()).append('/').append(maxStop).append(')');
                    }
                }
            }
            return sb.toString();
        }
    }

    public static final class CutData {
        public static int nextNum;
        public int visitNum;
        public int minVisitNum;

        public CutData() {
            this.visitNum = ++nextNum;
            this.minVisitNum = visitNum;
        }
    }

    public static final class TransitPairCount {
        private final Contig prevContig;
        private final Contig nextContig;
        private final TransitPairCount rc;
        private int count;

        public TransitPairCount( final Contig prevContig, final Contig nextContig ) {
            this.prevContig = prevContig;
            this.nextContig = nextContig;
            this.rc = new TransitPairCount(this);
            this.count = 1;
        }

        private TransitPairCount( final TransitPairCount rc ) {
            this.prevContig = rc.nextContig.rc();
            this.nextContig = rc.prevContig.rc();
            this.rc = rc;
            this.count = 1;
        }

        public Contig getPrevContig() { return prevContig; }
        public Contig getNextContig() { return nextContig; }
        public TransitPairCount getRC() { return rc; }
        public void observe() { count += 1; rc.count += 1; }
        public void resetCount() { count = 0; rc.count = 0; }
        public int getCount() { return count; }

        @Override public boolean equals( final Object obj ) {
            if ( !(obj instanceof TransitPairCount) ) return false;
            final TransitPairCount that = (TransitPairCount)obj;
            return this.prevContig == that.prevContig && this.nextContig == that.nextContig;
        }

        @Override public int hashCode() {
            return 47 * (47 * prevContig.hashCode() + nextContig.hashCode());
        }

        @Override public String toString() {
            return prevContig + "<-->" + nextContig + " " + count + "x";
        }
    }

    public static final class Traversal {
        private final List<Contig> contigs;
        private final boolean isInextensible;
        private final int minMaxObservations;

        public Traversal( final Collection<Contig> contigs ) {
            this(contigs, false);
        }

        public Traversal( final Collection<Contig> contigs, final boolean isInextensible ) {
            if ( contigs == null || contigs.isEmpty() ) {
                throw new GATKException("null or empty list of contigs in traversal");
            }
            this.contigs = new ArrayList<>(contigs);
            this.isInextensible = isInextensible;
            int minMaxObservations = Integer.MAX_VALUE;
            for ( Contig contig : contigs ) {
                minMaxObservations = Math.min(minMaxObservations, contig.getMaxObservations());
            }
            this.minMaxObservations = minMaxObservations;
        }

        // RC constructor
        private Traversal( final Traversal traversal ) {
            this.contigs = new ListRC(traversal.contigs);
            this.isInextensible = false;
            this.minMaxObservations = traversal.minMaxObservations;
        }

        public List<Contig> getContigs() { return Collections.unmodifiableList(contigs); }
        public Contig getFirstContig() { return contigs.get(0); }
        public Contig getLastContig() { return contigs.get(contigs.size() - 1); }
        public Traversal rc() { return new Traversal(this); }
        public boolean isRC() { return contigs instanceof ListRC; }
        public boolean isInextensible() { return isInextensible; }
        public int getMinMaxObservations() { return minMaxObservations; }

        public String getName() {
            final StringBuilder sb = new StringBuilder();
            String prefix = "";
            for ( final Contig contig : contigs ) {
                sb.append(prefix).append(contig.toString());
                prefix = "+";
            }
            return sb.toString();
        }

        public int getSequenceLength() {
            int len = 0;
            for ( final Contig contig : contigs ) {
                len += contig.getNKmers();
            }
            return len + Kmer.KSIZE - 1;
        }

        public String getSequence() {
            if ( contigs.size() == 0 ) return "";
            final StringBuilder sb =
                    new StringBuilder(contigs.get(0).getSequence().subSequence(0, Kmer.KSIZE - 1));
            for ( final Contig contig : contigs ) {
                final CharSequence seq = contig.getSequence();
                sb.append(seq.subSequence(Kmer.KSIZE - 1, seq.length()));
            }
            return sb.toString();
        }

        @Override public int hashCode() { return contigs.hashCode(); }
        @Override public boolean equals( final Object obj ) {
            if ( this == obj ) return true;
            if ( !(obj instanceof Traversal) ) return false;
            final Traversal that = (Traversal)obj;
            return contigs.equals(that.contigs);
        }

        public static Traversal combine( final Traversal trav1, final Traversal trav2 ) {
            return combineOverlappers(trav1, trav2, 1);
        }

        public static Traversal combineOverlappers( final Traversal trav1,
                                                    final Traversal trav2,
                                                    final int overlapLen ) {
            final int len1 = trav1.contigs.size();
            if ( !trav1.contigs.subList(len1 - overlapLen, len1)
                    .equals(trav2.contigs.subList(0, overlapLen)) ) {
                throw new GATKException("combining non-overlapping traversals");
            }
            final int len2 = trav2.contigs.size();
            final List<Contig> combinedList =
                    new ArrayList<>(trav1.contigs.size() + len2 - overlapLen);
            combinedList.addAll(trav1.contigs);
            combinedList.addAll(trav2.contigs.subList(overlapLen, len2));
            return new Traversal(combinedList);
        }
    }

    public final static class AssemblyTooComplexException extends RuntimeException {
        static final long serialVersionUID = -1L;
    }
}
