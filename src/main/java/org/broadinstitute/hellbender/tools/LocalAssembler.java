package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.PairWalker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.HopscotchSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Performs local assembly of small regions to discover structural variants.",
        oneLineSummary = "Local assembler for SVs",
        usageExample = "gatk LocalAssembler -L chr21:16187360-16187360 --ip 500 -R 38.fa.gz " +
                "-I NA19240.cram -I NA19240.distantmate.bam " +
                "--assembly-name chr21_16187360_16187360_INS --gfa-file test.gfa --fasta-file test.fa.gz",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class LocalAssembler extends PairWalker {
    @Argument(fullName="assembly-name", doc="Name of assembly used as a prefix for traversal names.")
    public String assemblyName;

    @Argument(fullName="gfa-file", doc="Path to assembly output in gfa format.", optional=true)
    public GATKPath gfaFile;

    @Argument(fullName="fasta-file", doc="Path to scaffolds in fasta format.", optional=true)
    public GATKPath fastaFile;

    public static final byte QMIN_DEFAULT = 25;
    @Argument(fullName="q-min", doc="Minimum base quality when kmerizing reads.", optional=true)
    private byte qMin = QMIN_DEFAULT;

    public static final int MIN_THIN_OBS_DEFAULT = 4;
    @Argument(fullName="min-thin-observations",
            doc="Minimum number of observations of some kmer within the contig required to " +
                    "retain the contig.", optional=true)
    private int minThinObs = MIN_THIN_OBS_DEFAULT;

    public static final int MIN_GAPFILL_COUNT_DEFAULT = 3;
    @Argument(fullName="min-gapfill-count",
            doc="Minimum number of observations of a sequence that patches a gap.", optional=true)
    private int minGapfillCount = MIN_GAPFILL_COUNT_DEFAULT;

    public static final int TOO_MANY_TRAVERSALS_DEFAULT = 100000;
    @Argument(fullName="too-many-traversals",
            doc="If the assembly graph produces this many traversals, just emit contigs instead.",
            optional=true)
    private int tooManyTraversals = TOO_MANY_TRAVERSALS_DEFAULT;

    public static final int TOO_MANY_SCAFFOLDS_DEFAULT = 50000;
    @Argument(fullName="too-many-scaffolds",
            doc="If the assembly graph produces this many scaffolds, just emit traversals instead.",
            optional=true)
    private int tooManyScaffolds = TOO_MANY_SCAFFOLDS_DEFAULT;

    public static final int MIN_SV_SIZE_DEFAULT = 50;
    @Argument(fullName="min-sv-size",
            doc="Smallest variation size to count as a structural variant.", optional=true)
    public int minSVSize = MIN_SV_SIZE_DEFAULT;

    @Argument(fullName="no-scaffolding", doc="turn off scaffolding -- write traversals instead", optional=true)
    private boolean noScaffolding = false;

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

        if ( gfaFile == null ) {
            gfaFile = new GATKPath(assemblyName + ".gfa.gz");
        }
        if ( fastaFile == null ) {
            fastaFile = new GATKPath(assemblyName + ".fa.gz");
        }

        final int regionSize = getTraversalIntervals().stream().mapToInt(SimpleInterval::size).sum();
        final KmerSet<KmerAdjacency> kmerAdjacencySet = new KmerSet<>(10 * regionSize);
        kmerizeReads(reads, qMin, kmerAdjacencySet);

        List<ContigImpl> contigs = createAssembly(kmerAdjacencySet, minThinObs);
        if ( fillGaps(kmerAdjacencySet, minGapfillCount, reads) ) {
            contigs = createAssembly(kmerAdjacencySet, minThinObs);
        }

        markCycles(contigs);

        final List<Path> readPaths = pathReads(kmerAdjacencySet, reads);
        final Map<Contig,List<TransitPairCount>> contigTransitsMap =
                collectTransitPairCounts(contigs, readPaths);
        try {
            final List<Traversal> allTraversals = new ArrayList<>(
                    traverseAllPaths(contigs, readPaths, tooManyTraversals, contigTransitsMap));
            contigs.sort(Comparator.comparingInt(ContigImpl::getId));
            writeGFA(gfaFile, contigs, allTraversals);
            if ( noScaffolding ) {
                writeTraversals(fastaFile, assemblyName, allTraversals);
                return null;
            }
            try {
                writeTraversals(fastaFile, assemblyName,
                        createScaffolds(allTraversals, tooManyScaffolds, minSVSize));
            } catch ( final AssemblyTooComplexException x ) {
                logger.warn("Assembly too complex for scaffolding. Writing traversals to fasta-file");
                writeTraversals(fastaFile, assemblyName, allTraversals);
            }
        } catch ( final AssemblyTooComplexException x ) {
            logger.warn("Assembly too complex to traverse.  Writing contigs as traversals to fasta-file");
            final Collection<Traversal> contigTraversals = new ArrayList<>(contigs.size());
            for ( final Contig contig : contigs ) {
                contigTraversals.add(new Traversal(Collections.singletonList(contig)));
            }
            writeTraversals(fastaFile, assemblyName, contigTraversals);
        }
        return null;
    }

    private static List<ContigImpl> createAssembly( final KmerSet<KmerAdjacency> kmerAdjacencySet,
                                                    final int minThinObs ) {
        final List<ContigImpl> contigs = buildContigs(kmerAdjacencySet);
        connectContigs(contigs);
        removeThinContigs(contigs, minThinObs, kmerAdjacencySet);
        weldPipes(contigs);
        return contigs;
    }

    /** trim read pairs of base calls that have gone past the end of a short fragment */
    @VisibleForTesting
    static void trimOverruns( final GATKRead read, final GATKRead mate ) {
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

    private static void trimClips( final GATKRead fwd, final GATKRead rev ) {
        final List<CigarElement> fwdElements = fwd.getCigarElements();
        final List<CigarElement> revElements = rev.getCigarElements();
        final int lastFwdElementIdx = fwdElements.size() - 1;
        final int lastRevElementIdx = revElements.size() - 1;
        final CigarElement fwdLastElement = fwdElements.get(lastFwdElementIdx);
        final CigarElement revLastElement = revElements.get(lastRevElementIdx);
        final CigarElement fwdFirstElement = fwdElements.get(0);
        final CigarElement revFirstElement = revElements.get(0);
        if ( fwdFirstElement.getOperator() == CigarOperator.M &&
                fwdLastElement.getOperator() == CigarOperator.S &&
                revFirstElement.getOperator() == CigarOperator.S &&
                revLastElement.getOperator() == CigarOperator.M ) {
            final byte[] fwdBases = fwd.getBasesNoCopy();
            final int lastElementLen = fwdLastElement.getLength();
            fwd.setBases(Arrays.copyOfRange(fwdBases, 0, fwdBases.length - lastElementLen));
            final byte[] fwdQuals = fwd.getBaseQualitiesNoCopy();
            if ( fwdQuals.length > 0 ) {
                final int qualsLen = fwdQuals.length - lastElementLen;
                fwd.setBaseQualities(Arrays.copyOfRange(fwdQuals, 0, qualsLen));
            }
            final List<CigarElement> newFwdElements = new ArrayList<>(fwdElements);
            newFwdElements.set(lastFwdElementIdx, new CigarElement(lastElementLen, CigarOperator.H));
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

    @VisibleForTesting
    static void kmerizeReads( final List<GATKRead> reads,
                              final byte qMin,
                              final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        for ( final GATKRead read : reads ) {
            final byte[] calls = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            KmerAdjacency.kmerize(calls, quals, qMin, kmerAdjacencySet);
        }
    }

    /** gather unbranched strings of kmers into contigs */
    @VisibleForTesting
    static List<ContigImpl> buildContigs( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        // This method identifies each KmerAdjacency that is a contig start or end, and then builds
        // a contig from that start or end. That's actually a lie: it builds contigs at the starts,
        // and builds contigs from the reverse complement of the ends (because that's the start of a
        // contig on the other strand) to economize on code.
        // A KmerAdjacency is a contig start if:
        // 1) it has more than 1 predecessor, or no predecessors
        // 2) it has a single predecessor, but that predecessor has multiple successors
        // 3) its predecessor is its own reverse complement (i.e., a palindromic hairpin)
        // Similarly, a KmerAdjacency is a contig end if:
        // 1) it has more than 1 successor, or no successors
        // 2) it has a single successor, but that successor has multiple predecessors
        // 3) its successor is its own reverse complement
        final List<ContigImpl> contigs = new ArrayList<>();
        int nContigs = 0;
        for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
            if ( kmerAdjacency.getContig() == null ) {
                ContigImpl contig = null;
                final KmerAdjacency predecessor = kmerAdjacency.getSolePredecessor();
                // if we should start a contig with this kmer
                if ( predecessor == null || // yes, no predecessors or more than 1
                        predecessor.getSuccessorCount() > 1 || // yes, predecessor has multiple successors
                        predecessor == kmerAdjacency.rc() ) { // yes, predecessor folds back as a palindrome
                    contig = new ContigImpl(++nContigs, kmerAdjacency);
                } else {
                    // if we should end a contig with this kmer (actually we'll start a contig with the RC)
                    final KmerAdjacency successor = kmerAdjacency.getSoleSuccessor();
                    if ( successor == null || // yes, no successors or more than 1
                            successor.getPredecessorCount() > 1 || // yes, successor has multiple predecessors
                            successor == kmerAdjacency.rc() ) { // yes, successor folds back as a palindrome
                        contig = new ContigImpl(++nContigs, kmerAdjacency.rc());
                    }
                }
                if ( contig != null ) {
                    setKmerContig(contig);
                    contigs.add(contig);
                }
            }
        }

        // if there are smooth circles like a plasmid, gather them together as a contig, too
        for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
            if ( kmerAdjacency.getContig() == null ) {
                final ContigImpl contig = new ContigImpl(++nContigs, kmerAdjacency);
                setKmerContig(contig);
                contigs.add(contig);
            }
        }

        return contigs;
    }

    /** connect contigs when the final kmer of one contig is adjacent to the inital contig of another */
    @VisibleForTesting
    static void connectContigs( final List<ContigImpl> contigs ) {
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
                        if ( contigEndKmer == null ) {
                            throw new GATKException("missing contig end kmer");
                        }
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
                        if ( contigEndKmer == null ) {
                            throw new GATKException("missing contig end kmer");
                        }
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

    /** remove contigs that have little evidence */
    @VisibleForTesting
    static void removeThinContigs( final List<ContigImpl> contigs,
                                   final int minThinObs,
                                   final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        contigs.sort(Comparator.comparingInt(ContigImpl::getMaxObservations));
        boolean contigRemoved;
        final WalkDataFactory walkDataFactory = new WalkDataFactory();
        do {
            // figure out which contigs are cut points
            // i.e., those contigs which, if removed, would result in a graph with more connected components
            final int nContigs = contigs.size();
            final Map<Contig, WalkData> cutDataMap = new HashMap<>(nContigs * 3);

            for ( final ContigImpl contig : contigs ) {
                if ( cutDataMap.containsKey(contig) ) {
                    continue;
                }

                cutDataMap.put(contig, walkDataFactory.create());
                int children = 0;
                for ( final Contig nextContig : contig.getSuccessors() ) {
                    if ( !cutDataMap.containsKey(nextContig) ) {
                        findCuts(nextContig, contig, walkDataFactory, cutDataMap);
                        children += 1;
                    }
                }
                for ( final Contig nextContig : contig.getPredecessors() ) {
                    if ( !cutDataMap.containsKey(nextContig) ) {
                        findCuts(nextContig, contig, walkDataFactory, cutDataMap);
                        children += 1;
                    }
                }
                if ( children >= 2 ) {
                    contig.setMarked(true);
                }
            }

            // remove poorly attested (low max observations) contigs, unless they are cut points
            contigRemoved = false;
            final Iterator<ContigImpl> itr = contigs.iterator();
            while ( itr.hasNext() ) {
                final Contig contig = itr.next();
                // TODO: Think about replacing the heuristic "minThinObs" with something that
                //       takes the observation depth of adjacent contigs into account.
                if ( contig.getMaxObservations() < minThinObs && !contig.isMarked() ) {
                    unlinkContig(contig, kmerAdjacencySet);
                    itr.remove();
                    contigRemoved = true;
                    break;
                }
            }
        } while ( contigRemoved );
        contigs.sort(Comparator.comparingInt(ContigImpl::getId));
    }

    private static WalkData findCuts( final Contig contig,
                                      final Contig parent,
                                      final WalkDataFactory walkDataFactory,
                                      final Map<Contig, WalkData> cutDataMap ) {
        final WalkData walkData = walkDataFactory.create();
        cutDataMap.put(contig, walkData);
        for ( final Contig nextContig : contig.getSuccessors() ) {
            if ( nextContig == parent ) {
                continue;
            }
            WalkData nextWalkData = cutDataMap.get(nextContig);
            if ( nextWalkData != null ) {
                walkData.minVisitNum = Math.min(walkData.minVisitNum, nextWalkData.visitNum);
            } else {
                nextWalkData = findCuts(nextContig, contig, walkDataFactory, cutDataMap);
                walkData.minVisitNum = Math.min(walkData.minVisitNum, nextWalkData.minVisitNum);
                if ( nextWalkData.minVisitNum >= walkData.visitNum ) {
                    contig.setMarked(true);
                }
            }
        }
        for ( final Contig nextContig : contig.getPredecessors() ) {
            if ( nextContig == parent ) {
                continue;
            }
            WalkData nextWalkData = cutDataMap.get(nextContig);
            if ( nextWalkData != null ) {
                walkData.minVisitNum = Math.min(walkData.minVisitNum, nextWalkData.visitNum);
            } else {
                nextWalkData = findCuts(nextContig, contig, walkDataFactory, cutDataMap);
                walkData.minVisitNum = Math.min(walkData.minVisitNum, nextWalkData.minVisitNum);
                if ( nextWalkData.minVisitNum >= walkData.visitNum ) {
                    contig.setMarked(true);
                }
            }
        }
        return walkData;
    }

    @VisibleForTesting
    static void unlinkContig( final Contig contig, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
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

    private static void clearKmerContig( final Contig contig ) {
        int count = 0;
        final KmerAdjacency firstKmer = contig.getFirstKmer();
        final KmerAdjacency lastKmer = contig.getLastKmer();
        for ( KmerAdjacency kmer = firstKmer; kmer != lastKmer; kmer = kmer.getSoleSuccessor() ) {
            if ( kmer == null ) {
                throw new GATKException("contig does not have a flat pipeline of kmers");
            }
            if ( kmer.getContig() == null ) {
                throw new GATKException("we've returned to a kmer we've already cleared");
            }
            kmer.clearContig();
            count += 1;
        }
        lastKmer.clearContig();
        if ( count + Kmer.KSIZE != contig.size() ) {
            throw new GATKException("kmer chain length does not equal contig size");
        }
    }

    /** Sets a pointer back to its unique containing contig for each kmer comprised by the contig */
    private static void setKmerContig( final Contig contig ) {
        int offset = 0;
        final KmerAdjacency firstKmer = contig.getFirstKmer();
        final KmerAdjacency lastKmer = contig.getLastKmer();
        for ( KmerAdjacency kmer = firstKmer; kmer != lastKmer; kmer = kmer.getSoleSuccessor() ) {
            if ( kmer == null ) {
                throw new GATKException("contig does not have a flat pipeline of kmers");
            }
            if ( kmer.getContig() != null ) {
                throw new GATKException("we've returned to a kmer we've already updated");
            }
            kmer.setContigOffset(contig, offset++);
        }
        lastKmer.setContigOffset(contig, offset);
        if ( offset + Kmer.KSIZE != contig.size() ) {
            throw new GATKException("kmer chain length does not equal contig size");
        }
    }

    /** replace adjacent contigs without branches with a single, larger contig */
    @VisibleForTesting
    static void weldPipes( final List<ContigImpl> contigs ) {
        for ( int contigIdx = 0; contigIdx != contigs.size(); ++contigIdx ) {
            final ContigImpl contig = contigs.get(contigIdx);
            if ( contig.getSuccessors().size() == 1 ) {
                final Contig successor = contig.getSuccessors().get(0);
                if ( successor != contig && successor != contig.rc() &&
                        successor.getPredecessors().size() == 1 ) {
                    contigs.set(contigIdx, join(contig.getId(), contig, successor));
                    if ( !contigs.remove(successor.canonical()) ) {
                        throw new GATKException("successor linkage is messed up");
                    }
                    contigIdx -= 1; // reconsider the new contig -- there might be more joining possible
                    continue;
                }
            }
            if ( contig.getPredecessors().size() == 1 ) {
                final Contig predecessor = contig.getPredecessors().get(0);
                if ( predecessor != contig && predecessor != contig.rc() &&
                        predecessor.getSuccessors().size() == 1 ) {
                    contigs.set(contigIdx, join(contig.getId(), predecessor, contig));
                    if ( !contigs.remove(predecessor.canonical()) ) {
                        throw new GATKException("predecessor linkage is messed up");
                    }
                    contigIdx -= 1; // reconsider
                }
            }
        }
    }

    private static ContigImpl join( final int id,
                                    final Contig predecessor,
                                    final Contig successor ) {
        final ContigImpl joinedContig = new ContigImpl(id, predecessor, successor);
        clearKmerContig(joinedContig);
        setKmerContig(joinedContig);
        return joinedContig;
    }

    @VisibleForTesting
    static void markCycles( final List<ContigImpl> contigs ) {
        for ( final Contig contig : contigs ) {
            contig.setIsCycleMember(false);
        }

        final int nContigs = contigs.size();
        final Deque<Contig> deque = new ArrayDeque<>(nContigs);
        final Map<Contig, WalkData> cutDataMap = new HashMap<>(nContigs * 3);
        final WalkDataFactory walkDataFactory = new WalkDataFactory();
        for ( final Contig contig : contigs ) {
            if ( !cutDataMap.containsKey(contig) ) {
                markCyclesRecursion(contig, deque, walkDataFactory, cutDataMap);
            }
        }
    }

    private static WalkData markCyclesRecursion( final Contig contig,
                                                 final Deque<Contig> deque,
                                                 final WalkDataFactory walkDataFactory,
                                                 final Map<Contig, WalkData> cutDataMap ) {
        final WalkData walkData = walkDataFactory.create();
        cutDataMap.put(contig, walkData);
        deque.addFirst(contig);

        for ( final Contig successor : contig.getSuccessors() ) {
            final WalkData successorWalkData = cutDataMap.get(successor);
            if ( successorWalkData == null ) {
                final int recursionVisitNum =
                        markCyclesRecursion(successor, deque, walkDataFactory, cutDataMap).minVisitNum;
                walkData.minVisitNum = Math.min(walkData.minVisitNum, recursionVisitNum);
            } else {
                walkData.minVisitNum = Math.min(walkData.minVisitNum, successorWalkData.visitNum);
            }
        }

        if ( walkData.visitNum == walkData.minVisitNum ) {
            Contig tig = deque.removeFirst();
            if ( tig == contig ) {
                cutDataMap.get(tig).visitNum = Integer.MAX_VALUE;

                // single-vertex component -- cyclic only if self-referential
                if ( tig.getSuccessors().contains(tig) ) {
                    tig.setIsCycleMember(true);
                }
            } else {
                while ( true ) {
                    // kill cross-links
                    cutDataMap.get(tig).visitNum = Integer.MAX_VALUE;
                    tig.setIsCycleMember(true);
                    if ( tig == contig ) break;
                    tig = deque.removeFirst();
                }
            }
        }
        return walkData;
    }

    @VisibleForTesting
    static boolean fillGaps( final KmerSet<KmerAdjacency> kmerAdjacencySet,
                             final int minGapfillCount,
                             final List<GATKRead> reads ) {
        final Map<String, Integer> gapFillCounts = new HashMap<>();
        final PathBuilder pathBuilder = new PathBuilder(kmerAdjacencySet);
        for ( final GATKRead read : reads ) {
            final Path path = new Path(read.getBasesNoCopy(), pathBuilder);
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
            if ( nObservations >= minGapfillCount ) {
                KmerAdjacency.kmerize(entry.getKey(), nObservations, kmerAdjacencySet);
                newKmers = true;
            }
        }

        if ( newKmers ) {
            for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
                kmerAdjacency.clearContig();
            }
        }
        return newKmers;
    }

    @VisibleForTesting
    static List<Path> pathReads( final KmerSet<KmerAdjacency> kmerAdjacencySet,
                          final List<GATKRead> reads ) {
        final List<Path> readPaths = new ArrayList<>(reads.size());
        final PathBuilder pathBuilder = new PathBuilder(kmerAdjacencySet);
        for ( final GATKRead read : reads ) {
            readPaths.add(new Path(read.getBasesNoCopy(), pathBuilder));
        }
        return readPaths;
    }

    @VisibleForTesting
    static Map<Contig,List<TransitPairCount>> collectTransitPairCounts(
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

    @VisibleForTesting
    static Set<Traversal> traverseAllPaths(
            final List<ContigImpl> contigs,
            final List<Path> readPaths,
            final int tooManyTraversals,
            final Map<Contig, List<TransitPairCount>> contigTransitsMap ) {
        final TraversalSet traversalSet = new TraversalSet(tooManyTraversals);
        final List<Contig> contigsList = new ArrayList<>();
        // build traversals from untransited contigs
        // untransited contigs are sources, sinks, or large contigs that can't be crossed by a read
        for ( final Contig contig : contigs ) {
            if ( !contigTransitsMap.containsKey(contig) ) {
                if ( contig.getSuccessors().isEmpty() && contig.getPredecessors().isEmpty() ) {
                    traversalSet.add(new Traversal(Collections.singletonList(contig)));
                } else {
                    for ( final Contig successor : contig.getSuccessors() ) {
                        traverse(successor, contig, contigsList,
                                    readPaths, contigTransitsMap, traversalSet);
                    }
                    for ( final Contig predecessor : contig.getPredecessors() ) {
                        traverse(predecessor.rc(), contig.rc(), contigsList,
                                    readPaths, contigTransitsMap, traversalSet);
                    }
                }
            }
        }

        // build traversals across acyclic contigs with transits that haven't been used yet
        for ( final Map.Entry<Contig, List<TransitPairCount>> entry :
                contigTransitsMap.entrySet() ) {
            final Contig contig = entry.getKey();
            if ( contig.isCycleMember() ) {
                continue;
            }
            for ( final TransitPairCount tpc : entry.getValue() ) {
                if ( tpc.getCount() > 0 ) {
                    tpc.resetCount();
                    final TraversalSet fwdTraversalSet = new TraversalSet(tooManyTraversals);
                    traverse(tpc.getNextContig(), contig, contigsList,
                            readPaths, contigTransitsMap, fwdTraversalSet);
                    final TraversalSet revTraversalSet = new TraversalSet(tooManyTraversals);
                    traverse(tpc.getPrevContig().rc(), contig.rc(), contigsList,
                            readPaths, contigTransitsMap, revTraversalSet);
                    for ( final Traversal revTraversal : revTraversalSet ) {
                        final Traversal revTraversalRC = revTraversal.rc();
                        for ( final Traversal fwdTraversal : fwdTraversalSet ) {
                            traversalSet.add(Traversal.combine(revTraversalRC, fwdTraversal));
                        }
                    }
                }
            }
        }

        // build traversals from any remaining unexplored transits
        for ( final Map.Entry<Contig, List<TransitPairCount>> entry :
                contigTransitsMap.entrySet() ) {
            final Contig contig = entry.getKey();
            for ( final TransitPairCount tpc : entry.getValue() ) {
                if ( tpc.getCount() > 0 ) {
                    tpc.resetCount();
                    final int lastElement = contigsList.size();
                    contigsList.add(tpc.getPrevContig());
                    traverse(tpc.getNextContig(), contig, contigsList,
                            readPaths, contigTransitsMap, traversalSet);
                    contigsList.set(lastElement, tpc.getNextContig().rc());
                    traverse(tpc.getPrevContig().rc(), contig.rc(), contigsList,
                            readPaths, contigTransitsMap, traversalSet);
                    contigsList.remove(lastElement);
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
                                  final TraversalSet traversalSet ) {
        contigsList.add(predecessor); // extend the contigsList
        if ( contig.isCycleMember() ) {
            traverseCycle(contig, contigsList, readPaths, contigTransitsMap, traversalSet);
            // about to return, remove element added 3 lines above
            contigsList.remove(contigsList.size() - 1);
            return;
        }
        final List<TransitPairCount> transits = contigTransitsMap.get(contig);
        boolean foundTransitInMap = false;
        if ( transits != null ) {
            for ( final TransitPairCount tpc : transits ) {
                if ( tpc.getPrevContig() == predecessor ) {
                    final Contig successor = tpc.getNextContig();
                    if ( predecessor == contig.rc() ) { // if palindromic
                        final int nContigs = contigsList.size();
                        if ( nContigs > 1 ) {
                            if ( successor.rc() == contigsList.get(nContigs - 2) ) {
                                continue;
                            }
                        }
                    }
                    tpc.resetCount();
                    traverse(successor, contig, contigsList,
                                readPaths, contigTransitsMap, traversalSet);
                    foundTransitInMap = true;
                }
            }
        }
        if ( !foundTransitInMap ) {
            contigsList.add(contig);
            traversalSet.add(new Traversal(contigsList));
            contigsList.remove(contigsList.size() - 1);
        }

        // undo what we did on the 1st line of this method
        contigsList.remove(contigsList.size() - 1);
    }

    private static void traverseCycle( final Contig contig,
                                       final List<Contig> contigsList,
                                       final List<Path> readPaths,
                                       final Map<Contig, List<TransitPairCount>> contigTransitsMap,
                                       final TraversalSet traversalSet ) {
        contigsList.add(contig);
        final int nContigs = contigsList.size();

        // We're here because the final element of the list is cyclic.
        // The previous element, if it exists, is a good place from which to start figuring out how
        // far the read paths lead us.
        if ( !contig.isCycleMember() ) {
            throw new GATKException("not called from a cyclic contig");
        }
        final List<Contig> contigsSubList =
                nContigs <= 2 ? contigsList : contigsList.subList(nContigs - 2, nContigs);
        final List<List<Contig>> longestPaths = findLongestPaths(contigsSubList, readPaths);
        // didn't get anywhere -- just complete the traversal
        if ( longestPaths.isEmpty() ) {
            traversalSet.add(new Traversal(contigsList));
        } else {
            // for each unique extension into the cycle
            for ( final List<Contig> path : longestPaths ) {
                // don't think this can happen, but still
                if ( path.isEmpty() ) {
                    traversalSet.add(new Traversal(contigsList));
                    continue;
                }
                final List<Contig> extendedContigsList =
                        new ArrayList<>(contigsList.size() + path.size());
                extendedContigsList.addAll(contigsList);
                // if we didn't get out of the cycle
                if ( path.get(path.size() - 1).isCycleMember() ) {
                    extendedContigsList.addAll(path);
                    traversalSet.add(new Traversal(extendedContigsList));
                } else {
                    // we found a cycle-exiting path, so extend that normally
                    for ( final Contig curContig : path ) {
                        if ( curContig.isCycleMember() ) {
                            extendedContigsList.add(curContig);
                        } else {
                            final Contig prevContig =
                                    extendedContigsList.remove(extendedContigsList.size() - 1);
                            traverse(curContig, prevContig, extendedContigsList, readPaths,
                                        contigTransitsMap, traversalSet);
                            extendedContigsList.add(prevContig);
                            break;
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

    private static List<List<Contig>> findLongestPaths( final List<Contig> toMatch,
                                                        final List<Path> readPaths ) {
        final List<List<Contig>> longestPaths = new ArrayList<>();
        for ( final Path path : readPaths ) {
            testPath(path, toMatch, longestPaths);
            testPath(path.rc(), toMatch, longestPaths);
        }
        return longestPaths;
    }

    private static void testPath( final Path path,
                                  final List<Contig> toMatch,
                                  final List<List<Contig>> longestPaths ) {
        final List<PathPart> pathParts = path.getParts();
        final int nPathParts = pathParts.size();
        final List<Contig> pathContigs =
                pathParts.stream()
                        .map(PathPart::getContig)
                        .collect(Collectors.toCollection(() -> new ArrayList<>(nPathParts)));
        final int matchIdx = Collections.indexOfSubList(pathContigs, toMatch);
        if ( matchIdx != -1 ) {
            final int suffixIdx = matchIdx + toMatch.size();
            if ( suffixIdx < nPathParts ) {
                addSubPathToLongestPaths(extractSubPath(pathContigs, suffixIdx), longestPaths);
            }
        }
    }

    private static List<Contig> extractSubPath( final List<Contig> pathContigs,
                                                final int suffixIdx ) {
        final int nPathContigs = pathContigs.size();
        Contig prev = pathContigs.get(suffixIdx - 1);
        final List<Contig> subPath = new ArrayList<>(nPathContigs - suffixIdx);
        for ( int idx = suffixIdx; idx != nPathContigs; ++idx ) {
            final Contig tig = pathContigs.get(idx);
            if ( tig == null || !prev.getSuccessors().contains(tig) ) break;
            subPath.add(tig);
            prev = tig;
        }
        return subPath;
    }

    private static void addSubPathToLongestPaths( final List<Contig> subPath,
                                                  final List<List<Contig>> longestPaths ) {
        final int nResults = longestPaths.size();
        for ( int idx = 0; idx != nResults; ++idx ) {
            final List<Contig> test = longestPaths.get(idx);
            if ( isPrefix(subPath, test) ) return;
            if ( isPrefix(test, subPath) ) {
                longestPaths.set(idx, subPath);
                return;
            }
        }
        longestPaths.add(subPath);
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

    @VisibleForTesting
    static Collection<Traversal> createScaffolds( final List<Traversal> allTraversals,
                                                  final int tooManyScaffolds,
                                                  final int minSVSize ) {
        removeTriviallyDifferentTraversals(allTraversals, minSVSize);

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
                expandTraversal(idx, touched, traversalsByFirstContig, allTraversals,
                        tooManyScaffolds, scaffolds);
            }
        }
        return scaffolds;
    }

    private static void expandTraversal( final int traversalIdx,
                                         final boolean[] touched,
                                         final Map<Contig, List<Integer>> traversalsByFirstContig,
                                         final List<Traversal> allTraversals,
                                         final int tooManyScaffolds,
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
                if ( scaffolds.size() >= tooManyScaffolds ) {
                    throw new AssemblyTooComplexException();
                }
                scaffolds.add(
                        Traversal.combineOverlappers(up.rc(), down, traversal.getContigs().size()));
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
                                            final Collection<Traversal> allTraversals,
                                            final int minSVSize ) {
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
            if ( isTriviallyDifferent(prevTraversal, curTraversal, minSVSize) ) {
                traversalIterator.remove();
            } else {
                prevTraversal = curTraversal;
            }
        }
        // remove duplicates where we have both strands surviving
        final Iterator<Traversal> traversalIterator2 = sortedTraversals.iterator();
        while ( traversalIterator2.hasNext() ) {
            final Traversal traversal = traversalIterator2.next();
            if ( sortedTraversals.contains(traversal.rc()) ) {
                traversalIterator2.remove();
            }
        }
        allTraversals.clear();
        allTraversals.addAll(sortedTraversals);
    }

    private static boolean isTriviallyDifferent( final Traversal traversal1,
                                                 final Traversal traversal2,
                                                 final int minSVSize ) {
        final Contig firstContig1 = traversal1.getFirstContig();
        final Contig lastContig1 = traversal1.getLastContig();
        final Contig firstContig2 = traversal2.getFirstContig();
        final Contig lastContig2 = traversal2.getLastContig();
        if ( firstContig1 != firstContig2 || lastContig1 != lastContig2 ) {
            return false;
        }
        final int interiorSize1 =
                traversal1.getSequenceLength() - firstContig1.size() - lastContig1.size();
        final int interiorSize2 =
                traversal2.getSequenceLength() - firstContig2.size() - lastContig2.size();

        // if the path lengths are so different that one could harbor an SV, they're not trivially different
        if ( Math.abs(interiorSize1 - interiorSize2) >= minSVSize ) {
            return false;
        }

        // if the paths are small enough that there can't be an SV's worth of differences, they're trivially different
        final int maxInteriorSize = Math.max(interiorSize1, interiorSize2);
        if ( maxInteriorSize < minSVSize ) {
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
        return (maxInteriorSize - commonLen) < minSVSize;
    }

    public static class TraversalEndpointComparator implements Comparator<Traversal> {
        @Override
        public int compare( final Traversal traversal1, final Traversal traversal2 ) {
            final List<Contig> contigs1 = traversal1.getContigs();
            final List<Contig> contigs2 = traversal2.getContigs();
            int cmp = Integer.compare(contigs1.get(0).getId(), contigs2.get(0).getId());
            if ( cmp != 0 ) {
                return cmp;
            }
            final int last1 = contigs1.size() - 1;
            final int last2 = contigs2.size() - 1;
            cmp = Integer.compare(contigs1.get(last1).getId(), contigs2.get(last2).getId());
            if ( cmp != 0 ) {
                return cmp;
            }
            // among those starting and ending at the same place, sort least observed last
            cmp = -Integer.compare(traversal1.getMinMaxObservations(),
                                    traversal2.getMinMaxObservations());
            if ( cmp != 0 ) {
                return cmp;
            }
            final int end = Math.min(last1, last2);
            for ( int idx = 1; idx < end; ++idx ) {
                cmp = Integer.compare(contigs1.get(idx).getId(), contigs2.get(idx).getId());
                if ( cmp != 0 ) {
                    return cmp;
                }
            }
            return Integer.compare(last1, last2);
        }
    }

    private static void writeGFA( final GATKPath gfaFile,
                                  final Collection<ContigImpl> contigs,
                                  final Collection<Traversal> traversals ) {
        for ( final ContigImpl contig : contigs ) {
            contig.setMarked(false);
        }

        try ( final BufferedWriter writer = createBufferedWriter(gfaFile) ) {
            writer.write("H\tVN:Z:2.0");
            writer.newLine();
            for ( final Contig contig : contigs ) {
                if ( !contig.isMarked() ) {
                    writeContig(contig, writer);
                }
            }
            for ( final Traversal traversal : traversals ) {
                writer.write(traversal.getContigs().stream()
                        .map(Contig::toRef)
                        .collect(Collectors.joining(" ", "O\t*\t", "")));
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Failed to write gfa-file " + gfaFile, ioe);
        }
    }

    private static void writeContig( final Contig contig,
                                     final BufferedWriter writer ) throws IOException {
        final Contig canonicalContig = contig.canonical();
        canonicalContig.setMarked(true);
        final CharSequence seq = canonicalContig.getSequence();
        writer.write("S\t" + canonicalContig + "\t" + seq.length() + "\t" + seq +
                "\tMO:i:" + canonicalContig.getMaxObservations() +
                "\tFO:i:" + canonicalContig.getFirstKmer().getNObservations() +
                "\tLO:i:" + canonicalContig.getLastKmer().getNObservations());
        writer.newLine();

        for ( final Contig successor : contig.getSuccessors() ) {
            if ( !successor.isMarked() ) {
                writeContig(successor, writer);
            }
            writeEdge(contig, successor, writer);
        }
        for ( final Contig predecessor : contig.getPredecessors() ) {
            if ( !predecessor.isMarked() ) {
                writeContig(predecessor, writer);
            }
        }
    }

    private static void writeEdge( final Contig contig,
                                   final Contig successor,
                                   final BufferedWriter writer ) throws IOException {
        final int contigLength = contig.getSequence().length();
        writer.write("E\t*\t" + contig.toRef() + "\t" + successor.toRef() + "\t" +
                (contigLength - Kmer.KSIZE + 1) + "\t" + contigLength + "$\t0\t" +
                (Kmer.KSIZE - 1) + "\t" + (Kmer.KSIZE - 1) + "M");
        writer.newLine();
    }

    private static void writeTraversals( final GATKPath fastaFile,
                                         final String assemblyName,
                                         final Collection<Traversal> traversals ) {
        try ( final BufferedWriter writer = createBufferedWriter(fastaFile) ) {
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
                writer.write(traversal.toString());
                writer.newLine();
                writer.write(traversal.getSequence());
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Failed to write fasta-file " + fastaFile, ioe);
        }
    }

    private static BufferedWriter createBufferedWriter( final GATKPath path ) throws IOException {
        final OutputStream os = path.getOutputStream();
        final String pathString = path.getRawInputString();
        if ( !pathString.endsWith(".gz") && !pathString.endsWith(".GZ") ) {
            return new BufferedWriter(new OutputStreamWriter(os));
        }
        return new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(os)));
    }

    /** fixed-size, immutable kmer.  usual 2-bit encoding: ACGT->0123.  low order bits are final call. */
    public static class Kmer {
        public static final int KSIZE = 31; // must be odd number less than 32
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
            return (int)(kVal ^ (kVal >>> 32));
        }
    }

    /** Set of Kmers.  Uses HopscotchSet, customized to find correct starting bin for Kmers and derivatives. */
    public static final class KmerSet<KMER extends Kmer> extends HopscotchSet<KMER> {
        public KmerSet( final int capacity ) { super(capacity); }

        @Override
        protected int hashToIndex( final Object kmer ) {
            final long positiveKval =
                    ((HopscotchSet.SPREADER * ((Kmer)kmer).getKVal()) & Long.MAX_VALUE);
            return (int)(positiveKval % capacity());
        }
    }

    /**
     *  A Kmer that remembers its predecessors and successors, and the number of times it's been observed
     *  in the assembly's input set of reads.
     *  The masks are bit-wise (1=A, 2=C, 4=G, 8=T) to show which predecessors or successors have been observed.
     *  The Kmer's position on a Contig is also tracked (in later phases of the assembly process).
     */
    public static abstract class KmerAdjacency extends Kmer {
        public KmerAdjacency( final long kVal ) { super(kVal); }

        public abstract KmerAdjacency getSolePredecessor(); // returns null if there's 0 or >1 predecessors
        public abstract int getPredecessorMask();
        public abstract int getPredecessorCount();
        public abstract void removePredecessor( final int callToRemove,
                                                final KmerSet<KmerAdjacency> kmerAdjacencySet );

        public abstract KmerAdjacency getSoleSuccessor(); // returns null if there's 0 or > 1 successors
        public abstract int getSuccessorMask();
        public abstract int getSuccessorCount();
        public abstract void removeSuccessor( final int callToRemove,
                                              final KmerSet<KmerAdjacency> kmerAdjacencySet );

        public abstract Contig getContig();
        public abstract int getContigOffset();
        // offset is 0-based measure on the contig sequence of the beginning of the kmer
        public abstract void setContigOffset( final Contig contig, final int contigOffset );
        public abstract void clearContig();

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
            sb.reverse(); // low order bits were loaded into sb first:  fix that now by reversing the sb.
            return sb.toString();
        }

        /**
         * Transform a read's calls into KmerAdjacencies, and add them to a KmerSet.
         * Skip kmers that include a call with a quality < qMin.
         * Skip kmers with non-ACGT calls.
         */
        public static void kmerize( final byte[] calls,
                                    final byte[] quals,
                                    final byte qMin,
                                    final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            int currentCount = 0; // number of calls loaded into currentKVal
            long currentKVal = 0;
            KmerAdjacency prevAdjacency = null;
            KmerAdjacency currentAdjacency = null;
            for ( int idx = 0; idx < calls.length; ++idx ) {
                if ( quals[idx] < qMin ) { // if we encounter a low-quality call
                    // take care of the most recent valid KmerAdjacency, if any
                    if ( currentAdjacency != null ) {
                        currentAdjacency.observe(prevAdjacency, null);
                    }
                    // ready ourselves to accumulate calls afresh
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
                if ( ++currentCount >= KSIZE ) { // if we've loaded enough calls to make a complete kmer
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

        /**
         * Kmerize a String.  This version is for gap fills.
         * The number of observations applies to all kmers except the 1st and last.
         */
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
                    default: throw new GATKException("unexpected base call in string to kmerize.");
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

        // Kmer lookup in KmerSet.
        // KmerSets holding KmerAdjacencies have only canonical Kmers, so RC non-canonical kmers before lookup.
        public static KmerAdjacency find( final long kVal,
                                          final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            if ( isCanonical(kVal) ) return kmerAdjacencySet.find(new Kmer(kVal & KMASK));
            final KmerAdjacency result = kmerAdjacencySet.find(new Kmer(reverseComplement(kVal)));
            return result == null ? null : result.rc();
        }

        // Kmer lookup in KmerSet.
        // KmerSets holding KmerAdjacencies have only canonical Kmers, so RC non-canonical kmers before lookup.
        // Add missing Kmers.
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

    /**
     * Class to implement KmerAdjacency for canonical Kmers.
     * In particular, a KmerSet created on KmerAdjacency contains only canonical Kmers.
     */
    public static final class KmerAdjacencyImpl extends KmerAdjacency {
        private final KmerAdjacencyRC rc; // the reverse-complement of this kmer

        // these fields won't change much after all reads have been kmerized, but they do change
        // during that process.  and maybe just a little during gap filling
        private KmerAdjacency solePredecessor; // set to null if there are no predecessors, or multiple predecessors
        private KmerAdjacency soleSuccessor; // set to null if there are no successors, or multiple successors
        private int predecessorMask; // bit mask of observed kmers preceding this one
        private int successorMask; // bit mask observed kmers following this one
        private int nObservations; // the reads in which this kmer was observed

        // these fields refer back to the graph:  they'll change as the graph structure is refined
        private Contig contig; // the contig that contains this Kmer
        private int contigOffset; // the offset within the contig where this kmer is found

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
        @Override public void setContigOffset( final Contig contig, final int contigOffset ) {
            if ( this.contig != null ) {
                throw new GATKException("Internal error: overwriting kmer contig and offset.");
            }
            this.contig = contig;
            this.contigOffset = contigOffset;
        }
        @Override public void clearContig() { contig = null; contigOffset = 0; }

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

    /**
     * Class to implement KmerAdjacency for Kmers that are the reverse-complement of a canonical Kmer.
     * In particular, a KmerSet created on KmerAdjacency contains only canonical Kmers.
     * A KmerAdjacencyRC represents the RC of each Kmer in the KmerSet.
     */
    public static final class KmerAdjacencyRC extends KmerAdjacency {
        private final KmerAdjacencyImpl rc;

        // lookup table to bit-reverse nibbles
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
        @Override public void setContigOffset( final Contig contig, final int contigOffset ) {
            rc.setContigOffset(contig.rc(), contig.size() - contigOffset - KSIZE);
        }
        @Override public void clearContig() { rc.clearContig(); }

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

    public enum ContigOrientation {
        FWD, // k-mer appears at the 5' end of the contig
        REV, // k-mer appears at the 5' end of the reverse-complemented contig
        BOTH // k-mer occurs on 5' end of the contig and its RC (can happen when the contig is a palindrome)
    }

    /** Initial or final Kmer in a Contig. */
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

    /**
     * An unbranched sequence of Kmers.
     * Each Kmer (except the last one) has a single successor, which allows enumerating the sequence
     * of Kmers in the Contig.  The sequence of base calls in the Contig is just the sequence of kmers
     * with the K-1 overlapping calls removed from adjacent kmers.
     */
    public interface Contig {
        int getId();
        String toRef(); // a GFA-format reference
        CharSequence getSequence();
        int getMaxObservations();
        KmerAdjacency getFirstKmer();
        KmerAdjacency getLastKmer();
        List<Contig> getPredecessors();
        List<Contig> getSuccessors();
        int size();
        default int getNKmers() { return size() - Kmer.KSIZE + 1; }
        Contig rc();
        boolean isCycleMember();
        void setIsCycleMember( final boolean isCycleMember );
        boolean isMarked();
        void setMarked( final boolean marked );
        boolean isCanonical();
        ContigImpl canonical();

        default boolean isPredecessor( final Contig contig ) {
            return findContig(getPredecessors(), contig);
        }
        default boolean isSuccessor( final Contig contig ) {
            return findContig(getSuccessors(), contig);
        }

        static boolean findContig( final List<Contig> contigs, final Contig contig ) {
            for ( final Contig tig : contigs ) {
                if ( contig == tig ) {
                    return true;
                }
            }
            return false;
        }
    }

    /** Simple implementation of Contig interface. */
    public static final class ContigImpl implements Contig {
        private final int id;
        private final CharSequence sequence;
        private final int maxObservations;
        private final KmerAdjacency firstKmer;
        private final KmerAdjacency lastKmer;
        private final List<Contig> predecessors;
        private final List<Contig> successors;
        private boolean cyclic;
        private boolean marked;
        private final Contig rc;

        public ContigImpl( final int id, final KmerAdjacency firstKmerAdjacency ) {
            this.id = id;
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
        }

        // create a new contig by joining two contigs
        public ContigImpl( final int id, final Contig predecessor, final Contig successor ) {
            if ( predecessor == successor || predecessor == successor.rc() ) {
                throw new GATKException("can't self-join");
            }
            if ( !checkOverlap(predecessor.getSequence(), successor.getSequence()) ) {
                throw new GATKException("sequences can't be joined");
            }
            this.id = id;
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

        @Override
        public int getId() {
            return id;
        }

        @Override
        public String toRef() { return toString() + "+"; }

        @Override
        public CharSequence getSequence() {
            return sequence;
        }

        @Override
        public int getMaxObservations() {
            return maxObservations;
        }

        @Override
        public KmerAdjacency getFirstKmer() {
            return firstKmer;
        }

        @Override
        public KmerAdjacency getLastKmer() {
            return lastKmer;
        }

        @Override
        public List<Contig> getPredecessors() {
            return predecessors;
        }

        @Override
        public List<Contig> getSuccessors() {
            return successors;
        }

        @Override
        public int size() {
            return sequence.length();
        }

        @Override
        public Contig rc() {
            return rc;
        }

        @Override
        public boolean isCycleMember() {
            return cyclic;
        }

        @Override
        public void setIsCycleMember( final boolean isCycleMember ) {
            this.cyclic = isCycleMember;
        }

        @Override
        public boolean isMarked() {
            return marked;
        }

        @Override
        public void setMarked( final boolean marked ) {
            this.marked = marked;
        }

        @Override
        public boolean isCanonical() {
            return true;
        }

        @Override
        public ContigImpl canonical() {
            return this;
        }

        @Override
        public String toString() {
            return "c" + id;
        }

        private static boolean checkOverlap( final CharSequence seq1, final CharSequence seq2 ) {
            final int seq1Len = seq1.length();
            final CharSequence seq1SubSeq = seq1.subSequence(seq1Len - Kmer.KSIZE + 1, seq1Len);
            final CharSequence seq2SubSeq = seq2.subSequence(0, Kmer.KSIZE - 1);
            return seq1SubSeq.equals(seq2SubSeq);
        }
    }

    /**
     * Implementation of Contig for the reverse-complement of some other Contig.
     * Which one is the "real" Contig, and which is the "RC" is completely arbitrary, since there
     * is no notion of canonical for Contigs.
     */
    public static final class ContigRCImpl implements Contig {
        private final CharSequence sequence;
        private final List<Contig> predecessors;
        private final List<Contig> successors;
        private final ContigImpl rc;

        public ContigRCImpl( final ContigImpl contig ) {
            this.sequence = new SequenceRC(contig.getSequence());
            this.predecessors = new ContigListRC(contig.getSuccessors());
            this.successors = new ContigListRC(contig.getPredecessors());
            this.rc = contig;
        }

        @Override public int getId() { return ~rc.getId(); }
        @Override public String toRef() { return rc.toString() + "-"; }
        @Override public CharSequence getSequence() { return sequence; }
        @Override public int getMaxObservations() { return rc.getMaxObservations(); }
        @Override public KmerAdjacency getFirstKmer() { return rc.getLastKmer().rc(); }
        @Override public KmerAdjacency getLastKmer() { return rc.getFirstKmer().rc(); }
        @Override public List<Contig> getPredecessors() { return predecessors; }
        @Override public List<Contig> getSuccessors() { return successors; }
        @Override public int size() { return sequence.length(); }
        @Override public Contig rc() { return rc; }
        @Override public boolean isCycleMember() { return rc.isCycleMember(); }
        @Override public void setIsCycleMember( final boolean isCycleMember ) { rc.setIsCycleMember(isCycleMember); }
        @Override public boolean isMarked() { return rc.isMarked(); }
        @Override public void setMarked( final boolean marked ) { rc.setMarked(marked); }
        @Override public boolean isCanonical() { return false; }
        @Override public ContigImpl canonical() { return rc; }
        @Override public String toString() { return rc.toString() + "RC"; }
    }

    /** A CharSequence that is a view of the reverse-complement of another sequence. */
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

    /** A list of Contigs that presents a reverse-complemented view of a List of Contigs. */
    public static final class ContigListRC extends AbstractList<Contig> {
        private final List<Contig> contigList;

        public ContigListRC( final List<Contig> contigList ) {
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

        public List<Contig> rc() { return contigList; }

        private int reflectIndex( final int index ) { return size() - 1 - index; }
    }

    /** A single-Contig portion of a path across the assembly graph. */
    public interface PathPart {
        Contig getContig(); // will be null for PathParts that depart from the graph (PathPartGap)
        CharSequence getSequence(); // will be null for PathParts on the graph (PathPartContig)
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

    /** A part of a path that isn't present in the graph. */
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

    /** A part of a path that is present as a sub-sequence of some Contig. */
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

    /** A helper class for Path building.
     * It transforms an array of base calls into a list of contigs (and gaps).
     */
    public static final class PathBuilder {
        private final KmerSet<KmerAdjacency> kmerAdjacencySet;
        private final List<PathPart> parts;
        private long kVal = 0;
        private int count = 0;
        private PathPart currentPathPart = null;

        public PathBuilder( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            this.kmerAdjacencySet = kmerAdjacencySet;
            this.parts = new ArrayList<>();
        }

        public List<PathPart> processCalls( final byte[] calls ) {
            parts.clear();
            kVal = 0;
            count = 0;
            currentPathPart = null;

            for ( int idx = 0; idx != calls.length; ++idx ) {
                processCall( (char)calls[idx] );
            }
            return new ArrayList<>(parts);
        }

        private void processCall( final char call ) {
            kVal <<= 2;
            switch ( call ) {
                case 'C': case 'c': kVal += 1; break;
                case 'G': case 'g': kVal += 2; break;
                case 'T': case 't': kVal += 3; break;
            }
            if ( ++count >= Kmer.KSIZE ) { // if we've seen enough calls to make a whole kmer
                processKval(call);
            }
        }

        private void processKval( final char call ) {
            final KmerAdjacency kmer = KmerAdjacencyImpl.find(kVal, kmerAdjacencySet);
            if ( kmer == null ) { // if the kmer isn't in the KmerSet
                processGap(call); // we'll make a "gap" PathPart
            } else {
                processKmer(call, kmer); // it's a part of some contig
            }
        }

        private void processGap( final char call ) {
            // if there's no current path part, or if the current part isn't a gap, create a PathPartGap
            if ( currentPathPart == null || !currentPathPart.isGap() ) {
                currentPathPart = new PathPartGap(new KmerAdjacencyImpl(kVal));
                parts.add(currentPathPart);
            } else {
                // the current path part is a PathPartGap, so just extend it
                currentPathPart.extend(call);
            }
        }

        private void processKmer( final char call, final KmerAdjacency kmer ) {
            final Contig contig = kmer.getContig();
            final int contigOffset = kmer.getContigOffset();
            if ( currentPathPart == null ) {
                // we've found a kmer, but don't have a current path part -- create one
                currentPathPart = new PathPartContig(contig, contigOffset);
                parts.add(currentPathPart);
            } else if ( contig == currentPathPart.getContig() ) {
                // our lookup is on the current path part's contig -- extend it
                extendCurrentPathPart(call, contig, contigOffset);
            } else {
                // we're moving from one contig to a new one, or from a gap onto a contig
                processNewContig(contig, contigOffset);
            }
        }

        private void extendCurrentPathPart( final char call,
                                            final Contig contig,
                                            final int contigOffset ) {
            final int curStop = currentPathPart.getStop();
            if ( contigOffset == curStop ) { // if this is just the next kmer on our current contig
                currentPathPart.extend(call);
            } else if ( contigOffset == 0 && contig.getNKmers() == curStop ) {
                // starting over on the same contig (i.e., a tight cycle)
                currentPathPart = new PathPartContig(contig, 0);
                parts.add(currentPathPart);
            } else {
                // kmer is non-contiguous on the current contig.
                // start a new path part after a zero-length gap.
                parts.add(zeroLengthGap(currentPathPart));
                currentPathPart = new PathPartContig(contig, contigOffset);
                parts.add(currentPathPart);
            }
        }

        private void processNewContig( final Contig contig, final int contigOffset ) {
            if ( currentPathPart.isGap() ) {
                // try to squash captured gaps caused by isolated, single-base sequencing errors
                if ( currentPathPart.getLength() == Kmer.KSIZE ) {
                    final int prevPartIdx = parts.size() - 2;
                    if ( prevPartIdx >= 0 ) {
                        if ( gapIsSquashed(prevPartIdx, contig, contigOffset) ) {
                            return;
                        }
                    }
                }
            } else if ( !currentPathPart.stopsAtEnd() || contigOffset != 0 ||
                    !currentPathPart.getContig().isSuccessor(contig) ) {
                // not an end-to-end join across connected contigs -- record a zero-length gap
                parts.add(zeroLengthGap(currentPathPart));
            }

            // we're jumping to a new contig.  start a new path part
            currentPathPart = new PathPartContig(contig, contigOffset);
            parts.add(currentPathPart);
        }

        private boolean gapIsSquashed( final int prevPartIdx,
                                       final Contig contig,
                                       final int contigOffset ) {
            final PathPart prevPart = parts.get(prevPartIdx);
            final Contig prevPartContig = prevPart.getContig();
            final int prevPartStart = prevPart.getStart();
            final int prevPartStop = prevPart.getStop();
            final int prevPartMaxStop =
                    prevPartContig.size() - Kmer.KSIZE + 1;
            final int newStop = contigOffset + 1;

            if ( prevPartContig == contig ) { // if the gap is internal to a single contig
                if ( contigOffset - prevPartStop == Kmer.KSIZE ) {
                    // smooth over a kmer-sized gap in a single contig by backfilling the gap
                    currentPathPart = new PathPartContig(prevPartContig, prevPartStart, newStop);
                    parts.set(prevPartIdx, currentPathPart);
                    parts.remove(prevPartIdx + 1);
                    return true;
                }
            } else if ( prevPartMaxStop - prevPartStop + contigOffset == Kmer.KSIZE ) {
                // kmer-size gap crosses from one contig to another:  backfill the gap
                parts.set(prevPartIdx,
                        new PathPartContig(prevPartContig, prevPartStart, prevPartMaxStop));
                currentPathPart = new PathPartContig(contig, 0, newStop);
                parts.set(prevPartIdx + 1, currentPathPart);
                return true;
            }
            return false;
        }

        private static PathPart zeroLengthGap( final PathPart currentPathPart ) {
            final int currentStop = currentPathPart.getStop();
            final CharSequence currentSequence = currentPathPart.getContig().getSequence();
            final CharSequence almostAKmer =
                    currentSequence.subSequence(currentStop, currentStop + Kmer.KSIZE - 1);
            return new PathPartGap(almostAKmer);
        }
    }

    /** A path through the assembly graph for something (probably a read). */
    public static final class Path {
        private final List<PathPart> parts;

        public Path( final byte[] calls, final PathBuilder pathBuilder ) {
            parts = pathBuilder.processCalls(calls);
        }

        // RCing constructor
        private Path( final Path that ) {
            final List<PathPart> thoseParts = that.parts;
            final int nParts = thoseParts.size();
            parts = new ArrayList<>(nParts);
            for ( int idx = nParts - 1; idx >= 0; --idx ) {
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

    public static final class WalkDataFactory {
        private int nextNum;

        public WalkData create() { return new WalkData(++nextNum); }
    }

    /** Per-Contig storage for depth-first graph walking. */
    public static final class WalkData {
        public int visitNum;
        public int minVisitNum;

        private WalkData( final int visitNum ) {
            this.visitNum = visitNum;
            this.minVisitNum = visitNum;
        }
    }

    /**
     * A count of the number of read Paths that cross through some Contig from some previous Contig
     * to some subsequent Contig.  This helps us with phasing, so that we limit ourselves to
     * graph traversals truly represented by the data.
     */
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

    /**
     * A list of Contigs through the assembly graph.
     * Differs from a Path in that there are no gaps, and all the Contigs are properly adjacent
     * in the graph.
     * In an earlier phase of assembly, all valid phased Traversals are produced.
     * Later, the same Traversal class is used to hook up Traversals across which we cannot phase.
     * (These are somewhat improperly called Scaffolds.)
     */
    public static final class Traversal {
        private final List<Contig> contigs;
        private final int minMaxObservations;
        private int hashCode; // Traversal is immutable, so we can stash the hashCode

        public Traversal( final Collection<Contig> contigs ) {
            if ( contigs == null || contigs.isEmpty() ) {
                throw new GATKException("null or empty list of contigs in traversal");
            }
            this.contigs = new ArrayList<>(contigs);
            int minMaxObservations = Integer.MAX_VALUE;
            for ( Contig contig : contigs ) {
                minMaxObservations = Math.min(minMaxObservations, contig.getMaxObservations());
            }
            this.minMaxObservations = minMaxObservations;
            this.hashCode = 0;
        }

        // RC constructor
        private Traversal( final Traversal traversal ) {
            final List<Contig> thoseContigs = traversal.contigs;
            this.contigs = thoseContigs instanceof ContigListRC ?
                            ((ContigListRC)thoseContigs).rc() : new ContigListRC(thoseContigs);
            this.minMaxObservations = traversal.minMaxObservations;
            this.hashCode = 0;
        }

        public List<Contig> getContigs() { return Collections.unmodifiableList(contigs); }
        public Contig getFirstContig() { return contigs.get(0); }
        public Contig getLastContig() { return contigs.get(contigs.size() - 1); }
        public Traversal rc() { return new Traversal(this); }
        public boolean isRC() { return contigs instanceof ContigListRC; }
        public int getMinMaxObservations() { return minMaxObservations; }
        public boolean isInextensible() { return getLastContig().isCycleMember(); }

        @Override
        public String toString() {
            return contigs.stream()
                    .map(Contig::toString)
                    .collect(Collectors.joining("+"));
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

        @Override public int hashCode() {
            if ( hashCode == 0 ) {
                hashCode = contigs.hashCode();
            }
            return hashCode;
        }

        @Override public boolean equals( final Object obj ) {
            if ( this == obj ) return true;
            if ( !(obj instanceof Traversal) ) return false;
            final Traversal that = (Traversal)obj;
            if ( hashCode != that.hashCode || contigs.size() != that.contigs.size() ) return false;
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

    /** Set of traversals.
     *  Rejects adding RC's of existing traversals.
     *  Explodes when it gets too big. */
    public static class TraversalSet extends HashSet<Traversal> {
        private static final long serialVersionUID = 1L;
        final int tooManyTraversals;

        public TraversalSet( final int tooManyTraversals ) {
            this.tooManyTraversals = tooManyTraversals;
        }

        @Override public boolean add( final Traversal traversal ) {
            if ( contains(traversal.rc()) ) {
                return false;
            }
            if ( size() >= tooManyTraversals ) {
                throw new AssemblyTooComplexException();
            }
            return super.add(traversal);
        }
    }

    /** Something to throw when we have too many Contigs or Traversals to proceed with assembly. */
    public final static class AssemblyTooComplexException extends RuntimeException {
        static final long serialVersionUID = -1L;
    }
}
