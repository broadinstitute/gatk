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
import org.broadinstitute.hellbender.tools.walkers.PairWalker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.Traversal;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.AssemblyTooComplexException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.UnalignedRead;

import java.util.*;

@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Performs local assembly of small regions to discover structural variants.",
        oneLineSummary = "Local assembler for SVs",
        usageExample = "gatk LocalSVAssembler -L chr21:16187360-16187360 --ip 500 -R 38.fa.gz " +
                "-I NA19240.cram -I NA19240.distantmate.bam " +
                "--assembly-name chr21_16187360_16187360_INS --gfa-file test.gfa --fasta-file test.fa.gz",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class LocalSVAssembler extends PairWalker {
    @Argument(fullName="assembly-name", doc="Name of assembly used as a prefix for traversal names.")
    public String assemblyName;

    @Argument(fullName="gfa-file", doc="Path to assembly output in gfa format.", optional=true)
    public GATKPath gfaFile;

    @Argument(fullName="fasta-file", doc="Path to scaffolds in fasta format.", optional=true)
    public GATKPath fastaFile;

    @Argument(fullName="q-min", doc="Minimum base quality when kmerizing reads.", optional=true)
    private byte qMin = LocalAssembler.QMIN_DEFAULT;

    @Argument(fullName="min-thin-observations",
            doc="Minimum number of observations of some kmer within the contig required to " +
                    "retain the contig.", optional=true)
    private int minThinObs = LocalAssembler.MIN_THIN_OBS_DEFAULT;

    @Argument(fullName="min-gapfill-count",
            doc="Minimum number of observations of a sequence that patches a gap.", optional=true)
    private int minGapfillCount = LocalAssembler.MIN_GAPFILL_COUNT_DEFAULT;

    @Argument(fullName="too-many-traversals",
            doc="If the assembly graph produces this many traversals, just emit contigs instead.",
            optional=true)
    private int tooManyTraversals = LocalAssembler.TOO_MANY_TRAVERSALS_DEFAULT;

    @Argument(fullName="too-many-scaffolds",
            doc="If the assembly graph produces this many scaffolds, just emit traversals instead.",
            optional=true)
    private int tooManyScaffolds = LocalAssembler.TOO_MANY_SCAFFOLDS_DEFAULT;

    @Argument(fullName="min-sv-size",
            doc="Smallest variation size to count as a structural variant.", optional=true)
    public int minSVSize = LocalAssembler.MIN_SV_SIZE_DEFAULT;

    @Argument(fullName="no-scaffolding", doc="turn off scaffolding -- write traversals instead", optional=true)
    private boolean noScaffolding = false;

    private final List<UnalignedRead> reads = new ArrayList<>();

    @Override public boolean requiresIntervals() { return true; }

    @Override public void apply( final GATKRead read, final GATKRead mate ) {
        trimOverruns(read, mate);
        reads.add(new UnalignedRead(read));
        reads.add(new UnalignedRead(mate));
    }

    @Override public void applyUnpaired( final GATKRead read ) {
        reads.add(new UnalignedRead(read));
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
        final LocalAssembler localAssembler =
                new LocalAssembler(regionSize, qMin, minThinObs, minGapfillCount, reads);

        try {
            final List<Traversal> traversals = localAssembler.getAllTraversals(tooManyTraversals);
            if ( noScaffolding ) {
                writeTraversals(localAssembler, traversals);
            } else {
                try {
                    writeTraversals(localAssembler, localAssembler.getScaffolds(tooManyScaffolds, minSVSize, traversals));
                } catch ( final AssemblyTooComplexException ex ) {
                    logger.warn("Assembly too complex to scaffold: writing traversals as scaffolds");
                    writeTraversals(localAssembler, traversals);
                }
            }
        } catch ( final AssemblyTooComplexException ex ) {
            logger.warn("Assembly too complex to traverse: writing contigs as scaffolds");
            writeTraversals(localAssembler, localAssembler.getContigsAsTraversals());
        }

        return null;
    }

    private void writeTraversals( final LocalAssembler localAssembler,
                                  final List<Traversal> traversals ) {
        localAssembler.writeGFA(gfaFile, traversals);
        localAssembler.writeTraversals(fastaFile, assemblyName, traversals);
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
}
