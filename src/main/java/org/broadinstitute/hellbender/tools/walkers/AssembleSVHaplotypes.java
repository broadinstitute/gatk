package org.broadinstitute.hellbender.tools.walkers;

import edu.mit.broad.tedsUtils.align.GlobalAligner;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.Traversal;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.tools.walkers.ErrorCorrectHiFi.CallIterator;
import org.broadinstitute.hellbender.utils.read.UnalignedRead;
import org.broadinstitute.hellbender.utils.read.ByteSequence;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Assemble long read fragments near SV breakpoint.",
        oneLineSummary = "Assemble long read fragments near SV breakpoint.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class AssembleSVHaplotypes extends VariantWalker {
    public static final int PADDING = 150;
    public static final int WINDOW_SIZE = 2 * PADDING;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext refContext,
                       final FeatureContext featureContext ) {
        final int paddedStart = Math.max(1, variant.getStart() + 1 - PADDING);
        final String contig = variant.getContig();
        final List<UnalignedRead> correctedReads = errorCorrect(contig, paddedStart, readsContext);
        final List<UnalignedRead> assembledReads = assemble(correctedReads);
        final ByteSequence refCalls = GenotypeSVs.getRefCalls(contig, paddedStart, WINDOW_SIZE, refContext);

        if ( GATKSVVCFConstants.SYMB_ALT_STRING_INS.equals(variant.getAttribute(GATKSVVCFConstants.SVTYPE)) ) {
            final ByteSequence altPiece = new ByteSequence(variant.getAlternateAllele(0).getBases()).subSequence(1);
            final ByteSequence altCalls = refCalls.subSequence(0, PADDING).append(altPiece).append(refCalls.subSequence(PADDING));
            System.out.println(variant.getID() + ": INS " + altPiece.length());
            for ( final UnalignedRead path : assembledReads ) {
                final ByteSequence pathCalls = path.getCalls();
                int refScore = Math.round(new GlobalAligner(pathCalls, refCalls).getScore().getScore());
                int altScore = Math.round(new GlobalAligner(pathCalls, altCalls).getScore().getScore());
                System.out.println(path.getName() + "\t" + refScore + "\t" + altScore);
            }
            System.out.println();
//        } else if ( GATKSVVCFConstants.SYMB_ALT_STRING_DEL.equals(variant.getAttribute(GATKSVVCFConstants.SVTYPE)) ) {
//            if ( variant.getLengthOnReference() < 50 ) {
//                return;
//            }
//            final int altStart = variant.getEnd() + 1;
//            final ByteSequence altPiece = GenotypeSVs.getRefCalls(contig, altStart, PADDING, refContext);
//            final ByteSequence altCalls = refCalls.subSequence(0, PADDING).append(altPiece);
//
//            System.out.println(variant.getID() + ": DEL " + variant.getLengthOnReference());
//            for ( final UnalignedRead path : assembledReads ) {
//                final ByteSequence pathCalls = path.getCalls();
//                int refScore = Math.round(new GlobalAligner(pathCalls, refCalls).getScore().getScore());
//                int altScore = Math.round(new GlobalAligner(pathCalls, altCalls).getScore().getScore());
//                System.out.println(path.getName() + "\t" + refScore + "\t" + altScore);
//            }
//            System.out.println();
        }
    }

    private List<UnalignedRead> errorCorrect( final String contig,
                                              final int windowStart,
                                              final ReadsContext readsContext ) {
        final SimpleInterval window = new SimpleInterval(contig, windowStart, windowStart);
        final Iterator<GATKRead> readsIterator = readsContext.iterator(window);
        final List<GATKRead> reads = new ArrayList<>();
        final int windowEnd = windowStart + 2 * PADDING - 1;
        while ( readsIterator.hasNext() ) {
            final GATKRead read = readsIterator.next();
            if ( read.getStart() <= windowStart && read.getEnd() >= windowEnd ) {
                reads.add(read);
            }
        }
        final List<CallIterator> callIteratorList =
                ErrorCorrectHiFi.errorCorrectWindow(reads, windowStart, windowEnd);
        final List<UnalignedRead> correctedReads = new ArrayList<>(callIteratorList.size());
        for ( final CallIterator callIterator : callIteratorList ) {
            correctedReads.add(callIterator.getUnalignedRead());
        }
        return correctedReads;
    }

    private List<UnalignedRead> assemble( final List<UnalignedRead> reads ) {
        final LocalAssembler assembler = new LocalAssembler(2*PADDING, (byte)0, 4, 3, reads);
        final List<Traversal> traversals = assembler.getPathedTraversals();
        assembler.writeGFA(new GATKPath("/tmp/assembly.gfa"), traversals);
        final List<UnalignedRead> assembledSequences = new ArrayList<>();
        int idx = 0;
        for ( final Traversal trav : traversals ) {
            final ByteSequence calls = new ByteSequence(trav.getSequence().getBytes());
            assembledSequences.add(new UnalignedRead("t" + ++idx, calls, null));
        }
        return assembledSequences;
    }
}
