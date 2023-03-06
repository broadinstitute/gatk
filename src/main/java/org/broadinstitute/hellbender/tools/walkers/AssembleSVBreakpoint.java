package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.tools.walkers.ErrorCorrectHiFi.CallIterator;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Assemble read fragments near SV breakpoint.",
        oneLineSummary = "Assemble read fragments near SV breakpoint.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class AssembleSVBreakpoint extends VariantWalker {
    public static final int PADDING = 70;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext refContext,
                       final FeatureContext featureContext ) {
        if ( variant.hasAttribute(GATKSVVCFConstants.SVTYPE) ) {
            final int start = variant.getStart();
            final int paddedStart = Math.max(1, start - PADDING);
            final int paddedEnd = start + PADDING - 1;
            final SimpleInterval interval =
                    new SimpleInterval(variant.getContig(), paddedStart, paddedEnd);
            final List<CallIterator> callIteratorList = new ArrayList<>();
            final Iterator<GATKRead> readsIterator = readsContext.iterator(interval);
            while ( readsIterator.hasNext() ) {
                callIteratorList.add(new CallIterator(readsIterator.next(), paddedStart));
            }
            ErrorCorrectHiFi.errorCorrectWindow(callIteratorList, paddedStart, paddedEnd);
        }
    }
}
