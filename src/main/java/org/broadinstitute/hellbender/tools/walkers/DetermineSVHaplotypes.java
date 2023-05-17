package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.walkers.ErrorCorrectHiFi.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ByteSequence;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;

@CommandLineProgramProperties(
        summary = "Assemble long read fragments near SV breakpoint.",
        oneLineSummary = "Assemble long read fragments near SV breakpoint.",
        programGroup = VariantEvaluationProgramGroup.class
)
public class DetermineSVHaplotypes extends VariantWalker {
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
        final SimpleInterval window = new SimpleInterval(variant.getContig(), paddedStart, paddedStart);
        final List<CallIterator> callIterators = new ArrayList<>();
        readsContext.iterator(window).forEachRemaining(read -> callIterators.add(new CallIterator(read, paddedStart)));
        final int paddedEnd = paddedStart + WINDOW_SIZE;
        for ( int refLoc = paddedStart; refLoc < paddedEnd; ++refLoc ) {
            final Map<ByteSequence, Integer> callMap = new HashMap<>(callIterators.size());
            for ( final CallIterator callIterator : callIterators ) {
                if ( callIterator.hasNext() ) {
                    final Call call = callIterator.next();
                    if ( call != null ) {
                        callMap.merge(call.getCalls(), call.getMeanQual(), Integer::sum);
                    }
                }
            }
        }
    }
}
