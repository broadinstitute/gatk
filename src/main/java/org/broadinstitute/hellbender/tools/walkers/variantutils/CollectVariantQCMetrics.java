package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;
import java.util.function.Function;

public final class CollectVariantQCMetrics extends VariantWalker{

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    public File out = null;

    private PrintStream outputStream = null;

    private List<Function<VariantContext, VariantQCMetric>> metricComputers;//Note: it's a list because the order is defined

    @Override
    public void onTraversalStart() {
        try {
            outputStream = new PrintStream(out);
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(out, e);
        }
    }

    interface VariantQCMetric{

    }
    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
         metricComputers.stream().map(computer -> computer.apply())
    }

    /**
     * Close out the new variants file.
     */
    @Override
    public Object onTraversalDone() {
        try {
            return null;
        } finally {
            outputStream.close();
        }
    }
}
