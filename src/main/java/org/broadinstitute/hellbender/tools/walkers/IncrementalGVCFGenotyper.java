package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantLocusWalker;

import java.io.File;
import java.util.List;

public class IncrementalGVCFGenotyper extends VariantLocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private File outputFile;


    @Override
    public void apply(Locatable loc, List<VariantContext> variants, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

    }
}
