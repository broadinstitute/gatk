package org.broadinstitute.hellbender.tools.walkers.longread;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

/**
 * Blah
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Blah ",
        oneLineSummary = "Blah",
        programGroup = LongReadAnalysisProgramGroup.class)
@BetaFeature
public class AnnotatePacBioSVVCFWithReadSupport extends VariantWalker {

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {

    }
}
