package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;

import java.util.Arrays;
import java.util.List;

abstract class STRAnalysisWalker extends MultiVariantWalker {

    @ArgumentCollection
    public STRModel model = new STRModel();

    class MyMultiVariantInputArgumentCollection extends MultiVariantInputArgumentCollection {

        @Argument(doc = "input VCF", shortName = "V", fullName = "variant")
        public FeatureInput<VariantContext> input;

        @Argument(doc = "truth VCF calls", shortName = "truth", fullName = "truthVariant", optional = true)
        public FeatureInput<VariantContext> truth;

        @Override
        public List<String> getDrivingVariantPaths() {
            return Arrays.asList(input.getFeaturePath(), truth.getFeaturePath());

        }

        public boolean isTruth(final VariantContext context) {
            return context.getSource().equals(truth.getFeaturePath());
        }

        public boolean isCall(final VariantContext context) {
            return context.getSource().equals(input.getFeaturePath());
        }
    }

    @Override
    public MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MyMultiVariantInputArgumentCollection();
    }

    protected boolean isTruth(final VariantContext contex) {
        return ((MyMultiVariantInputArgumentCollection) this.multiVariantInputArgumentCollection).isTruth(contex);
    }

    protected boolean isCall(final VariantContext contex) {
        return ((MyMultiVariantInputArgumentCollection) this.multiVariantInputArgumentCollection).isCall(contex);
    }

    protected List<VariantContext> getTruthFeatures(final FeatureContext context) {
        return context.getValues(((MyMultiVariantInputArgumentCollection) this.multiVariantInputArgumentCollection).truth);
    }

    protected List<VariantContext> getCallFeatures(final FeatureContext context) {
        return context.getValues(((MyMultiVariantInputArgumentCollection) this.multiVariantInputArgumentCollection).input);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }


    protected VariantContext findBestMatch(final VariantContext query, final List<VariantContext> targets) {
        if (query == null || targets.isEmpty()) {
            return null;
        }
        for (final VariantContext target : targets) {
            if (query.getStart() == target.getStart()) {
                return target;
            }
        }
        for (final VariantContext target : targets) {
            if (target.overlaps(query) && target.getAlternateAlleles().size() == 1 && target.getAlternateAllele(0).isNonRefAllele()) {
                return target;
            }
        }
        return null;
    }
}
