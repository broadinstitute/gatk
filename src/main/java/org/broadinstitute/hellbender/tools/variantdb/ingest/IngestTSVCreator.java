package org.broadinstitute.hellbender.tools.variantdb.ingest;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.List;

public abstract class IngestTSVCreator {

    public abstract List<String> createRow(final long start, final VariantContext variant, final String sampleId);
    public abstract void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext);
    public abstract void closeTool();
}
