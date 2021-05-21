package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A wrapper used internally by VariantEval and related classes to pass information related to the evaluation/stratification context, without exposing the entire walker to the consumer.
 */
public class VariantEvalContext {
    private final ReferenceContext referenceContext;
    private final Map<FeatureInput<VariantContext>, List<VariantContext>> variantMap;
    private final FeatureContext featureContext;
    private final VariantEvalEngine engine;

    public VariantEvalContext(ReferenceContext referenceContext, FeatureContext featureContext, Map<FeatureInput<VariantContext>, List<VariantContext>> variantMap, VariantEvalEngine engine) {
        this.referenceContext = referenceContext;
        this.variantMap = variantMap;
        this.featureContext = featureContext;
        this.engine = engine;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    //NOTE: this will return all variants from the same start site
    public List<VariantContext> getVariantsForFeature(FeatureInput<VariantContext> featureInput) {
        return variantMap.getOrDefault(featureInput, Collections.emptyList());
    }

    // This issues a query, which will return variants overlapping the interval, and therefore is not limiting to just those starting within the provided interval.
    public <T extends Feature> List<T> queryFeaturesIncludingOverlapping(FeatureInput<T> featureInput, SimpleInterval interval) {
        return featureContext.getValues(featureInput, interval);
    }

    public Set<String> getSampleNamesForEvaluation() { return engine.getSampleNamesForEvaluation(); }

    public SAMSequenceDictionary getSequenceDictionaryForDrivingVariants() {
        return engine.getSequenceDictionaryForDrivingVariants();
    }
}
