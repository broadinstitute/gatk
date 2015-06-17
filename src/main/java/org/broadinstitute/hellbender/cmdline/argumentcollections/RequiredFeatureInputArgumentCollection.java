package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureInput;

import java.util.List;

/**
 * An argument collection for use with tools that accept one or more input files containing Feature records
 * (eg., BED files, hapmap files, etc.), and require at least one such input.
 *
 * For variant inputs, tools should typically use one of the VariantInputArgumentCollections instead.
 */
public final class RequiredFeatureInputArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = StandardArgumentDefinitions.FEATURE_LONG_NAME, shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME, doc = "File containing features", optional = false)
    public List<FeatureInput<Feature>> featureFiles;

}
