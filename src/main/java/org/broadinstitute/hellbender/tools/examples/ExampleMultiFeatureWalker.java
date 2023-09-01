package org.broadinstitute.hellbender.tools.examples;

import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiFeatureWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.ArrayList;
import java.util.List;

/**
 * <p> Example subclass that shows how to use the MultiFeatureWalker class.</p>
 * <h3>Inputs</h3>
 * <ul>
 *     <li>
 *         One or more feature files.
 *         These must be locus-sorted, and must all contain the same type of feature.
 *     </li>
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Example of a MultiFeatureWalker subclass.",
        oneLineSummary = "Example of a MultiFeatureWalker subclass.",
        programGroup = ExampleProgramGroup.class
)
public class ExampleMultiFeatureWalker extends MultiFeatureWalker<Feature> {
    @Argument( fullName = StandardArgumentDefinitions.FEATURE_LONG_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME )
    private List<FeatureInput<Feature>> featureInputs;

    final List<Feature> features = new ArrayList<>();

    /** Do something with the merged stream of features. */
    @Override public void apply( final Feature feature,
                                 final Object header,
                                 final ReadsContext readsContext,
                                 final ReferenceContext referenceContext ) {
        // We'll just keep track of the Features we see, in the order that we see them.
        features.add(feature);
        // And print them
        System.out.println(feature);
    }
}
