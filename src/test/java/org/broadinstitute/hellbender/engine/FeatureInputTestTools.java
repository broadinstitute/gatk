package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;

/**
 * Test utilities involving {@link FeatureInput}s.
 * Created by jonn on 11/7/18.
 */
public class FeatureInputTestTools {

    /**
     * Create a feature input based on an input path and a name.
     * @param path A {@link String} containing the path to the backing data file for the resulting {@link FeatureInput}.
     * @param name A {@link String} containing the name of the feature input type.
     * @return A {@link FeatureInput} for the given {@code path} and {@code name}.
     */
    public static FeatureInput<? extends Feature> createFeatureInput(final String path, final String name) {
        return new FeatureInput<>(path, name);
    }

}
