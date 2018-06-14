package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Class that defines the variant arguments used for a MultiVariantWalker.  The default implementation below uses the standard --variant argument; however,
 * subclasses of {@link MultiVariantWalker} can override {@link MultiVariantWalker#getMultiVariantInputArgumentCollection} and provide their own argument pattern.
 */
public abstract class MultiVariantInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * @return List of paths to variants over which to iterate.  These will be merged and iterated as a single data source.
     */
    abstract public List<String> getDrivingVariantPaths();

    public static class DefaultMultiVariantInputArgumentCollection extends MultiVariantInputArgumentCollection {
        private static final long serialVersionUID = 1L;

        // NOTE: using List<String> rather than List<FeatureInput> here so that we can initialize the driving source of variants separately
        // from any other potential sources of Features.
        @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
                doc = "One or more VCF files containing variants", common = false, optional = false)
        public List<String> drivingVariantPaths = new ArrayList<>();

        @Override
        public List<String> getDrivingVariantPaths() {
            return drivingVariantPaths;
        }
    }
}
