package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;

import java.io.File;

/**
 * Super class for tools that traverse through targets.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class TargetWalker extends FeatureWalker<Target> {

    @Argument(
            doc = "target BED file.",
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            optional = false
    )
    protected File targetsFile;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(Target.class);
    }

    @Override
    public File getDrivingFeatureFile() { return targetsFile; }
}
