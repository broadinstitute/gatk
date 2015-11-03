package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
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

    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection();

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(Target.class);
    }

    @Override
    public File getDrivingFeatureFile() {
        final File targetsFile = targetArguments.getTargetsFile();
        if (targetsFile == null) {
            throw new UserException.BadArgumentValue(TargetArgumentCollection.TARGET_FILE_FULL_NAME, "you must specify a target file");
        }
        return targetsFile;
    }
}
