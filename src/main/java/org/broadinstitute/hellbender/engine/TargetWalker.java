package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.Target;

import java.io.File;

/**
 * Super class for tools that traverse through targets.
 *
 * Target files have the format described in {@link org.broadinstitute.hellbender.tools.exome.TargetWriter}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class TargetWalker extends FeatureWalker<Target> {

    @Argument(
            doc = "target file -- not a BED file.  Should be formatted as a tsv with at least the following header columns: contig, start, stop, name.",
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
