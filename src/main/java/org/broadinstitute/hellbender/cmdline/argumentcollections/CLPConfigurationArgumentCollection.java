package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.util.Log;

import java.io.Serializable;
import java.nio.file.Path;
import java.util.List;

/**
 * Arguments used to configure the {@link org.broadinstitute.hellbender.cmdline.CommandLineProgram}.
 *
 * <p>The purpose of this abstract class is to allow downstream projects to decide which configuration arguments they
 * want to expose to the CLP.
 *
 * <p>In addition, this class might contain arguments for documentation-only (e.g., arguments from a wrapper script).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public interface CLPConfigurationArgumentCollection extends Serializable {

    /**
     * Return a non-empty list of temp directories.
     */
    public abstract List<Path> getTmpDirectories();

    public abstract Log.LogLevel getVerbosity();

    public abstract boolean isQuiet();

    public abstract boolean useJdkDeflater();

    public abstract boolean useJdkInflater();

    public abstract int getNioMaxReopens();
}
