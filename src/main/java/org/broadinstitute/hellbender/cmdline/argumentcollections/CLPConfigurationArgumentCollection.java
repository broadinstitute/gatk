package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.SpecialArgumentsCollection;

import java.io.Serializable;
import java.nio.file.Path;

/**
 * Arguments used to configure the {@link org.broadinstitute.hellbender.cmdline.CommandLineProgram}.
 *
 * <p>The purpose of this abstract class is to allow downstream projects to decide which configuration
 * arguments they want to expose to the CLP.
 *
 * <p>In addition, this class might contain arguments for documentation-only (e.g., arguments from a wrapper script).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public interface CLPConfigurationArgumentCollection extends Serializable {

    /**
     * Returns {@link SpecialArgumentsCollection} for barclay.
     */
    public SpecialArgumentsCollection getSpecialArguments();

    /**
     * Return a non-empty list of temp directories.
     */
    public abstract Path getTmpDirectory();

    /**
     * Returns the verbosity used for the tool.
     */
    public abstract Log.LogLevel getVerbosity();

    /**
     * Returns {@code true} if no job summary shouldn't be output to {@link System#err}
     */
    public abstract boolean isQuiet();

    /**
     * Returns {@code true} if JDK deflater should be used.
     */
    public abstract boolean useJdkDeflater();

    /**
     * Returns {@code true} if JDK inflater should be used.
     */
    public abstract boolean useJdkInflater();

    /**
     * Returns the number of maximum re-opens for NIO.
     */
    public abstract int getNioMaxReopens();

    /**
     * Returns the project name to bill when accessing "requester pays" buckets or empty String if
     * those buckets should not be accessed.
     */
    public String getNioProjectForRequesterPays();
}
