package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Default arguments for GATK configuration.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class GATKDefaultCLPConfigurationArgumentCollection implements CLPConfigurationArgumentCollection {
    private static final long serialVersionUID = 1L;

    // TODO - TMP_DIR and QUIET are not in the standard kebab-case in the latest master
    // TODO - it was changed in the previous version of this commit/PR
    public static final String TMP_DIR_NAME = "TMP_DIR";
    public static final String QUIET_NAME = "QUIET";
    public static final String USE_JDK_DEFLATER_LONG_NAME = "use-jdk-deflater";
    public static final String USE_JDK_DEFLATER_SHORT_NAME = "jdk-deflater";
    public static final String USE_JDK_INFLATER_LONG_NAME = "use-jdk-inflater";
    public static final String USE_JDK_INFLATER_SHORT_NAME = "jdk-inflater";
    public static final String NIO_MAX_REOPENS_LONG_NAME = "gcs-max-retries";
    public static final String NIO_MAX_REOPENS_SHORT_NAME = "gcs-retries";

    @Argument(fullName = TMP_DIR_NAME, doc = "List of temp directories to use.", common=true, optional=true)
    public List<String> tmpDir = new ArrayList<>();

    @Argument(fullName = StandardArgumentDefinitions.VERBOSITY_NAME, shortName = StandardArgumentDefinitions.VERBOSITY_NAME, doc = "Control verbosity of logging.", common = true, optional = true)
    public Log.LogLevel verbosity = Log.LogLevel.INFO;

    @Argument(fullName = QUIET_NAME, doc = "Whether to suppress job-summary info on System.err.", common=true)
    public Boolean quiet = false;

    @Argument(fullName = USE_JDK_DEFLATER_LONG_NAME, shortName = USE_JDK_DEFLATER_SHORT_NAME, doc = "Whether to use the JdkDeflater (as opposed to IntelDeflater)", common=true)
    public boolean useJdkDeflater = false;

    @Argument(fullName = USE_JDK_INFLATER_LONG_NAME, shortName = USE_JDK_INFLATER_SHORT_NAME, doc = "Whether to use the JdkInflater (as opposed to IntelInflater)", common=true)
    public boolean useJdkInflater = false;

    @Argument(fullName = NIO_MAX_REOPENS_LONG_NAME, shortName = NIO_MAX_REOPENS_SHORT_NAME, doc = "If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection", optional = true)
    public int nioMaxReopens = ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries();

    // This option is here for documentation completeness.
    // This is actually parsed out in Main to initialize configuration files because
    // we need to have the configuration completely set up before we create our CommandLinePrograms.
    // (Some of the CommandLinePrograms have default values set to config values, and these are loaded
    // at class load time as static initializers).
    @Argument(fullName = StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION,
            doc = "A configuration file to use with the GATK.",
            common = true,
            optional = true)
    public String gatkConfigFile = null;

    private final List<Path> resolvedTmpDirectories = new ArrayList<>();

    @Override
    public List<Path> getTmpDirectories() {
        if (resolvedTmpDirectories.isEmpty()) {
            // Provide one temp directory if the caller didn't
            if (tmpDir == null || tmpDir.isEmpty()) {
                resolvedTmpDirectories.add(IOUtil.getDefaultTmpDirPath());
            } else {
                tmpDir.forEach(s -> resolvedTmpDirectories.add(IOUtils.getPath(s)));
            }
        }
        return Collections.unmodifiableList(resolvedTmpDirectories);
    }

    @Override
    public Log.LogLevel getVerbosity() {
        return verbosity;
    }

    @Override
    public boolean isQuiet() {
        return quiet;
    }

    @Override
    public boolean useJdkDeflater() {
        return useJdkDeflater;
    }

    @Override
    public boolean useJdkInflater() {
        return useJdkInflater;
    }

    @Override
    public int getNioMaxReopens() {
        return nioMaxReopens;
    }
}
