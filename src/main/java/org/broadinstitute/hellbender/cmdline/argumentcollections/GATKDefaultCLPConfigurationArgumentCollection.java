package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.SpecialArgumentsCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.nio.file.Path;

/**
 * Default arguments for GATK configuration.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class GATKDefaultCLPConfigurationArgumentCollection implements CLPConfigurationArgumentCollection {
    private static final long serialVersionUID = 1L;

    // CONSTANTS FOR THE DOCUMENTATION LINE
    public static final String TMP_DIR_DOC = "Temp directory to use.";
    public static final String VERBOSITY_DOC = "Control verbosity of logging.";
    public static final String QUIET_DOC = "Whether to suppress job-summary info on System.err.";
    public static final String USE_JDK_DEFLATER_DOC = "Whether to use the JdkDeflater (as opposed to IntelDeflater)";
    public static final String USE_JDK_INFLATER_DOC = "Whether to use the JdkInflater (as opposed to IntelInflater)";
    public static final String NIO_MAX_REOPENS_DOC = "If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection";
    public static final String NIO_PROJECT_FOR_REQUESTER_PAYS_LONG_DOC = "Project to bill when accessing \"requester pays\" buckets. If unset, these buckets cannot be accessed.";

    @ArgumentCollection(doc="Special Arguments that have meaning to the argument parsing system.  " +
            "It is unlikely these will ever need to be accessed by the command line program")
    public SpecialArgumentsCollection specialArgumentsCollection = new SpecialArgumentsCollection();

    @Argument(fullName = StandardArgumentDefinitions.TMP_DIR_NAME, common=true, optional=true, doc = TMP_DIR_DOC)
    public String tmpDir;

    @Argument(fullName = StandardArgumentDefinitions.VERBOSITY_NAME, shortName = StandardArgumentDefinitions.VERBOSITY_NAME, doc = VERBOSITY_DOC, common = true, optional = true)
    public Log.LogLevel verbosity = Log.LogLevel.INFO;

    @Argument(fullName = StandardArgumentDefinitions.QUIET_NAME, doc = QUIET_DOC, common=true)
    public Boolean quiet = false;

    @Argument(fullName = StandardArgumentDefinitions.USE_JDK_DEFLATER_LONG_NAME, shortName = StandardArgumentDefinitions.USE_JDK_DEFLATER_SHORT_NAME, doc = USE_JDK_DEFLATER_DOC, common=true)
    public boolean useJdkDeflater = false;

    @Argument(fullName = StandardArgumentDefinitions.USE_JDK_INFLATER_LONG_NAME, shortName = StandardArgumentDefinitions.USE_JDK_INFLATER_SHORT_NAME, doc = USE_JDK_INFLATER_DOC, common=true)
    public boolean useJdkInflater = false;

    @Argument(fullName = StandardArgumentDefinitions.NIO_MAX_REOPENS_LONG_NAME, shortName = StandardArgumentDefinitions.NIO_MAX_REOPENS_SHORT_NAME, doc = NIO_MAX_REOPENS_DOC, optional = true)
    public int nioMaxReopens = ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries();

    @Argument(fullName = StandardArgumentDefinitions.NIO_PROJECT_FOR_REQUESTER_PAYS_LONG_NAME, doc = NIO_PROJECT_FOR_REQUESTER_PAYS_LONG_DOC, optional = true)
    public String nioProjectForRequesterPays = ConfigFactory.getInstance().getGATKConfig().gcsProjectForRequesterPays();

    // This option is here for documentation completeness.
    // This is actually parsed out in Main to initialize configuration files because
    // we need to have the configuration completely set up before we create our CommandLinePrograms.
    // (Some of the CommandLinePrograms have default values set to config values, and these are loaded
    // at class load time as static initializers).
    @Argument(fullName = StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION,
            doc = "A configuration file to use with the GATK.",
            common = true,
            optional = true)
    public String GATK_CONFIG_FILE = null;

    @Override
    public SpecialArgumentsCollection getSpecialArguments() {
        return specialArgumentsCollection;
    }

    @Override
    public Path getTmpDirectory() {
        // Provide one temp directory if the caller didn't
        // TODO - this should use the HTSJDK IOUtil.getDefaultTmpDirPath, which is somehow broken in the current HTSJDK version
        if (tmpDir == null || tmpDir.isEmpty()) {
            return IOUtils.getPath(System.getProperty("java.io.tmpdir"));
        }
        return IOUtils.getPath(tmpDir);
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

    @Override
    public String getNioProjectForRequesterPays() {
        return nioProjectForRequesterPays;
    }
}
