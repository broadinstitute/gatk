package org.broadinstitute.hellbender.cmdline;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Strings;
import com.intel.gkl.compression.IntelDeflaterFactory;
import com.intel.gkl.compression.IntelInflaterFactory;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockGunzipper;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.RuntimeUtils;

import java.io.IOException;
import java.net.InetAddress;
import java.nio.file.*;
import java.text.DecimalFormat;
import java.time.Duration;
import java.time.ZonedDateTime;
import java.util.*;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import java.util.stream.Collectors;

/**
 * Abstract class to facilitate writing command-line programs.
 *
 * To use:
 *
 * 1. Extend this class with a concrete class that has data members annotated with @Argument, @PositionalArguments
 * and/or @Usage annotations.
 *
 * 2. If there is any custom command-line validation, override customCommandLineValidation().  When this method is
 * called, the command line has been parsed and set into the data members of the concrete class.
 *
 * 3. Implement a method doWork().  This is called after successful command-line processing.
 * The doWork() method may return null or a result object (they are not interpreted by the toolkit and passed onto the caller).
 * doWork() may throw unchecked exceptions, which are NOT caught and passed onto the VM.
 *
 */
public abstract class CommandLineProgram implements CommandLinePluginProvider {

    // Logger is a protected instance variable here to output the correct class name
    // with concrete sub-classes of CommandLineProgram.  Since CommandLineProgram is
    // abstract, this is fine (as long as no logging has to happen statically in this class).
    protected final Logger logger = LogManager.getLogger(this.getClass());

    private static final String DEFAULT_TOOLKIT_SHORT_NAME = "GATK";

    @Argument(fullName = StandardArgumentDefinitions.TMP_DIR_NAME, common=true, optional=true, doc = "Temp directory to use.")
    public GATKPath tmpDir;

    @ArgumentCollection(doc="Special Arguments that have meaning to the argument parsing system.  " +
            "It is unlikely these will ever need to be accessed by the command line program")
    public SpecialArgumentsCollection specialArgumentsCollection = new SpecialArgumentsCollection();

    @Argument(fullName = StandardArgumentDefinitions.VERBOSITY_NAME, shortName = StandardArgumentDefinitions.VERBOSITY_NAME, doc = "Control verbosity of logging.", common = true, optional = true)
    public Log.LogLevel VERBOSITY = Log.LogLevel.INFO;

    @Argument(fullName = StandardArgumentDefinitions.QUIET_NAME, doc = "Whether to suppress job-summary info on System.err.", common=true)
    public Boolean QUIET = false;

    @Argument(fullName = StandardArgumentDefinitions.USE_JDK_DEFLATER_LONG_NAME, shortName = StandardArgumentDefinitions.USE_JDK_DEFLATER_SHORT_NAME, doc = "Whether to use the JdkDeflater (as opposed to IntelDeflater)", common=true)
    public boolean useJdkDeflater = false;

    @Argument(fullName = StandardArgumentDefinitions.USE_JDK_INFLATER_LONG_NAME, shortName = StandardArgumentDefinitions.USE_JDK_INFLATER_SHORT_NAME, doc = "Whether to use the JdkInflater (as opposed to IntelInflater)", common=true)
    public boolean useJdkInflater = false;

    @Argument(fullName = StandardArgumentDefinitions.NIO_MAX_REOPENS_LONG_NAME, shortName = StandardArgumentDefinitions.NIO_MAX_REOPENS_SHORT_NAME, doc = "If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection", optional = true)
    public int NIO_MAX_REOPENS = ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries();

    @Argument(fullName = StandardArgumentDefinitions.NIO_PROJECT_FOR_REQUESTER_PAYS_LONG_NAME, doc = "Project to bill when accessing \"requester pays\" buckets. " +
            "If unset, these buckets cannot be accessed.  User must have storage.buckets.get permission on the bucket being accessed.", optional = true)
    public String NIO_PROJECT_FOR_REQUESTER_PAYS = ConfigFactory.getInstance().getGATKConfig().gcsProjectForRequesterPays();

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

    private CommandLineParser commandLineParser;

    private final List<Header> defaultHeaders = new ArrayList<>();

    /**
     * The reconstructed commandline used to run this program. Used for logging
     * and debugging.
     */
    private String commandLine;

    /**
     * Perform initialization/setup after command-line argument parsing but before doWork() is invoked.
     * Default implementation does nothing.
     * Subclasses can override to perform initialization.
     */
    protected void onStartup() {}

    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     * @return the return value or null is there is none.
     */
    protected abstract Object doWork();

    /**
     * Perform cleanup after doWork() is finished. Always executes even if an exception is thrown during the run.
     * Default implementation does nothing.
     * Subclasses can override to perform cleanup.
     */
    protected void onShutdown() {}

    /**
     * Template method that runs the startup hook, doWork and then the shutdown hook.
     */
    public final Object runTool(){
        try {
            logger.info("Initializing engine");
            onStartup();
            logger.info("Done initializing engine");
            return doWork();
        } finally {
            logger.info("Shutting down engine");
            onShutdown();
        }
    }

    public Object instanceMainPostParseArgs() {
        // Provide one temp directory if the caller didn't
        if (tmpDir == null) {
            tmpDir = new GATKPath(System.getProperty("java.io.tmpdir"));
        }

        // Build the default headers
        final ZonedDateTime startDateTime = ZonedDateTime.now();
        this.defaultHeaders.add(new StringHeader(commandLine));
        this.defaultHeaders.add(new StringHeader("Started on: " + Utils.getDateTimeForDisplay(startDateTime)));

        LoggingUtils.setLoggingLevel(VERBOSITY);  // propagate the VERBOSITY level to logging frameworks

        // set the temp directory as a java property, checking for existence and read/write access
        final Path p = tmpDir.toPath();
        try {
            p.getFileSystem().provider().checkAccess(p, AccessMode.READ, AccessMode.WRITE);
            System.setProperty("java.io.tmpdir", IOUtils.getAbsolutePathWithoutFileProtocol(p));
        } catch (final AccessDeniedException | NoSuchFileException e) {
            // TODO: it may be that the program does not need a tmp dir
            // TODO: if it fails, the problem can be discovered downstream
            // TODO: should log a warning instead?
            throw new UserException.BadTempDir(p, "should exist and have read/write access", e);
        } catch (final IOException e) {
            // other exceptions with the tmp directory
            throw new UserException.BadTempDir(p, e.getMessage(), e);
        }

        //Set defaults (note: setting them here means they are not controllable by the user)
        if (! useJdkDeflater) {
            BlockCompressedOutputStream.setDefaultDeflaterFactory(new IntelDeflaterFactory());
        }
        if (! useJdkInflater) {
            BlockGunzipper.setDefaultInflaterFactory(new IntelInflaterFactory());
        }

        BucketUtils.setGlobalNIODefaultOptions(NIO_MAX_REOPENS, NIO_PROJECT_FOR_REQUESTER_PAYS);

        if (!QUIET) {
            printStartupMessage(startDateTime);
        }

        warnOnToolStatus();

        try {
            return runTool();
        } finally {
            // Emit the time even if program throws
            if (!QUIET) {
                final ZonedDateTime endDateTime = ZonedDateTime.now();
                final double elapsedMinutes = (Duration.between(startDateTime, endDateTime).toMillis()) / (1000d * 60d);
                final String elapsedString  = new DecimalFormat("#,##0.00").format(elapsedMinutes);
                System.err.println("[" + Utils.getDateTimeForDisplay(endDateTime) + "] " +
                        getClass().getName() + " done. Elapsed time: " + elapsedString + " minutes.");
                System.err.println("Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
            }
        }
    }

    public Object instanceMain(final String[] argv) {
        if (!parseArgs(argv)) {
            //an information only argument like help or version was specified, just exit
            return 0;
        }
        return instanceMainPostParseArgs();
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any arguments set by command-line parser can be validated.
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     * to be written to the appropriate place.
     * @throws CommandLineException if command line is invalid and handling as exception is preferred.
     */
    protected String[] customCommandLineValidation() {
        return null;
    }

    /**
     * Parse arguments and initialize any values annotated with {@link Argument}
     * @return true if program should be executed, false if an information only argument like {@link SpecialArgumentsCollection#HELP_FULLNAME} was specified
     * @throws CommandLineException if command line validation fails
     */
    protected final boolean parseArgs(final String[] argv) {

        final boolean ret = getCommandLineParser().parseArguments(System.err, argv);
        commandLine = getCommandLineParser().getCommandLine();
        if (!ret) {
            return false;
        }
        final String[] customErrorMessages = customCommandLineValidation();
        if (customErrorMessages != null) {
            throw new CommandLineException("Command Line Validation failed:" + Arrays.stream(customErrorMessages).collect(
                    Collectors.joining(", ")));
        }
        return true;

    }

    /**
     * Return the list of GATKCommandLinePluginDescriptors to be used for this CLP.
     * Default implementation returns null. Subclasses can override this to return a custom list.
     */
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() { return new ArrayList<>(); }

    /** Gets a MetricsFile with default headers already written into it. */
    protected <A extends MetricBase,B extends Comparable<?>> MetricsFile<A,B> getMetricsFile() {
        final MetricsFile<A,B> file = new MetricsFile<>();
        for (final Header h : this.defaultHeaders) {
            file.addHeader(h);
        }

        return file;
    }

    /**
     * Prints a user-friendly message on startup with some information about who we are and the
     * runtime environment.
     *
     * May be overridden by subclasses to provide a custom implementation if desired.
     *
     * @param startDateTime Startup date/time
     */
    protected void printStartupMessage(final ZonedDateTime startDateTime) {
        try {
            logger.info(Utils.dupChar('-', 60));
            logger.info(String.format("%s v%s", getToolkitName(), getVersion()));
            logger.info(getSupportInformation());
            logger.info(String.format("Executing as %s@%s on %s v%s %s",
                    System.getProperty("user.name"), InetAddress.getLocalHost().getHostName(),
                    System.getProperty("os.name"), System.getProperty("os.version"), System.getProperty("os.arch")));
            logger.info(String.format("Java runtime: %s v%s",
                    System.getProperty("java.vm.name"), System.getProperty("java.runtime.version")));
            logger.info("Start Date/Time: " + Utils.getDateTimeForDisplay(startDateTime));
            logger.info(Utils.dupChar('-', 60));
            logger.info(Utils.dupChar('-', 60));

            // Print versions of important dependencies
            printLibraryVersions();

            // Print important settings to the logger:
            printSettings();
        }
        catch (final Exception e) { /* Unpossible! */ }
    }

    /**
     * If this tool is either Experimental or Beta, return a warning message advising against use in production
     * envirogetnment.
     * @param useTerminalColor true if the message should include highlighting terminal colorization
     * @return a warning message if the tool is Beta or Experimental, otherwise null
     */
    protected String getToolStatusWarning(final boolean useTerminalColor) {
        final String KNRM = "\u001B[0m"; // reset
        final String BOLDRED = "\u001B[1m\u001B[31m";
        final int BORDER_LENGTH = 60;

        String warningMessage = null;
        if (isBetaFeature()) {
            warningMessage = String.format(
                    "\n\n%s   %s\n\n   Warning: %s is a BETA tool and is not yet ready for use in production\n\n   %s%s\n\n",
                    useTerminalColor ? BOLDRED : "",
                    Utils.dupChar('!', BORDER_LENGTH),
                    this.getClass().getSimpleName(),
                    Utils.dupChar('!', BORDER_LENGTH),
                    useTerminalColor ? KNRM : ""
            );
        }
        else if (isExperimentalFeature()) {
            warningMessage = String.format(
                    "\n\n%s   %s\n\n   Warning: %s is an EXPERIMENTAL tool and should not be used for production\n\n   %s%s\n\n",
                    useTerminalColor ? BOLDRED : "",
                    Utils.dupChar('!', BORDER_LENGTH),
                    this.getClass().getSimpleName(),
                    Utils.dupChar('!', BORDER_LENGTH),
                    useTerminalColor ? KNRM : ""
            );
        }
        return warningMessage;
    }

    /**
     * If a tool is either Experimental or Beta, log a warning against use in production a environment.
     */
    protected void warnOnToolStatus() {
        final String warningMessage = getToolStatusWarning(true);
        if (warningMessage != null) {
            logger.warn(warningMessage);
        }
    }

    /**
     * @return true if this tool has {@code BetaFeature} status.
     */
    public boolean isBetaFeature() { return this.getClass().getAnnotation(BetaFeature.class) != null; }

    /**
     * @return true if this tool has {@code ExperimentalFeature} status.
     */
    public boolean isExperimentalFeature() { return this.getClass().getAnnotation(ExperimentalFeature.class) != null; }

    /**
     * @return The name of this toolkit. The default implementation uses "Implementation-Title" from the
     *         jar manifest, or (if that's not available) the package name.
     *
     * May be overridden by subclasses to provide a custom implementation if desired.
     */
    protected String getToolkitName() {
        return RuntimeUtils.getToolkitName(this.getClass());
    }

    /**
     * @return An abbreviated name of the toolkit for this tool. Looks for "Tool-Short-Name" in the manifest by default.
     *         Uses {@link #DEFAULT_TOOLKIT_SHORT_NAME} if the manifest is unavailable.
     * Subclasses may override to do something different.
     */
    protected String getToolkitShortName() {
        final Manifest manifest = RuntimeUtils.getManifest(this.getClass());
        if( manifest != null){
            return manifest.getMainAttributes().getValue("Toolkit-Short-Name");
        } else {
            return DEFAULT_TOOLKIT_SHORT_NAME;
        }
    }

    /**
     * @return the version of this tool. It is the version stored in the manifest of the jarfile
     *          by default, or "Unavailable" if that's not available.
     *
     * May be overridden by subclasses to provide a custom implementation if desired.
     */
    public String getVersion() {
       return RuntimeUtils.getVersion(this.getClass());
    }

    /**
     * @return A String containing information about how to get support for this toolkit.
     *
     * May be overridden by subclasses to provide a custom implementation if desired.
     */
    protected String getSupportInformation() {
        return "For support and documentation go to " + HelpConstants.GATK_MAIN_SITE;
    }

    /**
     * Output versions of important dependencies to the logger.
     *
     * May be overridden by subclasses to provide a custom implementation if desired.
     */
    protected void printLibraryVersions() {
        final Manifest manifest = RuntimeUtils.getManifest(this.getClass());
        if( manifest != null ){
                final Attributes manifestAttributes = manifest.getMainAttributes();
                final String htsjdkVersion = manifestAttributes.getValue("htsjdk-Version");
                final String picardVersion = manifestAttributes.getValue("Picard-Version");
                final String sparkVersion = manifestAttributes.getValue("Spark-Version");
                logger.info("HTSJDK Version: " + (htsjdkVersion != null ? htsjdkVersion : "unknown"));
                logger.info("Picard Version: " + (picardVersion != null ? picardVersion : "unknown"));
                logger.info("Built for Spark Version: " + (sparkVersion != null ? sparkVersion : "unknown"));
        }
    }

    /**
     * Output a curated set of important settings to the logger.
     *
     * May be overridden by subclasses to specify a different set of settings to output.
     */
    protected void printSettings() {
        if ( VERBOSITY != Log.LogLevel.DEBUG ) {
            logger.info("HTSJDK Defaults.COMPRESSION_LEVEL : " + Defaults.COMPRESSION_LEVEL);
            logger.info("HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : " + Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS);
            logger.info("HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : " + Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS);
            logger.info("HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : " + Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE);
        }
        else {
            // At DEBUG verbosity, print all the HTSJDK defaults:
            Defaults.allDefaults()
                    .forEach((key, value) -> logger.info("HTSJDK " + Defaults.class.getSimpleName() + "." + key + " : " + value));
        }

        // Log the configuration options:
        ConfigFactory.logConfigFields(ConfigFactory.getInstance().getGATKConfig(), Log.LogLevel.DEBUG);

        final boolean usingIntelDeflater = (BlockCompressedOutputStream.getDefaultDeflaterFactory() instanceof IntelDeflaterFactory && ((IntelDeflaterFactory)BlockCompressedOutputStream.getDefaultDeflaterFactory()).usingIntelDeflater());
        logger.info("Deflater: " + (usingIntelDeflater ? "IntelDeflater": "JdkDeflater"));
        final boolean usingIntelInflater = (BlockGunzipper.getDefaultInflaterFactory() instanceof IntelInflaterFactory && ((IntelInflaterFactory)BlockGunzipper.getDefaultInflaterFactory()).usingIntelInflater());
        logger.info("Inflater: " + (usingIntelInflater ? "IntelInflater": "JdkInflater"));

        logger.info("GCS max retries/reopens: " + BucketUtils.getCloudStorageConfiguration(NIO_MAX_REOPENS, "").maxChannelReopens());
        if (Strings.isNullOrEmpty(NIO_PROJECT_FOR_REQUESTER_PAYS)) {
            logger.info("Requester pays: disabled");
        } else {
            logger.info("Requester pays: enabled. Billed to: " + NIO_PROJECT_FOR_REQUESTER_PAYS);
        }
    }

    /**
     * @return the commandline used to run this program, will be null if arguments have not yet been parsed
     */
    public final String getCommandLine() {
        return commandLine;
    }

    /**
     * @return get usage and help information for this command line program if it is available
     *
     */
    public final String getUsage(){
        return getCommandLineParser().usage(true, specialArgumentsCollection.SHOW_HIDDEN);
    }

    /**
     * Replaces the set of default metrics headers by the given argument.
     * The given list is copied.
     */
    public final void setDefaultHeaders(final List<Header> headers) {
        Utils.nonNull(headers);
        this.defaultHeaders.clear();
        this.defaultHeaders.addAll(headers);
    }

    /**
     * Returns the (live) list of default metrics headers used by this tool.
     */
    public final List<Header> getDefaultHeaders() {
        return this.defaultHeaders;
    }

    /**
     * @return this programs CommandLineParser.  If one is not initialized yet this will initialize it.
     */
    @VisibleForTesting
    public final CommandLineParser getCommandLineParser() {
        if( commandLineParser == null) {
            commandLineParser = new CommandLineArgumentParser(this, getPluginDescriptors(), Collections.emptySet());
        }
        return commandLineParser;
    }
}
