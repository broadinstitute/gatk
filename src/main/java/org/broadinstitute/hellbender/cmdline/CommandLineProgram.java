package org.broadinstitute.hellbender.cmdline;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.net.InetAddress;
import java.text.DecimalFormat;
import java.time.Duration;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.time.format.FormatStyle;
import java.util.ArrayList;
import java.util.List;

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
public abstract class CommandLineProgram {
    protected final Logger logger = LogManager.getLogger(this.getClass());

    @Argument(common=true, optional=true)
    public List<File> TMP_DIR = new ArrayList<>();

    @ArgumentCollection(doc="Special Arguments that have meaning to the argument parsing system.  " +
            "It is unlikely these will ever need to be accessed by the command line program")
    public SpecialArgumentsCollection specialArgumentsCollection = new SpecialArgumentsCollection();

    @Argument(fullName = StandardArgumentDefinitions.VERBOSITY_NAME, shortName = StandardArgumentDefinitions.VERBOSITY_NAME, doc = "Control verbosity of logging.", common = true, optional = true)
    public Log.LogLevel VERBOSITY = Log.LogLevel.INFO;

    @Argument(doc = "Whether to suppress job-summary info on System.err.", common=true)
    public Boolean QUIET = false;
    private final String standardUsagePreamble = CommandLineParser.getStandardUsagePreamble(getClass());

    /**
    * Initialized in parseArgs.  Subclasses may want to access this to do their
    * own validation, and then print usage using commandLineParser.
    */
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
        if (this.TMP_DIR == null) this.TMP_DIR = new ArrayList<>();
        if (this.TMP_DIR.isEmpty()) TMP_DIR.add(IOUtil.getDefaultTmpDir());

        // Build the default headers
        final ZonedDateTime startDateTime = ZonedDateTime.now();
        this.defaultHeaders.add(new StringHeader(commandLine));
        this.defaultHeaders.add(new StringHeader("Started on: " +
                                startDateTime.format(DateTimeFormatter.ofLocalizedDateTime(FormatStyle.LONG))));

        LoggingUtils.setLoggingLevel(VERBOSITY);  // propagate the VERBOSITY level to logging frameworks

        for (final File f : TMP_DIR) {
            // Intentionally not checking the return values, because it may be that the program does not
            // need a tmp_dir. If this fails, the problem will be discovered downstream.
            if (!f.exists()) f.mkdirs();
            f.setReadable(true, false);
            f.setWritable(true, false);
            System.setProperty("java.io.tmpdir", f.getAbsolutePath()); // in loop so that last one takes effect
        }

        if (!QUIET) {
            System.err.println("[" + ZonedDateTime.now().format(DateTimeFormatter.ofLocalizedDateTime(FormatStyle.LONG)) +
                                "] " + commandLine);

            // Output a one liner about who/where and what software/os we're running on
            try {
                System.err.println("[" + ZonedDateTime.now().format(DateTimeFormatter.ofLocalizedDateTime(FormatStyle.LONG)) +
                        "] Executing as " +
                        System.getProperty("user.name") + "@" + InetAddress.getLocalHost().getHostName() +
                        " on " + System.getProperty("os.name") + " " + System.getProperty("os.version") +
                        " " + System.getProperty("os.arch") + "; " + System.getProperty("java.vm.name") +
                        " " + System.getProperty("java.runtime.version") +
                        "; Version: " + commandLineParser.getVersion());

                Defaults.allDefaults().entrySet().stream().forEach(e->
                        logger.info(Defaults.class.getSimpleName() + "." + e.getKey() + " : " + e.getValue())
                );

                logger.info("Deflater " + (DeflaterFactory.usingIntelDeflater()? "IntelDeflater": "JdkDeflater"));
            }
            catch (final Exception e) { /* Unpossible! */ }
        }

        try {
            return runTool();
        } finally {
            // Emit the time even if program throws
            if (!QUIET) {
                final ZonedDateTime endDateTime = ZonedDateTime.now();
                final double elapsedMinutes = (Duration.between(startDateTime, endDateTime).toMillis()) / (1000d * 60d);
                final String elapsedString  = new DecimalFormat("#,##0.00").format(elapsedMinutes);
                System.err.println("[" + endDateTime.format(DateTimeFormatter.ofLocalizedDateTime(FormatStyle.LONG)) +
                                    "] " + getClass().getName() + " done. Elapsed time: " + elapsedString + " minutes.");
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
    */
    protected String[] customCommandLineValidation() {
        return null;
    }

    /**
    *
    * @return true if command line is valid
    */
    protected boolean parseArgs(final String[] argv) {

        commandLineParser = new CommandLineParser(this);
        final boolean ret = commandLineParser.parseArguments(System.err, argv);
        commandLine = commandLineParser.getCommandLine();
        if (!ret) {
            return false;
        }
        final String[] customErrorMessages = customCommandLineValidation();
        if (customErrorMessages != null) {
            for (final String msg : customErrorMessages) {
                System.err.println(msg);
            }
            commandLineParser.usage(System.err, false);
            return false;
        }
        return true;
    }

    /** Gets a MetricsFile with default headers already written into it. */
    protected <A extends MetricBase,B extends Comparable<?>> MetricsFile<A,B> getMetricsFile() {
        final MetricsFile<A,B> file = new MetricsFile<>();
        for (final Header h : this.defaultHeaders) {
            file.addHeader(h);
        }

        return file;
    }

    /**
     * Returns the version of this tool. It is the version stored in the manifest of the jarfile.
     */
    public final String getVersion() {
        return commandLineParser == null ? "SNAPSHOT" : commandLineParser.getVersion();
    }

    /**
     * Returns the commandline used to run this program.
     */
    public final String getCommandLine() {
        return commandLine;
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

}
