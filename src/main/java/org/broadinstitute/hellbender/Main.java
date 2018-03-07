package org.broadinstitute.hellbender;

import com.google.cloud.storage.StorageException;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.DeprecatedToolsRegistry;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgramExecutor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.PicardNonZeroExitException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.runtime.RuntimeUtils;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.*;

/**
 * This is the main class of Hellbender and is the way of executing individual command line programs.
 *
 * CommandLinePrograms are listed in a single command line interface based on the java package specified to instanceMain.
 *
 * If you want your own single command line program, extend this class and override if required:
 *
 * - {@link #getPackageList()} to return a list of java packages in which to search for classes that extend CommandLineProgram.
 * - {@link #getClassList()} to return a list of single classes to include (e.g. required input pre-processing tools).
 * - {@link #getCommandLineName()} for the name of the toolkit.
 * - {@link #handleResult(Object)} for handle the result of the tool.
 * - {@link #handleNonUserException(Exception)} for handle non {@link UserException}.
 * - {@link #parseArgsForConfigSetup(String[])} for pulling command-line configuration options out and initializing the {@link org.broadinstitute.hellbender.utils.config.GATKConfig}
 *
 * Note: If any of the previous methods was overrided, {@link #main(String[])} should be implemented to instantiate your class
 * and call {@link #mainEntry(String[])} to make the changes effective.
 */
public class Main {

    static {
        /**
         * The very first thing that any GATK application does is forces the JVM locale into US English, so that we don't have
         * to think about number formatting issues.
         */
        Utils.forceJVMLocaleToUSEnglish();

        // Turn off the Picard legacy parser and opt in to Barclay syntax for Picard tools. This should be replaced
        // with a config setting once PR https://github.com/broadinstitute/gatk/pull/3447 is merged.
        System.setProperty("picard.useLegacyParser", "false");
    }

    /**
     * Provides ANSI colors for the terminal output *
     */
    private static final String KNRM = "\u001B[0m"; // reset
    private static final String RED = "\u001B[31m";
    private static final String GREEN = "\u001B[32m";
    private static final String CYAN = "\u001B[36m";
    private static final String WHITE = "\u001B[37m";
    private static final String BOLDRED = "\u001B[1m\u001B[31m";

    /**
     * exit value when an issue with the commandline is detected, ie CommandLineException.
     * This is the same value as picard uses.
     */
    private static final int COMMANDLINE_EXCEPTION_EXIT_VALUE = 1;

    /**
     * Exit value when an unrecoverable {@link UserException} occurs.
     */
    public static final int USER_EXCEPTION_EXIT_VALUE = 2;

    /**
     * Exit value used when a Picard tool returns a non-zero exit code (the actual value is displayed on the command line)
     */
    public static final int PICARD_TOOL_EXCEPTION = 4;

    /**
     * exit value when any unrecoverable exception other than {@link UserException} occurs
     */
    private static final int ANY_OTHER_EXCEPTION_EXIT_VALUE = 3;
    private static final String STACK_TRACE_ON_USER_EXCEPTION_PROPERTY = "GATK_STACKTRACE_ON_USER_EXCEPTION";

    /**
     * Prints the given message (may be null) to the provided stream, adding adornments and formatting.
     */
    protected static void printDecoratedExceptionMessage(final PrintStream ps, final Exception e, String prefix){
        Utils.nonNull(ps, "stream");
        Utils.nonNull(e, "exception");
        ps.println("***********************************************************************");
        ps.println();
        ps.println(prefix + e.getMessage());
        ps.println();
        ps.println("***********************************************************************") ;
    }

    /**
     * The packages we wish to include in our command line.
     */
    protected List<String> getPackageList() {
        final List<String> packageList = new ArrayList<>();
        packageList.addAll(Arrays.asList("org.broadinstitute.hellbender"));
        packageList.addAll(Arrays.asList("picard"));
        return packageList;
    }


    /**
     * Reads from the given command-line arguments, pulls out configuration options,
     * and initializes the configuration for this instance of Main.
     *
     * Suggested use for this is to handle downstream project configuration options and overrides.
     * For example this would allow:
     *
     *      Custom command-line arguments for use in tools
     *      Custom config file loading and initialization
     */
    protected void parseArgsForConfigSetup(final String[] args) {
        ConfigFactory.getInstance().initializeConfigurationsFromCommandLineArgs(args, "--" + StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION);
    }

    /**
     * The single classes we wish to include in our command line.
     */
    protected List<Class<? extends CommandLineProgram>> getClassList() {
        return Collections.emptyList();
    }

    /** Returns the command line that will appear in the usage. */
    protected String getCommandLineName() {
        return "";
    }

    /**
     * The main method.
     * <p/>
     * Give a list of java packages in which to search for classes that extend CommandLineProgram and a list of single CommandLineProgram classes.
     * Those will be included on the command line.
     *
     * This method is not intended to be used outside of the GATK framework and tests.
     *
     */
    public Object instanceMain(final String[] args, final List<String> packageList, final List<Class<? extends CommandLineProgram>> classList, final String commandLineName) {

        final CommandLineProgram program = setupConfigAndExtractProgram(args, packageList, classList, commandLineName);
        return runCommandLineProgram(program, args);
    }

    /**
     * Run the given command line program with the raw arguments from the command line
     * @param rawArgs these are the raw arguments from the command line, the first will be stripped off
     * @return the result of running  {program} with the given args, possibly null
     */
    protected static Object runCommandLineProgram(final CommandLineProgram program, final String[] rawArgs) {

        if (null == program) return null; // no program found!  This will happen if help was specified with no other arguments
        // we can lop off the first two arguments but it requires an array copy or alternatively we could update CLP to remove them
        // in the constructor do the former in this implementation.
        final String[] mainArgs = Arrays.copyOfRange(rawArgs, 1, rawArgs.length);
        return program.instanceMain(mainArgs);
    }

    /**
     * Set up the configuration file store and create the {@link CommandLineProgram} to run.
     * @param args Argument array passed into this invocation of {@link Main}.
     * @param packageList List of packages to include in the command-line.
     * @param classList List of single classes to include in the command-line.
     * @param commandLineName The command-line name as it appears in the usage.
     * @return The {@link CommandLineProgram} to run from this invocation of {@link Main}.
     */
    protected CommandLineProgram setupConfigAndExtractProgram(final String[] args,
                                                              final List<String> packageList,
                                                              final List<Class<? extends CommandLineProgram>> classList,
                                                              final String commandLineName ){
        // Parse our config file path from our arguments and initialize the configuration file.
        // Note: this must be here because the command-line invocation inserts into here:
        parseArgsForConfigSetup(args);

        // Get our command-line program:
        return extractCommandLineProgram(args, packageList, classList, commandLineName);
    }

    /**
     * This method is not intended to be used outside of the GATK framework and tests.
     */
    public Object instanceMain(final String[] args) {
        return instanceMain(args, getPackageList(), getClassList(), getCommandLineName());
    }

    /**
     * The entry point to the toolkit from commandline: it uses {@link #instanceMain(String[])} to run the command line
     * program and handle the returned object with {@link #handleResult(Object)}, and exit with 0.
     * If any error occurs, it handles the exception (if non-user exception, through {@link #handleNonUserException(Exception)})
     * and exit with the concrete error exit value.
     *
     * Note: this is the only method that is allowed to call System.exit (because gatk tools may be run from test harness etc)
     */
    protected final void mainEntry(final String[] args) {

        CommandLineProgram program = null;
        try {
            program = setupConfigAndExtractProgram(args, getPackageList(), getClassList(), getCommandLineName());
            final Object result = runCommandLineProgram(program, args);
            handleResult(result);
            //no explicit System.exit(0) since that causes issues when running in Yarn containers
        } catch (final CommandLineException e){
            if (program != null) {
                System.err.println(program.getUsage());
            }
            handleUserException(e);
            System.exit(COMMANDLINE_EXCEPTION_EXIT_VALUE);
        } catch (final PicardNonZeroExitException e) {
            // a Picard tool returned a non-zero exit code
            handleResult(e.getToolReturnCode());
            System.exit(PICARD_TOOL_EXCEPTION);
        } catch (final UserException e){
            handleUserException(e);
            System.exit(USER_EXCEPTION_EXIT_VALUE);
        } catch (final StorageException e) {
            handleStorageException(e);
            System.exit(ANY_OTHER_EXCEPTION_EXIT_VALUE);
        } catch (final Exception e){
            handleNonUserException(e);
            System.exit(ANY_OTHER_EXCEPTION_EXIT_VALUE);
        }
    }


    /**
     * Handle the result returned for a tool. Default implementation prints a message with the string value of the object if it is not null.
     * @param result the result of the tool (may be null)
     */
    protected void handleResult(final Object result) {
        if (result != null) {
            System.out.println("Tool returned:\n" + result);
        }
    }

    /**
     * Handle an exception that was likely caused by user error.
     * This includes {@link UserException} and {@link CommandLineException}
     *
     * Default implementation produces a pretty error message
     * and a stack trace iff {@link #printStackTraceOnUserExceptions()}
     *
     * @param e the exception to handle
     */
    protected void handleUserException(Exception e) {
        printDecoratedExceptionMessage(System.err, e, "A USER ERROR has occurred: ");

        if(printStackTraceOnUserExceptions()) {
            e.printStackTrace();
        } else {
            System.err.println(String.format(
                    "Set the system property %s (--java-options '-D%s=true') to print the stack trace.",
                    STACK_TRACE_ON_USER_EXCEPTION_PROPERTY,
                    STACK_TRACE_ON_USER_EXCEPTION_PROPERTY));
        }
    }

    /**
     * Handle any exception that does not come from the user. Default implementation prints the stack trace.
     * @param exception the exception to handle (never an {@link UserException}).
     */
    protected void handleNonUserException(final Exception exception) {
        exception.printStackTrace();
    }

    /**
     * Handle any exception that does not come from the user. Default implementation prints the stack trace.
     * @param exception the exception to handle (never an {@link UserException}).
     */
    protected void handleStorageException(final StorageException exception) {
        // HTTP error code
        System.err.println("code:      " + exception.getCode());
        // user-friendly message
        System.err.println("message:   " + exception.getMessage());
        // short reason code, eg. "invalidArgument"
        System.err.println("reason:    " + exception.getReason());
        // eg. the name of the argument that was invalid
        System.err.println("location:  " + exception.getLocation());
        // true indicates the server thinks the same request may succeed later
        System.err.println("retryable: " + exception.isRetryable());
        exception.printStackTrace();
    }

    /** The entry point to GATK from commandline. It calls {@link #mainEntry(String[])} from this instance. */
    public static void main(final String[] args) {
        new Main().mainEntry(args);
    }

    private static boolean printStackTraceOnUserExceptions() {
        return "true".equals(System.getenv(STACK_TRACE_ON_USER_EXCEPTION_PROPERTY)) || Boolean.getBoolean(STACK_TRACE_ON_USER_EXCEPTION_PROPERTY);
    }

    /**
     * Returns the command line program specified, or prints the usage and exits with exit code 1 *
     */
    private CommandLineProgram extractCommandLineProgram( final String[] args,
                                                          final List<String> packageList,
                                                          final List<Class<? extends CommandLineProgram>> classList,
                                                          final String commandLineName ) {
        /** Get the set of classes that are our command line programs **/
        final ClassFinder classFinder = new ClassFinder();
        for (final String pkg : packageList) {
            classFinder.find(pkg, picard.cmdline.CommandLineProgram.class);
            classFinder.find(pkg, CommandLineProgram.class);
        }
        String missingAnnotationClasses = "";
        final Set<Class<?>> toCheck = classFinder.getClasses();
        toCheck.addAll(classList);
        final Map<String, Class<?>> simpleNameToClass = new LinkedHashMap<>();
        for (final Class<?> clazz : toCheck) {
            if (clazz.equals(PicardCommandLineProgramExecutor.class)) {
                continue;
            }
            // No interfaces, synthetic, primitive, local, or abstract classes.
            if (ClassUtils.canMakeInstances(clazz)) {
                final CommandLineProgramProperties property = getProgramProperty(clazz);
                // Check for missing annotations
                if (null == property) {
                    if (missingAnnotationClasses.isEmpty()) missingAnnotationClasses += clazz.getSimpleName();
                    else missingAnnotationClasses += ", " + clazz.getSimpleName();
                } else { /** We should check for missing annotations later **/
                    if (simpleNameToClass.containsKey(clazz.getSimpleName())) {
                        throw new RuntimeException("Simple class name collision: " + clazz.getName());
                    }
                    simpleNameToClass.put(clazz.getSimpleName(), clazz);
                }
            }
        }
        if (!missingAnnotationClasses.isEmpty()) {
            throw new RuntimeException("The following classes are missing the required CommandLineProgramProperties annotation: " + missingAnnotationClasses);
        }

        final Set<Class<?>> classes = new LinkedHashSet<>();
        classes.addAll(simpleNameToClass.values());

        if (args.length < 1 || args[0].equals("-h") || args[0].equals("--help")) {
            printUsage(System.out, classes, commandLineName);
        } else {
            if (simpleNameToClass.containsKey(args[0])) {
                final Class<?> clazz = simpleNameToClass.get(args[0]);
                try {
                    final Object commandLineProgram = clazz.newInstance();
                    // wrap Picard CommandLinePrograms in a PicardCommandLineProgramExecutor
                    return commandLineProgram instanceof picard.cmdline.CommandLineProgram ?
                            new PicardCommandLineProgramExecutor((picard.cmdline.CommandLineProgram) commandLineProgram) :
                            (CommandLineProgram) commandLineProgram;
                } catch (final InstantiationException | IllegalAccessException e) {
                    throw new RuntimeException(e);
                }
            }
            printUsage(System.err, classes, commandLineName);
            throw new UserException(getUnknownCommandMessage(classes, args[0]));
        }
        return null;
    }

    public static CommandLineProgramProperties getProgramProperty(Class<?> clazz) {
        return clazz.getAnnotation(CommandLineProgramProperties.class);
    }

    private static class SimpleNameComparator implements Comparator<Class<?>>, Serializable {
        private static final long serialVersionUID = 1096632824580028876L;

        @Override
        public int compare(final Class<?> aClass, final Class<?> bClass) {
            return RuntimeUtils.toolDisplayName(aClass).compareTo(RuntimeUtils.toolDisplayName(bClass));
        }
    }

    private void printUsage(final PrintStream destinationStream, final Set<Class<?>> classes, final String commandLineName) {
        final StringBuilder builder = new StringBuilder();
        builder.append(BOLDRED + "USAGE: " + commandLineName + " " + GREEN + "<program name>" + BOLDRED + " [-h]\n\n" + KNRM)
                .append(BOLDRED + "Available Programs:\n" + KNRM);

        /** Group CommandLinePrograms by CommandLineProgramGroup **/
        final Map<Class<? extends CommandLineProgramGroup>, CommandLineProgramGroup> programGroupClassToProgramGroupInstance = new LinkedHashMap<>();
        final Map<CommandLineProgramGroup, List<Class<?>>> programsByGroup = new TreeMap<>(CommandLineProgramGroup.comparator);
        final Map<Class<?>, CommandLineProgramProperties> programsToProperty = new LinkedHashMap<>();
        for (final Class<?> clazz : classes) {
            // Get the command line property for this command line program
            final CommandLineProgramProperties property = getProgramProperty(clazz);
            if (null == property) {
                throw new RuntimeException(String.format("The class '%s' is missing the required CommandLineProgramProperties annotation.", clazz.getSimpleName()));
            } else if (!property.omitFromCommandLine()) { // only if they are not omit from the command line
                programsToProperty.put(clazz, property);
                // Get the command line program group for the command line property
                // NB: we want to minimize the number of times we make a new instance, hence programGroupClassToProgramGroupInstance
                CommandLineProgramGroup programGroup = programGroupClassToProgramGroupInstance.get(property.programGroup());
                if (null == programGroup) {
                    try {
                        programGroup = property.programGroup().newInstance();
                    } catch (final InstantiationException | IllegalAccessException e) {
                        throw new RuntimeException(e);
                    }
                    programGroupClassToProgramGroupInstance.put(property.programGroup(), programGroup);
                }
                List<Class<?>> programs = programsByGroup.get(programGroup);
                if (null == programs) {
                    programsByGroup.put(programGroup, programs = new ArrayList<>());
                }
                programs.add(clazz);
            }
        }

        /** Print out the programs in each group **/
        for (final Map.Entry<CommandLineProgramGroup, List<Class<?>>> entry : programsByGroup.entrySet()) {
            final CommandLineProgramGroup programGroup = entry.getKey();

            builder.append(WHITE + "--------------------------------------------------------------------------------------\n" + KNRM);
            builder.append(String.format("%s%-48s %-45s%s\n", RED, programGroup.getName() + ":", programGroup.getDescription(), KNRM));

            final List<Class<?>> sortedClasses = new ArrayList<>();
            sortedClasses.addAll(entry.getValue());
            Collections.sort(sortedClasses, new SimpleNameComparator());

            for (final Class<?> clazz : sortedClasses) {
                final CommandLineProgramProperties clpProperties = programsToProperty.get(clazz);
                if (null == clpProperties) {
                    throw new RuntimeException(String.format("Unexpected error: did not find the CommandLineProgramProperties annotation for '%s'", clazz.getSimpleName()));
                }
                builder.append(getDisplaySummaryForTool(clazz, clpProperties));
            }
            builder.append(String.format("\n"));
        }
        builder.append(WHITE + "--------------------------------------------------------------------------------------\n" + KNRM);
        destinationStream.println(builder.toString());
    }

    /**
     * Return a summary string for a command line tool suitable for display.
     * @param toolClass tool class
     * @param clpProperties {@CommandLineProgramProperties} for the tool
     * @return
     */
    protected String getDisplaySummaryForTool(final Class<?> toolClass, final CommandLineProgramProperties clpProperties) {
        final BetaFeature betaFeature = toolClass.getAnnotation(BetaFeature.class);
        final ExperimentalFeature experimentalFeature = toolClass.getAnnotation(ExperimentalFeature.class);

        final StringBuilder builder = new StringBuilder();
        final String summaryLine;

        if (experimentalFeature != null) {
            summaryLine = String.format("%s%s %s%s", RED, "(EXPERIMENTAL Tool)", CYAN, clpProperties.oneLineSummary());
        } else if (betaFeature != null) {
            summaryLine = String.format("%s%s %s%s", RED, "(BETA Tool)", CYAN, clpProperties.oneLineSummary());
        } else {
            summaryLine = String.format("%s%s", CYAN, clpProperties.oneLineSummary());
        }
        final String annotatedToolName = getDisplayNameForToolClass(toolClass);
        if (toolClass.getSimpleName().length() >= 45) {
            builder.append(String.format("%s    %s    %s%s\n", GREEN, annotatedToolName, summaryLine, KNRM));
        } else {
            builder.append(String.format("%s    %-45s%s%s\n", GREEN, annotatedToolName, summaryLine, KNRM));
        }
        return builder.toString();
    }

    /**
     * @return A display name to be used for the tool who's class is {@code clazz}.
     */
    protected String getDisplayNameForToolClass(final Class<?> clazz) {
        return RuntimeUtils.toolDisplayName(clazz);
    }

    /**
     * Get deprecation message for a tool
     * @param toolName command specified by the user
     * @return deprecation message string, or null if none
     */
    public String getToolDeprecationMessage(final String toolName) {
        return DeprecatedToolsRegistry.getToolDeprecationInfo(toolName);
    }

    /**
     * When a command does not match any known command, searches for a deprecation message, if any, or for similar
     * commands.
     * @return returns an error message including the closes match if relevant.
     */
    public String getUnknownCommandMessage(final Set<Class<?>> classes, final String command) {
        final String deprecationMessage = getToolDeprecationMessage(command);
        if (deprecationMessage != null) {
            return deprecationMessage;
        }
        return getSuggestedAlternateCommand(classes, command);
    }

    /**
     * similarity floor for matching in getUnknownCommandMessage *
     */
    private static final int HELP_SIMILARITY_FLOOR = 7;
    private static final int MINIMUM_SUBSTRING_LENGTH = 5;

    /**
     * When a command does not match any known command, searches for similar commands, using the same method as GIT *
     * @return returns an error message including the closes match if relevant.
     */
    public String getSuggestedAlternateCommand(final Set<Class<?>> classes, final String command) {
        final Map<Class<?>, Integer> distances = new LinkedHashMap<>();

        int bestDistance = Integer.MAX_VALUE;
        int bestN = 0;

        // Score against all classes
        for (final Class<?> clazz : classes) {
            final String name = clazz.getSimpleName();
            final int distance;
            if (name.equals(command)) {
                throw new RuntimeException("Command matches: " + command);
            }
            if (name.startsWith(command) || (MINIMUM_SUBSTRING_LENGTH <= command.length() && name.contains(command))) {
                distance = 0;
            } else {
                distance = StringUtil.levenshteinDistance(command, name, 0, 2, 1, 4);
            }
            distances.put(clazz, distance);

            if (distance < bestDistance) {
                bestDistance = distance;
                bestN = 1;
            } else if (distance == bestDistance) {
                bestN++;
            }
        }

        // Upper bound on the similarity score
        if (0 == bestDistance && bestN == classes.size()) {
            bestDistance = HELP_SIMILARITY_FLOOR + 1;
        }

        StringBuilder message = new StringBuilder();
        // Output similar matches
        message.append(String.format("'%s' is not a valid command.", command));
        message.append(System.lineSeparator());
        if (bestDistance < HELP_SIMILARITY_FLOOR) {
            message.append(String.format("Did you mean %s?", (bestN < 2) ? "this" : "one of these"));
            message.append(System.lineSeparator());
            for (final Class<?> clazz : classes) {
                if (bestDistance == distances.get(clazz)) {
                    message.append(String.format("        %s", clazz.getSimpleName()));
                }
            }
        }
        return message.toString();
    }
}
