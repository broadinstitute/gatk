package org.broadinstitute.hellbender;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.ClassFinder;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.conf.GATKConf;
import org.broadinstitute.hellbender.conf.GATKConfBuilder;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.ClassUtils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This is the main class of Hellbender and is the way of executing individual command line programs.
 *
 * CommandLinePrograms are listed in a single command line interface based on the java package specified to instanceMain.
 *
 * If you want your own single command line program, extend this class and override if required:
 *
 * - {@link #getConfiguration()} to setup the configuration on startup.
 * - {@link #getCommandLineName()} for the name of the toolkit.
 * - {@link #handleResult(Object)} for handle the result of the tool.
 * - {@link #handleNonUserException(Exception)} for handle non {@link UserException}.
 *
 * Note: If any of the previous methods was overrided, {@link #main(String[])} should be implemented to instantiate your class
 * and call {@link #mainEntry(String[])} to make the changes effective.
 */
public class Main {

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
     * exit value when an unrecoverable {@link UserException} occurs
     */
    private static final int USER_EXCEPTION_EXIT_VALUE = 2;

    /**
     * exit value when any unrecoverable exception other than {@link UserException} occurs
     */
    private static final int ANY_OTHER_EXCEPTION_EXIT_VALUE = 3;
    private static final String STACK_TRACE_ON_USER_EXCEPTION_PROPERTY = "GATK_STACKTRACE_ON_USER_EXCEPTION";

    protected GATKConf getConfiguration() {
        return GATKConfBuilder.getDefaultConfiguration();
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
    public Object instanceMain(final String[] args, final GATKConf configuration, final String commandLineName) {
        // first set up the codecs
        FeatureManager.setCodecPackages(getConfiguration().getGATKProperty(GATKConf.CODEC_PACKAGE_LIST_KEY));
        final List<Class<? extends CommandLineProgram>> classList = configuration.getGATKProperty(GATKConf.TOOL_CLASS_LIST_KEY).stream()
                .map(c -> {
                    try {
                        final Class<? extends CommandLineProgram> clazz = Class.forName(c).asSubclass(CommandLineProgram.class);
                        return clazz;
                    } catch (ClassNotFoundException e) {
                        throw new GATKException.ConfigurationException(GATKConf.TOOL_CLASS_LIST_KEY, "class " + c + " not found", e);
                    } catch (ClassCastException e) {
                        throw new GATKException.ConfigurationException(GATKConf.TOOL_CLASS_LIST_KEY, "class " +  c + "is not a CommandLineProgram", e);
                    }
                }).collect(Collectors.toList());
        final CommandLineProgram program = extractCommandLineProgram(args, configuration.getGATKProperty(GATKConf.TOOL_PACKAGE_LIST_KEY), classList, commandLineName);
        if (null == program) return null; // no program found!
        // we can lop off the first two arguments but it requires an array copy or alternatively we could update CLP to remove them
        // in the constructor do the former in this implementation.
        final String[] mainArgs = Arrays.copyOfRange(args, 1, args.length);
        // we set the configuration to the clp
        program.setConfiguration(configuration);
        return program.instanceMain(mainArgs);
    }

    /**
     * This method is not intended to be used outside of the GATK framework and tests.
     */
    public Object instanceMain(final String[] args) {
        final GATKConf configuration = getConfiguration();
        return instanceMain(args, configuration, getCommandLineName());
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
        try {
            final Object result = instanceMain(args);
            handleResult(result);
            System.exit(0);
        } catch (final CommandLineException e){
            //the usage has already been printed so don't print it here.
            if(printStackTraceOnUserExceptions()) {
                e.printStackTrace();
            }
            System.exit(COMMANDLINE_EXCEPTION_EXIT_VALUE);
        } catch (final UserException e){
            CommandLineProgram.printDecoratedUserExceptionMessage(System.err, e);

            if(printStackTraceOnUserExceptions()) {
                e.printStackTrace();
            }
            System.exit(USER_EXCEPTION_EXIT_VALUE);
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
     * Handle any exception that does not come from the user. Default implementation prints the stack trace.
     * @param exception the exception to handle (never an {@link UserException}).
     */
    protected void handleNonUserException(final Exception exception) {
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
    private static CommandLineProgram extractCommandLineProgram(final String[] args, final List<String> packageList, final List<Class<? extends CommandLineProgram>> classList, final String commandLineName) {
        /** Get the set of classes that are our command line programs **/
        final ClassFinder classFinder = new ClassFinder();
        for (final String pkg : packageList) {
            classFinder.find(pkg, CommandLineProgram.class);
        }
        String missingAnnotationClasses = "";
        final Set<Class<?>> toCheck = classFinder.getClasses();
        toCheck.addAll(classList);
        final Map<String, Class<?>> simpleNameToClass = new LinkedHashMap<>();
        for (final Class<?> clazz : toCheck) {
            // No interfaces, synthetic, primitive, local, or abstract classes.
            if (ClassUtils.canMakeInstances(clazz)) {
                final CommandLineProgramProperties property = getProgramProperty(clazz);
                // Check for missing annotations
                if (null == property) {
                    if (missingAnnotationClasses.isEmpty()) missingAnnotationClasses += clazz.getSimpleName();
                    else missingAnnotationClasses += ", " + clazz.getSimpleName();
                } else if (!property.omitFromCommandLine()) { /** We should check for missing annotations later **/
                    if (simpleNameToClass.containsKey(clazz.getSimpleName())) {
                        throw new RuntimeException("Simple class name collision: " + clazz.getSimpleName());
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

        if (args.length < 1) {
            printUsage(classes, commandLineName);
        } else {
            if (args[0].equals("-h") || args[0].equals("--help")) {
                printUsage(classes, commandLineName);
            } else {
                if (simpleNameToClass.containsKey(args[0])) {
                    final Class<?> clazz = simpleNameToClass.get(args[0]);
                    try {
                        return (CommandLineProgram) clazz.newInstance();
                    } catch (final InstantiationException | IllegalAccessException e) {
                        throw new RuntimeException(e);
                    }
                }
                printUsage(classes, commandLineName);
                throw new UserException(getUnknownCommandMessage(classes, args[0]));
            }
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
            return aClass.getSimpleName().compareTo(bClass.getSimpleName());
        }
    }

    private static void printUsage(final Set<Class<?>> classes, final String commandLineName) {
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
            }
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

        /** Print out the programs in each group **/
        for (final Map.Entry<CommandLineProgramGroup, List<Class<?>>> entry : programsByGroup.entrySet()) {
            final CommandLineProgramGroup programGroup = entry.getKey();

            builder.append(WHITE + "--------------------------------------------------------------------------------------\n" + KNRM);
            builder.append(String.format("%s%-48s %-45s%s\n", RED, programGroup.getName() + ":", programGroup.getDescription(), KNRM));

            final List<Class<?>> sortedClasses = new ArrayList<>();
            sortedClasses.addAll(entry.getValue());
            Collections.sort(sortedClasses, new SimpleNameComparator());

            for (final Class<?> clazz : sortedClasses) {
                final CommandLineProgramProperties property = programsToProperty.get(clazz);
                if (null == property) {
                    throw new RuntimeException(String.format("Unexpected error: did not find the CommandLineProgramProperties annotation for '%s'", clazz.getSimpleName()));
                }
                if (clazz.getSimpleName().length() >= 45) {
                    builder.append(String.format("%s    %s    %s%s%s\n", GREEN, clazz.getSimpleName(), CYAN, property.oneLineSummary(), KNRM));
                } else {
                    builder.append(String.format("%s    %-45s%s%s%s\n", GREEN, clazz.getSimpleName(), CYAN, property.oneLineSummary(), KNRM));
                }
            }
            builder.append(String.format("\n"));
        }
        builder.append(WHITE + "--------------------------------------------------------------------------------------\n" + KNRM);
        System.err.println(builder.toString());
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
    public static String getUnknownCommandMessage(final Set<Class<?>> classes, final String command) {
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
