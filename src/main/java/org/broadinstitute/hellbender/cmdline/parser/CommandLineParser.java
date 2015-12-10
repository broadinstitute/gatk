package org.broadinstitute.hellbender.cmdline.parser;

import joptsimple.OptionException;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpecBuilder;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PositionalArguments;
import org.broadinstitute.hellbender.cmdline.SpecialArgumentsCollection;
import org.broadinstitute.hellbender.cmdline.StrictBooleanConverter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Annotation-driven utility for parsing command-line arguments, checking for errors, and producing usage message.
 * <p/>
 * This class supports arguments of the form KEY=VALUE, plus positional arguments.  Positional arguments must not contain
 * an equal sign lest they be mistaken for a KEY=VALUE pair.
 * <p/>
 * The caller must supply an object that both defines the command line and has the parsed arguments set into it.
 * For each possible KEY=VALUE argument, there must be a public data member annotated with @Argument.  The KEY name is
 * the name of the data member.  An abbreviated name may also be specified with the shortName attribute of @Argument.
 * If the data member is a List<T>, then the argument may be specified multiple times.  The type of the data member,
 * or the type of the List element must either have a ctor T(String), or must be an Enum.  List arguments must
 * be initialized by the caller with some kind of list.  Any other argument that is non-null is assumed to have the given
 * value as a default.  If an argument has no default value, and does not have the optional attribute of @Argument set,
 * is required.  For List arguments, minimum and maximum number of elements may be specified in the @Argument annotation.
 * <p/>
 * A single List data member may be annotated with the @PositionalArguments.  This behaves similarly to a Argument
 * with List data member: the caller must initialize the data member, the type must be constructable from String, and
 * min and max number of elements may be specified.  If no @PositionalArguments annotation appears in the object,
 * then it is an error for the command line to contain positional arguments.
 * <p/>
 * A single String public data member may be annotated with @Usage.  This string, if present, is used to
 * construct the usage message.  Details about the possible arguments are automatically appended to this string.
 * If @Usage does not appear, a boilerplate usage message is used.
 */
public final class CommandLineParser {
    // For formatting argument section of usage message.
    static final int ARGUMENT_COLUMN_WIDTH = 30;
    static final int DESCRIPTION_COLUMN_WIDTH = 90;

    private static final String ENUM_OPTION_DOC_PREFIX = "Possible values: {";
    private static final String ENUM_OPTION_DOC_SUFFIX = "} ";

    // Use these if no @Usage annotation
    private static final String defaultUsagePreamble = "Usage: program [arguments...]\n";
    private static final String defaultUsagePreambleWithPositionalArguments =
            "Usage: program [arguments...] [positional-arguments...]\n";
    static final String NULL_STRING = "null";
    public static final String COMMENT = "#";
    public static final String POSITIONAL_ARGUMENTS_NAME = "Positional Argument";


    /**
     * A typical command line program will call this to get the beginning of the usage message,
     * and then append a description of the program, like this:
     * <p/>
     * \@Usage
     * public String USAGE = CommandLineParser.getStandardUsagePreamble(getClass()) + "Frobnicates the freebozzle."
     */
    public static String getStandardUsagePreamble(final Class<?> mainClass) {
        return "USAGE: " + mainClass.getSimpleName() + " [arguments]\n\n";
    }


    private void putInArgumentMap(ArgumentDefinition arg){
        for (String key: arg.getNames()){
            argumentMap.put(key, arg);
        }
    }

    private boolean inArgumentMap(ArgumentDefinition arg){
        for (String key: arg.getNames()){
            if(argumentMap.containsKey(key)){
                return true;
            }
        }
        return false;
    }

    // This is the object that the caller has provided that contains annotations,
    // and into which the values will be assigned.
    private final Object callerArguments;


    // null if no @PositionalArguments annotation
    private Field positionalArguments;
    private int minPositionalArguments;
    private int maxPositionalArguments;
    private Object positionalArgumentsParent;

    // List of all the data members with @Argument annotation
    private final List<ArgumentDefinition> argumentDefinitions = new ArrayList<>();

    // Maps long name, and short name, if present, to an argument definition that is
    // also in the argumentDefinitions list.
    private final Map<String, ArgumentDefinition> argumentMap = new HashMap<>();

    // In case implementation wants to get at arg for some reason.
    private String[] argv;


    // The associated program properties using the CommandLineProgramProperties annotation
    private final CommandLineProgramProperties programProperties;

    private String getUsagePreamble() {
        String usagePreamble = "";
        if (null != programProperties) {
            usagePreamble += programProperties.summary();
        } else if (positionalArguments == null) {
            usagePreamble += defaultUsagePreamble;
        } else {
            usagePreamble += defaultUsagePreambleWithPositionalArguments;
        }
        return usagePreamble;
    }


    public CommandLineParser(final Object callerArguments) {
        this.callerArguments = callerArguments;

        createArgumentDefinitions(callerArguments);

        this.programProperties = this.callerArguments.getClass().getAnnotation(CommandLineProgramProperties.class);
    }

    private void createArgumentDefinitions(final Object callerArguments) {
        for (final Field field : getAllFields(callerArguments.getClass())) {
            if (field.getAnnotation(Argument.class) != null && field.getAnnotation(ArgumentCollection.class) != null){
                throw new GATKException.CommandLineParserInternalException("An Argument cannot be an argument collection: "
                        +field.getName() + " in " + callerArguments + " is annotated as both.");
            }
            if (field.getAnnotation(PositionalArguments.class) != null) {
                handlePositionalArgumentAnnotation(field, callerArguments);
            }
            if (field.getAnnotation(Argument.class) != null) {
                handleArgumentAnnotation(field, callerArguments);
            }
            if (field.getAnnotation(ArgumentCollection.class) != null) {
                try {
                    field.setAccessible(true);
                    createArgumentDefinitions(field.get(callerArguments));
                } catch (final IllegalAccessException e) {
                    throw new GATKException.ShouldNeverReachHereException("should never reach here because we setAccessible(true)", e);
                }
            }
        }
    }

    private static List<Field> getAllFields(Class<?> clazz) {
        final List<Field> ret = new ArrayList<>();
        do {
            ret.addAll(Arrays.asList(clazz.getDeclaredFields()));
            clazz = clazz.getSuperclass();
        } while (clazz != null);
        return ret;
    }

    public String getVersion() {
        return "Version:" + this.callerArguments.getClass().getPackage().getImplementationVersion();
    }

    /**
     * Print a usage message based on the arguments object passed to the ctor.
     *
     * @param stream Where to write the usage message.
     */
    public void usage(final PrintStream stream) {
        stream.print(getStandardUsagePreamble(callerArguments.getClass()) + getUsagePreamble());
        stream.println('\n' + getVersion());

        final Map<Boolean, List<ArgumentDefinition>> argMap = argumentDefinitions.stream()
                .collect(Collectors.partitioningBy(ArgumentDefinition::isOptional));
        printArgumentGroup(stream, argMap.get(false), "\n\nRequired Arguments:\n");
        printArgumentGroup(stream, argMap.get(true), "\nOptional Arguments:\n");
    }

    private void printArgumentGroup(PrintStream stream, List<ArgumentDefinition> args, String groupName) {
        if (args != null && !args.isEmpty()) {
            stream.println(groupName);
            args.stream().sorted(Comparator.comparing(ArgumentDefinition::getLongName))
                    .forEachOrdered(argumentDefinition -> printArgumentUsage(stream, argumentDefinition));
        }
    }

    private String generateUsageFromRequiredArguments(){
        return argumentDefinitions.stream()
                .filter(arg -> !arg.isOptional())
                .map(ArgumentDefinition::getPrettyNameAndType)
                .collect(Collectors.joining(" ", "Usage: " + callerArguments.getClass().getSimpleName(), " [options]"));
    }

    /**
     * Parse command-line arguments, and store values in callerArguments object passed to ctor.
     * @param messageStream Where to write error messages.
     * @param args          Command line tokens.
     * @return true if command line is valid and the program should run, false if help or version was requested
     * @throws UserException.CommandLineException if there is an invalid command line
     */
    @SuppressWarnings("unchecked")
    public boolean parseArguments(final PrintStream messageStream, final String[] args) {
        this.argv = args;
        try {
            final OptionParser optionParser = setupOptionParser(argumentDefinitions, positionalArguments != null);
            final OptionSet parsedArguments = parseArgumentsIntoOptionSet(args, optionParser);

            //check if special short circuiting arguments are set
            if (isSpecialFlagSet(parsedArguments, SpecialArgumentsCollection.HELP_FULLNAME)) {
                usage(messageStream);
                return false;
            } else if (isSpecialFlagSet(parsedArguments, SpecialArgumentsCollection.VERSION_FULLNAME)) {
                messageStream.println(getVersion());
                return false;
            }

            parsedArguments.asMap().keySet().stream()
                    .filter(optSpec -> parsedArguments.has(optSpec))
                    .forEach(optSpec -> {
                        //get the ArgumentDefinitions that corresponds to the name of the OptionParser Opt Spec
                        final ArgumentDefinition argDef = argumentMap.get(optSpec.options().get(0));
                        argDef.setArgument((List<String>) optSpec.values(parsedArguments));
                    });

            for (Object arg : parsedArguments.nonOptionArguments()) {
                setPositionalArgument((String) arg);
            }

            assertArgumentsAreValid();

            return true;
        }catch (UserException.CommandLineException e){
                throw new UserException.CommandLineException(e.getMessage(), getCommandLineAsInput(), e);
        }
    }

    /**
     * Parse the arguments from the command line into an {@link OptionSet}
     * If argument files are specified, this will be recursively invoked until all files are loaded.
     * @param args arguments from the command line
     * @param optionParser a parser pre configured with all of the valid OptionSpecs
     * @return an OptionSet containing a valid parsing of the argument values
     * @throws UserException.CommandLineException if the the arguments cannot be parsed
     */
    private static OptionSet parseArgumentsIntoOptionSet(String[] args, OptionParser optionParser){
        return parseArgumentsIntoOptionSet(args, optionParser, new HashSet<>());
    }

    private static OptionSet parseArgumentsIntoOptionSet(String[] args, OptionParser optionParser, Set<String> argumentsFilesLoadedAlready){
        final OptionSet parsedArguments;
        try {
            parsedArguments = optionParser.parse(args);
        } catch (final OptionException e) {
            throw new UserException.CommandLineException(e.getMessage(), String.join(" ", Arrays.asList(args)));
        }
        //Check for the special arguments file flag
        //if it's seen, read arguments from that file and recursively call parseArguments()
        if (parsedArguments.has(SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME)) {
            final List<String> argfiles = parsedArguments.valuesOf(SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME).stream()
                    .map(f -> (String)f)
                    .collect(Collectors.toList());

            final List<String> newargs = argfiles.stream()
                    .distinct()
                    .filter(file -> !argumentsFilesLoadedAlready.contains(file))
                    .flatMap(file -> loadArgumentsFile(file).stream())
                    .collect(Collectors.toList());

            argumentsFilesLoadedAlready.addAll(argfiles);

            if (!newargs.isEmpty()) {
                newargs.addAll(Arrays.asList(args));
                return parseArgumentsIntoOptionSet(newargs.toArray(new String[newargs.size()]), optionParser, argumentsFilesLoadedAlready);
            }
        }
        return parsedArguments;
    }

    /**
     * Setup a new {@link OptionParser} by configuring it with the given argumentDefinitions
     */
    private static OptionParser setupOptionParser(Collection<ArgumentDefinition> argumentDefinitions, boolean hasPositionalArguments) {
        final OptionParser parser = new OptionParser();

        final StrictBooleanConverter converter = new StrictBooleanConverter();
        for (ArgumentDefinition arg : argumentDefinitions){
            final OptionSpecBuilder bld = parser.acceptsAll(arg.getNames(), arg.getDoc());
            if (arg.isFlag()) {
                bld.withOptionalArg().withValuesConvertedBy(converter);
            } else {
                bld.withRequiredArg();
            }
        }

        if(hasPositionalArguments){
            parser.nonOptions();
        }
        return parser;
    }

    /**
     *  helper to deal with the case of special flags that are evaluated before the options are properly set
     */
    private static boolean isSpecialFlagSet(OptionSet parsedArguments, String flagName){
        if (parsedArguments.has(flagName)){
            final Object value = parsedArguments.valueOf(flagName);
            return  (value == null || !value.equals("false"));
        } else{
            return false;
        }
    }

    /**
     * After command line has been parsed, make sure that all required arguments have values, and that
     * lists with minimum # of elements have sufficient.
     *
     * @throws UserException.CommandLineException if arguments requirements are not satisfied.
     */
    private void assertArgumentsAreValid() {
        try {
            final Set<ArgumentDefinition> missingRequiredArguments = new HashSet<>();

            for (final ArgumentDefinition argumentDefinition : argumentDefinitions) {
                final String name = argumentDefinition.getLongName();
                final Set<ArgumentDefinition> mutextArguments = getConflictingMutallyExclusiveArguments(argumentDefinition, argumentMap);
                if (argumentDefinition.hasBeenSet() && !mutextArguments.isEmpty()) {
                    throw new UserException.ConflictingMutuallyExclusiveArguments("Argument '" + name +
                            "' cannot be used in conjunction with argument(s) " +
                            mutextArguments.stream().map(ArgumentDefinition::getLongName).collect(Collectors.joining(" ")), getCommandLineAsInput());
                }
                if (!argumentDefinition.isOptional() && !argumentDefinition.hasBeenSet() && mutextArguments.isEmpty()) {
                    missingRequiredArguments.add(argumentDefinition);
                }
            }

            if (positionalArguments != null) {
                @SuppressWarnings("rawtypes")
                final Collection c = (Collection) positionalArguments.get(positionalArgumentsParent);
                if (c.size() < minPositionalArguments) {
                    throw new UserException.MissingArgument(POSITIONAL_ARGUMENTS_NAME + " was missing. At least " + minPositionalArguments +
                            " positional arguments must be specified.", getCommandLineAsInput());
                }
            }
            if(!missingRequiredArguments.isEmpty()){
                throw new UserException.MissingArgument(composeMissingArgumentsMessage(missingRequiredArguments), getCommandLineAsInput());
            }
        } catch (final IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("Should never happen",e);
        }


    }

    private String composeMissingArgumentsMessage(Set<ArgumentDefinition> missingRequiredArguments) {
        final Set<ArgumentDefinition> alreadyOutput = new HashSet<>();

        //Sort the missing arguments so that single arguments are sorted before mutually exclusive groups, and then alphabetically within group
        final List<ArgumentDefinition> sortedArgs = missingRequiredArguments.stream()
                .sorted(Comparator.comparing((ArgumentDefinition arg) -> !arg.getMutuallyExclusive().isEmpty())
                        .thenComparing(ArgumentDefinition::getLongName))
                .collect(Collectors.toList());

        final List<String> messages = new ArrayList<>();

        for( ArgumentDefinition arg: sortedArgs) {
            if(!alreadyOutput.contains(arg)) {
                messages.add(arg.composeMissingArgumentMessage(argumentMap));
                alreadyOutput.addAll(arg.getMutuallyExclusiveArguments(argumentMap));
            }
        }

        return String.join("\n", messages);
    }

    /**
     * Get a Set of other {@link ArgumentDefinition}s that have been set that are mutually exclusive to arg
     * @param arg any argument, set or unset
     * @param argumentMap a mapping from argument name to argument
     * @return a Set of other argument definitions that have been set, that are mutually exclusive with arg
     */
    private static Set<ArgumentDefinition> getConflictingMutallyExclusiveArguments(ArgumentDefinition arg, Map<String, ArgumentDefinition> argumentMap) {
        return arg.getMutuallyExclusive().stream()
                .sorted()
                .map(argumentMap::get)
                .filter(mutexArg -> mutexArg != null && mutexArg.hasBeenSet())
                .collect(Collectors.toSet());
    }

    @SuppressWarnings("unchecked")
    private void setPositionalArgument(final String stringValue) {
        if (positionalArguments == null) {
            throw new UserException.CommandLineException("Invalid argument '" + stringValue + "'.", getCommandLineAsInput());
        }
        final Object value = constructFromString(getUnderlyingType(positionalArguments), stringValue, POSITIONAL_ARGUMENTS_NAME);
        @SuppressWarnings("rawtypes")
        final Collection c;
        try {
            c = (Collection) positionalArguments.get(callerArguments);
        } catch (final IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException(e);
        }
        if (c.size() >= maxPositionalArguments) {  //we're checking if there is space to add another argument
            throw new UserException.CommandLineException("No more than " + maxPositionalArguments +
                    " positional arguments may be specified on the command line.", getCommandLineAsInput());
        }
        c.add(value);
    }

    /**
     * Read an argument file and return a list of the args contained in it
     * A line that starts with {@link #COMMENT}  is ignored.
     *
     * @param argumentsFile a text file containing args
     * @return false if a fatal error occurred
     */
    private static List<String> loadArgumentsFile(final String argumentsFile) {
        final List<String> args = new ArrayList<>();
            try (BufferedReader reader = new BufferedReader(new FileReader(argumentsFile))){
                String line;
                while ((line = reader.readLine()) != null) {
                    if (!line.startsWith(COMMENT) && !line.trim().isEmpty()) {
                        args.addAll(Arrays.asList(StringUtils.split(line)));
                    }
                }
            } catch (final IOException e) {
                throw new UserException("I/O error loading arguments file:" + argumentsFile, e);
            }
        return args;
    }

    private void printArgumentUsage(final PrintStream stream, final ArgumentDefinition argumentDefinition) {
        stream.print(argumentDefinition.getArgumentParamUsage(argumentMap));
        stream.println();
    }

    /**
     * Generates the option help string for a {@code boolean} or {@link Boolean} typed argument.
     * @return never {@code null}.
     */
    private static String getBooleanOptions() {
        return String.format("%s%s, %s%s", ENUM_OPTION_DOC_PREFIX, Boolean.TRUE, Boolean.FALSE, ENUM_OPTION_DOC_SUFFIX);
    }

    /**
     * Composes the help string on the possible options an {@link Enum} typed argument can take.
     *
     * @param clazz target enum class. Assumed no to be {@code null}.
     * @param <T> enum class type.
     * @param <U> ClpEnum implementing version of <code>&lt;T&gt</code>;.
     * @throws GATKException if {@code &lt;T&gt;} has no constants.
     * @return never {@code null}.
     */
    private static <T extends Enum<T>, U extends Enum<U> & ClpEnum> String getEnumOptions(final Class<T> clazz) {
        // We assume that clazz is guaranteed to be a Class<? extends Enum>, thus
        // getEnumConstants() won't ever return a null.
        final T[] enumConstants = clazz.getEnumConstants();
        if (enumConstants.length == 0) {
            throw new GATKException(String.format("Bad argument enum type '%s' with no options", clazz.getName()));
        }

        if (ClpEnum.class.isAssignableFrom(clazz)) {
            @SuppressWarnings("unchecked")
            final U[] clpEnumCastedConstants = (U[]) enumConstants;
            return getEnumOptionsWithDescription(clpEnumCastedConstants);
        } else {
            return getEnumOptionsWithoutDescription(enumConstants);
        }
    }

    /**
     * Composes the help string for enum options that do not provide additional help documentation.
     * @param enumConstants the enum constants. Assumed non-null.
     * @param <T> the enum type.
     * @return never {@code null}.
     */
    private static <T extends Enum<T>> String getEnumOptionsWithoutDescription(final T[] enumConstants) {
        return Stream.of(enumConstants)
                .map(T::name)
                .collect(Collectors.joining(", ",ENUM_OPTION_DOC_PREFIX,ENUM_OPTION_DOC_SUFFIX));
    }

    /**
     * Composes the help string for enum options that provide additional documentation.
     * @param enumConstants the enum constants. Assumed non-null.
     * @param <T> the enum type.
     * @return never {@code null}.
     */
    private static <T extends Enum<T> & ClpEnum> String getEnumOptionsWithDescription(final T[] enumConstants) {
        final String optionsString = Stream.of(enumConstants)
                .map(c -> String.format("%s (%s)",c.name(),c.getHelpDoc()))
                .collect(Collectors.joining("\n"));
        return String.join("\n",ENUM_OPTION_DOC_PREFIX,optionsString,ENUM_OPTION_DOC_SUFFIX);
    }

    /**
     * Returns the help string with details about valid options for the given argument class.
     *
     * <p>
     *     Currently this only make sense with {@link Boolean} and {@link Enum}. Any other class
     *     will result in an empty string.
     * </p>
     *
     * @param clazz the target argument's class.
     * @return never {@code null}.
     */
    @SuppressWarnings({"unchecked","rawtypes"})
    static String getOptions(final Class<?> clazz) {
        if (clazz == Boolean.class) {
            return getBooleanOptions();
        } else if (clazz.isEnum()) {
            final Class<? extends Enum> enumClass = (Class<? extends Enum>)clazz;
            return getEnumOptions(enumClass);
        } else {
            return "";
        }
    }

    private void handleArgumentAnnotation(final Field field, final Object parent) {
        try {
            field.setAccessible(true);
            final Argument argumentAnnotation = field.getAnnotation(Argument.class);
            final boolean isCollection = isCollectionField(field);
            if (isCollection) {
                field.setAccessible(true);
                if (field.get(parent) == null) {
                    createCollection(field, parent, "@Argument");
                }
            }
            if (!canBeMadeFromString(getUnderlyingType(field))) {
                throw new GATKException.CommandLineParserInternalException("@Argument member \"" + field.getName() +
                        "\" must have a String constructor or be an enum");
            }

            final ArgumentDefinition argumentDefinition = new ArgumentDefinition(field, argumentAnnotation, parent);

            for (final String argumentName : argumentAnnotation.mutex()) {
                final ArgumentDefinition mutuallyExclusiveArgument = argumentMap.get(argumentName);
                if (mutuallyExclusiveArgument != null) {
                    mutuallyExclusiveArgument.getMutuallyExclusive().add(field.getName());
                }
            }
            if (inArgumentMap(argumentDefinition)) {
                throw new GATKException.CommandLineParserInternalException(argumentDefinition.getNames() + " has already been used.");
            } else {
                putInArgumentMap(argumentDefinition);
                argumentDefinitions.add(argumentDefinition);
            }
        } catch (final IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("We should not have reached here because we set accessible to true", e);
        }
    }

    private void handlePositionalArgumentAnnotation(final Field field, Object parent) {
        if (positionalArguments != null) {
            throw new GATKException.CommandLineParserInternalException
                    ("@PositionalArguments cannot be used more than once in an argument class.");
        }
        field.setAccessible(true);
        positionalArguments = field;
        positionalArgumentsParent = parent;
        if (!isCollectionField(field)) {
            throw new GATKException.CommandLineParserInternalException("@PositionalArguments must be applied to a Collection");
        }

        if (!canBeMadeFromString(getUnderlyingType(field))) {
            throw new GATKException.CommandLineParserInternalException("@PositionalParameters member " + field.getName() +
                    "does not have a String ctor");
        }

        final PositionalArguments positionalArgumentsAnnotation = field.getAnnotation(PositionalArguments.class);
        minPositionalArguments = positionalArgumentsAnnotation.minElements();
        maxPositionalArguments = positionalArgumentsAnnotation.maxElements();
        if (minPositionalArguments > maxPositionalArguments) {
            throw new GATKException.CommandLineParserInternalException("In @PositionalArguments, minElements cannot be > maxElements");
        }
        try {
            field.setAccessible(true);
            if (field.get(parent) == null) {
                createCollection(field, parent, "@PositionalParameters");
            }
        } catch (final IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("We should not have reached here because we set accessible to true", e);

        }
    }


    public static boolean isCollectionField(final Field field) {
        try {
            field.getType().asSubclass(Collection.class);
            return true;
        } catch (final ClassCastException e) {
            return false;
        }
    }

    private static void createCollection(final Field field, final Object callerArguments, final String annotationType)
            throws IllegalAccessException {
        try {
            field.set(callerArguments, field.getType().newInstance());
        } catch (final Exception ex) {
            try {
                field.set(callerArguments, new ArrayList<>());
            } catch (final IllegalArgumentException e) {
                throw new GATKException.CommandLineParserInternalException("In collection " + annotationType +
                        " member " + field.getName() +
                        " cannot be constructed or auto-initialized with ArrayList, so collection must be initialized explicitly.");
            }

        }

    }

    /**
     * Returns the type that each instance of the argument needs to be converted to. In
     * the case of primitive fields it will return the wrapper type so that String
     * constructors can be found.
     */
    static Class<?> getUnderlyingType(final Field field) {
        if (isCollectionField(field)) {
            final ParameterizedType clazz = (ParameterizedType) (field.getGenericType());
            final Type[] genericTypes = clazz.getActualTypeArguments();
            if (genericTypes.length != 1) {
                throw new GATKException.CommandLineParserInternalException("Strange collection type for field " +
                        field.getName());
            }

            // If the Collection's parametrized type is itself parametrized (eg., List<Foo<Bar>>),
            // return the raw type of the outer parameter (Foo.class, in this example) to avoid a
            // ClassCastException. Otherwise, return the Collection's type parameter directly as a Class.
            return (Class<?>) (genericTypes[0] instanceof ParameterizedType ?
                               ((ParameterizedType)genericTypes[0]).getRawType() :
                               genericTypes[0]);

        } else {
            final Class<?> type = field.getType();
            if (type == Byte.TYPE) { return Byte.class; }
            if (type == Short.TYPE) { return Short.class; }
            if (type == Integer.TYPE) { return Integer.class; }
            if (type == Long.TYPE) { return Long.class; }
            if (type == Float.TYPE) { return Float.class; }
            if (type == Double.TYPE) { return Double.class; }
            if (type == Boolean.TYPE) { return Boolean.class; }

            return type;
        }
    }

    // True if clazz is an enum, or if it has a ctor that takes a single String argument.
    private static boolean canBeMadeFromString(final Class<?> clazz) {
        if (clazz.isEnum()) {
            return true;
        }
        try {
            // Need to use getDeclaredConstructor() instead of getConstructor() in case the constructor
            // is non-public
            clazz.getDeclaredConstructor(String.class);
            return true;
        } catch (final NoSuchMethodException e) {
            return false;
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    static Object constructFromString(final Class clazz, final String s, final String argumentName) {
        try {
            if (clazz.isEnum()) {
                try {
                    return Enum.valueOf(clazz, s);
                } catch (final IllegalArgumentException e) {
                    throw new UserException.BadArgumentValue(argumentName, s, '\'' + s + "' is not a valid value for " +
                            clazz.getSimpleName() + ". "+ getEnumOptions(clazz) );
                }
            }
            // Need to use getDeclaredConstructor() instead of getConstructor() in case the constructor
            // is non-public. Set it to be accessible if it isn't already.
            final Constructor<?> ctor = clazz.getDeclaredConstructor(String.class);
            ctor.setAccessible(true);
            return ctor.newInstance(s);
        } catch (final NoSuchMethodException e) {
            // Shouldn't happen because we've checked for presence of ctor
            throw new GATKException.ShouldNeverReachHereException("Cannot find string ctor for " + clazz.getName(), e);
        } catch (final InstantiationException e) {
            throw new GATKException.CommandLineParserInternalException("Abstract class '" + clazz.getSimpleName() +
                    "'cannot be used for an argument value type.", e);
        } catch (final IllegalAccessException e) {
            throw new GATKException.CommandLineParserInternalException("String constructor for argument value type '" + clazz.getSimpleName() +
                    "' must be public.", e);
        } catch (final InvocationTargetException e) {
            throw new UserException.BadArgumentValue(argumentName, s, "Problem constructing " + clazz.getSimpleName() +
                    " from the string '" + s + "'.");
        }
    }

    public interface ClpEnum {
        String getHelpDoc();
    }

    /**
     * The commandline used to run this program, including any default args that
     * weren't necessarily specified. This is used for logging and debugging.
     * <p/>
     * NOTE: {@link #parseArguments(PrintStream, String[])} must be called before
     * calling this method.
     *
     * @return The commandline, or null if {@link #parseArguments(PrintStream, String[])}
     * hasn't yet been called, or didn't complete successfully.
     */
    @SuppressWarnings("unchecked")
    public String getFullySpecifiedCommandLine() {
        final String toolName = callerArguments.getClass().getName();
        final StringBuilder commandLineString = new StringBuilder();

        final List<Object> positionalArgs;
        if( positionalArguments != null) {
            try {
                positionalArguments.setAccessible(true);
                positionalArgs = (List<Object>) positionalArguments.get(positionalArgumentsParent);
            } catch (IllegalAccessException e) {
                throw new GATKException.ShouldNeverReachHereException("Should never reach here because we setAccessible(true)", e);
            }
            for (final Object posArg : positionalArgs) {
                commandLineString.append(' ').append(posArg);
            }
        }

        //first, append args that were explicitly set
        commandLineString.append(argumentDefinitions.stream()
                .filter(ArgumentDefinition::hasBeenSet)
                .map(ArgumentDefinition::toCommandLineString)
                .collect(Collectors.joining(" "," ","  ")));

        //next, append args that weren't explicitly set, but have a default value
        commandLineString.append(argumentDefinitions.stream()
                .filter(argumentDefinition -> !argumentDefinition.hasBeenSet() && !argumentDefinition.getDefaultValue().equals(NULL_STRING))
                .map(ArgumentDefinition::toCommandLineString)
                .collect(Collectors.joining(" ")));

        return toolName + ' ' + commandLineString;
    }

    /**
     *
     * This must be called after {@link CommandLineParser#argv} is initialized by {@link CommandLineParser#parseArguments(PrintStream, String[])} or it throw {@link GATKException}
     * @return what was actually entered as the command line joined into a String
     */
    private String getCommandLineAsInput(){
        if( argv == null){
            throw new GATKException("No commandline was specified for parsing yet");
        } else{
            return String.join(" ", argv);
        }
    }

    /**
     * Locates and returns the VALUES of all Argument-annotated fields of a specified type in a given object,
     * pairing each field value with its corresponding Field object.
     *
     * Must be called AFTER argument parsing and value injection into argumentSource is complete (otherwise there
     * will be no values to gather!). As a result, this is implemented as a static utility method into which
     * the fully-initialized tool instance must be passed.
     *
     * Locates Argument-annotated fields of the target type, subtypes of the target type, and Collections of
     * the target type or one of its subtypes. Unpacks Collection fields, returning a separate Pair for each
     * value in each Collection.
     *
     * Searches argumentSource itself, as well as ancestor classes, and also recurses into any ArgumentCollections
     * found.
     *
     * Will return Pairs containing a null second element for fields having no value, including empty Collection fields
     * (these represent arguments of the target type that were not specified on the command line and so never initialized).
     *
     * @param type Target type. Search for Argument-annotated fields that are either of this type, subtypes of this type, or Collections of this type or one of its subtypes.
     * @param argumentSource Object whose fields to search. Must have already undergone argument parsing and argument value injection.
     * @param <T> Type parameter representing the type to search for and return
     * @return A List of Pairs containing all Argument-annotated field values found of the target type. First element in each Pair
     *         is the Field object itself, and the second element is the actual value of the argument field. The second
     *         element will be null for uninitialized fields.
     */
    public static <T> List<Pair<Field, T>> gatherArgumentValuesOfType( final Class<T> type, final Object argumentSource ) {
        final List<Pair<Field, T>> argumentValues = new ArrayList<>();

        // Examine all fields in argumentSource (including superclasses)
        for ( Field field : getAllFields(argumentSource.getClass()) ) {
            field.setAccessible(true);

            try {
                // Consider only fields that have Argument annotations and are either of the target type,
                // subtypes of the target type, or Collections of the target type or one of its subtypes:
                if ( field.getAnnotation(Argument.class) != null && type.isAssignableFrom(getUnderlyingType(field)) ) {

                    if ( isCollectionField(field) ) {
                        // Collection arguments are guaranteed by the parsing system to be non-null (at worst, empty)
                        final Collection<?> argumentContainer = (Collection<?>)field.get(argumentSource);

                        // Emit a Pair with an explicit null value for empty Collection arguments
                        if ( argumentContainer.isEmpty() ) {
                            argumentValues.add(Pair.of(field, null));
                        }
                        // Unpack non-empty Collections of the target type into individual values,
                        // each paired with the same Field object.
                        else {
                            for ( Object argumentValue : argumentContainer ) {
                                argumentValues.add(Pair.of(field, type.cast(argumentValue)));
                            }
                        }
                    }
                    else {
                        // Add values for non-Collection arguments of the target type directly
                        argumentValues.add(Pair.of(field, type.cast(field.get(argumentSource))));
                    }
                }
                else if ( field.getAnnotation(ArgumentCollection.class) != null ) {
                    // Recurse into ArgumentCollections for more potential matches.
                    argumentValues.addAll(gatherArgumentValuesOfType(type, field.get(argumentSource)));
                }
            }
            catch ( IllegalAccessException e ) {
                throw new GATKException.ShouldNeverReachHereException("field access failed after setAccessible(true)");
            }
        }

        return argumentValues;
    }
}
