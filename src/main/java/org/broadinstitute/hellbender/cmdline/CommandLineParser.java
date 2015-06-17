package org.broadinstitute.hellbender.cmdline;

import htsjdk.samtools.util.StringUtil;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import joptsimple.OptionSpecBuilder;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.*;
import java.util.*;
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
    private static final int ARGUMENT_COLUMN_WIDTH = 30;
    private static final int DESCRIPTION_COLUMN_WIDTH = 90;

    private static final String ENUM_OPTION_DOC_PREFIX = "Possible values: {";
    private static final String ENUM_OPTION_DOC_SUFFIX = "} ";

    // Use these if no @Usage annotation
    private static final String defaultUsagePreamble = "Usage: program [arguments...]\n";
    private static final String defaultUsagePreambleWithPositionalArguments =
            "Usage: program [arguments...] [positional-arguments...]\n";
    private static final String NULL_STRING = "null";
    public static final String COMMENT = "#";
    public static final String POSITIONAL_ARGUMENTS_NAME = "Positional Argument";


    private Set<String> argumentsFilesLoadedAlready = new HashSet<>();

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
            usagePreamble += programProperties.usage();
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

    private List<ArgumentDefinition> createArgumentDefinitions(final Object callerArguments) {
        for (final Field field : getAllFields(callerArguments.getClass())) {
            if (field.getAnnotation(Argument.class) != null && field.getAnnotation(ArgumentCollection.class) != null){
                throw new GATKException.CommandLineParserInternalException("An Argument cannot be an argument collection: "
                        +field.getName() + " in " + callerArguments.toString() + " is annotated as both.");
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
        return null;
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
    public void usage(final PrintStream stream, final boolean printCommon) {
        stream.print(getStandardUsagePreamble(callerArguments.getClass()) + getUsagePreamble());
        stream.println("\n" + getVersion());
        stream.println("\n\nArguments:\n");

            argumentDefinitions.stream()
                    .filter(argumentDefinition -> printCommon || !argumentDefinition.isCommon)
                    .forEach(argumentDefinition -> printArgumentUsage(stream, argumentDefinition));
    }

    /**
     * Parse command-line arguments, and store values in callerArguments object passed to ctor.
     * @param messageStream Where to write error messages.
     * @param args          Command line tokens.
     * @return true if command line is valid and the program should run, false if help or version was requested
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CommandLineException if there is an invalid command line
     */
    @SuppressWarnings("unchecked")
    public boolean parseArguments(final PrintStream messageStream, final String[] args) {
        this.argv = args;

        OptionParser parser = new OptionParser();

        for (ArgumentDefinition arg : argumentDefinitions){
            OptionSpecBuilder bld = parser.acceptsAll(arg.getNames(), arg.doc);
            if (arg.isFlag()) {
                bld.withOptionalArg().withValuesConvertedBy(new StrictBooleanConverter());
            } else {
                bld.withRequiredArg();
            }
        }
        if(positionalArguments != null){
            parser.nonOptions();
        }

        OptionSet parsedArguments;
        try {
            parsedArguments = parser.parse(args);
        } catch (final joptsimple.OptionException e) {
            throw new UserException.CommandLineException(e.getMessage());
        }
        //Check for the special arguments file flag
        //if it's seen, read arguments from that file and recursively call parseArguments()
        if (parsedArguments.has(SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME)) {
            List<String> argfiles = parsedArguments.valuesOf(SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME).stream()
                    .map(f -> (String)f)
                    .collect(Collectors.toList());

            List<String> newargs = argfiles.stream()
                    .distinct()
                    .filter(file -> !argumentsFilesLoadedAlready.contains(file))
                    .flatMap(file -> loadArgumentsFile(file).stream())
                    .collect(Collectors.toList());
            argumentsFilesLoadedAlready.addAll(argfiles);

            if (!newargs.isEmpty()) {
                newargs.addAll(Arrays.asList(args));
                return parseArguments(messageStream, newargs.toArray(new String[newargs.size()]));
            }
        }

        //check if special short circuiting arguments are set
        if (isSpecialFlagSet(parsedArguments, SpecialArgumentsCollection.HELP_FULLNAME)) {
            usage(messageStream, true);
            return false;
        } else if (isSpecialFlagSet(parsedArguments, SpecialArgumentsCollection.VERSION_FULLNAME)) {
            messageStream.println(getVersion());
            return false;
        }

        for (OptionSpec<?> optSpec : parsedArguments.asMap().keySet()) {
            if (parsedArguments.has(optSpec)) {
                ArgumentDefinition argDef = argumentMap.get(optSpec.options().get(0));
                setArgument(argDef, (List<String>) optSpec.values(parsedArguments));
            }
        }

        for (Object arg : parsedArguments.nonOptionArguments()) {
            setPositionalArgument((String) arg);
        }

        assertArgumentsAreValid();

        return true;
    }

    /**
     *  helper to deal with the case of special flags that are evaluated before the options are properly set
     */
    private boolean isSpecialFlagSet(OptionSet parsedArguments, String flagName){
        if (parsedArguments.has(flagName)){
            Object value = parsedArguments.valueOf(flagName);
            return  (value == null || !((String)value).equals("false"));
        } else{
            return false;
        }

    }

    /**
     * After command line has been parsed, make sure that all required arguments have values, and that
     * lists with minimum # of elements have sufficient.
     *
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CommandLineException if arguments requirements are not satisfied.
     */
    private void assertArgumentsAreValid() {
        try {
            for (final ArgumentDefinition argumentDefinition : argumentDefinitions) {
                final String fullName = argumentDefinition.getLongName();
                final StringBuilder mutextArgumentNames = new StringBuilder();
                for (final String mutexArgument : argumentDefinition.mutuallyExclusive) {
                    final ArgumentDefinition mutextArgumentDef = argumentMap.get(mutexArgument);
                    if (mutextArgumentDef != null && mutextArgumentDef.hasBeenSet) {
                        mutextArgumentNames.append(" ").append(mutextArgumentDef.getLongName());
                    }
                }
                if (argumentDefinition.hasBeenSet && mutextArgumentNames.length() > 0) {
                    throw new UserException.CommandLineException("Argument '" + fullName +
                            "' cannot be used in conjunction with argument(s)" +
                            mutextArgumentNames.toString());
                }
                if (argumentDefinition.isCollection && !argumentDefinition.optional) {
                    @SuppressWarnings("rawtypes")
                    final Collection c = (Collection) argumentDefinition.getFieldValue();
                    if (c.isEmpty()) {
                        throw new UserException.MissingArgument(fullName, "Argument '" + fullName + "' must be specified at least once.");
                    }
                } else if (!argumentDefinition.optional && !argumentDefinition.hasBeenSet && mutextArgumentNames.length() == 0) {
                    throw new UserException.MissingArgument(fullName, "Argument '" + fullName + "' is required" +
                            (argumentDefinition.mutuallyExclusive.isEmpty() ? "." : " unless any of " + argumentDefinition.mutuallyExclusive +
                                    " are specified."));
                }

            }
            if (positionalArguments != null) {
                @SuppressWarnings("rawtypes")
                final Collection c = (Collection) positionalArguments.get(positionalArgumentsParent);
                if (c.size() < minPositionalArguments) {
                    throw new UserException.MissingArgument(POSITIONAL_ARGUMENTS_NAME,"At least " + minPositionalArguments +
                            " positional arguments must be specified.");
                }
            }
        } catch (final IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("Should never happen",e);
        }


    }

    @SuppressWarnings("unchecked")
    private void setPositionalArgument(final String stringValue) {
        if (positionalArguments == null) {
            throw new UserException.CommandLineException("Invalid argument '" + stringValue + "'.");
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
                    " positional arguments may be specified on the command line.");
        }
        c.add(value);
    }

    @SuppressWarnings("unchecked")
    private void setArgument(ArgumentDefinition argumentDefinition, final List<String> values) {
        //special treatment for flags
        if (argumentDefinition.isFlag() && values.isEmpty()){
            argumentDefinition.hasBeenSet = true;
            argumentDefinition.setFieldValue(true);
            return;
        }

        if (!argumentDefinition.isCollection && (argumentDefinition.hasBeenSet || values.size() > 1)) {
                throw new UserException.CommandLineException("Argument '" + argumentDefinition.getNames() + "' cannot be specified more than once.");
        }

        for (String stringValue: values) {
            final Object value;
            if (stringValue.equals(NULL_STRING)) {
                //"null" is a special value that allows the user to override any default
                //value set for this arg
                if (argumentDefinition.optional) {
                    value = null;
                } else {
                    throw new UserException.CommandLineException("Non \"null\" value must be provided for '" + argumentDefinition.getNames() + "'.");
                }
            } else {
                value = constructFromString(getUnderlyingType(argumentDefinition.field), stringValue, argumentDefinition.getLongName());
            }

            if (argumentDefinition.isCollection) {
                @SuppressWarnings("rawtypes")
                final Collection c = (Collection) argumentDefinition.getFieldValue();
                if (value == null) {
                    //user specified this arg=null which is interpreted as empty list
                    c.clear();
                } else {
                    c.add(value);
                }
                argumentDefinition.hasBeenSet = true;
            } else {
                argumentDefinition.setFieldValue(value);
                argumentDefinition.hasBeenSet = true;
            }
        }
    }

    /**
     * Read an argument file and return a list of the args contained in it
     * A line that starts with {@link #COMMENT}  is ignored.
     *
     * @param argumentsFile a text file containing args
     * @return false if a fatal error occurred
     */
    private List<String> loadArgumentsFile(final String argumentsFile) {
        List<String> args = new ArrayList<>();
            try (BufferedReader reader = new BufferedReader(new FileReader(argumentsFile))){
                String line;
                while ((line = reader.readLine()) != null) {
                    if (!line.startsWith(COMMENT) && line.trim().length() != 0) {
                        args.addAll(Arrays.asList(StringUtils.split(line)));
                    }
                }
            } catch (final IOException e) {
                throw new UserException("I/O error loading arguments file:" + argumentsFile, e);
            }
        return args;
    }

    private void printArgumentUsage(final PrintStream stream, final ArgumentDefinition argumentDefinition) {
        printArgumentParamUsage(stream, argumentDefinition.getLongName(), argumentDefinition.shortName,
                getUnderlyingType(argumentDefinition.field).getSimpleName(),
                makeArgumentDescription(argumentDefinition));
    }


    private void printArgumentParamUsage(final PrintStream stream, final String name, final String shortName,
                                       final String type, final String argumentDescription) {
        String argumentLabel = name;
        if (type != null) argumentLabel = "--"+ argumentLabel;

        if (shortName.length() > 0) {
            argumentLabel+=",-" + shortName;
        }
        argumentLabel += ":" + type;
        stream.print(argumentLabel);

        int numSpaces = ARGUMENT_COLUMN_WIDTH - argumentLabel.length();
        if (argumentLabel.length() > ARGUMENT_COLUMN_WIDTH) {
            stream.println();
            numSpaces = ARGUMENT_COLUMN_WIDTH;
        }
        printSpaces(stream, numSpaces);
        final String wrappedDescription = StringUtil.wordWrap(argumentDescription, DESCRIPTION_COLUMN_WIDTH);
        final String[] descriptionLines = wrappedDescription.split("\n");
        for (int i = 0; i < descriptionLines.length; ++i) {
            if (i > 0) {
                printSpaces(stream, ARGUMENT_COLUMN_WIDTH);
            }
            stream.println(descriptionLines[i]);
        }
        stream.println();
    }

    private String makeArgumentDescription(final ArgumentDefinition argumentDefinition) {
        final StringBuilder sb = new StringBuilder();
        if (argumentDefinition.doc.length() > 0) {
            sb.append(argumentDefinition.doc);
            sb.append("  ");
        }
        if (argumentDefinition.isCollection) {
            if (argumentDefinition.optional) {
                sb.append("This argument may be specified 0 or more times. ");
            } else {
                sb.append("This argument must be specified at least once. ");
            }
        }
        if (argumentDefinition.optional) {
            sb.append("Default value: ");
            sb.append(argumentDefinition.defaultValue);
            sb.append(". ");
        } else {
            sb.append("Required. ");
        }
        sb.append(getOptions(getUnderlyingType(argumentDefinition.field)));
        if (!argumentDefinition.mutuallyExclusive.isEmpty()) {
            sb.append(" Cannot be used in conjuction with argument(s)");
            for (final String argument : argumentDefinition.mutuallyExclusive) {
                final ArgumentDefinition mutextArgumentDefinition = argumentMap.get(argument);

                if (mutextArgumentDefinition == null) {
                    throw new GATKException("Invalid argument definition in source code.  " + argument +
                                                  " doesn't match any known argument.");
                }

                sb.append(" ").append(mutextArgumentDefinition.fieldName);
                if (mutextArgumentDefinition.shortName.length() > 0) {
                    sb.append(" (").append(mutextArgumentDefinition.shortName).append(")");
                }
            }
        }
        return sb.toString();
    }

    /**
     * Generates the option help string for a {@code boolean} or {@link Boolean} typed argument.
     * @return never {@code null}.
     */
    private String getBooleanOptions() {
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
    private <T extends Enum<T>,U extends Enum<U> & ClpEnum> String getEnumOptions(final Class<T> clazz) {
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
    private <T extends Enum<T>> String getEnumOptionsWithoutDescription(final T[] enumConstants) {
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
    private <T extends Enum<T> & ClpEnum> String getEnumOptionsWithDescription(final T[] enumConstants) {
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
    private String getOptions(final Class<?> clazz) {
        if (clazz == Boolean.class) {
            return getBooleanOptions();
        } else if (clazz.isEnum()) {
            final Class<? extends Enum> enumClass = (Class<? extends Enum>)clazz;
            return getEnumOptions(enumClass);
        } else {
            return "";
        }
    }

    private void printSpaces(final PrintStream stream, final int numSpaces) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numSpaces; ++i) {
            sb.append(" ");
        }
        stream.print(sb);
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
                throw new GATKException.CommandLineParserInternalException("@Argument member " + field.getName() +
                        " must have a String ctor or be an enum");
            }

            final ArgumentDefinition argumentDefinition = new ArgumentDefinition(field, argumentAnnotation, parent);

            for (final String argument : argumentAnnotation.mutex()) {
                final ArgumentDefinition mutextArgumentDef = argumentMap.get(argument);
                if (mutextArgumentDef != null) {
                    mutextArgumentDef.mutuallyExclusive.add(field.getName());
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

    private void createCollection(final Field field, final Object callerArguments, final String annotationType)
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
    private static Class<?> getUnderlyingType(final Field field) {
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
            if (type == Byte.TYPE) return Byte.class;
            if (type == Short.TYPE) return Short.class;
            if (type == Integer.TYPE) return Integer.class;
            if (type == Long.TYPE) return Long.class;
            if (type == Float.TYPE) return Float.class;
            if (type == Double.TYPE) return Double.class;
            if (type == Boolean.TYPE) return Boolean.class;

            return type;
        }
    }

    // True if clazz is an enum, or if it has a ctor that takes a single String argument.
    private boolean canBeMadeFromString(final Class<?> clazz) {
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
    private Object constructFromString(final Class clazz, final String s, final String argumentName) {
        try {
            if (clazz.isEnum()) {
                try {
                    return Enum.valueOf(clazz, s);
                } catch (final IllegalArgumentException e) {
                    throw new UserException.BadArgumentValue(argumentName, s, "'" + s + "' is not a valid value for " +
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

    protected static class ArgumentDefinition {
        final Field field;
        final String fieldName;
        final String fullName;
        final String shortName;
        final String doc;
        final boolean optional;
        final boolean isCollection;
        final String defaultValue;
        final boolean isCommon;
        boolean hasBeenSet = false;
        final Set<String> mutuallyExclusive;
        final Object parent;
        final boolean isSpecial;

        public ArgumentDefinition(final Field field, final Argument annotation, final Object parent){
            this.field = field;
            this.fieldName = field.getName();
            this.parent = parent;
            this.fullName = annotation.fullName();
            this.shortName = annotation.shortName();
            this.doc = annotation.doc();
            this.isCollection = isCollectionField(field);

            this.isCommon = annotation.common();
            this.isSpecial = annotation.special();

            this.mutuallyExclusive = new HashSet<>(Arrays.asList(annotation.mutex()));

            Object tmpDefault = getFieldValue();
            if (tmpDefault != null) {
                if (isCollection && ((Collection) tmpDefault).isEmpty()) {
                    //treat empty collections the same as uninitialized primitive types
                    this.defaultValue = NULL_STRING;
                } else {
                    //this is an initialized primitive type or a non-empty collection
                    this.defaultValue = tmpDefault.toString();
                }
            } else {
                this.defaultValue = NULL_STRING;
            }

            //null collections have been initialized by createCollection which is called in handleArgumentAnnotation
            //this is optional if it's specified as being optional or if there is a default value specified
            this.optional = annotation.optional() || ! this.defaultValue.equals(NULL_STRING);
        }


        public Object getFieldValue(){
            try {
                field.setAccessible(true);
                return field.get(parent);
            } catch (IllegalAccessException e) {
                throw new GATKException.ShouldNeverReachHereException("This shouldn't happen since we setAccessible(true).", e);
            }
        }

        public void setFieldValue(final Object value){
            try {
                field.setAccessible(true);
                field.set(parent, value);
            } catch (IllegalAccessException e) {
                throw new GATKException.ShouldNeverReachHereException("BUG: couldn't set field value. For "
                        + fieldName +" in " + parent.toString() + " with value " + value.toString()
                        + " This shouldn't happen since we setAccessible(true)", e);
            }
        }

        public boolean isFlag(){
            return field.getType().equals(boolean.class) || field.getType().equals(Boolean.class);
        }

        public List<String> getNames(){
            List<String> names = new ArrayList<>();
            if (!shortName.isEmpty()){
                names.add(shortName);
            }
            if (!fullName.isEmpty()){
                names.add(fullName);
            } else {
                names.add(fieldName);
            }
            return names;
        }

        public String getLongName(){
            return !fullName.isEmpty() ? fullName : fieldName;
        }

        /**
         * Helper for pretty printing this option.
         * @param value A value this argument was given
         * @return a string
         *
         */
        private String prettyNameValue(Object value) {
            if(value != null){
                return String.format("--%s %s", getLongName(), value);
            }
            return "";
        }

        /**
         * @return A string representation of this argument and it's value(s) which would be valid if copied and pasted
         * back as a command line argument
         */
        public String toCommandLineString(){
            Object value = getFieldValue();
            if (this.isCollection){
                Collection<?> collect = (Collection<?>)value;
                return collect.stream()
                        .map(this::prettyNameValue)
                        .collect(Collectors.joining(" "));

            } else {
                return prettyNameValue(value);
            }
        }

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
    public String getCommandLine() {
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
                commandLineString.append(" ").append(posArg.toString());
            }
        }

        //first, append args that were explicitly set
        commandLineString.append(argumentDefinitions.stream()
                .filter(argumentDefinition -> argumentDefinition.hasBeenSet)
                .map(ArgumentDefinition::toCommandLineString)
                .collect(Collectors.joining(" "," ","  ")));

        //next, append args that weren't explicitly set, but have a default value
        commandLineString.append(argumentDefinitions.stream()
                .filter(argumentDefinition -> !argumentDefinition.hasBeenSet && !argumentDefinition.defaultValue.equals(NULL_STRING))
                .map(ArgumentDefinition::toCommandLineString)
                .collect(Collectors.joining(" ")));

        return toolName + " " + commandLineString.toString();
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
        List<Pair<Field, T>> argumentValues = new ArrayList<>();

        // Examine all fields in argumentSource (including superclasses)
        for ( Field field : getAllFields(argumentSource.getClass()) ) {
            field.setAccessible(true);

            try {
                // Consider only fields that have Argument annotations and are either of the target type,
                // subtypes of the target type, or Collections of the target type or one of its subtypes:
                if ( field.getAnnotation(Argument.class) != null && type.isAssignableFrom(getUnderlyingType(field)) ) {

                    if ( isCollectionField(field) ) {
                        // Collection arguments are guaranteed by the parsing system to be non-null (at worst, empty)
                        Collection<?> argumentContainer = (Collection<?>)field.get(argumentSource);

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
