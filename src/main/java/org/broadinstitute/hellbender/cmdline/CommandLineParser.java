/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.hellbender.cmdline;

import htsjdk.samtools.util.StringUtil;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import joptsimple.OptionSpecBuilder;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.*;
import java.util.*;
import java.util.stream.Collectors;

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
public class CommandLineParser {
    // For formatting argument section of usage message.
    private static final int ARGUMENT_COLUMN_WIDTH = 30;
    private static final int DESCRIPTION_COLUMN_WIDTH = 90;

    private static final Boolean[] TRUE_FALSE_VALUES = {Boolean.TRUE, Boolean.FALSE};

    // Use these if no @Usage annotation
    private static final String defaultUsagePreamble = "Usage: program [arguments...]\n";
    private static final String defaultUsagePreambleWithPositionalArguments =
            "Usage: program [arguments...] [positional-arguments...]\n";
    private static final String NULL_STRING = "null";
    public static final String COMMENT = "#";


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

    // For printing error messages when parsing command line.
    private PrintStream messageStream;

    // In case implementation wants to get at arg for some reason.
    private String[] argv;

    private String programVersion = null;

    // The command line used to launch this program, including non-null default arguments that
    // weren't explicitly specified. This is used for logging and debugging.
    private String commandLine = "";


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

        if (null != this.programVersion && 0 < this.programVersion.length()) {
            usagePreamble += "Version: " + getVersion() + "\n";
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
        return this.callerArguments.getClass().getPackage().getImplementationVersion();
    }

    /**
     * Print a usage message based on the arguments object passed to the ctor.
     *
     * @param stream Where to write the usage message.
     */
    public void usage(final PrintStream stream, final boolean printCommon) {
        stream.print(getStandardUsagePreamble(callerArguments.getClass()) + getUsagePreamble());
        stream.println("\nVersion: " + getVersion());
        stream.println("\n\nArguments:\n");

            argumentDefinitions.stream()
                    .filter(argumentDefinition -> printCommon || !argumentDefinition.isCommon)
                    .forEach(argumentDefinition -> printArgumentUsage(stream, argumentDefinition));
    }


    /**
     * Parse command-line arguments, and store values in callerArguments object passed to ctor.
     *
     * @param messageStream Where to write error messages.
     * @param args          Command line tokens.
     * @return true if command line is valid and the program should run, false if only help was requested.
     */
    @SuppressWarnings("unchecked")
    public boolean parseArguments(final PrintStream messageStream, final String[] args) {
        this.argv = args;
        this.messageStream = messageStream;

        OptionParser parser = new OptionParser();
        for (ArgumentDefinition arg : argumentDefinitions){
            OptionSpecBuilder bld = parser.acceptsAll(arg.getNames(), arg.doc);
            if (arg.isFlag()) {
                bld.withOptionalArg();
            } else {
                bld.withRequiredArg();
            }
        }
        if(positionalArguments != null){
            parser.nonOptions();
        }

        OptionSet parsedArguments = null;
        try {
            parsedArguments = parser.parse(args);
        } catch (final joptsimple.OptionException e) {
            usage(messageStream, true);
            return false;
        }
        //Check for the special arguments file flag
        //if it's seen, read arguments from that file and recursively call parseArguments()
        if(parsedArguments.has(SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME)){
            List<String> argfiles = parsedArguments.valuesOf(SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME).stream()
                    .map(f -> (String)f)
                    .collect(Collectors.toList());

            List<String> newargs = argfiles.stream()
                    .distinct()
                    .filter( file -> !argumentsFilesLoadedAlready.contains(file) )
                    .flatMap(file -> loadArgumentsFile(file).stream())
                    .collect(Collectors.toList());
            argumentsFilesLoadedAlready.addAll(argfiles);

            if ( !newargs.isEmpty()) {
                newargs.addAll(Arrays.asList(args));
                return parseArguments(messageStream, newargs.toArray(new String[newargs.size()]));
            }
        }

        if(parsedArguments.has(SpecialArgumentsCollection.HELP_FULLNAME)){
            usage(messageStream, true);
            return true;
        } else if (parsedArguments.has(SpecialArgumentsCollection.VERSION_FULLNAME)) {
            messageStream.println(getVersion());
            return true;
        }

        for (OptionSpec<?> optSpec : parsedArguments.asMap().keySet()){
            if(parsedArguments.has(optSpec)) {
                ArgumentDefinition argDef = argumentMap.get(optSpec.options().get(0));
                if (!setArgument(argDef, (List<String>) optSpec.values(parsedArguments))) {
                    messageStream.println();
                    usage(messageStream, true);
                    return false;
                }
            }
        }
        for (Object arg: parsedArguments.nonOptionArguments()) {
            if (!parsePositionalArgument((String)arg)) {
                messageStream.println();
                usage(messageStream, false);
                return false;
            }
        }

        if (!checkNumArguments()) {
            messageStream.println();
            usage(messageStream, false);
            return false;
        }

        return true;
    }

    /**
     * After command line has been parsed, make sure that all required arguments have values, and that
     * lists with minimum # of elements have sufficient.
     *
     * @return true if valid
     */
    private boolean checkNumArguments() {
        try {
            for (final ArgumentDefinition argumentDefinition : argumentDefinitions) {
                final String fullName = argumentDefinition.fieldName;
                final StringBuilder mutextArgumentNames = new StringBuilder();
                for (final String mutexArgument : argumentDefinition.mutuallyExclusive) {
                    final ArgumentDefinition mutextArgumentDef = argumentMap.get(mutexArgument);
                    if (mutextArgumentDef != null && mutextArgumentDef.hasBeenSet) {
                        mutextArgumentNames.append(" ").append(mutextArgumentDef.fieldName);
                    }
                }
                if (argumentDefinition.hasBeenSet && mutextArgumentNames.length() > 0) {
                    messageStream.println("ERROR: Argument '" + fullName +
                            "' cannot be used in conjunction with argument(s)" +
                            mutextArgumentNames.toString());
                    return false;
                }
                if (argumentDefinition.isCollection) {
                    @SuppressWarnings("rawtypes")
                    final Collection c = (Collection) argumentDefinition.getFieldValue();
                    if (c.size() < argumentDefinition.minElements) {
                        messageStream.println("ERROR: Argument '" + fullName + "' must be specified at least " +
                                argumentDefinition.minElements + " times.");
                        return false;
                    }
                } else if (!argumentDefinition.optional && !argumentDefinition.hasBeenSet && mutextArgumentNames.length() == 0) {
                    messageStream.print("ERROR: Argument '" + fullName + "' is required");
                    if (argumentDefinition.mutuallyExclusive.isEmpty()) {
                        messageStream.println(".");
                    } else {
                        messageStream.println(" unless any of " + argumentDefinition.mutuallyExclusive +
                                " are specified.");
                    }
                    return false;
                }

            }
            if (positionalArguments != null) {
                @SuppressWarnings("rawtypes")
                final Collection c = (Collection) positionalArguments.get(positionalArgumentsParent);
                if (c.size() < minPositionalArguments) {
                    messageStream.println("ERROR: At least " + minPositionalArguments +
                            " positional arguments must be specified.");
                    return false;
                }

            }

            return true;
        } catch (final IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("Should never happen",e);
        }


    }

    @SuppressWarnings("unchecked")
    private boolean parsePositionalArgument(final String stringValue) {
        if (positionalArguments == null) {
            messageStream.println("ERROR: Invalid argument '" + stringValue + "'.");
            return false;
        }
        final Object value;
        try {
            value = constructFromString(getUnderlyingType(positionalArguments), stringValue);
        } catch (final GATKException.CommandLineParserInternalException e) {
            messageStream.println("ERROR: " + e.getMessage());
            return false;
        }
        @SuppressWarnings("rawtypes")
        final Collection c;
        try {
            c = (Collection) positionalArguments.get(callerArguments);
        } catch (final IllegalAccessException e) {
            throw new RuntimeException(e);
        }
        if (c.size() >= maxPositionalArguments) {
            messageStream.println("ERROR: No more than " + maxPositionalArguments +
                    " positional arguments may be specified on the command line.");
            return false;
        }

        c.add(value);
        return true;
    }

    @SuppressWarnings("unchecked")
    private boolean setArgument(ArgumentDefinition argumentDefinition, final List<String> values) {
        if (argumentDefinition.isFlag() && values.isEmpty()){
            argumentDefinition.hasBeenSet = true;
            argumentDefinition.setFieldValue(true);
            return true;
        }

        if (!argumentDefinition.isCollection) {
            if (argumentDefinition.hasBeenSet || values.size() > 1) {
                messageStream.println("ERROR: Argument '" + argumentDefinition.getNames() + "' cannot be specified more than once.");
                return false;
            }
        }


        for (String stringValue: values) {
            final Object value;
            try {
                if (stringValue.equals(NULL_STRING)) {
                    //"null" is a special value that allows the user to override any default
                    //value set for this arg
                    if (argumentDefinition.optional) {
                        value = null;
                    } else {
                        messageStream.println("ERROR: non-null value must be provided for '" + argumentDefinition.getNames() + "'.");
                        return false;
                    }
                } else {
                    value = constructFromString(getUnderlyingType(argumentDefinition.field), stringValue);
                }

            } catch (final GATKException.CommandLineParserInternalException e) {
                messageStream.println("ERROR: " + e.getMessage());
                return false;
            }
            if (argumentDefinition.isCollection) {
                @SuppressWarnings("rawtypes")
                final Collection c = (Collection) argumentDefinition.getFieldValue();
                if (value == null) {
                    //user specified this arg=null which is interpreted as empty list
                    c.clear();
                } else if (c.size() >= argumentDefinition.maxElements) {
                    messageStream.println("ERROR: Argument '" + argumentDefinition.getNames() + "' cannot be used more than " +
                            argumentDefinition.maxElements + " times.");
                    return false;
                } else {
                    c.add(value);
                }
                argumentDefinition.hasBeenSet = true;
            } else {
                argumentDefinition.setFieldValue(value);
                argumentDefinition.hasBeenSet = true;
            }
        }
        return true;
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
        if (argumentDefinition.optional && !argumentDefinition.isCollection) {
            sb.append("Default value: ");
            sb.append(argumentDefinition.defaultValue);
            sb.append(". ");
            if (!argumentDefinition.defaultValue.equals(NULL_STRING)) {
                sb.append("This argument can be set to 'null' to clear the default value. ");
            }
        } else if (!argumentDefinition.isCollection) {
            sb.append("Required. ");
        }
        Object[] enumConstants = getUnderlyingType(argumentDefinition.field).getEnumConstants();
        if (enumConstants == null && getUnderlyingType(argumentDefinition.field) == Boolean.class) {
            enumConstants = TRUE_FALSE_VALUES;
        }

        if (enumConstants != null) {
            final Boolean isClpEnum = enumConstants.length > 0 && (enumConstants[0] instanceof ClpEnum);

            sb.append("Possible values: {");
            if (isClpEnum) sb.append("\n");

            for (int i = 0; i < enumConstants.length; ++i) {
                if (i > 0 && !isClpEnum) {
                    sb.append(", ");
                }
                sb.append(enumConstants[i].toString());

                if (isClpEnum) {
                    sb.append(" (").append(((ClpEnum) enumConstants[i]).getHelpDoc()).append(")\n");
                }
            }
            sb.append("} ");
        }
        if (argumentDefinition.isCollection) {
            if (argumentDefinition.minElements == 0) {
                if (argumentDefinition.maxElements == Integer.MAX_VALUE) {
                    sb.append("This argument may be specified 0 or more times. ");
                } else {
                    sb.append("This argument must be specified no more than ").append(argumentDefinition.maxElements).append(
                            " times. ");
                }
            } else if (argumentDefinition.maxElements == Integer.MAX_VALUE) {
                sb.append("This argument must be specified at least ").append(argumentDefinition.minElements).append(" times. ");
            } else {
                sb.append("This argument may be specified between ").append(argumentDefinition.minElements).append(
                        " and ").append(argumentDefinition.maxElements).append(" times. ");
            }

            if (!argumentDefinition.defaultValue.equals(NULL_STRING)) {
                sb.append("This argument can be set to 'null' to clear the default list. ");
            }

        }
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
                if (argumentAnnotation.maxElements() == 0) {
                    throw new GATKException.CommandLineParserInternalException("@Argument member " + field.getName() +
                            "has maxElements = 0");
                }
                if (argumentAnnotation.minElements() > argumentAnnotation.maxElements()) {
                    throw new GATKException.CommandLineParserInternalException("In @Argument member " + field.getName() +
                            ", minElements cannot be > maxElements");
                }
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


    private static boolean isCollectionField(final Field field) {
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
    private Class<?> getUnderlyingType(final Field field) {
        if (isCollectionField(field)) {
            final ParameterizedType clazz = (ParameterizedType) (field.getGenericType());
            final Type[] genericTypes = clazz.getActualTypeArguments();
            if (genericTypes.length != 1) {
                throw new GATKException.CommandLineParserInternalException("Strange collection type for field " +
                        field.getName());
            }
            return (Class<?>) genericTypes[0];

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
            clazz.getConstructor(String.class);
            return true;
        } catch (final NoSuchMethodException e) {
            return false;
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private Object constructFromString(final Class clazz, final String s) {
        try {
            if (clazz.isEnum()) {
                try {
                    return Enum.valueOf(clazz, s);
                } catch (final IllegalArgumentException e) {
                    throw new GATKException.CommandLineParserInternalException("'" + s + "' is not a valid value for " +
                            clazz.getSimpleName() + ".", e);
                }
            }
            final Constructor<?> ctor = clazz.getConstructor(String.class);
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
            throw new GATKException.CommandLineParserInternalException("Problem constructing " + clazz.getSimpleName() +
                    " from the string '" + s + "'.", e.getCause());
        }
    }

    public String[] getArgv() {
        return argv;
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
        final int minElements;
        final int maxElements;
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
            this.optional = annotation.optional() || getFieldValue() != null ;
            this.isCollection = isCollectionField(field);

            if ( this.isFlag()){
                this.minElements = 0;
                this.maxElements = 1;
            } else {
                this.minElements = annotation.minElements();
                this.maxElements = annotation.maxElements();
            }
            this.isCommon = annotation.common();
            this.isSpecial = annotation.special();


            this.mutuallyExclusive = new HashSet<>(Arrays.asList(annotation.mutex()));

            Object tmpDefault = getFieldValue();
            if (tmpDefault != null) {
                if (isCollection && ((Collection) tmpDefault).isEmpty()) {
                    //treat empty collections the same as uninitialized primitive types
                    this.defaultValue = NULL_STRING;
                } else {
                    //this is an intialized primitive type or a non-empty collection
                    this.defaultValue = tmpDefault.toString();
                }
            } else {
                this.defaultValue = NULL_STRING;
            }
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
        argumentDefinitions.stream()
                .filter(argumentDefinition -> argumentDefinition.hasBeenSet)
                .forEach(argumentDefinition ->
                        commandLineString.append(" --")
                                .append(argumentDefinition.getLongName())
                                .append(" ")
                                .append(argumentDefinition.getFieldValue()));

        commandLineString.append("   "); //separator to tell the 2 apart
        //next, append args that weren't explicitly set, but have a default value
        argumentDefinitions.stream()
                .filter(argumentDefinition -> !argumentDefinition.hasBeenSet && !argumentDefinition.defaultValue.equals(NULL_STRING))
                .forEach(argumentDefinition ->
                        commandLineString.append(" --").append(argumentDefinition.getLongName()).append(" ")
                                .append(argumentDefinition.defaultValue));

        return toolName + " " + commandLineString.toString();
    }

    /**
     * This method is only needed when calling one of the public methods that doesn't take a messageStream argument.
     */
    public void setMessageStream(final PrintStream messageStream) {
        this.messageStream = messageStream;
    }
}
