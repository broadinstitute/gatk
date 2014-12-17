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
package picard.cmdline;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil.MultiMap;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Annotation-driven utility for parsing command-line arguments, checking for errors, and producing usage message.
 * <p/>
 * This class supports options of the form KEY=VALUE, plus positional arguments.  Positional arguments must not contain
 * an equal sign lest they be mistaken for a KEY=VALUE pair.
 * <p/>
 * The caller must supply an object that both defines the command line and has the parsed options set into it.
 * For each possible KEY=VALUE option, there must be a public data member annotated with @Option.  The KEY name is
 * the name of the data member.  An abbreviated name may also be specified with the shortName attribute of @Option.
 * If the data member is a List<T>, then the option may be specified multiple times.  The type of the data member,
 * or the type of the List element must either have a ctor T(String), or must be an Enum.  List options must
 * be initialized by the caller with some kind of list.  Any other option that is non-null is assumed to have the given
 * value as a default.  If an option has no default value, and does not have the optional attribute of @Option set,
 * is required.  For List options, minimum and maximum number of elements may be specified in the @Option annotation.
 * <p/>
 * A single List data member may be annotated with the @PositionalArguments.  This behaves similarly to a Option
 * with List data member: the caller must initialize the data member, the type must be constructable from String, and
 * min and max number of elements may be specified.  If no @PositionalArguments annotation appears in the object,
 * then it is an error for the command line to contain positional arguments.
 * <p/>
 * A single String public data member may be annotated with @Usage.  This string, if present, is used to
 * construct the usage message.  Details about the possible options are automatically appended to this string.
 * If @Usage does not appear, a boilerplate usage message is used.
 */
public class CommandLineParser {
    // For formatting option section of usage message.
    private static final int OPTION_COLUMN_WIDTH = 30;
    private static final int DESCRIPTION_COLUMN_WIDTH = 90;

    private static final Boolean[] TRUE_FALSE_VALUES = {Boolean.TRUE, Boolean.FALSE};

    private static final String[] PACKAGES_WITH_WEB_DOCUMENTATION = {"picard"};

    // Use these if no @Usage annotation
    private static final String defaultUsagePreamble = "Usage: program [options...]\n";
    private static final String defaultUsagePreambleWithPositionalArguments =
            "Usage: program [options...] [positional-arguments...]\n";
    private static final String OPTIONS_FILE = "OPTIONS_FILE";

    private static final String PRECEDENCE_SYMBOL = "++";

    /** name, shortName, description for options built in to framework */
    private static final String[][] FRAMEWORK_OPTION_DOC = {
            {"--help", "-h", "Displays options specific to this tool."},
            {"--stdhelp", "-H", "Displays options specific to this tool AND " +
                    "options common to all Picard command line tools."},
            {"--version", null, "Displays program version."}
    };

    private final Set<String> optionsThatCannotBeOverridden = new HashSet<String>();

    /**
     * A typical command line program will call this to get the beginning of the usage message,
     * and then append a description of the program, like this:
     * <p/>
     * \@Usage
     * public String USAGE = CommandLineParser.getStandardUsagePreamble(getClass()) + "Frobnicates the freebozzle."
     */
    public static String getStandardUsagePreamble(final Class mainClass) {
        return "USAGE: " + mainClass.getSimpleName() + " [options]\n\n" +
                (hasWebDocumentation(mainClass) ?
                        "Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#" +
                                mainClass.getSimpleName() + "\n\n"
                        : "");
    }

    /**
     * Determines if a class has web documentation based on its package name
     *
     * @param clazz
     * @return true if the class has web documentation, false otherwise
     */
    public static boolean hasWebDocumentation(final Class clazz) {
        for (final String pkg : PACKAGES_WITH_WEB_DOCUMENTATION) {
            if (clazz.getPackage().getName().startsWith(pkg)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the link to a FAQ
     */
    public static String getFaqLink() {
        return "To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp";
    }

    /**
     * Find all of the members annotated with @NestedOptions.
     * This is package scope and static so that CommandLineProgram can use it to provide default implementation
     * of its own getNestedOptions() method.
     */
    static Map<String, Object> getNestedOptions(final Object callerOptions) {
        // LinkedHashMap so usage message is generated in order of declaration
        final Map<String, Object> ret = new LinkedHashMap<String, Object>();
        final Class<?> clazz = callerOptions.getClass();
        for (final Field field : getAllFields(clazz)) {
            if (field.getAnnotation(NestedOptions.class) != null) {
                field.setAccessible(true);
                try {
                    ret.put(field.getName(), field.get(callerOptions));
                } catch (final IllegalAccessException e) {
                    throw new RuntimeException("Should never happen.", e);
                }
            }
        }
        return ret;
    }

    // This is the object that the caller has provided that contains annotations,
    // and into which the values will be assigned.
    private final Object callerOptions;

    // For child CommandLineParser, this contains the prefix for the option names, which is needed for generating
    // the command line.  For non-nested, this is the empty string.
    private final String prefix;
    // For non-nested, empty string.  For nested, prefix + "."
    private final String prefixDot;

    // null if no @PositionalArguments annotation
    private Field positionalArguments;
    private int minPositionalArguments;
    private int maxPositionalArguments;

    // List of all the data members with @Option annotation
    private final List<OptionDefinition> optionDefinitions = new ArrayList<OptionDefinition>();

    // Maps long name, and short name, if present, to an option definition that is
    // also in the optionDefinitions list.
    private final Map<String, OptionDefinition> optionMap = new HashMap<String, OptionDefinition>();

    // Maps child options prefix to CommandLineParser for the child object.
    // Key: option prefix.
    private final Map<String, CommandLineParser> childOptionsMap = new LinkedHashMap<String, CommandLineParser>();

    // Holds the command-line arguments for a child option parser.
    // Key: option prefix.  Value: List of arguments for child corresponding to that prefix (with prefix stripped).
    private final MultiMap<String, ChildOptionArg> childOptionArguments = new MultiMap<String, ChildOptionArg>();

    // For printing error messages when parsing command line.
    private PrintStream messageStream;

    // In case implementation wants to get at arg for some reason.
    private String[] argv;

    private String programVersion = null;

    // The command line used to launch this program, including non-null default options that
    // weren't explicitly specified. This is used for logging and debugging.
    private String commandLine = "";

    /**
     * This attribute is here just to facilitate printing usage for OPTIONS_FILE
     */
    public File IGNORE_THIS_PROPERTY;

    // The associated program properties using the CommandLineProgramProperties annotation
    private final CommandLineProgramProperties programProperties;

    /**
     * Prepare for parsing command line arguments, by validating annotations.
     *
     * @param callerOptions This object contains annotations that define the acceptable command-line options,
     *                      and ultimately will receive the settings when a command line is parsed.
     */
    public CommandLineParser(final Object callerOptions) {
        this(callerOptions, "");
    }

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

    /**
     * @param prefix Non-empty for child options object.
     */
    private CommandLineParser(final Object callerOptions, final String prefix) {
        this.callerOptions = callerOptions;

        this.prefix = prefix;
        if (prefix.isEmpty()) {
            prefixDot = "";
        } else {
            prefixDot = prefix + ".";
        }


        for (final Field field : getAllFields(this.callerOptions.getClass())) {
            if (field.getAnnotation(PositionalArguments.class) != null) {
                handlePositionalArgumentAnnotation(field);
            }
            if (field.getAnnotation(Option.class) != null) {
                handleOptionAnnotation(field);
            } else if (!isCommandLineProgram() && field.getAnnotation(NestedOptions.class) != null) {
                // If callerOptions is an instance of CommandLineProgram, defer creation of child
                // CommandLineParsers until after parsing options for this parser, in case CommandLineProgram
                // wants to do something dynamic based on values for this parser.
                handleNestedOptionsAnnotation(field);
            }
        }

        this.programProperties = this.callerOptions.getClass().getAnnotation(CommandLineProgramProperties.class);
    }

    private boolean isCommandLineProgram() {
        return callerOptions instanceof CommandLineProgram;
    }

    private static List<Field> getAllFields(Class clazz) {
        final List<Field> ret = new ArrayList<Field>();
        do {
            ret.addAll(Arrays.asList(clazz.getDeclaredFields()));
            clazz = clazz.getSuperclass();
        } while (clazz != null);
        return ret;
    }

    public String getVersion() {
        return this.callerOptions.getClass().getPackage().getImplementationVersion();
    }

    /**
     * Print a usage message based on the options object passed to the ctor.
     *
     * @param stream Where to write the usage message.
     */
    public void usage(final PrintStream stream, final boolean printCommon) {
        if (prefix.isEmpty()) {
            stream.print(getStandardUsagePreamble(callerOptions.getClass()) + getUsagePreamble());
            stream.println("\nVersion: " + getVersion());
            stream.println("\n\nOptions:\n");

            for (final String[] optionDoc : FRAMEWORK_OPTION_DOC) {
                printOptionParamUsage(stream, optionDoc[0], optionDoc[1], null, optionDoc[2]);
            }
        }

        if (!optionDefinitions.isEmpty()) {
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                if (printCommon || !optionDefinition.isCommon) printOptionUsage(stream, optionDefinition);
            }
        }

        if (printCommon) {
            final Field fileField;
            try {
                fileField = getClass().getField("IGNORE_THIS_PROPERTY");
            } catch (final NoSuchFieldException e) {
                throw new PicardException("Should never happen", e);
            }
            final OptionDefinition optionsFileOptionDefinition =
                    new OptionDefinition(fileField, OPTIONS_FILE, "",
                            "File of OPTION_NAME=value pairs.  No positional parameters allowed.  Unlike command-line options, " +
                                    "unrecognized options are ignored.  " + "A single-valued option set in an options file may be overridden " +
                                    "by a subsequent command-line option.  " +
                                    "A line starting with '#' is considered a comment.",
                            false, true, false, 0, Integer.MAX_VALUE, null, true, new String[0]);
            printOptionUsage(stream, optionsFileOptionDefinition);
        }

        // Generate usage for child parsers.
        final Collection<CommandLineParser> childClps;
        childClps = getChildParsersForHelp();
        for (final CommandLineParser childClp : childClps) {
            childClp.usage(stream, printCommon);
        }
    }

    private Collection<CommandLineParser> getChildParsersForHelp() {
        final Collection<CommandLineParser> childClps;
        if (isCommandLineProgram()) {
            childClps = new ArrayList<CommandLineParser>();
            for (final Map.Entry<String, Object> entry :
                    ((CommandLineProgram) callerOptions).getNestedOptionsForHelp().entrySet()) {
                if (entry.getKey().contains(".")) {
                    throw new IllegalArgumentException("Prefix for nested options should not contain period: " + entry.getKey());
                }
                childClps.add(new CommandLineParser(entry.getValue(), prefixDot + entry.getKey()));
            }
        } else {
            childClps = childOptionsMap.values();
        }
        return childClps;
    }


    public void htmlUsage(final PrintStream stream, final String programName, final boolean printCommon) {
        // TODO: Should HTML escape usage preamble and option usage, including line breaks
        stream.println("<a id=\"" + programName + "\"/>");
        stream.println("<h3>" + programName + "</h3>");
        stream.println("<section>");
        stream.println("<p>" + htmlEscape(getUsagePreamble()) + "</p>");
        boolean hasOptions = false;
        for (final OptionDefinition optionDefinition : optionDefinitions) {
            if (!optionDefinition.isCommon || printCommon) {
                hasOptions = true;
                break;
            }
        }
        if (hasOptions) {
            htmlPrintOptions(stream, printCommon);
        }
        stream.println("</section>");
    }

    public void htmlPrintOptions(final PrintStream stream, final boolean printCommon) {
        stream.println("<table>");
        stream.println("<tr><th>Option</th><th>Description</th></tr>");
        if (printCommon) {
            for (final String[] optionDoc : FRAMEWORK_OPTION_DOC) {
                stream.println("<tr><td>" + optionDoc[0] + "</td><td>" +
                        htmlEscape(optionDoc[2]) + "</td></tr>");
            }
        }
        htmlPrintOptionTableRows(stream, printCommon);
        stream.println("</table>");
    }

    /**
     * Prints options as rows in an HTML table.
     *
     * @param stream
     * @param printCommon
     */
    private void htmlPrintOptionTableRows(final PrintStream stream, final boolean printCommon) {
        for (final OptionDefinition optionDefinition : optionDefinitions) {
            if (!optionDefinition.isCommon || printCommon) {
                printHtmlOptionUsage(stream, optionDefinition);
            }
        }
        for (final CommandLineParser childParser : getChildParsersForHelp()) {
            childParser.htmlPrintOptionTableRows(stream, false);
        }

    }

    private static String htmlEscape(String str) {
        // May need more here
        str = str.replaceAll("<", "&lt;");
        str = str.replaceAll("\n", "\n<p>");
        return str;
    }

    /**
     * Parse command-line options, and store values in callerOptions object passed to ctor.
     *
     * @param messageStream Where to write error messages.
     * @param args          Command line tokens.
     * @return true if command line is valid.
     */
    public boolean parseOptions(final PrintStream messageStream, final String[] args) {
        this.argv = args;
        this.messageStream = messageStream;
        if (prefix.isEmpty()) {
            commandLine = callerOptions.getClass().getName();
        }
        for (int i = 0; i < args.length; ++i) {
            final String arg = args[i];
            if (arg.equals("-h") || arg.equals("--help")) {
                usage(messageStream, false);
                return false;
            }
            if (arg.equals("-H") || arg.equals("--stdhelp")) {
                usage(messageStream, true);
                return false;
            }

            if (arg.equals("--version")) {
                messageStream.println(getVersion());
                return false;
            }


            final String[] pair = arg.split("=", 2);
            if (pair.length == 2 && pair[1].length() == 0) {

                if (i < args.length - 1) {
                    pair[1] = args[++i];
                }
            }
            if (pair.length == 2) {
                if (!parseOption(pair[0], pair[1], false)) {
                    messageStream.println();
                    usage(messageStream, true);
                    return false;
                }
            } else if (!parsePositionalArgument(arg)) {
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

        if (!parseChildOptions()) {
            messageStream.println();
            usage(messageStream, false);
            return false;
        }

        return true;
    }

    /**
     * After command line has been parsed, make sure that all required options have values, and that
     * lists with minimum # of elements have sufficient.
     *
     * @return true if valid
     */
    private boolean checkNumArguments() {
        //Also, since we're iterating over all options and args, use this opportunity to recreate the commandLineString
        final StringBuilder commandLineString = new StringBuilder();
        try {
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                final String fullName = prefixDot + optionDefinition.name;
                final StringBuilder mutextOptionNames = new StringBuilder();
                for (final String mutexOption : optionDefinition.mutuallyExclusive) {
                    final OptionDefinition mutextOptionDef = optionMap.get(mutexOption);
                    if (mutextOptionDef != null && mutextOptionDef.hasBeenSet) {
                        mutextOptionNames.append(" ").append(prefixDot).append(mutextOptionDef.name);
                    }
                }
                if (optionDefinition.hasBeenSet && mutextOptionNames.length() > 0) {
                    messageStream.println("ERROR: Option '" + fullName +
                            "' cannot be used in conjunction with option(s)" +
                            mutextOptionNames.toString());
                    return false;
                }
                if (optionDefinition.isCollection) {
                    final Collection c = (Collection) optionDefinition.field.get(callerOptions);
                    if (c.size() < optionDefinition.minElements) {
                        messageStream.println("ERROR: Option '" + fullName + "' must be specified at least " +
                                optionDefinition.minElements + " times.");
                        return false;
                    }
                } else if (!optionDefinition.optional && !optionDefinition.hasBeenSet &&
                        !optionDefinition.hasBeenSetFromParent && mutextOptionNames.length() == 0) {
                    messageStream.print("ERROR: Option '" + fullName + "' is required");
                    if (optionDefinition.mutuallyExclusive.isEmpty()) {
                        messageStream.println(".");
                    } else {
                        messageStream.println(" unless any of " + optionDefinition.mutuallyExclusive +
                                " are specified.");
                    }
                    return false;
                }

            }
            if (positionalArguments != null) {
                final Collection c = (Collection) positionalArguments.get(callerOptions);
                if (c.size() < minPositionalArguments) {
                    messageStream.println("ERROR: At least " + minPositionalArguments +
                            " positional arguments must be specified.");
                    return false;
                }
                for (final Object posArg : c) {
                    commandLineString.append(" ").append(posArg.toString());
                }
            }
            //first, append args that were explicitly set
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                if (optionDefinition.hasBeenSet) {
                    commandLineString.append(" ").append(prefixDot).append(optionDefinition.name).append("=").append(
                            optionDefinition.field.get(callerOptions));
                }
            }
            commandLineString.append("   "); //separator to tell the 2 apart
            //next, append args that weren't explicitly set, but have a default value
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                if (!optionDefinition.hasBeenSet && !optionDefinition.defaultValue.equals("null")) {
                    commandLineString.append(" ").append(prefixDot).append(optionDefinition.name).append("=").append(
                            optionDefinition.defaultValue);
                }
            }
            this.commandLine += commandLineString.toString();
            return true;
        } catch (final IllegalAccessException e) {
            // Should never happen because lack of publicness has already been checked.
            throw new RuntimeException(e);
        }


    }

    private boolean parsePositionalArgument(final String stringValue) {
        if (positionalArguments == null) {
            messageStream.println("ERROR: Invalid argument '" + stringValue + "'.");
            return false;
        }
        final Object value;
        try {
            value = constructFromString(getUnderlyingType(positionalArguments), stringValue);
        } catch (final CommandLineParseException e) {
            messageStream.println("ERROR: " + e.getMessage());
            return false;
        }
        final Collection c;
        try {
            c = (Collection) positionalArguments.get(callerOptions);
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

    private boolean parseOption(final String key, final String stringValue, final boolean optionsFile) {
        return parseOption(key, stringValue, optionsFile, false);
    }

    private boolean parseOption(String key, final String stringValue, final boolean optionsFile,
                                boolean precedenceSet) {
        key = key.toUpperCase();
        if (key.equals(OPTIONS_FILE)) {
            commandLine += " " + prefix + OPTIONS_FILE + "=" + stringValue;
            return parseOptionsFile(stringValue);
        }

        // Check to see if the precedence symbol was used
        if (key.startsWith(PRECEDENCE_SYMBOL)) {
            key = key.substring(PRECEDENCE_SYMBOL.length());
            precedenceSet = true;
        }

        // Save child options for later processing.
        final Integer prefixIndex = key.indexOf('.');
        if (prefixIndex != -1) {
            final String prefix = key.substring(0, prefixIndex);
            final String subKey = key.substring(prefixIndex + 1);
            if (!subKey.isEmpty()) {
                childOptionArguments.append(prefix, new ChildOptionArg(subKey, stringValue, optionsFile,
                        precedenceSet));
                return true;
            } else {
                messageStream.println("ERROR: Unrecognized option: " + key);
                return false;
            }
        }

        final OptionDefinition optionDefinition = optionMap.get(key);
        if (optionDefinition == null) {
            if (optionsFile) {
                // Silently ignore unrecognized option from options file
                return true;
            }
            messageStream.println("ERROR: Unrecognized option: " + key);
            return false;
        }

        // Check to see if the option has been "fixed" already
        if (this.optionsThatCannotBeOverridden.contains(optionDefinition.name)) {
            return true;
        } else if (precedenceSet) {
            this.optionsThatCannotBeOverridden.add(optionDefinition.name);
        }

        if (!optionDefinition.isCollection) {
            if (optionDefinition.hasBeenSet && !optionDefinition.hasBeenSetFromOptionsFile) {
                messageStream.println("ERROR: Option '" + key + "' cannot be specified more than once.");
                return false;
            }
        }
        final Object value;
        try {
            if (stringValue.equals("null")) {
                //"null" is a special value that allows the user to override any default
                //value set for this arg. It can only be used for optional args. When
                //used for a list arg, it will clear the list.
                if (optionDefinition.optional) {
                    value = null;
                } else {
                    messageStream.println("ERROR: non-null value must be provided for '" + key + "'.");
                    return false;
                }
            } else {
                value = constructFromString(getUnderlyingType(optionDefinition.field), stringValue);
            }

        } catch (final CommandLineParseException e) {
            messageStream.println("ERROR: " + e.getMessage());
            return false;
        }
        try {
            if (optionDefinition.isCollection) {
                final Collection c = (Collection) optionDefinition.field.get(callerOptions);
                if (value == null) {
                    //user specified this arg=null which is interpreted as empty list
                    c.clear();
                } else if (c.size() >= optionDefinition.maxElements) {
                    messageStream.println("ERROR: Option '" + key + "' cannot be used more than " +
                            optionDefinition.maxElements + " times.");
                    return false;
                } else {
                    c.add(value);
                }
                optionDefinition.hasBeenSet = true;
                optionDefinition.hasBeenSetFromOptionsFile = optionsFile;
            } else {
                optionDefinition.field.set(callerOptions, value);
                optionDefinition.hasBeenSet = true;
                optionDefinition.hasBeenSetFromOptionsFile = optionsFile;
            }
        } catch (final IllegalAccessException e) {
            // Should never happen because we only iterate through public fields.
            throw new RuntimeException(e);
        }
        return true;
    }

    /**
     * Parsing of options from file is looser than normal.  Any unrecognized options are
     * ignored, and a single-valued option that is set in a file may be overridden by a
     * subsequent appearance of that option.
     * A line that starts with '#' is ignored.
     *
     * @param optionsFile
     * @return false if a fatal error occurred
     */
    private boolean parseOptionsFile(final String optionsFile) {
        return parseOptionsFile(optionsFile, true);
    }

    /**
     * @param optionFileStyleValidation true: unrecognized options are silently ignored; and a single-valued option may be overridden.
     *                                  false: standard rules as if the options in the file were on the command line directly.
     * @return
     */
    public boolean parseOptionsFile(final String optionsFile, final boolean optionFileStyleValidation) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(optionsFile));
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#") || line.trim().length() == 0) {
                    continue;
                }
                final String[] pair = line.split("=", 2);
                if (pair.length == 2) {
                    if (!parseOption(pair[0], pair[1], optionFileStyleValidation)) {
                        messageStream.println();
                        usage(messageStream, true);
                        return false;
                    }
                } else {
                    messageStream.println("Strange line in OPTIONS_FILE " + optionsFile + ": " + line);
                    usage(messageStream, true);
                    return false;
                }
            }
            reader.close();
            return true;

        } catch (final IOException e) {
            throw new PicardException("I/O error loading OPTIONS_FILE=" + optionsFile, e);
        } finally {
            CloserUtil.close(reader);
        }
    }

    private void printHtmlOptionUsage(final PrintStream stream, final OptionDefinition optionDefinition) {
        final String type = getUnderlyingType(optionDefinition.field).getSimpleName();
        final String optionLabel = prefixDot + optionDefinition.name + " (" + type + ")";
        stream.println("<tr><td>" + optionLabel + "</td><td>" +
                htmlEscape(makeOptionDescription(optionDefinition)) + "</td></tr>");
    }

    private void printOptionUsage(final PrintStream stream, final OptionDefinition optionDefinition) {
        printOptionParamUsage(stream, optionDefinition.name, optionDefinition.shortName,
                getUnderlyingType(optionDefinition.field).getSimpleName(),
                makeOptionDescription(optionDefinition));
    }


    private void printOptionParamUsage(final PrintStream stream, final String name, final String shortName,
                                       final String type, final String optionDescription) {
        String optionLabel = prefixDot + name;
        if (type != null) optionLabel += "=" + type;

        stream.print(optionLabel);
        if (shortName != null && shortName.length() > 0) {
            stream.println();
            optionLabel = prefixDot + shortName;
            if (type != null) optionLabel += "=" + type;
            stream.print(optionLabel);
        }

        int numSpaces = OPTION_COLUMN_WIDTH - optionLabel.length();
        if (optionLabel.length() > OPTION_COLUMN_WIDTH) {
            stream.println();
            numSpaces = OPTION_COLUMN_WIDTH;
        }
        printSpaces(stream, numSpaces);
        final String wrappedDescription = StringUtil.wordWrap(optionDescription, DESCRIPTION_COLUMN_WIDTH);
        final String[] descriptionLines = wrappedDescription.split("\n");
        for (int i = 0; i < descriptionLines.length; ++i) {
            if (i > 0) {
                printSpaces(stream, OPTION_COLUMN_WIDTH);
            }
            stream.println(descriptionLines[i]);
        }
        stream.println();
    }

    private String makeOptionDescription(final OptionDefinition optionDefinition) {
        final StringBuilder sb = new StringBuilder();
        if (optionDefinition.doc.length() > 0) {
            sb.append(optionDefinition.doc);
            sb.append("  ");
        }
        if (optionDefinition.optional && !optionDefinition.isCollection) {
            sb.append("Default value: ");
            sb.append(optionDefinition.defaultValue);
            sb.append(". ");
            if (!optionDefinition.defaultValue.equals("null")) {
                sb.append("This option can be set to 'null' to clear the default value. ");
            }
        } else if (!optionDefinition.isCollection) {
            sb.append("Required. ");
        }
        Object[] enumConstants = getUnderlyingType(optionDefinition.field).getEnumConstants();
        if (enumConstants == null && getUnderlyingType(optionDefinition.field) == Boolean.class) {
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
        if (optionDefinition.isCollection) {
            if (optionDefinition.minElements == 0) {
                if (optionDefinition.maxElements == Integer.MAX_VALUE) {
                    sb.append("This option may be specified 0 or more times. ");
                } else {
                    sb.append("This option must be specified no more than ").append(optionDefinition.maxElements).append(
                            " times. ");
                }
            } else if (optionDefinition.maxElements == Integer.MAX_VALUE) {
                sb.append("This option must be specified at least ").append(optionDefinition.minElements).append(" times. ");
            } else {
                sb.append("This option may be specified between ").append(optionDefinition.minElements).append(
                        " and ").append(optionDefinition.maxElements).append(" times. ");
            }

            if (!optionDefinition.defaultValue.equals("null")) {
                sb.append("This option can be set to 'null' to clear the default list. ");
            }

        }
        if (!optionDefinition.mutuallyExclusive.isEmpty()) {
            sb.append(" Cannot be used in conjuction with option(s)");
            for (final String option : optionDefinition.mutuallyExclusive) {
                final OptionDefinition mutextOptionDefinition = optionMap.get(option);

                if (mutextOptionDefinition == null) {
                    throw new PicardException("Invalid option definition in source code.  " + option +
                            " doesn't match any known option.");
                }

                sb.append(" ").append(mutextOptionDefinition.name);
                if (mutextOptionDefinition.shortName.length() > 0) {
                    sb.append(" (").append(mutextOptionDefinition.shortName).append(")");
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

    private void handleOptionAnnotation(final Field field) {
        try {
            field.setAccessible(true);
            final Option optionAnnotation = field.getAnnotation(Option.class);
            final boolean isCollection = isCollectionField(field);
            if (isCollection) {
                if (optionAnnotation.maxElements() == 0) {
                    throw new CommandLineParserDefinitionException("@Option member " + field.getName() +
                            "has maxElements = 0");
                }
                if (optionAnnotation.minElements() > optionAnnotation.maxElements()) {
                    throw new CommandLineParserDefinitionException("In @Option member " + field.getName() +
                            ", minElements cannot be > maxElements");
                }
                if (field.get(callerOptions) == null) {
                    createCollection(field, callerOptions, "@Option");
                }
            }
            if (!canBeMadeFromString(getUnderlyingType(field))) {
                throw new CommandLineParserDefinitionException("@Option member " + field.getName() +
                        " must have a String ctor or be an enum");
            }

            final OptionDefinition optionDefinition = new OptionDefinition(field,
                    field.getName(),
                    optionAnnotation.shortName(),
                    optionAnnotation.doc(), optionAnnotation.optional() || (field.get(callerOptions) != null),
                    optionAnnotation.overridable(), isCollection, optionAnnotation.minElements(),
                    optionAnnotation.maxElements(), field.get(callerOptions), optionAnnotation.common(),
                    optionAnnotation.mutex());

            for (final String option : optionAnnotation.mutex()) {
                final OptionDefinition mutextOptionDef = optionMap.get(option);
                if (mutextOptionDef != null) {
                    mutextOptionDef.mutuallyExclusive.add(field.getName());
                }
            }
            if (!optionDefinition.overridable && optionMap.containsKey(optionDefinition.name)) {
                throw new CommandLineParserDefinitionException(optionDefinition.name + " has already been used.");
            }
            if (optionDefinition.shortName.length() > 0) {
                if (optionMap.containsKey(optionDefinition.shortName)) {
                    if (!optionDefinition.overridable) {
                        throw new CommandLineParserDefinitionException(optionDefinition.shortName +
                                " has already been used");
                    }
                } else {
                    optionMap.put(optionDefinition.shortName, optionDefinition);
                }
            }
            //if we are overridable and we already exist don't add again to the option defs
            if (!(optionDefinition.overridable && optionMap.containsKey(optionDefinition.name))) {
                optionDefinitions.add(optionDefinition);
                optionMap.put(optionDefinition.name, optionDefinition);
            }
        } catch (final IllegalAccessException e) {
            throw new CommandLineParserDefinitionException(field.getName() +
                    " must have public visibility to have @Option annotation");
        }
    }

    private void handlePositionalArgumentAnnotation(final Field field) {
        if (positionalArguments != null) {
            throw new CommandLineParserDefinitionException
                    ("@PositionalArguments cannot be used more than once in an option class.");
        }
        field.setAccessible(true);
        positionalArguments = field;
        if (!isCollectionField(field)) {
            throw new CommandLineParserDefinitionException("@PositionalArguments must be applied to a Collection");
        }

        if (!canBeMadeFromString(getUnderlyingType(field))) {
            throw new CommandLineParserDefinitionException("@PositionalParameters member " + field.getName() +
                    "does not have a String ctor");
        }

        final PositionalArguments positionalArgumentsAnnotation = field.getAnnotation(PositionalArguments.class);
        minPositionalArguments = positionalArgumentsAnnotation.minElements();
        maxPositionalArguments = positionalArgumentsAnnotation.maxElements();
        if (minPositionalArguments > maxPositionalArguments) {
            throw new CommandLineParserDefinitionException("In @PositionalArguments, minElements cannot be > maxElements");
        }
        try {
            if (field.get(callerOptions) == null) {
                createCollection(field, callerOptions, "@PositionalParameters");
            }
        } catch (final IllegalAccessException e) {
            throw new CommandLineParserDefinitionException(field.getName() +
                    " must have public visibility to have @PositionalParameters annotation");

        }
    }

    private void handleNestedOptionsAnnotation(final Field field) {
        field.setAccessible(true);
        try {
            childOptionsMap.put(field.getName(),
                    new CommandLineParser(field.get(this.callerOptions), prefixDot + field.getName()));
        } catch (final IllegalAccessException e) {
            throw new CommandLineParserDefinitionException("Should never happen.", e);
        }
    }

    private boolean isCollectionField(final Field field) {
        try {
            field.getType().asSubclass(Collection.class);
            return true;
        } catch (final ClassCastException e) {
            return false;
        }
    }

    private void createCollection(final Field field, final Object callerOptions, final String annotationType)
            throws IllegalAccessException {
        try {
            field.set(callerOptions, field.getType().newInstance());
        } catch (final Exception ex) {
            try {
                field.set(callerOptions, new ArrayList());
            } catch (final IllegalArgumentException e) {
                throw new CommandLineParserDefinitionException("In collection " + annotationType +
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
    private Class getUnderlyingType(final Field field) {
        if (isCollectionField(field)) {
            final ParameterizedType clazz = (ParameterizedType) (field.getGenericType());
            final Type[] genericTypes = clazz.getActualTypeArguments();
            if (genericTypes.length != 1) {
                throw new CommandLineParserDefinitionException("Strange collection type for field " +
                        field.getName());
            }
            return (Class) genericTypes[0];

        } else {
            final Class type = field.getType();
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
    private boolean canBeMadeFromString(final Class clazz) {
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

    private Object constructFromString(final Class clazz, final String s) {
        try {
            if (clazz.isEnum()) {
                try {
                    return Enum.valueOf(clazz, s);
                } catch (final IllegalArgumentException e) {
                    throw new CommandLineParseException("'" + s + "' is not a valid value for " +
                            clazz.getSimpleName() + ".", e);
                }
            }
            final Constructor ctor = clazz.getConstructor(String.class);
            return ctor.newInstance(s);
        } catch (final NoSuchMethodException e) {
            // Shouldn't happen because we've checked for presence of ctor
            throw new CommandLineParseException("Cannot find string ctor for " + clazz.getName(), e);
        } catch (final InstantiationException e) {
            throw new CommandLineParseException("Abstract class '" + clazz.getSimpleName() +
                    "'cannot be used for an option value type.", e);
        } catch (final IllegalAccessException e) {
            throw new CommandLineParseException("String constructor for option value type '" + clazz.getSimpleName() +
                    "' must be public.", e);
        } catch (final InvocationTargetException e) {
            throw new CommandLineParseException("Problem constructing " + clazz.getSimpleName() +
                    " from the string '" + s + "'.", e.getCause());
        }
    }

    public String[] getArgv() {
        return argv;
    }

    public interface ClpEnum {
        String getHelpDoc();
    }

    protected static class OptionDefinition {
        final Field field;
        final String name;
        final String shortName;
        final String doc;
        final boolean optional;
        final boolean overridable;
        final boolean isCollection;
        final int minElements;
        final int maxElements;
        final String defaultValue;
        final boolean isCommon;
        boolean hasBeenSet = false;
        boolean hasBeenSetFromOptionsFile = false;
        boolean hasBeenSetFromParent = false;
        final Set<String> mutuallyExclusive;

        private OptionDefinition(final Field field, final String name, final String shortName, final String doc,
                                 final boolean optional, final boolean overridable, boolean collection, final int minElements,
                                 final int maxElements, final Object defaultValue, final boolean isCommon,
                                 final String[] mutuallyExclusive) {
            this.field = field;
            this.name = name.toUpperCase();
            this.shortName = shortName.toUpperCase();
            this.doc = doc;
            this.optional = optional;
            this.overridable = overridable;
            isCollection = collection;
            this.minElements = minElements;
            this.maxElements = maxElements;
            if (defaultValue != null) {
                if (isCollection && ((Collection) defaultValue).isEmpty()) {
                    //treat empty collections the same as uninitialized primitive types
                    this.defaultValue = "null";
                } else {
                    //this is an intialized primitive type or a non-empty collection
                    this.defaultValue = defaultValue.toString();
                }
            } else {
                this.defaultValue = "null";
            }
            this.isCommon = isCommon;
            this.mutuallyExclusive = new HashSet<String>(Arrays.asList(mutuallyExclusive));
        }
    }

    /**
     * Holds a command-line argument that is destined for a child parser.  Prefix has been stripped from name.
     */
    private static class ChildOptionArg {
        final String name;
        final String value;
        final boolean fromFile;
        final boolean precedenceSet;

        private ChildOptionArg(final String name, final String value, final boolean fromFile, final boolean precedenceSet) {
            this.name = name;
            this.value = value;
            this.fromFile = fromFile;
            this.precedenceSet = precedenceSet;
        }
    }

    /**
     * Propagate options from parent to children as appropriate, parse command line options for
     * children, and then validate that children have been properly initialized.  This is done recursively
     * for any child that itself has a child.
     *
     * @return true if parsing is successful.  Writes any errors to the message stream.
     */
    private boolean parseChildOptions() {

        // If callerOptions is not an instance of CommandLineProgram, then the child options are populated
        // when annotations are processed.
        if (isCommandLineProgram()) {
            final CommandLineProgram commandLineProgram = (CommandLineProgram) callerOptions;
            for (final Map.Entry<String, Object> entry : commandLineProgram.getNestedOptions().entrySet()) {
                if (entry.getKey().contains(".")) {
                    throw new IllegalArgumentException("Prefix for nested options should not contain period: " + entry.getKey());
                }
                childOptionsMap.put(entry.getKey(),
                        new CommandLineParser(entry.getValue(), prefixDot + entry.getKey()));
            }
        }
        boolean retval = true;

        // Check for child options for which there is no parser
        for (final String prefix : childOptionArguments.keySet()) {
            if (!childOptionsMap.containsKey(prefix)) {
                messageStream.println("ERROR: Option prefix '" + prefix + "' is not valid.");
                retval = false;
            }
        }

        try {
            // Propagate options from this parser to child parsers
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                // Handling collection value propagation is confusing, so just don't do it.
                if (optionDefinition.isCollection) continue;
                final Object value = optionDefinition.field.get(callerOptions);
                if (value == null) continue;
                for (final CommandLineParser childParser : childOptionsMap.values()) {
                    maybePropagateValueToChild(childParser, optionDefinition, value);
                }
            }
        } catch (final IllegalAccessException e) {
            throw new RuntimeException("Should never happen", e);
        }

        for (final Map.Entry<String, CommandLineParser> entry : childOptionsMap.entrySet()) {
            final String prefix = entry.getKey();
            final CommandLineParser childParser = entry.getValue();
            childParser.messageStream = this.messageStream;
            final Collection<ChildOptionArg> childOptionArgs = this.childOptionArguments.get(prefix);
            if (childOptionArgs != null) {
                for (final ChildOptionArg arg : childOptionArgs) {
                    childParser.parseOption(arg.name, arg.value, arg.fromFile, arg.precedenceSet);
                }
            }
            if (!childParser.checkNumArguments()) {
                retval = false;
            }
            if (!childParser.parseChildOptions()) {
                retval = false;
            }
            this.commandLine += " " + childParser.getCommandLine();
        }

        return retval;
    }

    /**
     * Propagate value from parent to child if appropriate to do so.
     */
    private void maybePropagateValueToChild(final CommandLineParser childParser,
                                            final OptionDefinition optionDefinition,
                                            final Object value) {
        try {
            final OptionDefinition childOptionDefinition = childParser.optionMap.get(optionDefinition.name);
            if (childOptionDefinition != null) {
                final Object childValue = childOptionDefinition.field.get(childParser.callerOptions);
                if (childValue == null || optionDefinition.hasBeenSet) {
                    childOptionDefinition.field.set(childParser.callerOptions, value);
                    childOptionDefinition.hasBeenSetFromParent = true;
                    childOptionDefinition.hasBeenSetFromOptionsFile = optionDefinition.hasBeenSetFromOptionsFile;
                }
            }
        } catch (final IllegalAccessException e) {
            throw new RuntimeException("Should never happen", e);
        }
    }

    /**
     * The commandline used to run this program, including any default args that
     * weren't necessarily specified. This is used for logging and debugging.
     * <p/>
     * NOTE: {@link #parseOptions(PrintStream, String[])} must be called before
     * calling this method.
     *
     * @return The commandline, or null if {@link #parseOptions(PrintStream, String[])}
     * hasn't yet been called, or didn't complete successfully.
     */
    public String getCommandLine() { return commandLine; }

    /**
     * This method is only needed when calling one of the public methods that doesn't take a messageStream argument.
     */
    public void setMessageStream(final PrintStream messageStream) {
        this.messageStream = messageStream;
    }
}
