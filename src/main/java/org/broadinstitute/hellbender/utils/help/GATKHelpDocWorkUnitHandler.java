package org.broadinstitute.hellbender.utils.help;

import com.netflix.servo.util.VisibleForTesting;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.help.DefaultDocWorkUnitHandler;
import org.broadinstitute.barclay.help.DocWorkUnit;

import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The GATK Documentation work unit handler class that is the companion to GATKHelpDoclet.
 *
 * NOTE: Methods in this class are intended to be called by Gradle/Javadoc only, and should not be called
 * by methods that are used by the GATK runtime, as this class assumes a dependency on com.sun.javadoc classes
 * which may not be present.
 */
public class GATKHelpDocWorkUnitHandler extends DefaultDocWorkUnitHandler {

    private final static String GATK_JAVADOC_TAG_PREFIX = "GATK"; // prefix for custom javadoc tags used by GATK

    private final static String GATK_FREEMARKER_TEMPLATE_NAME = "generic.template.html";

    private final static String GATK_COMMAND_PROGRAM_NAME = "gatk";
    
    private final static String GATK_JAVA_OPTIONS_ARGNAME = "--javaOptions";
    
    /**
     * Pattern to detect the presence of code example blocks that might contain 
     * Picard command examples. 
     * <p>
     * These would consists {@code <pre> ... </pre>} blocks. 
     * </p>
     */
    private final Pattern CODE_EXAMPLE_BLOCK = Pattern.compile("<pre\\s*>(.|\\s)*?<\\/pre\\s*>",
            Pattern.CASE_INSENSITIVE);

    /**
     * Pattern to detect the presence of a Picard command start.
     * <p>
     * This would look like {@code java -jar picard.jar ToolName}.
     * </p>
     */
    private final Pattern PICARD_CODE_EXAMPLE_COMMAND_START = Pattern.compile("java\\s*(.*?)\\s*-jar\\s*\\S*?picard\\S*?\\.jar",
            Pattern.CASE_INSENSITIVE | Pattern.MULTILINE);

    /**
     * Pattern that matches a Picard tool name.
     * <p>
     * Any sequence of non-space characters does.
     * </p>
     */
    private final Pattern PICARD_CODE_EXAMPLE_COMMAND_TOOLNAME = Pattern.compile("\\S+");

    /**
     * Pattern to match a Picard command argument value pair.
     * <p>
     * These would look like {@code arg-name=value}.
     * </p>
     * NOTE: Currently we don't support spaces in neither the argument name nor the value; I guess these are possible using 
     * quotes. 
     */
    private final Pattern PICARD_CODE_EXAMPLE_COMMAND_ARGUMENT_VALUE = Pattern.compile("(\\S+?)\\s*=\\s*(\\S+)");

    /**
     * Pattern to match {@code "true"} argument value for flag-typed arguments.
     * <p>
     * Suppored alternatives: "t", "T", "true", "TRUE".
     * </p> 
     */
    private final Pattern PICARD_FLAG_TRUE_VALUE = Pattern.compile("\\s*t(rue)?\\s*", Pattern.CASE_INSENSITIVE);

    /**
     * Detects the end of the Picard command line.
     * <p>
     * This is any of: {@code </pre>} or an empty line or a non-empty line where
     * the last non-space character is not the slash "\".
     * </p>
     */
    private final Pattern PICARD_CODE_EXAMPLE_COMMAND_END = Pattern.compile("(<\\/pre\\s*>)|(^\\s*$)|((?<=[^ \\t\\\\])\\s*$)",
            Pattern.MULTILINE | Pattern.CASE_INSENSITIVE);

    /**
     * {@inheritDoc}
     */
    public GATKHelpDocWorkUnitHandler(final GATKHelpDoclet doclet) {
        super(doclet);
    }
    
    /**
     * @return Prefix for custom GATK tags that should be lifted from the javadoc and stored in the
     * FreeMarker map. These will be available in the template returned by {@link #getTemplateName}.
     */
    @Override
    protected String getTagFilterPrefix() { return GATK_JAVADOC_TAG_PREFIX; }

    /**
     * @param workUnit the classdoc object being processed
     * @return the name of a the freemarker template to be used for the class being documented.
     * Must reside in the folder passed to the Barclay Doclet via the "-settings-dir" parameter to
     * Javadoc.
     */
    @Override
    public String getTemplateName(final DocWorkUnit workUnit) { return GATK_FREEMARKER_TEMPLATE_NAME; }

    /**
     * {@inheritDoc}
     * <p>
     * Additionally, it convert any Picard command code example look-a-like into the equivalent GATK command call:
     * <p>
     * For example:
     * <pre>java -jar picard.jar ToolName ARG1=VAL1 ARG2=VAL2 ...</pre>
     * becomes
     * <pre>gatk ToolName --ARG1 VAL1 --ARG2 VAL2</pre>
     * </p>
     */
    @Override
    public void processWorkUnit(
            final DocWorkUnit workUnit,
            final List<Map<String, String>> featureMaps,
            final List<Map<String, String>> groupMaps) {
        super.processWorkUnit(workUnit, featureMaps, groupMaps);
        final CharSequence description = new StringBuilder((CharSequence) workUnit.getProperty("description"));
        workUnit.setProperty("description", translatePicardCodeBlocks(description));
    }

    /**
     * Translates Picard calls into GATK calls in all code example blocks in the input character
     * sequence. 
     * @param input the input character-sequence.
     * @return the input sequence with the command translations if any.
     */
    @VisibleForTesting
    String translatePicardCodeBlocks(final CharSequence input) {
    		final Matcher matcher = CODE_EXAMPLE_BLOCK.matcher(input);
        if (!matcher.find()) { // we cannot find a Picard code example.
            return input.toString();
        } else {
            return translatePicardCodeBlocks(input, matcher);     
        }
    }
    
    private String translatePicardCodeBlocks(final CharSequence description, final Matcher matcher) {
        final StringBuilder result = new StringBuilder(description.length() << 1);
        int lastMatchedOffset = 0;
        do {
            result.append(description.subSequence(lastMatchedOffset, matcher.start()));
            translatePicardCodeBlock(description, matcher.start(), matcher.end(), result);
            lastMatchedOffset = matcher.end();
        } while(matcher.find());
        result.append(description.subSequence(lastMatchedOffset, description.length()));
        return result.toString();
    }

    /**
     * Translates the content of a block given the character sequence that contains it,
     * the start and end of its content.
     * @param block the character sequence that contains the block.
     * @param start position of the block content to translate (inclusive).
     * @param end position of the block content to translate (exclusive).
     * @param result the result of translating the block content.
     */
    @VisibleForTesting
    void translatePicardCodeBlock(final CharSequence block, final int start, final int end, final StringBuilder result) {
        final Matcher startMatcher = PICARD_CODE_EXAMPLE_COMMAND_START.matcher(block).region(start, end);
        final Matcher endMatcher = PICARD_CODE_EXAMPLE_COMMAND_END.matcher(block).region(start, end);
        int lastMatcherOffset = start;
        while (startMatcher.find()) {
            result.append(block.subSequence(lastMatcherOffset, startMatcher.start()));
            endMatcher.region(startMatcher.start(), end);
            if (!endMatcher.find()) { // we cannot find and end!? (in practice this hould never happen though.
                lastMatcherOffset = startMatcher.start();
            } else {
                final String javaOptions = startMatcher.group(1).trim();
                final CharSequence toolAndArguments = block.subSequence(startMatcher.end(), endMatcher.start());
                result.append(GATK_COMMAND_PROGRAM_NAME);
                if (!javaOptions.isEmpty()) {
                    result.append(' ').append(GATK_JAVA_OPTIONS_ARGNAME)
                          .append(" '").append(javaOptions.replace("'","\\'")).append('\'');
                }
                translatePicardToolAndArguments(toolAndArguments, result);
                result.append(endMatcher.group());
                lastMatcherOffset = endMatcher.end();
                startMatcher.region(endMatcher.end(), end);
            }
        }
        result.append(block.subSequence(lastMatcherOffset, end));
    }

    private void translatePicardToolAndArguments(final CharSequence toolAndArguments, final StringBuilder buffer) {
        final Matcher toolNameMatcher = PICARD_CODE_EXAMPLE_COMMAND_TOOLNAME.matcher(toolAndArguments);
        if (!toolNameMatcher.find()) {
            buffer.append(toolAndArguments);
        } else {
            buffer.append(toolAndArguments.subSequence(0, toolNameMatcher.start()));
            buffer.append(toolNameMatcher.group());
            final String toolName = toolNameMatcher.group();
            final CommandLineProgramIntrospector introspector = CommandLineProgramIntrospector.of(toolName);
            if (introspector == null) {
                printWarning("unknown picard tool-name: " + toolName);
            } else {
                printNotice("found picard code command example on tool: " + toolName);
            }
            final Matcher argumentMatcher = PICARD_CODE_EXAMPLE_COMMAND_ARGUMENT_VALUE.matcher(
                    toolAndArguments).region(toolNameMatcher.end(), toolAndArguments.length());
            int lastMatcherOffset = toolNameMatcher.end();
            while (argumentMatcher.find()) {
                buffer.append(toolAndArguments.subSequence(lastMatcherOffset, argumentMatcher.start()));
                lastMatcherOffset = argumentMatcher.end();
                final String name = argumentMatcher.group(1);
                final String value = argumentMatcher.group(2);
                final CommandLineArgumentParser.ArgumentDefinition definition = introspector != null ? introspector.getArgumentDefinition(name) : null;
                if (definition == null && introspector != null) {
                		printWarning("unknown argument named '" + name + "' for picard tool '" + toolName + "'");
                }
                final String prefix;
                if (definition != null) { // we can find out whether is long or short name using the descriptor.
                    prefix = definition.getLongName().equals(name) ? "--" : "-";
                } else if (name.contains("_")) { // heuristic, if it has a _ then is long.
                    prefix = "--";
                } else { // we apply a heuristic where 3 or less letter long names are short, 4 or more are long.
                    prefix = name.length() <= 3 ? "-" : "--";
                }
                if (definition != null && definition.isFlag() && PICARD_FLAG_TRUE_VALUE.matcher(value).matches()) {
                    buffer.append(prefix).append(name); // we compress ARG=true to -ARG or --ARG for flags.
                } else {
                    buffer.append(prefix).append(name).append(' ').append(value);
                }
            }
            buffer.append(toolAndArguments.subSequence(lastMatcherOffset, toolAndArguments.length()));
        }
    }

    @SuppressWarnings("unused")
    private void printError(final String message) {
        ((GATKHelpDoclet)getDoclet()).printError(message);
    }

    private void printWarning(final String message) {
        ((GATKHelpDoclet)getDoclet()).printWarning(message);
    }

    private void printNotice(final String message) {
        ((GATKHelpDoclet)getDoclet()).printNotice(message);
    }
}
