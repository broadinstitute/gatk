package org.broadinstitute.hellbender.cmdline;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * This class is a static helper for implementing 'mode arguments' by tools. A mode argument is defined as
 * an argument which controls the setting of multiple (other) arguments to preset values, thereby resulting
 * in the tool being used in a specific 'mode'. In this sense, a mode is simply a macro of sorts. A mode
 * argument might be binary (in which case the associated argument and their values will be set when it is true)
 * or be an enum (in which case a different set of arguments and value be set according to the selected value).
 *
 * The detection of the mode (i.e. the implementation of the mode argument) is the tool's responsibility. Following
 * that, a public method of this class, setArgValues, can be used to install the associated argument values
 * in the command parser, which will propagate them into the associate tools instance fields.
 *
 * The setting of the arguments is done such that argument already specified on the command line are not
 * overwritten. This allows for mode refinement - i.e. the setting of the mode while turing some of its arguments
 * to specific values.
 *
 * NOTE: This utility relies on the argument engine having already been executed. Consequently, the best place in the
 *       standard GATKTool would be in the #customCommandLineArgumentValidation() or at a similar stage since that gets
 *       executed after Barclay parsing is complete but before the tool has done anything with the filled arguments.
 */

public final class ModeArgumentUtils {

    protected static final Logger logger = LogManager.getLogger(ModeArgumentUtils.class);

    /**
     * Set a group of  arguments (and associated fields) to the specified values - thereby implementing a mode
     *
     * The method does not overwrite arguments already set on the command line (as indicated by the parser)
     * Additionally, the method emits a warning message indicating which argument were set by the mode and
     * to which values.
     *
     * NOTE: This does not validate that the specified mode was actually set, that responsibility is left to tool authors.
     *
     * @param parser - commandline parser associated with the tool. This should be an instance of CommandLineArgumentParser
     * @param argValues - an array of arguments and values to potentially set. The order of elements in the array
     *                  is arg0,value0,arg1,value1 ...
     * @param modeName - the name of the mode being set. This is used in textual message, such as the warning
     *                 issued to notify the user of the changed argument and not for validation that the argument was set.
     */
    public static void setArgValues(final CommandLineParser parser, final String[] argValues, final String modeName) {
        final Map<String, String>  modifiedArgs = new LinkedHashMap<>();

        for ( int i = 0 ; i < argValues.length ; i += 2 ) {
            if ( !hasBeenSet(parser, argValues[i]) ) {
                String parserMessage = setValue(parser, argValues[i], argValues[i+1]);

                if ( StringUtils.isEmpty(parserMessage) ) {
                    modifiedArgs.put(argValues[i], argValues[i + 1]);
                } else {
                    modifiedArgs.put(argValues[i], argValues[i + 1] + " (" + parserMessage + ")");
                }
            } else {
                logger.info("parameter not set by the '" + modeName + "' argument mode, as it was already set on the command line: " + argValues[i]);
            }
        }

        logModeNotice(modifiedArgs, modeName);
    }

    private static boolean hasBeenSet(final CommandLineParser parser, final String alias) {

        if ( parser instanceof CommandLineArgumentParser ) {
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(alias);

            return (namedArg != null) ? namedArg.getHasBeenSet() : false;
        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }
    }

    private static String setValue(final CommandLineParser parser, final String alias, final String value) {
        if ( parser instanceof CommandLineArgumentParser ) {
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(alias);
            if ( namedArg == null ) {
                throw new IllegalArgumentException("alias not found: " + alias);
            }

            PrintStream         ps = new PrintStream(new ByteArrayOutputStream());
            List<String>        values = Arrays.asList(value);
            namedArg.setArgumentValues((CommandLineArgumentParser)parser, ps, values);
            return ps.toString();
        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }
    }

    private static void logModeNotice(final Map<String, String> modifiedArgs, final String modeName) {
        logger.warn("*************************************************************************");
        logger.warn(String.format("* %-69s *", "--" + modeName + " was enabled"));
        logger.warn("* The following arguments have had their inputs overwritten:            *");
        modifiedArgs.forEach((name, value) -> {
            logger.warn(String.format("* %-69s *", "--" + name + " " + value));
        });
        logger.warn("*                                                                       *");
        logger.warn("* If you would like to run this mode with different inputs for any      *");
        logger.warn("* of the above arguments please manually construct the command or       *");
        logger.warn("* add your specific inputs after the mode argument. This mode           *");
        logger.warn("* will not override inputs explicitly provided.                         *");
        logger.warn("*************************************************************************");
    }
}
