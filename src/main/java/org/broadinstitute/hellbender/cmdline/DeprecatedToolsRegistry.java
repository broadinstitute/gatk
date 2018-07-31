package org.broadinstitute.hellbender.cmdline;

import org.apache.commons.lang3.tuple.Pair;

import java.util.HashMap;
import java.util.Map;

/**
 * When a tool is removed from GATK (after having been @Deprecated for a suitable period), an entry should
 * be added to this list to issue a message when the user tries to run that tool.
 *
 * NOTE: Picard tools should be listed here as well, since by definition such tools will not be found in
 * the Picard jar.
 */
public class DeprecatedToolsRegistry {

    // Mapping from tool name to string describing the major version number where the tool first disappeared and
    // optional recommended alternatives
    private static Map<String, Pair<String, String>> deprecatedTools = new HashMap<>();

    static {
        // Indicate version in which the tool disappeared, and recommended replacement in parentheses if applicable
        deprecatedTools.put("IndelRealigner", Pair.of("4.0.0.0", "Please use GATK3 to run this tool"));
        deprecatedTools.put("RealignerTargetCreator", Pair.of("4.0.0.0", "Please use GATK3 to run this tool"));
    }

    /**
     * Utility method to pull up the version number at which a tool was deprecated and the suggested replacement, if any
     *
     * @param toolName   the tool class name (not the full package) to check
     */
    public static String getToolDeprecationInfo(final String toolName) {
        return deprecatedTools.containsKey(toolName) ?
                String.format("%s is no longer included in GATK as of version %s. %s",
                        toolName,
                        deprecatedTools.get(toolName).getLeft(),
                        deprecatedTools.get(toolName).getRight()
                ) :
                null;
    }

}
