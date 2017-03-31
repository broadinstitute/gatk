package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.util.*;

/**
 * Default {@link GATKReadFilterArgumentCollection} applied in GATK for optional read filters in the command line.
 * It contains arguments that allow the user to:
 *
 * - Provide a list of read filters to apply.
 * - Disable some and/or all read filters.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DefaultGATKReadFilterArgumentCollection extends GATKReadFilterArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.READ_FILTER_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_FILTER_SHORT_NAME,
            doc="Read filters to be applied before analysis", optional=true, common = true)
    public final List<String> userEnabledReadFilterNames = new ArrayList<>(); // preserve order

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME,
            shortName = StandardArgumentDefinitions.DISABLE_READ_FILTER_SHORT_NAME,
            doc="Read filters to be disabled before analysis", optional=true, common = true)
    public final List<String> userDisabledReadFilterNames = new ArrayList<>();

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
            shortName = StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
            doc = "Disable all tool default read filters", common = true, optional = true)
    public boolean disableToolDefaultReadFilters = false;

    /** Returns the list with the read filters provided by the user, preserving the order. */
    @Override
    public List<String> getUserEnabledReadFilterNames() {
        return userEnabledReadFilterNames;
    }

    /** Returns the set of filters disabled by the user. */
    @Override
    public List<String> getUserDisabledReadFilterNames() {
        return userDisabledReadFilterNames;
    }

    /** {@inheritDoc}. */
    @Override
    public boolean getDisableToolDefaultReadFilters() {
        return disableToolDefaultReadFilters;
    }

}
