package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.util.*;

/**
 * {@link GATKReadFilterArgumentCollection} for optional read filters in the command line. It allows:
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

    /** Returns a list with the read filters provided by the user, preserving the order. */
    @Override
    public Collection<String> getUserEnabledReadFilterNames() {
        return userEnabledReadFilterNames;
    }

    /** Returns a set of filters disabled by the user. */
    @Override
    public Collection<String> getUserDisabledReadFilterNames() {
        return userDisabledReadFilterNames;
    }

    /** {@inheritDoc}. */
    @Override
    public boolean disableToolDefaultReadFilters() {
        return disableToolDefaultReadFilters;
    }

}
