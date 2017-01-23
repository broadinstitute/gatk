package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;

import java.util.*;

/**
 * A CommandLinePluginDescriptor for ReadFilter plugins
 */
public final class GATKReadFilterPluginDescriptor extends AbstractReadFilterPluginDescriptor {

    private static final String pluginPackageName = "org.broadinstitute.hellbender.engine.filters";

    @Argument(fullName = StandardArgumentDefinitions.READ_FILTER_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_FILTER_SHORT_NAME,
            doc="Read filters to be applied before analysis", optional=true)
    public final List<String> userReadFilterNames = new ArrayList<>(); // preserve order

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME,
            shortName = StandardArgumentDefinitions.DISABLE_READ_FILTER_SHORT_NAME,
            doc="Read filters to be disabled before analysis", optional=true)
    public final Set<String> disableFilters = new HashSet<>();

    @Argument(fullName = "disableAllReadFilters",
            shortName = "disableAllReadFilters",
            doc = "Disable all read filters", common = false, optional = true)
    public boolean disableAllReadFilters = false;

    /**
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKReadFilterPluginDescriptor(final List<ReadFilter> toolDefaultFilters) {
        super(toolDefaultFilters);
    }

    /** Only ReadFilters in {@link #pluginPackageName} are allowed. */
    @Override
    public List<String> getPackageNames() {return Collections.singletonList(pluginPackageName);};

    @Override
    public List<String> getCommandLineReadFilters() {
        return userReadFilterNames;
    }

    @Override
    public Set<String> getDisableFilters() {
        return disableFilters;
    }

    @Override
    public boolean shouldDisableAllReadFilters() {
        return disableAllReadFilters;
    }
}
