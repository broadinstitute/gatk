package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * A CommandLinePluginDescriptor for ReadFilter plugins
 */
public class GATKReadFilterPluginDescriptor extends CommandLinePluginDescriptor<ReadFilter> {

    protected transient Logger logger = LogManager.getLogger(this.getClass());

    private void readObject(java.io.ObjectInputStream in)
            throws IOException, ClassNotFoundException {
        in.defaultReadObject();
        logger = LogManager.getLogger(this.getClass()); // Logger is not serializable (even by Kryo)
    }

    private static final String pluginPackageName = "org.broadinstitute.hellbender.engine.filters";
    private static final Class<?> pluginBaseClass = org.broadinstitute.hellbender.engine.filters.ReadFilter.class;

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

    // Map of read filter (simple) class names to the corresponding discovered plugin instance
    private Map<String, ReadFilter> allDiscoveredReadFilters = new HashMap<>();

    // List of default filters in the order they were specified by the tool
    private List<String> toolDefaultReadFilterNamesInOrder = new ArrayList<>();

    // Map of read filter (simple) class names to the corresponding default plugin instance
    private Map<String, ReadFilter> toolDefaultReadFilters = new HashMap<>();

    // Set of dependent args for which we've seen values (requires predecessor)
    private Set<String> requiredPredecessors = new HashSet<>();

    /**
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKReadFilterPluginDescriptor(final List<ReadFilter> toolDefaultFilters) {
        if (null != toolDefaultFilters) {
            toolDefaultFilters.forEach(f -> {
                final Class<? extends ReadFilter> rfClass = f.getClass();
                // anonymous classes have a 0-length simple name, and thus cannot be accessed or
                // controlled by the user via the command line, but they should still be valid
                // as default filters, so use the full name to ensure that their map entries
                // don't clobber each other
                String className = rfClass.getSimpleName();
                if (className.length() == 0) {
                    className = rfClass.getName();
                }
                // we want to remember the order in which these were provided, and also keep the
                // actual instances in case they have any additional state provided by the tool
                // when they were created
                toolDefaultReadFilterNamesInOrder.add(className);
                toolDefaultReadFilters.put(className, f);
            });
        }
    }

    /////////////////////////////////////////////////////////
    // GATKCommandLinePluginDescriptor implementation methods

    /**
     * Return a display name to identify this plugin to the user
     * @return A short user-friendly name for this plugin.
     */
    @Override
    public String getDisplayName() { return StandardArgumentDefinitions.READ_FILTER_LONG_NAME; }

    /**
     * @return the class object for the base class of all plugins managed by this descriptor
     */
    @Override
    public Class<?> getPluginClass() {return pluginBaseClass;}

    /**
     * A list of package names which will be searched for plugins managed by the descriptor.
     * @return
     */
    @Override
    public List<String> getPackageNames() {return Collections.singletonList(pluginPackageName);};

    @Override
    public Predicate<Class<?>> getClassFilter() {
        return c -> {
            // don't use the ReadFilter base class, it's inner classes, the CountingReadFilter,
            // or the unit tests
            return !c.getName().equals(this.getPluginClass().getName()) &&
                    !c.getName().startsWith(CountingReadFilter.class.getName()) &&
                    !c.getName().startsWith(this.getPluginClass().getName() + "$") &&
                    !c.getName().contains("UnitTest$");
        };
    }

    // Instantiate a new ReadFilter derived object and save it in the list
    @Override
    public Object getInstance(final Class<?> pluggableClass) throws IllegalAccessException, InstantiationException {
        ReadFilter readFilter = null;
        final String simpleName = pluggableClass.getSimpleName();

        if (allDiscoveredReadFilters.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                            "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            allDiscoveredReadFilters.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultReadFilters.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default filters;
            // use the default instance as the target for command line argument values
            // rather than creating a new one, in case it has state provided by the tool
            readFilter = toolDefaultReadFilters.get(simpleName);
        } else {
            readFilter = (ReadFilter) pluggableClass.newInstance();
        }

        // Add all filters to the allDiscoveredReadFilters list, even if the instance came from the
        // tool defaults list (we want the actual instances to be shared to preserve state)
        allDiscoveredReadFilters.put(simpleName, readFilter);
        return readFilter;
    }

    @Override
    public boolean isDependentArgumentAllowed(final Class<?> predecessorClass) {
        // Make sure the predecessor for a dependent argument was either specified on the command line or
        // is a tool default, otherwise reject it.
        // NOTE: This method is called by the CLP during parsing at the time the dependet argument is seen
        // on the command line. Even if this check passes at the time this method is called, its possible
        // for the user to subsequently disable the required predecessor. That case is caught during final
        // validation done by the validateArguments method.
        String predecessorName = predecessorClass.getSimpleName();
        boolean isAllowed = userEnabledReadFilterNames.contains(predecessorName)
                || (toolDefaultReadFilters.get(predecessorName) != null);
        if (isAllowed) {
            // Keep track of the ones we allow so we can validate later that they weren't subsequently disabled
            requiredPredecessors.add(predecessorName);
        }
        return isAllowed;
    }

    /**
     * Pass back the list of ReadFilter instances that were actually seen on the command line in the same
     * order they were specified. Its possible for this to return a filter that was originally included
     * in the list of tool defaults only in the case where the user also specifies it on the command line.
     *
     * NOTE: this method is somewhat misnamed in that it doesn't return ALL instances since it leaves out
     * default filters (Except as noted above). The refactored interface in Barclay renames this method and
     * changes it's contract. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<ReadFilter> getAllInstances() {
        final ArrayList<ReadFilter> filters = new ArrayList<>(userEnabledReadFilterNames.size());
        userEnabledReadFilterNames.forEach(s -> {
            ReadFilter rf = allDiscoveredReadFilters.get(s);
            filters.add(rf);
        });
        return filters;
    }

    /**
     * Get the list of default plugins used for this instance of this descriptor. Used for help/doc generation.
     *
     * NOTE: this method does not account for disabled default filters and just return ALL default instances.
     * The refactored interface in Barclay changes it's contract to allows returning a list with only 'enabled' default
     * instances. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<Object> getDefaultInstances() { return new ArrayList<>(toolDefaultReadFilters.values()); }

    /**
     * Return the class representing the instance of the plugin specified by {@code pluginName}
     * @param pluginName Name of the plugin requested
     * @return Class object for the plugin instance requested
     */
    @Override
    public Class<?> getClassForInstance(final String pluginName) {
        return allDiscoveredReadFilters.get(pluginName).getClass();
    }

    /**
     * Return the allowed values for readFilterNames/disableReadFilter for use by the help system.
     * @param longArgName long name of the argument for which help is requested
     * @return
     */
    @Override
    public Set<String> getAllowedValuesForDescriptorArgument(final String longArgName) {
        if (longArgName.equals(StandardArgumentDefinitions.READ_FILTER_LONG_NAME)) {
            return allDiscoveredReadFilters.keySet();
        }
        if (longArgName.equals(StandardArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME)) {
            return toolDefaultReadFilters.keySet();
        }
        throw new IllegalArgumentException("Allowed values request for unrecognized string argument: " + longArgName);
    }

    /**
     * Validate the list of arguments and reduce the list of read filters to those
     * actually seen on the command line. This is called by the command line parser
     * after all arguments have been parsed.
     */
    @Override
    public void validateArguments() {
        // throw if a filter is *enabled* more than once by the user
        final Set<String> duplicateUserEnabledFilterNames = Utils.getDuplicatedItems(userEnabledReadFilterNames);
        if (!duplicateUserEnabledFilterNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The read filter(s) are enabled more than once: %s",
                            Utils.join(", ", duplicateUserEnabledFilterNames)));
        }

        // throw if a filter is *disabled* more than once by the user
        final Set<String> duplicateDisabledUserFilterNames = Utils.getDuplicatedItems(userDisabledReadFilterNames);
        if (!duplicateDisabledUserFilterNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The read filter(s) are disabled more than once: %s",
                            Utils.join(", ", duplicateDisabledUserFilterNames)));
        }

        // throw if a filter is both enabled *and* disabled by the user
        final Set<String> enabledAndDisabled = new HashSet<>(userEnabledReadFilterNames);
        enabledAndDisabled.retainAll(userDisabledReadFilterNames);
        if (!enabledAndDisabled.isEmpty()) {
            final String badFiltersList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The read filter(s): %s are both enabled and disabled", badFiltersList));
        }

        // throw if a disabled filter doesn't exist; warn if it wasn't enabled by the tool in the first place
        userDisabledReadFilterNames.forEach(s -> {
            if (!allDiscoveredReadFilters.containsKey(s)) {
                throw new CommandLineException.BadArgumentValue(String.format("Disabled filter (%s) does not exist", s));
            } else if (!toolDefaultReadFilters.containsKey(s)) {
                logger.warn(String.format("Disabled filter (%s) is not enabled by this tool", s));
            }
        });

        // warn if a filter is both default and enabled by the user
        final Set<String> redundant = new HashSet<>(toolDefaultReadFilters.keySet());
        redundant.retainAll(userEnabledReadFilterNames);
        redundant.forEach(
            s -> {
                logger.warn(String.format("Redundant enabled filter (%s) is enabled for this tool by default", s));
            });

        // Throw if args were specified for a filter that was also disabled, or that was not enabled by the
        // tool by default.
        //
        // Note that this is also checked during command line argument parsing, but needs to be checked again
        // here. Whenever the command line parser sees a dependent argument on the command line, it delegates
        // back to the descriptor's isDependentArgumentAllowed method to allow it to validate that the predecessor
        // for that dependent argument has been supplied, either by a default read filter, or by an explicitly
        // enabled read filter. However, its possible for the user to subsequently try to disable that
        // predecessor, which is what we want to catch here.
        //
        userDisabledReadFilterNames.forEach(s -> {
            if (requiredPredecessors.contains(s)) {
                String message = String.format("Values were supplied for (%s) that is also disabled", s);
                if (toolDefaultReadFilters.containsKey(s)) {
                    // NOTE: https://github.com/broadinstitute/barclay/issues/23
                    // This is a special case to work around the issue where we can't really tell if the
                    // predecessor was added as a result of a user-provided value, or a default value. The
                    // CLP doesn't distinguish, so we only warn here for now.
                    logger.warn(message);
                } else {
                    throw new CommandLineException(message);
                }
            }
        });

        // throw if a filter name was specified that has no corresponding instance
        final Map<String, ReadFilter> requestedReadFilters = new HashMap<>();
        userEnabledReadFilterNames.forEach(s -> {
            ReadFilter trf = allDiscoveredReadFilters.get(s);
            if (null == trf) {
                if (!toolDefaultReadFilters.containsKey(s)) {
                    throw new CommandLineException("Unrecognized read filter name: " + s);
                }
            } else {
                requestedReadFilters.put(s, trf);
            }
        });
    }

    /////////////////////////////////////////////////////////
    // ReadFilter plugin-specific helper methods

    /**
     * Determine if a particular ReadFilter was disabled on the command line, either directly of by disabling all
     * tool defaults.
     * @param filterName name of the filter to query.
     * @return {@code true} if the name appears in the list of disabled filters, or is a tool default not provided by
     * the user and all tool defaults are disabled; {@code false} otherwise.
     */
    public boolean isDisabledFilter(final String filterName) {
        return userDisabledReadFilterNames.contains(filterName)
                || (disableToolDefaultReadFilters && !userEnabledReadFilterNames.contains(filterName));
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read filter instances
     * @return Single merged read filter.
     */
    public ReadFilter getMergedReadFilter(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadFilter(
                samHeader,
                ReadFilter::fromList
        );
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read filter instances
     * @return Single merged counting read filter.
     */
    public CountingReadFilter getMergedCountingReadFilter(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadFilter(
                samHeader,
                CountingReadFilter::fromList
        );
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader a SAMFileHeader to initialize read filter instances. May not be null.
     * @param aggregateFunction function to use to merge ReadFilters, usually ReadFilter::fromList. The function
     *                          must return the ALLOW_ALL_READS filter wrapped in the appropriate type when passed
     *                          a null or empty list.
     * @param <T> extends ReadFilter, type returned by the wrapperFunction
     * @return Single merged read filter.
     */
    public <T extends ReadFilter> T getMergedReadFilter(
            final SAMFileHeader samHeader,
            final BiFunction<List<ReadFilter>, SAMFileHeader, T> aggregateFunction) {

        Utils.nonNull(samHeader);
        Utils.nonNull(aggregateFunction);

        // start with the tool's default filters in the order they were specified, and remove any that were disabled
        // on the command line
        // if --disableToolDefaultReadFilters is specified, just initialize an empty list with initial capacity of user filters
        final List<ReadFilter> finalFilters =
                disableToolDefaultReadFilters ?
                        new ArrayList<>(userEnabledReadFilterNames.size()) :
                        toolDefaultReadFilterNamesInOrder
                                .stream()
                                .filter(s -> !isDisabledFilter(s))
                                .map(s -> toolDefaultReadFilters.get(s))
                                .collect(Collectors.toList());

        // now add in any additional filters enabled on the command line (preserving order)
        final List<ReadFilter> clFilters = getAllInstances();
        if (clFilters != null) {
            clFilters.stream()
                    .filter(f -> !finalFilters.contains(f)) // remove redundant filters
                    .forEach(f -> finalFilters.add(f));
        }

        return aggregateFunction.apply(finalFilters, samHeader);
    }

}
