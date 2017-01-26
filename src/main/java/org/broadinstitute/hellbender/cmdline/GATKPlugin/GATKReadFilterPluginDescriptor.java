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

    // Map of read filter (simple) class names to the corresponding discovered plugin instance
    private Map<String, ReadFilter> readFilters = new HashMap<>();

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

        if (readFilters.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                            "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            readFilters.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultReadFilters.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default filters;
            // use the default instance as the target for command line argument values
            // rather than creating a new one in case it has state provided by the tool
            readFilter = toolDefaultReadFilters.get(simpleName);
        } else {
            readFilter = (ReadFilter) pluggableClass.newInstance();
            readFilters.put(simpleName, readFilter);
        }
        return readFilter;
    }

    @Override
    public boolean isDependentArgumentAllowed(final Class<?> dependentClass) {
        // make sure the predecessor for this dependent class was either specified
        // on the command line or is a tool default, otherwise reject it
        String predecessorName = dependentClass.getSimpleName();
        boolean isAllowed = userReadFilterNames.contains(predecessorName)
                || (toolDefaultReadFilters.get(predecessorName) != null);
        if (isAllowed) {
            // keep track of the ones we allow so we can validate later that they
            // weren't subsequently disabled
            requiredPredecessors.add(predecessorName);
        }
        return isAllowed;
    }

    /**
     * Pass back the list of ReadFilter instances that were actually seen on the
     * command line in the same order they were specified. This list does not
     * include the tool defaults.
     */
    @Override
    public List<ReadFilter> getAllInstances() {
        // Add the instances in the order they were specified on the command line
        // (use the order of userReadFilterNames list).
        //
        // NOTE: it's possible for the userReadFilterNames list to contain one or more
        // names for which there are no corresponding instances in the readFilters list.
        // This happens when the user specifies a filter name on the command line that's
        // already included in the toolDefault list, since in that case the descriptor
        // uses the tool-supplied instance and doesn't add a separate one to the
        // readFilters list, but the name from the command line still appears in
        // userReadFilterNames. In that case, we don't include the tool's instance in the
        // list returned by this method since it will be merged in later by the merge method.
        final ArrayList<ReadFilter> filters = new ArrayList<>(userReadFilterNames.size());
        userReadFilterNames.forEach(s -> {
            ReadFilter rf = readFilters.get(s);
            if (rf != null) {
                filters.add(rf);
            }
        });
        return filters;
    }

    // Return the allowable values for readFilterNames/disableReadFilter
    @Override
    public Set<String> getAllowedValuesForDescriptorArgument(final String longArgName) {
        if (longArgName.equals(StandardArgumentDefinitions.READ_FILTER_LONG_NAME)) {
            return readFilters.keySet();
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
        // throw if any filter is both enabled *and* disabled by the user
        final Set<String> enabledAndDisabled = new HashSet<>(userReadFilterNames);
        enabledAndDisabled.retainAll(disableFilters);
        if (!enabledAndDisabled.isEmpty()) {
            final String badFiltersList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The read filter(s): %s are both enabled and disabled", badFiltersList));
        }

        // warn if a disabled filter wasn't enabled by the tool in the first place
        disableFilters.forEach(s -> {
            if (!toolDefaultReadFilters.containsKey(s)) {
                logger.warn(String.format("Disabled filter (%s) is not enabled by this tool", s));
            }
        });

        // warn on redundant enabling of filters already enabled by default
        final Set<String> redundant = new HashSet<>(toolDefaultReadFilters.keySet());
        redundant.retainAll(userReadFilterNames);
        redundant.forEach(
            s -> {
                logger.warn(String.format("Redundant enabled filter (%s) is enabled for this tool by default", s));
            });

        // throw if args were specified for a filter that was also disabled
        disableFilters.forEach(s -> {
            if (requiredPredecessors.contains(s)) {
                throw new CommandLineException(
                        String.format("Values were supplied for (%s) that is also disabled", s));
            }
        });

        // throw if a filter name was specified that has no corresponding instance
        final Map<String, ReadFilter> requestedReadFilters = new HashMap<>();
        userReadFilterNames.forEach(s -> {
            ReadFilter trf = readFilters.get(s);
            if (null == trf) {
                if (!toolDefaultReadFilters.containsKey(s)) {
                    throw new CommandLineException("Unrecognized read filter name: " + s);
                }
            } else {
                requestedReadFilters.put(s, trf);
            }
        });

        // update the readFilters list with the final list of filters specified on the
        // command line; do not include tool defaults as these will be merged in at merge
        // time if they were not disabled
        readFilters = requestedReadFilters;
    }

    /////////////////////////////////////////////////////////
    // ReadFilter plugin-specific helper methods

    /**
     * Determine if a particular ReadFilter was disabled on the command line.
     * @param filterName name of the filter to query
     * @return true if the name appears in the list of disabled filters
     */
    public boolean isDisabledFilter(final String filterName) {
        return disableFilters.contains(filterName);
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
     *                          must return the ALLOW_ALL_READS filter wrapped in the appropriate  type when passed
     *                          a null list.
     * @param <T> extends ReadFilter, type returned by the wrapperFunction
     * @return Single merged read filter.
     */
    public <T extends ReadFilter> T getMergedReadFilter(
            final SAMFileHeader samHeader,
            final BiFunction<List<ReadFilter>, SAMFileHeader, T> aggregateFunction) {

        Utils.nonNull(samHeader);
        Utils.nonNull(aggregateFunction);

        if (disableAllReadFilters) {
            return aggregateFunction.apply(null, samHeader);
        }

        // start with the tool's default filters in the order they were specified, and remove any that were disabled
        // on the command line
        final List<ReadFilter> finalFilters = toolDefaultReadFilterNamesInOrder
                .stream()
                .filter(s -> !isDisabledFilter(s))
                .map(s -> toolDefaultReadFilters.get(s))
                .collect(Collectors.toList());

        // now add in any additional filters enabled on the command line (preserving order)
        final List<ReadFilter> clFilters = getAllInstances();
        if (clFilters != null) {
            clFilters.forEach(f -> finalFilters.add(f));
        }

        return aggregateFunction.apply(finalFilters, samHeader);
    }

}
