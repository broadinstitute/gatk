package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.ReadTransformerArgumentDefinitions;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.config.GATKConfig;

import java.io.IOException;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 * A CommandLinePluginDescriptor for ReadTransformer plugins
 */
public class GATKReadTransformerPluginDescriptor extends CommandLinePluginDescriptor<ReadTransformer> {

    protected transient Logger logger = LogManager.getLogger(this.getClass());

    private void readObject(java.io.ObjectInputStream in)
            throws IOException, ClassNotFoundException {
        in.defaultReadObject();
        logger = LogManager.getLogger(this.getClass()); // Logger is not serializable (even by Kryo)
    }

    /**
     * At startup, set the plugin package name to the one(s) in the configuration file.
     */
    private static final List<String> PLUGIN_PACKAGE_NAMES;
    static {
        // Get our configuration:
        final GATKConfig config = ConfigFactory.getInstance().getGATKConfig();
        // Exclude abstract classes and interfaces from the list of discovered codec classes
        PLUGIN_PACKAGE_NAMES = Collections.unmodifiableList(config.read_transformer_packages());
    }

    private static final Class<ReadTransformer> pluginBaseClass = ReadTransformer.class;
    private static final String READ_PLUGIN_DISPLAY_NAME = "readTransformer";

    // the purpose of this argument collection is to allow the caller to control the exposure of the command line arguments
    @VisibleForTesting
    @ArgumentCollection
    final GATKReadTransformerArgumentCollection userArgs;

    // Map of read transformer (simple) class names to the corresponding discovered plugin instance
    private final Map<String, ReadTransformer> allDiscoveredReadTransformers = new HashMap<>();

    // Map of read transformer (simple) class names to the corresponding default plugin instance
    // it is a LinkedHashMap because we want to remember the order in which these were provided, and also keep the
    // actual instances in case they have any additional state provided by the tool
    // when they were created
    private final Map<String, ReadTransformer> toolDefaultReadTransformers = new LinkedHashMap<>();

    // Set of predecessor readTransformers for which we've seen arguments that must exist either as a tool default or be supplied by the user
    // (eg. ReadLengthReadTransformer if we see "--maxReadLength" on the command line)
    private final Set<String> requiredPredecessors = new HashSet<>();

    /**
     * @param userArgs           Argument collection to control the exposure of the command line arguments.
     * @param toolDefaultTransformers Default transformers that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKReadTransformerPluginDescriptor(final GATKReadTransformerArgumentCollection userArgs, final List<ReadTransformer> toolDefaultTransformers) {
        this.userArgs = userArgs;
        if (null != toolDefaultTransformers) {
            toolDefaultTransformers.forEach(f -> {
                final Class<? extends ReadTransformer> rfClass = f.getClass();
                // anonymous classes have a 0-length simple name, and thus cannot be accessed or
                // controlled by the user via the command line, but they should still be valid
                // as default transformers, so use the full name to ensure that their map entries
                // don't clobber each other
                String className = rfClass.getSimpleName();
                if (className.length() == 0) {
                    className = rfClass.getName();
                }
                toolDefaultReadTransformers.put(className, f);
            });
        }
    }

    /**
     * @param toolDefaultTransformers Default transformers that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKReadTransformerPluginDescriptor(final List<ReadTransformer> toolDefaultTransformers) {
        this(new DefaultGATKReadTransformerArgumentCollection(), toolDefaultTransformers);
    }

    /////////////////////////////////////////////////////////
    // GATKCommandLinePluginDescriptor implementation methods

    /**
     * Return a display name to identify this plugin to the user
     * @return A short user-friendly name for this plugin.
     */
    @Override
    public String getDisplayName() {
        // The value returned by this method is placed into the freemarker property map by the docgen system,
        // and becomes a variable in the freemarker template language. We can't use the actual argument name
        // here because the kebabified version ("read-transformer") isn't a valid freemarker variable (and is
        // parsed as an expression with an embedded '-'). So use the pre-kebabified version, which matches
        // what is used in the template).
        return READ_PLUGIN_DISPLAY_NAME;
    }

    /**
     * @return the class object for the base class of all plugins managed by this descriptor
     */
    @Override
    public Class<ReadTransformer> getPluginBaseClass() {return pluginBaseClass;}

    /**
     * A list of package names which will be searched for plugins managed by the descriptor.
     * @return
     */
    @Override
    public List<String> getPackageNames() {return PLUGIN_PACKAGE_NAMES;}

    @Override
    public boolean includePluginClass(Class<?> c) {
        // don't use the ReadTransformer base class, it's inner classes, the CountingReadTransformer,
        // or the unit tests
        return !c.getName().equals(this.getPluginBaseClass().getName()) &&
                !c.getName().startsWith(this.getPluginBaseClass().getName() + "$") &&
                !c.getName().contains("UnitTest$");
    }

    // Instantiate a new ReadTransformer derived object and save it in the list
    @Override
    @SuppressWarnings("deprecation")
    public ReadTransformer createInstanceForPlugin(final Class<?> pluggableClass) throws IllegalAccessException, InstantiationException {
        ReadTransformer readTransformer = null;
        final String simpleName = pluggableClass.getSimpleName();

        if (allDiscoveredReadTransformers.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                            "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            allDiscoveredReadTransformers.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultReadTransformers.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default transformers;
            // use the default instance as the target for command line argument values
            // rather than creating a new one, in case it has state provided by the tool
            readTransformer = toolDefaultReadTransformers.get(simpleName);
        } else {
            readTransformer = (ReadTransformer) pluggableClass.newInstance();
        }

        // Add all transformers to the allDiscoveredReadTransformers list, even if the instance came from the
        // tool defaults list (we want the actual instances to be shared to preserve state)
        allDiscoveredReadTransformers.put(simpleName, readTransformer);
        return readTransformer;
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
        boolean isAllowed = userArgs.getUserEnabledReadTransformerNames().contains(predecessorName)
                || (toolDefaultReadTransformers.get(predecessorName) != null);
        if (isAllowed) {
            // Keep track of the ones we allow so we can validate later that they weren't subsequently disabled
            requiredPredecessors.add(predecessorName);
        }
        return isAllowed;
    }

    /**
     * @return an ordered list of ReadTransformer instances that result from resolving all command line
     * arguments with any default plugins that have been provided to this descriptor. This list
     * represents the ReadTransformers that will actually be used by the consumer.
     *
     * NOTE: calling this method before argument parsing (and thus before {@link #validateAndResolvePlugins}
     * has been called) may return a different list than calling it after parsing, since the resolution
     * policy is implementation-dependent and may depend on the actual arguments specified by the user
     * on the command line.
     */
    @Override
    public List<ReadTransformer> getResolvedInstances() {
        // start with the tool's default transformers in the order they were specified, and remove any that
        // were disabled on the command line (if --disableToolDefaultReadTransformers is specified, just initialize
        // an empty list with initial capacity of user transformers)
        final List<ReadTransformer> finalTransformers =
                userArgs.getDisableToolDefaultReadTransformers() ?
                        new ArrayList<>(userArgs.getUserEnabledReadTransformerNames().size()) :
                        toolDefaultReadTransformers.entrySet()
                                .stream()
                                .filter(e -> !isDisabledTransformer(e.getKey()))
                                .map(e -> e.getValue())
                                .collect(Collectors.toList());

        // now add in any additional transformers enabled on the command line (preserving order)
        final List<ReadTransformer> userEnabledTransformers = getUserEnabledInstances();
        if (userEnabledTransformers != null) {
            userEnabledTransformers.stream()
                    .filter(f -> !finalTransformers.contains(f)) // remove redundant transformers
                    .forEach(f -> finalTransformers.add(f));
        }
        return finalTransformers;
    }

    // Return a list with a read transformer instance for each read transformer enabled by the user.
    private List<ReadTransformer> getUserEnabledInstances() {
        final ArrayList<ReadTransformer> transformers = new ArrayList<>(userArgs.getUserEnabledReadTransformerNames().size());
        userArgs.getUserEnabledReadTransformerNames().forEach(s -> {
            ReadTransformer rf = allDiscoveredReadTransformers.get(s);
            transformers.add(rf);
        });
        return transformers;
    }

    /**
     * Get the list of default plugins used for this instance of this descriptor. Used for help/doc generation.
     *
     * NOTE: this method does not account for disabled default transformers and just return ALL default instances.
     * The refactored interface in Barclay changes it's contract to allows returning a list with only 'enabled' default
     * instances. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<ReadTransformer> getDefaultInstances() { return new ArrayList<>(toolDefaultReadTransformers.values()); }

    /**
     * Return the allowed values for read-transformer/disable-read-transformer names for use by the help system.
     * @param longArgName long name of the argument for which help is requested
     * @return
     */
    @Override
    public Set<String> getAllowedValuesForDescriptorHelp(final String longArgName) {
        if (longArgName.equals(ReadTransformerArgumentDefinitions.READ_TRANSFORMER_LONG_NAME)) {
            return allDiscoveredReadTransformers.keySet();
        }
        if (longArgName.equals(ReadTransformerArgumentDefinitions.DISABLE_READ_TRANSFORMER_LONG_NAME)) {
            return toolDefaultReadTransformers.keySet();
        }
        return null;
    }

    /**
     * Return the class object for the plugin with simple class name {@code pluginName}
     * Used for help/usage and documentation generation.
     *
     * @param pluginName Name of the plugin requested
     * @return Class object for the plugin instance requested, upper bounded by type {@code T}
     */
    @Override
    public Class<?> getClassForPluginHelp(final String pluginName) {
        if (allDiscoveredReadTransformers.get(pluginName) != null) {
            return allDiscoveredReadTransformers.get(pluginName).getClass();
        } else {
            throw new IllegalArgumentException(String.format("Can't resolve ReadTransformer plugin for name %s", pluginName));
        }
    }

    /**
     * Validate the list of arguments and reduce the list of read transformers to those
     * actually seen on the command line. This is called by the command line parser
     * after all arguments have been parsed.
     */
    @Override
    public void validateAndResolvePlugins() {
        // throw if a transformer is *enabled* more than once by the user
        final Set<String> duplicateUserEnabledTransformerNames = Utils.getDuplicatedItems(userArgs.getUserEnabledReadTransformerNames());
        if (!duplicateUserEnabledTransformerNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The read transformer(s) are enabled more than once: %s",
                            Utils.join(", ", duplicateUserEnabledTransformerNames)));
        }

        // throw if a transformer is *disabled* more than once by the user
        final Set<String> duplicateDisabledUserTransformerNames = Utils.getDuplicatedItems(userArgs.getUserDisabledReadTransformerNames());
        if (!duplicateDisabledUserTransformerNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The read transformer(s) are disabled more than once: %s",
                            Utils.join(", ", duplicateDisabledUserTransformerNames)));
        }

        // throw if a transformer is both enabled *and* disabled by the user
        final Set<String> enabledAndDisabled = new HashSet<>(userArgs.getUserEnabledReadTransformerNames());
        enabledAndDisabled.retainAll(userArgs.getUserDisabledReadTransformerNames());
        if (!enabledAndDisabled.isEmpty()) {
            final String badTransformersList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The read transformer(s): %s are both enabled and disabled", badTransformersList));
        }

        // throw if a disabled transformer doesn't exist; warn if it wasn't enabled by the tool in the first place
        userArgs.getUserDisabledReadTransformerNames().forEach(s -> {
            if (!allDiscoveredReadTransformers.containsKey(s)) {
                throw new CommandLineException.BadArgumentValue(String.format("Disabled transformer (%s) does not exist", s));
            } else if (!toolDefaultReadTransformers.containsKey(s)) {
                logger.warn(String.format("Disabled transformer (%s) is not enabled by this tool", s));
            }
        });

        // warn if a transformer is both default and enabled by the user
        final Set<String> redundant = new HashSet<>(toolDefaultReadTransformers.keySet());
        redundant.retainAll(userArgs.getUserEnabledReadTransformerNames());
        redundant.forEach(
            s -> {
                logger.warn(String.format("Redundant enabled transformer (%s) is enabled for this tool by default", s));
            });

        // Throw if args were specified for a transformer that was also disabled, or that was not enabled by the
        // tool by default.
        //
        // Note that this is also checked during command line argument parsing, but needs to be checked again
        // here. Whenever the command line parser sees a dependent argument on the command line, it delegates
        // back to the descriptor's isDependentArgumentAllowed method to allow it to validate that the predecessor
        // for that dependent argument has been supplied, either by a default read transformer, or by an explicitly
        // enabled read transformer. However, its possible for the user to subsequently try to disable that
        // predecessor, which is what we want to catch here.
        //
        userArgs.getUserDisabledReadTransformerNames().forEach(s -> {
            if (requiredPredecessors.contains(s)) {
                String message = String.format("Values were supplied for (%s) that is also disabled", s);
                if (toolDefaultReadTransformers.containsKey(s)) {
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

        // throw if a transformer name was specified that has no corresponding instance
        final Map<String, ReadTransformer> requestedReadTransformers = new HashMap<>();
        userArgs.getUserEnabledReadTransformerNames().forEach(s -> {
            ReadTransformer trf = allDiscoveredReadTransformers.get(s);
            if (null == trf) {
                if (!toolDefaultReadTransformers.containsKey(s)) {
                    throw new CommandLineException("Unrecognized read transformer name: " + s);
                }
            } else {
                requestedReadTransformers.put(s, trf);
            }
        });
    }

    /////////////////////////////////////////////////////////
    // ReadTransformer plugin-specific helper methods

    /**
     * Determine if a particular ReadTransformer was disabled on the command line, either directly of by disabling all
     * tool defaults.
     * @param transformerName name of the transformer to query.
     * @return {@code true} if the name appears in the list of disabled transformers, or is a tool default not provided by
     * the user and all tool defaults are disabled; {@code false} otherwise.
     */
    public boolean isDisabledTransformer(final String transformerName) {
        return userArgs.getUserDisabledReadTransformerNames().contains(transformerName)
                || (userArgs.getDisableToolDefaultReadTransformers() && !userArgs.getUserEnabledReadTransformerNames().contains(transformerName));
    }

    /**
     * Merge the default transformers with the users's command line read transformer requests, then initialize
     * the resulting transformers.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read transformer instances
     * @return Single merged read transformer.
     */
    public final ReadTransformer getMergedReadTransformer(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadTransformer(
                samHeader,
                ReadTransformer::fromList
        );
    }

    /**
     * Merge the default transformers with the users's command line read transformer requests, then initialize
     * the resulting transformers.
     *
     * @param samHeader a SAMFileHeader to initialize read transformer instances. May not be null.
     * @param aggregateFunction function to use to merge ReadTransformers, usually ReadTransformer::fromList. The function
     *                          must return the ALLOW_ALL_READS transformer wrapped in the appropriate type when passed
     *                          a null or empty list.
     * @param <T> extends ReadTransformer, type returned by the wrapperFunction
     * @return Single merged read transformer.
     */
    public <T extends ReadTransformer> T getMergedReadTransformer(
            final SAMFileHeader samHeader,
            final BiFunction<List<ReadTransformer>, SAMFileHeader, T> aggregateFunction) {

        Utils.nonNull(samHeader);
        Utils.nonNull(aggregateFunction);

        final List<ReadTransformer> finalTransformers = getResolvedInstances();
        return aggregateFunction.apply(finalTransformers, samHeader);
    }

}
