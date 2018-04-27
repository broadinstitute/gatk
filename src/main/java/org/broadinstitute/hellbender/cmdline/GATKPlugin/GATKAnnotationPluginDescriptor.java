package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.PedigreeAnnotation;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.lang.reflect.Modifier;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A plugin descriptor for managing the dynamic discovery of both {@link org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation} and {@link org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation} objects
 * within the packages defined by the method getPackageNames() (default {@link org.broadinstitute.hellbender.tools.walkers.annotator}).
 * Also handles integrating annotation specific arguments from the command line with tool specified defaults.
 *
 * Unlike {@link GATKReadFilterPluginDescriptor} annotation order is not important and thus argument order is not guaranteed to be
 * preserved in all cases, especially when group annotations are involved.
 *
 * An alternative method for discovering annotations is ClassUtils.knownSubInterfaceSimpleNames(Annotation.class), which can
 * be invoked in absence of command line inputs.
 *
 * NOTE: this class enforces that annotations with required arguments must see their arguments, yet this is not currently tested
 *       as no such annotations exist in the GATK.
 */
public class GATKAnnotationPluginDescriptor  extends CommandLinePluginDescriptor<Annotation> {
    //TODO this should be a configurable option or otherwise exposed to the user when configurations are fully supported.
    private static final String pluginPackageName = "org.broadinstitute.hellbender.tools.walkers.annotator";
    private static final Class<?> pluginBaseClass = org.broadinstitute.hellbender.tools.walkers.annotator.Annotation.class;

    protected transient Logger logger = LogManager.getLogger(this.getClass());

    @ArgumentCollection
    private final GATKAnnotationArgumentCollection userArgs;

    // Map of Annotation (simple) class names to the corresponding discovered plugin instance
    private final Map<String, Annotation> allDiscoveredAnnotations = new HashMap<>();

    // Map of Annotation (simple) class names to the corresponding default plugin instance.
    // We keep the actual instances in case they have any additional state provided by the tool
    // when they were created
    private final Map<String, Annotation> toolDefaultAnnotations = new HashMap<>();
    private final Set<String> toolDefaultGroups = new HashSet<>();

    // Set of predecessor annotations for which we've seen arguments exist either as a tool default or be supplied by the user
    // (eg. InbreedingCoeff if we see "--founderID" on the command line)
    private final Set<String> requiredPredecessors = new HashSet<>();

    // Map of annotation group name to list of annotations flagged with that group. The key here the simple name of the interface
    // that describes the annotation group.
    private final Map<String, Map<String, Annotation>> discoveredGroups = new HashMap<>();

    // Resolved instances calculated from the user input and tool defaults together
    private List<Annotation> resolvedInstances;

    // Annotation arguments that are shared and thus must be part of the plugin descriptor
    @Argument(fullName = "founder-id", shortName = "founder-id", doc="Samples representing the population \"founders\"", optional=true)
    private List<String> founderIds;

    @Argument(fullName = StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName = StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for determining the population \"founders\"", optional=true)
    private File pedigreeFile;

    /**
     * @return the class object for the base class of all plugins managed by this descriptor
     */
    @Override
    public Class<?> getPluginBaseClass() {
        return pluginBaseClass;
    }

    /**
     * A list of package names which will be searched for plugins managed by the descriptor.
     *
     * @return Strings of package names over which the plugin descriptor will search for Annotation Classes
     */
    @Override
    public List<String> getPackageNames() {
        return Collections.singletonList(pluginPackageName);
    }

    /**
     * Constructor that allows client tools to specify what annotations (optionally with parameters specified) to use as their defaults
     * before discovery of user specified annotations.
     *
     * @param userArgs           Argument collection to control the exposure of the command line arguments.
     * @param toolDefaultAnnotations Default annotations that may be supplied with arguments
     *                               on the command line. May be null.
     * @param toolDefaultGroups List of tool specified default annotation group names. Annotations specified this way
     *                          will be instantiated with default arguments. may be null.
     */
    public GATKAnnotationPluginDescriptor(final List<Annotation> toolDefaultAnnotations, final List<Class<? extends Annotation>> toolDefaultGroups) {
        this(new DefaultGATKVariantAnnotationArgumentCollection(), toolDefaultAnnotations, toolDefaultGroups);
    }

    public GATKAnnotationPluginDescriptor(final GATKAnnotationArgumentCollection userArgs, final List<Annotation> toolDefaultAnnotations, final List<Class<? extends Annotation>> toolDefaultGroups) {
        this.userArgs = userArgs;
        if (null != toolDefaultAnnotations) {
            toolDefaultAnnotations.forEach(f -> {
                final Class<? extends Annotation> annotClass = f.getClass();
                // anonymous classes have a 0-length simple name, and thus cannot be accessed or
                // controlled by the user via the command line, but they should still be valid
                // as default annotations, so use the full name to ensure that their map entries
                // don't clobber each other
                String className = annotClass.getSimpleName();
                if (className.length() == 0) {
                    className = annotClass.getName();
                }
                populateAnnotationGroups(className, f);
                this.toolDefaultAnnotations.put(className, f);
            });
        }
        if (null != toolDefaultGroups) {
            toolDefaultGroups.forEach(a -> {if (a.isInterface() && a!=Annotation.class) {
                this.toolDefaultGroups.add(a.getSimpleName());
            } else {
                throw new GATKException(String.format("Tool specified annotation group %s is not a valid annotation group, must be an interface extending Annotation", a.getSimpleName()));
            }
            });
        }
    }

    @Override
    public boolean includePluginClass(Class<?> c) {
        return !c.getName().equals(this.getPluginBaseClass().getName()) &&
                    !Modifier.isAbstract(c.getModifiers()) &&
                    !c.getName().contains("UnitTest$");
    }

    /**
     * Return a display name to identify this plugin to the user
     *
     * @return A short user-friendly name for this plugin.
     */
    @Override
    public String getDisplayName() {
        return StandardArgumentDefinitions.ANNOTATION_LONG_NAME;
    }

    /**
     * Return an instance of the specified pluggable class. The descriptor should
     * instantiate or otherwise obtain (possibly by having been provided an instance
     * through the descriptor's constructor) an instance of this plugin class.
     * The descriptor should maintain a list of these instances so they can later
     * be retrieved by {@link #getResolvedInstances()} ()}.
     * <p>
     * In addition, this method should recognize and reject any attempt to instantiate
     * a second instance of any specific Annotation since the user has no way to disambiguate
     * these on the command line).
     *
     * @param pluggableClass a plugin class discovered by the command line parser that
     *                       was not rejected by {@link #includePluginClass(Class)} ()}
     * @return the instantiated object that will be used by the command line parser
     * as an argument source
     * @throws IllegalAccessException if thrown when calling the {@code pluginClass} constructor
     * @throws InstantiationException if thrown when calling the {@code pluginClass} constructor
     */
    @Override
    public Annotation createInstanceForPlugin(Class<?> pluggableClass) throws IllegalAccessException, InstantiationException {
        Annotation annot = null;
        final String simpleName = pluggableClass.getSimpleName();

        if (allDiscoveredAnnotations.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                                    "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            allDiscoveredAnnotations.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultAnnotations.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default annotations;
            // use the default instance as the target for command line argument values
            // rather than creating a new one, in case it has state provided by the tool
            annot = toolDefaultAnnotations.get(simpleName);
        } else {
            annot = (Annotation) pluggableClass.newInstance();
        }

        // Add all annotations to the allDiscoveredAnnotations list, even if the instance came from the
        // tool defaults list (we want the actual instances to be shared to preserve state)
        allDiscoveredAnnotations.put(simpleName, annot);

        populateAnnotationGroups(simpleName, annot);

        return annot;
    }

    // Dynamic discovery of annotation groups
    // We must discover annotation groups and store them for each instance so we can resolve group membership
    // for command line including of groups based on their simple name.
    private void populateAnnotationGroups(final String simpleName, final Annotation annot) {
        Queue<Class<?>> interfaces = new LinkedList<>();
        Collections.addAll(interfaces, annot.getClass().getInterfaces());
        while (!interfaces.isEmpty()) {
            Class<?> inter = interfaces.poll();
            // Following with how groups are currently defined and discovered, namely they are interfaces that
            // extend Annotation, groups are discovered by interrogating annotations for their interfaces and
            // associating the discovered annotations with their defined groups.
            // If a duplicate annotation is added, the group will opt to keep the old instantiation around
            if ((inter != pluginBaseClass) && (pluginBaseClass.isAssignableFrom(inter))) {
                Map<String, Annotation> groupIdentity = (discoveredGroups.containsKey(inter.getSimpleName()) ? discoveredGroups.get(inter.getSimpleName()) : new HashMap<>());
                groupIdentity.putIfAbsent(simpleName, annot);
                discoveredGroups.put(inter.getSimpleName(), groupIdentity);
                Collections.addAll(interfaces, inter.getInterfaces());
            }
        }
    }

    /**
     * Return the allowed values for annotationNames/disableAnnotations/annotationGroups for use by the help system.
     *
     * @param longArgName long name of the argument for which help is requested
     * @return
     */
    @Override
    public Set<String> getAllowedValuesForDescriptorHelp(String longArgName) {
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATION_LONG_NAME)) {
            return allDiscoveredAnnotations.keySet();
        }
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME)) {
            Set<String> annotations = toolDefaultGroups.stream().map(k -> discoveredGroups.get(k).keySet()).flatMap(Collection::stream).collect(Collectors.toSet());
            annotations.addAll(toolDefaultAnnotations.keySet());
            return annotations;
        }
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATION_GROUP_LONG_NAME)) {
            return discoveredGroups.keySet();
        }
        throw new IllegalArgumentException("Allowed values request for unrecognized string argument: " + longArgName);
    }

    @Override
    public boolean isDependentArgumentAllowed(final Class<?> predecessorClass) {
        // Make sure the predecessor for a dependent argument was either specified on the command line or
        // is a tool default, otherwise reject it.
        // NOTE: This method is called by the CLP during parsing at the time the depended argument is seen
        // on the command line. Even if this check passes at the time this method is called, its possible
        // for the user to subsequently disable the required predecessor. That case is caught during final
        // validation done by the validateArguments method.
        String predecessorName = predecessorClass.getSimpleName();
        boolean isAllowed = (userArgs.getUserEnabledAnnotationNames().contains(predecessorName))
                || (toolDefaultAnnotations.get(predecessorName) != null);
        if (!isAllowed) {
            // Check whether any of the annotations have been added via groups (either tool default or user enabled)
            isAllowed = Stream.of(userArgs.getUserEnabledAnnotationGroups(), toolDefaultGroups)
                    .flatMap(Collection::stream)
                    .anyMatch(group ->
                            discoveredGroups.containsKey(group) && discoveredGroups.get(group).keySet().stream().anyMatch(s -> s.equals(predecessorName)));
        }
        if (isAllowed) {
            // Keep track of the ones we allow so we can validate later that they weren't subsequently disabled
            requiredPredecessors.add(predecessorName);
        }
        return isAllowed;
    }

    /**
     * Validate the list of arguments and reduce the list of annotations to those
     * actually seen on the command line. This is called by the command line parser
     * after all arguments have been parsed. Tries to catch most cases where the user
     * provides potentially confusing input.
     */
    @Override
    public void validateAndResolvePlugins() throws CommandLineException {
        // throw if an annotation group is *enabled* more than once by the user
        final Set<String> duplicateUserEnabledAnnotationNames = Utils.getDuplicatedItems(userArgs.getUserEnabledAnnotationNames());
        if (!duplicateUserEnabledAnnotationNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The annotation(s) are enabled more than once: %s",
                            Utils.join(", ", duplicateUserEnabledAnnotationNames)));
        }

        // throw if an annotation is *disabled* more than once by the user
        final Set<String> duplicateDisabledUserAnnotationNames = Utils.getDuplicatedItems(userArgs.getUserDisabledAnnotationNames());
        if (!duplicateDisabledUserAnnotationNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The annotation(s) are disabled more than once: %s",
                            Utils.join(", ", duplicateDisabledUserAnnotationNames)));
        }

        // throw if an annotation is both enabled *and* disabled by the user
        final Set<String> enabledAndDisabled = new HashSet<>(userArgs.getUserEnabledAnnotationNames());
        enabledAndDisabled.retainAll(userArgs.getUserDisabledAnnotationNames());
        if (!enabledAndDisabled.isEmpty()) {
            final String badAnnotationList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The annotation(s): %s are both enabled and disabled", badAnnotationList));
        }

        // throw if a disabled annotation doesn't exist; warn if it wasn't enabled by the tool in the first place
        userArgs.getUserDisabledAnnotationNames().forEach(s -> {
            if (!allDiscoveredAnnotations.containsKey(s)) {
                throw new CommandLineException.BadArgumentValue(String.format("Disabled annotation (%s) does not exist", s));
            } else if (!toolDefaultAnnotations.containsKey(s)) {
                logger.warn(String.format("Disabled annotation (%s) is not enabled by this tool", s));
            }
        });

        // warn if an annotation is both default and enabled by the user
        final Set<String> redundantAnnots = new HashSet<>(toolDefaultAnnotations.keySet());
        redundantAnnots.retainAll(userArgs.getUserEnabledAnnotationNames());
        redundantAnnots.forEach(
                s -> {
                    logger.warn(String.format("Redundant enabled annotation (%s) is enabled for this tool by default", s));
                });

        // warn if an annotation group is both default and enabled by the user
        final Set<String> redundantGroups = new HashSet<>(toolDefaultGroups);
        redundantGroups.retainAll(userArgs.getUserEnabledAnnotationGroups());
        redundantGroups.forEach(
                s -> {
                    logger.warn(String.format("Redundant enabled annotation group (%s) is enabled for this tool by default", s));
                });

        // Throw if args were specified for an annotation that was also disabled, or that was not enabled by the
        // tool by default.
        //
        // Note that this is also checked during command line argument parsing, but needs to be checked again
        // here. Whenever the command line parser sees a dependent argument on the command line, it delegates
        // back to the descriptor's isDependentArgumentAllowed method to allow it to validate that the predecessor
        // for that dependent argument has been supplied, either by a default annotation, or by an explicitly
        // enabled annotation. However, its possible for the user to subsequently try to disable that
        // predecessor, which is what we want to catch here.
        //
        userArgs.getUserDisabledAnnotationNames().forEach(s -> {
            if (requiredPredecessors.contains(s)) {
                String message = String.format("Values were supplied for (%s) that is also disabled", s);
                if (toolDefaultAnnotations.containsKey(s)) {
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

        // throw if an annotation name was specified that has no corresponding instance
        userArgs.getUserEnabledAnnotationNames().forEach(s -> {
            Annotation ta = allDiscoveredAnnotations.get(s);
            if (null == ta) {
                if (!toolDefaultAnnotations.containsKey(s)) {
                    throw new CommandLineException("Unrecognized annotation name: " + s);
                }
            }
        });

        // throw if an annotation group was specified that has no corresponding instance
        userArgs.getUserEnabledAnnotationGroups().forEach(s -> {
            if (!discoveredGroups.containsKey(s)) {
                throw new CommandLineException("Unrecognized annotation group name: " + s);
            }
        });

        // Populating the tool default annotations with the ones requested by groups
        for (String group : toolDefaultGroups ) {
            for (Annotation annot : discoveredGroups.get(group).values()) {
                toolDefaultAnnotations.put(annot.getClass().getSimpleName(), annot);
            }
        }

        // Populating any discovered pedigree annotations with the pedigree arguments from the command line.
        if (((founderIds!=null && !founderIds.isEmpty()) || (pedigreeFile!=null)) && getResolvedInstances().stream()
                    .filter(PedigreeAnnotation.class::isInstance)
                    .map(a -> (PedigreeAnnotation) a)
                    .peek(a -> {
                        if (!founderIds.isEmpty()) a.setFounderIds(founderIds);
                        if (pedigreeFile != null) a.setPedigreeFile(pedigreeFile);
                    })
                    .count() == 0) {
            // Throwing an exception if no pedigree annotations were found
            throw new CommandLineException(
                    String.format(
                            "Pedigree argument \"%s\" or \"%s\" was specified without a pedigree annotation being requested, (eg: %s))",
                            StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME,
                            "founder-id",
                            allDiscoveredAnnotations.entrySet().stream().filter(PedigreeAnnotation.class::isInstance).map(a -> a.getClass().getSimpleName()).collect(Collectors.joining(", "))));
        }
    }

    /**
     * Get the list of default plugins used for this instance of this descriptor. Used for help/doc generation.
     * <p>
     * NOTE: this method does not account for disabled default annotation and just return ALL default instances.
     * The refactored interface in Barclay changes it's contract to allows returning a list with only 'enabled' default
     * instances. We'll change the implementation when we integrate the updated interface.
     *
     * @return A list of Annotation objects that were enabled by the tool by default either by toolDefaultGroups or toolDefaultAnnotations
     */
    @Override
    public List<Annotation> getDefaultInstances() {
        return new ArrayList<>(toolDefaultAnnotations.values());
    }

    /**
     * Merge the default annotations with the users's command line annotation requests, then initialize
     * the resulting annotations. Specifically, unless the user disables all tool default annotations it will
     * first add all the tool enabled annotations which were not individually blocked by the user and then
     * adds in annotations defined by the users specified groups, then individual annotations.
     *
     * NOTE: calling this method before argument parsing (and thus before {@link #validateAndResolvePlugins}
     * has been called) may return a different list than calling it after parsing, because annotations associated
     * with pedigree files will not have had their arguments input.
     *
     * @return An unordered Collection of annotations.
     */
    @Override
    public List<Annotation> getResolvedInstances() {
        if (resolvedInstances == null) {
            final SortedSet<Annotation> annotations = new TreeSet<>(Comparator.comparing(t -> t.getClass().getSimpleName()));

            if (!userArgs.getDisableToolDefaultAnnotations()) {
                annotations.addAll(toolDefaultAnnotations.values());
            }
            for (String group : userArgs.getUserEnabledAnnotationGroups()) {
                annotations.addAll(discoveredGroups.get(group).values());
            }
            if (userArgs.getEnableAllAnnotations()) {
                annotations.addAll(allDiscoveredAnnotations.values());
            } else {
                for (String annotation : userArgs.getUserEnabledAnnotationNames()) {
                    annotations.add(allDiscoveredAnnotations.get(annotation));
                }
            }
            resolvedInstances = annotations.stream().filter(t -> !userArgs.getUserDisabledAnnotationNames().contains(t.getClass().getSimpleName())).collect(Collectors.toList());
        }
        return resolvedInstances;
    }


    /**
     * Return the class representing the instance of the plugin specified by {@code pluginName}
     *
     * @param pluginName Name of the plugin requested
     * @return Class object for the plugin instance requested
     */
    @Override
    public Class<?> getClassForPluginHelp(final String pluginName) {
        return allDiscoveredAnnotations.containsKey(pluginName) ? allDiscoveredAnnotations.get(pluginName).getClass() : null;
    }

}
