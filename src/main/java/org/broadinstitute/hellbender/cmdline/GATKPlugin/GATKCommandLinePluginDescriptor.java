package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

/**
 * A base class for descriptors for plugins that can be dynamically discovered by the
 * command line parser and specified as command line arguments. An instance of each
 * plugin descriptor to be used should be passed to the command line parser, and will
 * be queried to find the class and package names to search for all plugin classes
 * that should be discovered dynamically. The command line parser will find all such
 * classes, and delegate to the descriptor to obtain the corresponding plugin instance;
 * the object returned to the parser is then added to the parser's list of argument sources.
 *
 * Descriptors (sub)classes:
 *
 * - must live in the org.broadinstitute.hellbender.cmdline.GATKPlugin package
 * - should have at least one @Argument used to accumulate the user-specified instances
 *   of the plugin seen on the command line. Allowed values for this argument are the
 *   simple class names of the discovered plugin subclasses.
 *
 * Plugin (sub)classes:
 *
 * - should subclass a common base class (the name of which is returned by the descriptor)
 * - may live in any one of the packages returned by the descriptor {@Link #getPackageNames},
 *   but must have a unique simple name to avoid command line name collisions.
 * - should contain @Arguments for any values they wish to collect. @Arguments may be
 *   optional or required. If required, the arguments are in effect "provisionally
 *   required" in that they are contingent on the specific plugin being specified on
 *   the command line; they will only be marked by the command line parser as missing
 *   if the they have not been specified on the command line, and the plugin class
 *   containing the plugin argument *has* been specified on the command line (as
 *   determined by the command line parser via a call to isDependentArgumentAllowed).
 *
 * NOTE: plugin class @Arguments that are marked "optional=false" should be not have a primitive
 * type, and should not have an initial value, as the command line parser will interpret these as
 * having been set even if they have not been specified on the command line. Conversely, @Arguments
 * that are optional=true should have an initial value, since they parser will not require them
 * to be set in the command line.
 *
 * The methods for each descriptor are called in the following order:
 *
 *  getPluginClass()/getPackageNames() - once when argument parsing begins (if the descriptor
 *  has been passed to the command line parser as a target descriptor)
 *
 *  getClassFilter() - once for each plugin subclass found
 *  getInstance() - once for each plugin subclass that isn't filtered out by getClassFilter
 *  validateDependentArgumentAllowed  - once for each plugin argument value that has been
 *  specified on the command line for a plugin that is controlled by this descriptor
 *
 *  validateArguments() - once when argument parsing is complete
 *  getAllInstances() - whenever the pluggable class consumer wants the resulting plugin instances
 *
 *  getAllowedValuesForDescriptorArgument is only called when the command line parser is constructing
 *  a help/usage message.
 */
public abstract class GATKCommandLinePluginDescriptor<T> {

    /**
     * Return a display name to identify this plugin to the user
     * @return A short user-friendly name for this plugin.
     */
    public String getDisplayName() { return getPluginClass().getSimpleName(); }

    /**
     * Base class for all command line plugin classes managed by this descriptor. Subclasses of
     * this class in any of the packages returned by {@link #getPackageNames} will be command line
     * accessible.
     */
    public abstract Class<?> getPluginClass();

    /**
     * List of package names from which to load command line plugin classes.
     *
     * Note that the simple name of each class must be unique, even across packages.
     * @return List of package names.
     */
    public abstract List<String> getPackageNames();

    /**
     * Give this descriptor a chance to filter out any classes it doesn't want to be
     * dynamically discoverable.
     * @return false if the class shouldn't be used; otherwise true
     */
    public Predicate<Class<?>> getClassFilter() { return c -> true;}

    /**
     * Return an instance of the specified pluggable class. The descriptor should
     * instantiate or otherwise obtain (possibly by having been provided an instance
     * through the descriptor's constructor) an instance of this plugin class.
     * The descriptor should maintain a list of these instances so they can later
     * be retrieved by {@link #getAllInstances}.
     *
     * @param pluggableClass a plugin class discovered by the command line parser that
     *                       was not rejected by {@link #getClassFilter}
     * @return the instantiated object that will be used by the command line parser
     * as an argument source
     * @throws IllegalAccessException
     * @throws InstantiationException
     */
    public abstract Object getInstance(Class<?> pluggableClass)
            throws IllegalAccessException, InstantiationException;

    /**
     * Return the allowable values for the String argument of this plugin descriptor
     * that is specified by longArgName. Called by the command line parser to generate
     * a usage string. If the value is unrecognized, the implementation should throw
     * IllegalArgumentException.
     *
     * @param longArgName
     * @return Set<String> of allowable values, or empty set if any value is allowed
     */
    public abstract Set<String> getAllowedValuesForDescriptorArgument(String longArgName);

    /**
     * Called by the command line parser when an argument value from the class specified
     * by dependentClass has been seen on the command line.
     *
     * Return true if the argument is allowed (i.e., this name of this class was specified
     * as a predecessor on the command line) otherwise false.
     *
     * This method can be used by both the command line parser and the descriptor class for
     * determining when to issue error messages for "dangling" arguments (dependent arguments
     * for which a value has been supplied on the command line, but for which the predecessor
     * argument was not supplied).
     *
     * When this method returns "false", the parser will issue an error message if an argument
     * value in this class has been set on the command line.
     *
     * @param dependentClass
     * @return true if the plugin for this class was specified on the command line, or the
     * values in this class may be set byt he user, otherwise false
     */
    public abstract boolean isDependentArgumentAllowed(Class<?> dependentClass);

    /**
     * This method is called after all command line arguments have been processed to allow
     * the descriptor to validate the plugin arguments that have been specified.
     *
     * It is the descriptor's job to contain an argument list which will be populated
     * by the command line parser with the name of each plugin specified on the command line,
     * and for each such plugin, to maintain a list of the corresponding instance. This
     * method gives the descriptor a chance to reduce that list to include only those
     * instances actually seen on the command line.
     *
     * Implementations of this method should minimally validate that all of values that have
     * been specified on the command line have a corresponding plugin instance (this will
     * detect a user-specified value for which there is no corresponding plugin class).
     *
     * @throws UserException.CommandLineException if a plugin value has been specified that
     * has no corresponding plugin instance (i.e., the plugin class corresponding to the name
     * was not discovered)
     */
    public abstract void validateArguments() throws UserException.CommandLineException;

    /**
     * @return an ordered List of actual plugin instances that have been specified on the command
     * line, in the same order they were obtained/created by {@line #getInstance}).
     */
    public abstract List<T> getAllInstances();

}
