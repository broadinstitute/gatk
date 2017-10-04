package org.broadinstitute.hellbender.utils.config;

import com.google.common.annotations.VisibleForTesting;
import org.aeonbits.owner.*;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A static class to contain configuration utility methods.
 * Created by jonn on 7/19/17.
 */
public final class ConfigUtils {

    static private final Logger logger = LogManager.getLogger(ConfigUtils.class);

    // This class just has static methods to help with configuration, so no need for this constructor.
    private ConfigUtils() {}

    /**
     * A regex to use to look for variables in the Sources annotation
     */
    private static final Pattern sourcesAnnotationPathVariablePattern = Pattern.compile("\\$\\{(.*)}");

    /**
     * A set to keep track of the classes we've already resolved for configuration path purposes:
     */
    private static final Set<Class<? extends Config>> alreadyResolvedPathVariables = new HashSet<>();

    /**
     * Value to set each variable for configuration file paths when the variable
     * has not been set in either Java System properties or environment properties.
     */
    @VisibleForTesting
    static final String NO_PATH_VARIABLE_VALUE = "/dev/null";

    // =================================================================================================================

    /**
     * Checks each of the given {@code filenameProperties} for if they are defined in system {@link System#getProperties()}
     * or environment {@link System#getenv()} properties.  If they are not, this method will set them in the
     * {@link ConfigFactory} to a non-existent path so the {@link ConfigFactory} will know to try to resolve them as
     * variables at load-time (and not as raw paths).
     * @param filenameProperties A {@link List} of filename properties as specified in {@link Config} {@link org.aeonbits.owner.Config.Sources} annotations to check for existence in system and environment properties.
     */
    @VisibleForTesting
    static void checkFileNamePropertyExistenceAndSetConfigFactoryProperties(final List<String> filenameProperties) {
        // Grab the system properties:
        final Properties systemProperties = System.getProperties();

        // Grab the environment properties:
        final Map<String, String> environmentProperties = System.getenv();

        // Make sure that if our property isn't in the system, environment, and ConfigFactory
        // properties, that we set it to a neutral value that will not contain
        // anything (so that the property will fall back into the next value).
        for (final String property : filenameProperties) {

            if ( environmentProperties.keySet().contains(property) ) {
                logger.info("Config path variable found in Environment Properties: " + property + "=" + environmentProperties.get(property) + " - will search for config here.");
            }
            else if ( systemProperties.keySet().contains(property) ) {
                logger.info("Config path variable found in System Properties: " + property + "=" + systemProperties.get(property) + " - will search for config here.");
            }
            else if ( ConfigFactory.getProperties().keySet().contains(property) ) {
                logger.info("Config path variable found in Config Factory Properties(probably from the command-line): " + property + "=" + ConfigFactory.getProperty(property) + " - will search for config here.");
            }
            else {
                logger.warn("Config path variable not found: " + property +
                        " - setting value to default empty variable: " + NO_PATH_VARIABLE_VALUE);
                ConfigFactory.setProperty(property, NO_PATH_VARIABLE_VALUE);
            }
        }
    }

    /**
     * Get a list of the config file variables from the given {@link Config} classes.
     * @param configurationClasses A list of configuration classes from which to extract variable names in their {@link org.aeonbits.owner.Config.Sources}.
     * @return A list of variables in the {@link org.aeonbits.owner.Config.Sources} of the given {@code configurationClasses}
     */
    @VisibleForTesting
    static List<String> getConfigPathVariableNamesFromConfigClasses(final List<Class<?>> configurationClasses) {

        final List<String> configPathVariableNames = new ArrayList<>();

        // Loop through our classes and grab any sources with variables in there:
        for ( final Class<?> clazz : configurationClasses ) {

            // Make sure that we get config classes here:
            if ( Config.class.isAssignableFrom(clazz) ) {
                final Config.Sources annotation = clazz.getAnnotation(Config.Sources.class);

                if ( annotation != null ) {
                    final String[] annotationValues = annotation.value();

                    for (final String val : annotationValues) {

                        final Matcher m = sourcesAnnotationPathVariablePattern.matcher(val);
                        if (m.find()) {
                            configPathVariableNames.add(m.group(1));
                        }
                    }
                }
            }
        }

        return configPathVariableNames;
    }

    /**
     * Get a list of the config file variables from the given {@link Config} class.
     * @param configClass A lconfiguration class from which to extract variable names in its {@link org.aeonbits.owner.Config.Sources}.
     * @return A list of variables in the {@link org.aeonbits.owner.Config.Sources} of the given {@code configClass}
     */
    @VisibleForTesting
    static <T extends Config> List<String> getSourcesAnnotationPathVariables(final Class<? extends T> configClass) {

        final List<String> configPathVariableNames = new ArrayList<>();

        final Config.Sources annotation = configClass.getAnnotation(Config.Sources.class);

        if ( annotation != null ) {
            for (final String val : annotation.value()) {

                final Matcher m = sourcesAnnotationPathVariablePattern.matcher(val);
                if (m.find()) {
                    configPathVariableNames.add(m.group(1));
                }
            }
        }

        return configPathVariableNames;
    }

    /**
     * Gets all system properties from the given {@link Config}-derived object.
     * System properties are denoted via the presence of the {@link SystemProperty} annotation.
     * @param config A {@link Config}-derived object from which to read system properties.
     * @param <T> A {@link Config}-derived type from which to read System Properties.
     * @return A properties {@link Map} representing all System Properties in the given {@code config}.
     */
    @VisibleForTesting
    static <T extends Config> Map<String, String> getSystemPropertiesFromConfig(final T config) {

        final Map<String, String> properties = new HashMap<>();

        // This is gross and uses reflection to get all methods in the given config class
        // and then interrogates those methods for internal data on the config parameters.

        // We have to match our interfaces to the config interface that we're actually using.
        // It's not as simple as using getDeclaredMethods on the Class object because we'll get
        // a LOT of extraneous stuff that we don't care about.
        for ( final Class<?> classInterface : config.getClass().getInterfaces() ) {

            // If we have an interface that is a child of the OWNER Config interface, then
            // we must have the right object.
            if (Config.class.isAssignableFrom(classInterface)) {

                // Now we cycle through our interface methods, resolve parameter names,
                // and set the values in the system.
                for (final Method propertyMethod : classInterface.getDeclaredMethods()) {

                    // Get the property name:
                    String propertyName = propertyMethod.getName();

                    final Config.Key key = propertyMethod.getAnnotation(Config.Key.class);
                    if (key != null) {
                        propertyName = key.value();
                    }

                    // Only set properties that have the SystemProperty annotation:
                    if (propertyMethod.getAnnotation(SystemProperty.class) == null)  {
                        continue;
                    }

                    // Get the value of the property into a string:
                    final String propertyValue;
                    try {
                        propertyValue = String.valueOf( propertyMethod.invoke(config, new Object[]{}) );
                    } catch (final IllegalAccessException ex) {
                        throw new GATKException("Could not access the config getter: " +
                                config.getClass().getSimpleName() + "." +
                                propertyMethod.getName(), ex);

                    } catch (final InvocationTargetException ex) {
                        throw new GATKException("Could not invoke the config getter: " +
                                config.getClass().getSimpleName() + "." +
                                propertyMethod.getName(), ex);
                    }

                    // Add our property:
                    properties.put(propertyName, propertyValue);
                }
            }
        }
        return properties;
    }

    /**
     * Injects the given property to the System Properties.
     * Validates that this property was set after setting it.
     * @param properties A {@link Map} of key, value pairs of properties to add to the System Properties.
     */
    @VisibleForTesting
    static void injectToSystemProperties(final Map<String, String> properties) {

        for ( final String propertyName : properties.keySet() ) {

            final String propertyValue = properties.get(propertyName);

            System.setProperty(propertyName, propertyValue);

            // Test for validation that it worked:
            final String propertyValueThatWasSet = System.getProperty(propertyName);
            if (propertyValueThatWasSet == null) {
                throw new GATKException("Unable to set System Property (" + propertyName + "=" + propertyValue + ")!");
            }

            if (!propertyValueThatWasSet.equals(propertyValue)) {
                throw new GATKException("System Property corrupted (" + propertyName + "!=" + propertyValue + " -> " + propertyValueThatWasSet + ")!");
            }
        }
    }

    // =================================================================================================================

    /**
     * Override for {@link ConfigFactory#create(Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * Creates a {@link Config} instance from the specified interface
     *
     * @param clazz   the interface extending from {@link Config} that you want to instantiate.
     * @param imports additional variables to be used to resolve the properties.
     * @param <T>     type of the interface.
     * @return an object implementing the given interface, which maps methods to property values.
     */
    public static <T extends Config> T create(final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(clazz);

        if ( !alreadyResolvedPathVariables.contains(clazz) ) {
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(
                    getSourcesAnnotationPathVariables(clazz)
            );

            alreadyResolvedPathVariables.add(clazz);
        }

        return ConfigFactory.create(clazz, imports);
    }

    /**
     * Override for {@link ConfigCache#getOrCreate(Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * Gets from the cache or create, an instance of the given class using the given imports.
     * The factory used to create new instances is the static {@link ConfigFactory#INSTANCE}.
     *
     * @param clazz     the interface extending from {@link Config} that you want to instantiate.
     * @param imports   additional variables to be used to resolve the properties.
     * @param <T>       type of the interface.
     * @return          an object implementing the given interface, that can be taken from the cache,
     *                  which maps methods to property values.
     */
    public static <T extends Config> T getOrCreate(final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(clazz);

        if ( !alreadyResolvedPathVariables.contains(clazz) ) {
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(
                    getSourcesAnnotationPathVariables(clazz)
            );

            alreadyResolvedPathVariables.add(clazz);
        }

        return ConfigCache.getOrCreate(clazz, imports);
    }

    /**
     * Override for {@link ConfigCache#getOrCreate(Factory, Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * Gets from the cache or create, an instance of the given class using the given imports.
     *
     * @param factory   the factory to use to eventually create the instance.
     * @param clazz     the interface extending from {@link Config} that you want to instantiate.
     * @param imports   additional variables to be used to resolve the properties.
     * @param <T>       type of the interface.
     * @return          an object implementing the given interface, that can be taken from the cache,
     *                  which maps methods to property values.
     */
    public static <T extends Config> T getOrCreate(final Factory factory, final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(factory);
        Utils.nonNull(clazz);

        if ( !alreadyResolvedPathVariables.contains(clazz) ) {
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(
                    getSourcesAnnotationPathVariables(clazz)
            );

            alreadyResolvedPathVariables.add(clazz);
        }

        return ConfigCache.getOrCreate(factory, clazz, imports);
    }

    /**
     * Override for {@link ConfigCache#getOrCreate(Object, Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * Gets from the cache or create, an instance of the given class using the given imports.
     * The factory used to create new instances is the static {@link ConfigFactory#INSTANCE}.
     *
     * @param key       the key object to be used to identify the instance in the cache.
     * @param clazz     the interface extending from {@link Config} that you want to instantiate.
     * @param imports   additional variables to be used to resolve the properties.
     * @param <T>       type of the interface.
     * @return          an object implementing the given interface, that can be taken from the cache,
     *                  which maps methods to property values.
     */
    public static <T extends Config> T getOrCreate(final Object key, final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(key);
        Utils.nonNull(clazz);

        if ( !alreadyResolvedPathVariables.contains(clazz) ) {
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(
                    getSourcesAnnotationPathVariables(clazz)
            );

            alreadyResolvedPathVariables.add(clazz);
        }

        return ConfigCache.getOrCreate(key, clazz, imports);
    }

    /**
     * Override for {@link ConfigCache#getOrCreate(Factory, Object, Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * @param factory   the factory to use to eventually create the instance.
     * @param key       the key object to be used to identify the instance in the cache.
     * @param clazz     the interface extending from {@link Config} that you want to instantiate.
     * @param imports   additional variables to be used to resolve the properties.
     * @param <T>       type of the interface.
     * @return          an object implementing the given interface, that can be taken from the cache,
     *                  which maps methods to property values.
     */
    public static <T extends Config> T getOrCreate(final Factory factory, final Object key,
                                                   final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(factory);
        Utils.nonNull(key);
        Utils.nonNull(clazz);

        if ( !alreadyResolvedPathVariables.contains(clazz) ) {
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(
                    getSourcesAnnotationPathVariables(clazz)
            );

            alreadyResolvedPathVariables.add(clazz);
        }

        return ConfigCache.getOrCreate(key, clazz, imports);
    }

    /**
     * Override for {@link ConfigCache#get(Object)}.
     * This method is here to complete the interface for getting {@link Config} objects.
     * It is basically a wrapper.
     *
     * Gets from the cache the {@link Config} instance identified by the given key.
     *
     * @param key       the key object to be used to identify the instance in the cache.
     * @param <T>       type of the interface.
     * @return          the {@link Config} object from the cache if exists, or <tt>null</tt> if it doesn't.
     */
    public static <T extends Config> T get(final Object key) {
        Utils.nonNull(key);
        return ConfigCache.get(key);
    }

    /**
     * Get the configuration file name from the given arguments.
     * Modifies the given arguments to remove both the configuration file specification string
     * and the configuration file name from the args.
     *
     * NOTE: Does NOT validate that the resulting string is a valid configuration file.
     *
     * @param args Command-line arguments passed to this program.
     * @param configFileOption The command-line option indicating that the config file is next
     * @return The name of the configuration file for this program or {@code null}.
     */
    public static String getConfigFilenameFromArgs( final String[] args, final String configFileOption ) {

        Utils.nonNull(args);
        Utils.nonNull(configFileOption);

        String configFileName = null;

        for ( int i = 0 ; i < args.length ; ++i ) {
            if (args[i].equals(configFileOption)) {

                // Get the config file name:
                if ( (i+1) < args.length ) {
                    configFileName = args[i+1];
                    break;
                }
                else {
                    // Option was provided, but no file was specified.
                    // We cannot work under these conditions:
                    throw new UserException.BadInput("ERROR: Configuration file not given after config file option specified: " + configFileOption);
                }
            }
        }

        return configFileName;
    }

    /**
     * Removes the configuration file option, {@link StandardArgumentDefinitions#GATK_CONFIG_FILE_OPTION}, and the
     * corresponding config file path from the given {@code args}.
     * If the config file option is given and the file is not specified, throws a {@link UserException.BadInput}.
     * Assumes that if {@code configFileOption} is given in {@code args}, the following item is the config file path.
     * @param args The argument list from which to remove the config file option and config file path.
     * @return An {@code Array} of {@link String} containing all arguments without the config file option and path.
     */
    public static String[] removeConfigOptionAndFileFromArgs( final String[] args ) {
        return removeConfigOptionAndFileFromArgs(args, StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION);
    }

    /**
     * Removes the given {@code configFileOption}, and the corresponding config file path from the given {@code args}.
     * If the config file option is given and the file is not specified, throws a {@link UserException.BadInput}.
     * Assumes that if {@code configFileOption} is given in {@code args}, the following item is the config file path.
     * @param args The argument list from which to remove the config file option and config file path.
     * @param configFileOption The name of the option specifying a config file.
     * @return An {@code Array} of {@link String} containing all arguments without the config file option and path.
     */
    public static String[] removeConfigOptionAndFileFromArgs( final String[] args, final String configFileOption ) {
        final List<String> cleanArgs = new ArrayList<>(((args.length - 2) >= 0) ? args.length-2 : args.length);

        for ( int i = 0; i < args.length; ++i) {
            if ( args[i].equals(configFileOption) ) {
                if ( (i + 1) < args.length ) {
                    // Skip the config file:
                    ++i;
                }
                else {
                    // Option was provided, but no file was specified.
                    // We cannot work under these conditions:
                    throw new UserException.BadInput("ERROR: Configuration file not given after config file option specified: " + configFileOption);
                }
            }
            else {
                // Add our argument to the list:
                cleanArgs.add( args[i] );
            }
        }

        return cleanArgs.toArray(new String[cleanArgs.size()]);
    }

    /**
     * Get the configuration filename from the command-line (if it exists) and create a configuration for it.
     * Configuration type defaults to {@link GATKConfig}
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     */
    public static void initializeConfigurationsFromCommandLineArgs(final String[] argList,
                                                                   final String configFileOption) {
        initializeConfigurationsFromCommandLineArgs(
                argList,
                configFileOption,
                GATKConfig.class
        );
    }

    /**
     * Get the configuration from filename the command-line (if it exists) and create a configuration for it of the given type.
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     * @param configClass The class of the configuration file to instantiate.
     */
    public static <T extends Config> void initializeConfigurationsFromCommandLineArgs(final String[] argList,
                                                                                      final String configFileOption,
                                                                                      final Class<? extends T> configClass) {
        Utils.nonNull(argList);
        Utils.nonNull(configFileOption);
        Utils.nonNull(configClass);

        // Get main config from args:
        final String configFileName = getConfigFilenameFromArgs( argList, configFileOption );

        // Set the config path if we've specified it:
        if ( configFileName != null ){
            ConfigFactory.setProperty( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, configFileName );
        }

        // Set the config file to be the one we want to use from the command-line:
        final T gatkConfig = ConfigUtils.getOrCreate(configClass);

        // To start with we inject our system properties to ensure they are defined for downstream components:
        injectSystemPropertiesFromConfig( gatkConfig );
    }

    /**
     * Initializes and returns the configuration as specified by {@code configFileName}
     * Also caches this configuration in the {@link ConfigCache} for use elsewhere.
     * @param configFileName The name of the file from which to initialize the configuration
     * @param configClass The type of configuration in which to interpret the given {@code configFileName}
     * @return The configuration instance implementing {@link GATKConfig} containing any overrides in the given file.
     */
    public static <T extends Config> T initializeConfigurationAsProperties(final String configFileName,
                                                                           final Class<? extends T> configClass) {

        Utils.nonNull(configClass);

        // Get a place to store our properties:
        final Properties userConfigFileProperties = new Properties();

        // Try to get the config from the specified file:
        if ( configFileName != null ) {

            try (final FileInputStream userConfigFileInputStream = new FileInputStream(configFileName)) {
                userConfigFileProperties.load(userConfigFileInputStream);

                logger.info("Found " + configClass.getSimpleName() + " Configuration File: " + configFileName);

            } catch (final FileNotFoundException e) {
                throw new GATKException("Unable to find specified " + configClass.getSimpleName() + " configuration file: "
                        + configFileName + " - defaulting to built-in config settings.", e);
            }
            catch (final IOException e) {
                throw new GATKException("Unable to load specified " + configClass.getSimpleName() + " configuration file: "
                        + configFileName + " - defaulting to built-in config settings.", e);
            }
        }

        // Cache and return our configuration:
        // NOTE: The configuration will be stored in the ConfigCache under the key GATKConfig.class.
        //       This means that any future call to getOrCreate for this GATKConfig.class will return
        //       Not only the configuration itself, but also the overrides as specified in userConfigFileProperties
        return ConfigUtils.getOrCreate(configClass, userConfigFileProperties);
    }

    /**
     * Injects system properties from the given configuration file.
     * System properties are specified by the presence of the {@link SystemProperty} annotation.
     * @param config The {@link GATKConfig} object from which to inject system properties.
     */
    public static <T extends Config> void injectSystemPropertiesFromConfig(final T config) {
        
        Utils.nonNull(config);

        // Get our system properties:
        final Map<String, String> properties = getSystemPropertiesFromConfig(config);

        // Set our properties:
        injectToSystemProperties(properties);
    }

    /**
     * Logs all the parameters in the given {@link Config} object at {@link Level#DEBUG}
     * @param config A {@link Config} object from which to log all parameters and values.
     * @param <T> any {@link Config} type to use to log all configuration information.
     */
    public static <T extends Config> void logConfigFields(final T config) {
        logConfigFields(config, Level.DEBUG);
    }

    /**
     * Logs all the parameters in the given {@link Config} object at the given {@link Level}
     * @param config A {@link Config} object from which to log all parameters and values.
     * @param level The log {@link Level} at which to log the data in {@code config}
     * @param <T> any {@link Config} type to use to log all configuration information.
     */
    public static <T extends Config> void logConfigFields(final T config, final Level level) {

        Utils.nonNull(config);
        Utils.nonNull(level);

        // Only continue in this method here if we would log the given level:
        if ( !logger.isEnabled(level) ) {
            return;
        }

        // This is gross and uses reflection to get all methods in the given config class
        // and then interrogates those methods for internal data on the config parameters.

        // We have to match our interfaces to the config interface that we're actually using.
        // It's not as simple as using getDeclaredMethods on the Class object because we'll get
        // a LOT of extraneous stuff that we don't care about.
        for ( final Class<?> classInterface : config.getClass().getInterfaces() ) {

            // If we have an interface that is a child of the OWNER Config interface, then
            // we must have the right object.
            if (Config.class.isAssignableFrom(classInterface)) {

                logger.log(level, "Configuration file values: ");

                // Now we cycle through our interface methods, resolve parameter names,
                // and log their values at the given level:
                for (final Method propertyMethod : classInterface.getDeclaredMethods()) {

                    // Get the property name:
                    String propertyName = propertyMethod.getName();

                    final Config.Key key = propertyMethod.getAnnotation(Config.Key.class);
                    if (key != null) {
                        propertyName = key.value();
                    }

                    final StringBuilder sb = new StringBuilder();
                    sb.append('\t');
                    sb.append(propertyName);
                    sb.append(" = ");

                    try {
                        sb.append(propertyMethod.invoke(config, new Object[]{}));
                    } catch (final IllegalAccessException ex) {
                        throw new GATKException("Could not access the config getter: " +
                                config.getClass().getSimpleName() + "." +
                                propertyMethod.getName(), ex);

                    } catch (final InvocationTargetException ex) {
                        throw new GATKException("Could not invoke the config getter: " +
                                config.getClass().getSimpleName() + "." +
                                propertyMethod.getName(), ex);
                    }

                    logger.log(level, sb.toString());
                }
            }
        }
    }
}
