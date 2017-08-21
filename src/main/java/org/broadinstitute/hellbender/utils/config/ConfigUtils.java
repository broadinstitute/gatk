package org.broadinstitute.hellbender.utils.config;

import com.google.common.annotations.VisibleForTesting;
import org.aeonbits.owner.*;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
     * Whether we have already set the config factory variable defaults.
     */
    private static boolean hasSetConfigFactoryVariableDefaults = false;

    /**
     * Value to set each variable for configuration file paths when the variable
     * has not been set in either Java System properties or environment properties.
     */
    @VisibleForTesting
    static final String NO_PATH_VARIABLE_VALUE = "/dev/null";

    /**
     * Sets the {@link org.aeonbits.owner.ConfigFactory} variables so that it knows about
     * the variable paths for config files.
     */
    public static synchronized void setConfigFactoryVariableDefaults() {

        // You only need to do this once.

        if ( !hasSetConfigFactoryVariableDefaults ) {

            // Get the classes from which we need to look for sources:
            // We need to enumerate all of our config sources here:
            final List<Class<?>> configurationClasses = Arrays.asList( GATKConfig.class );

            // Get our property names from the configuration classes:
            final List<String> filenameProperties = getConfigPathVariableNamesFromConfigClasses(configurationClasses);

            // Check if the filename properties exist and if they do not, add them
            // to the ConfigFactory.
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(filenameProperties);

            hasSetConfigFactoryVariableDefaults = true;
        }
        else {
            logger.error("Attempted to setConfigFactoryVariableDefaults more than once!");
        }
    }

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

        // Make sure that if our property isn't in the system and environment
        // properties, that we set it to a neutral value that will not contain
        // anything (so that the property will fall back into the next value).
        for (final String property : filenameProperties) {

            if ((!environmentProperties.keySet().contains(property)) &&
                    (!systemProperties.containsKey(property))) {

                logger.warn("Config path variable not found: " + property +
                            " - setting value to default empty variable: " + NO_PATH_VARIABLE_VALUE);
                ConfigFactory.setProperty(property, NO_PATH_VARIABLE_VALUE);
            }
            else {
                logger.info("Config path variable found: " + property + " - will search for config here.");
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

        final List<String> propertyNames = new ArrayList<>();

        // Create a regex to use to look for variables in the Sources annotation:
        final Pattern p = Pattern.compile("\\$\\{(.*)}");

        // Loop through our classes and grab any sources with variables in there:
        for ( final Class<?> clazz : configurationClasses ) {

            // Make sure that we get config classes here:
            if ( Config.class.isAssignableFrom(clazz) ) {
                final Config.Sources annotation = clazz.getAnnotation(Config.Sources.class);

                if ( annotation != null ) {
                    final String[] annotationValues = annotation.value();

                    for (final String val : annotationValues) {

                        final Matcher m = p.matcher(val);
                        if (m.find()) {
                            propertyNames.add(m.group(1));
                        }
                    }
                }
            }
        }

        return propertyNames;
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
    public static String getConfigFilenameFromArgs( final ArrayList<String> args, final String configFileOption ) {

        Utils.nonNull(args);
        Utils.nonNull(configFileOption);

        String configFileName = null;

        for ( int i = 0 ; i < args.size() ; ++i ) {
            if (args.get(i).equals(configFileOption)) {

                // Get rid of the command-line argument name:
                args.remove(i);

                if ( i < args.size() ) {

                    // Get and remove the specified config file:
                    configFileName = args.remove(i);
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
     * Get the configuration filename from the command-line (if it exists) and create a configuration for it.
     * Configuration type defaults to {@link GATKConfig}
     * Removes the configuration filenames and configuration file options from the given {@code argList}.
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     */
    public static void initializeConfigurationsFromCommandLineArgs(final ArrayList<String> argList,
                                                                   final String configFileOption) {
        initializeConfigurationsFromCommandLineArgs(
                argList,
                configFileOption,
                GATKConfig.class
        );
    }

    /**
     * Get the configuration from filename the command-line (if it exists) and create a configuration for it of the given type.
     * Removes the configuration filenames and configuration file options from the given {@code argList}.
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     * @param configClass The class of the configuration file to instantiate.
     */
    public static <T extends Config> void initializeConfigurationsFromCommandLineArgs(final ArrayList<String> argList,
                                                                                      final String configFileOption,
                                                                                      final Class<? extends T> configClass) {
        Utils.nonNull(argList);
        Utils.nonNull(configFileOption);

        // Get main config from args:
        final String configFileName = getConfigFilenameFromArgs( argList, configFileOption );

        // Set the config path if we've specified it:
        if ( configFileName != null ){
            ConfigFactory.setProperty( GATKConfig.CONFIG_FILE_VARIABLE_NAME, configFileName );
        }

        // Set the config file to be the one we want to use from the command-line:
        final T gatkConfig = ConfigCache.getOrCreate(configClass);

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
        return ConfigCache.getOrCreate(configClass, userConfigFileProperties);
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
