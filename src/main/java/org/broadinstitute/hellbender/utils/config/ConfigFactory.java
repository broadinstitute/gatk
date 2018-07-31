package org.broadinstitute.hellbender.utils.config;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Log;
import org.aeonbits.owner.Config;
import org.aeonbits.owner.ConfigCache;
import org.aeonbits.owner.Factory;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A singleton class to act as a user interface for loading configuration files from {@link org.aeonbits.owner}.
 * This class wraps functionality in the {@link org.aeonbits.owner} configuration utilities to be a more GATK-specific
 * interface.
 * Created by jonn on 7/19/17.
 */
public final class ConfigFactory {

    private static final Logger logger = LogManager.getLogger(ConfigFactory.class);

    //=======================================
    // Singleton members / methods:
    private static final ConfigFactory instance;

    static {
        instance = new ConfigFactory();
    }

    /**
     * @return An instance of this {@link ConfigFactory}, which can be used to create a configuration.
     */
    public static ConfigFactory getInstance() {
        return instance;
    }

    // This class is a singleton, so no public construction.
    private ConfigFactory() {}

    //=======================================

    /**
     * A regex to use to look for variables in the Sources annotation
     */
    private static final Pattern sourcesAnnotationPathVariablePattern = Pattern.compile("\\$\\{(.*)}");

    /**
     * Value to set each variable for configuration file paths when the variable
     * has not been set in either Java System properties or environment properties.
     */
    @VisibleForTesting
    static final String NO_PATH_VARIABLE_VALUE = "/dev/null";

    /**
     * A set to keep track of the classes we've already resolved for configuration path purposes:
     */
    private final Set<Class<? extends Config>> alreadyResolvedPathVariables = new HashSet<>();

    // =================================================================================================================

    /**
     * Checks each of the given {@code filenameProperties} for if they are defined in system {@link System#getProperties()}
     * or environment {@link System#getenv()} properties.  If they are not, this method will set them in the
     * {@link org.aeonbits.owner.ConfigFactory} to an empty file path so the {@link org.aeonbits.owner.ConfigFactory} will know to try to resolve them as
     * variables at load-time (and not as raw paths).
     * @param filenameProperties A {@link List} of filename properties as specified in {@link Config} {@link org.aeonbits.owner.Config.Sources} annotations to check for existence in system and environment properties.
     */
    @VisibleForTesting
    void checkFileNamePropertyExistenceAndSetConfigFactoryProperties(final List<String> filenameProperties) {
        // Grab the system properties:
        final Properties systemProperties = System.getProperties();

        // Grab the environment properties:
        final Map<String, String> environmentProperties = System.getenv();

        // Make sure that if our property isn't in the system, environment, and ConfigFactory
        // properties, that we set it to a neutral value that will not contain
        // anything (so that the property will fall back into the next value).
        for (final String property : filenameProperties) {

            if ( environmentProperties.containsKey(property) ) {
                logger.debug("Config path variable found in Environment Properties: " + property + "=" + environmentProperties.get(property) + " - will search for config here.");
            }
            else if ( systemProperties.containsKey(property) ) {
                logger.debug("Config path variable found in System Properties: " + property + "=" + systemProperties.get(property) + " - will search for config here.");
            }
            else if ( org.aeonbits.owner.ConfigFactory.getProperties().containsKey(property) ) {
                logger.debug("Config path variable found in Config Factory Properties(probably from the command-line): " + property + "=" + org.aeonbits.owner.ConfigFactory.getProperty(property) + " - will search for config here.");
            }
            else {
                logger.debug("Config path variable not found: " + property +
                        " - setting value to default empty variable: " + (NO_PATH_VARIABLE_VALUE == null ? "null" : String.valueOf(NO_PATH_VARIABLE_VALUE)) );
                org.aeonbits.owner.ConfigFactory.setProperty(property, NO_PATH_VARIABLE_VALUE);
            }
        }
    }

    /**
     * Get a list of the config file variables from the given {@link Config} classes.
     * @param configurationClasses A list of configuration classes from which to extract variable names in their {@link org.aeonbits.owner.Config.Sources}.
     * @return A list of variables in the {@link org.aeonbits.owner.Config.Sources} of the given {@code configurationClasses}
     */
    @VisibleForTesting
    List<String> getConfigPathVariableNamesFromConfigClasses(final List<Class<?>> configurationClasses) {

        final List<String> configPathVariableNames = new ArrayList<>();

        // Loop through our classes and grab any sources with variables in there:
        for ( final Class<?> clazz : ClassUtils.getClassesOfType(Config.class, configurationClasses) ) {
            @SuppressWarnings("unchecked")
            final Class<? extends Config> castedClass = (Class<? extends Config>) clazz;
            configPathVariableNames.addAll( getSourcesAnnotationPathVariables(castedClass));
        }

        return configPathVariableNames;
    }

    /**
     * Get a list of the config file variables from the given {@link Config} class.
     * @param configClass A configuration class from which to extract variable names in its {@link org.aeonbits.owner.Config.Sources}.
     * @return A list of variables in the {@link org.aeonbits.owner.Config.Sources} of the given {@code configClass}
     */
    @VisibleForTesting
    <T extends Config> List<String> getSourcesAnnotationPathVariables(final Class<? extends T> configClass) {

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
     * Injects the given property to the System Properties.
     * Validates that this property was set after setting it.
     * This will NOT override properties that already exist in the system.
     * @param properties A {@link Map} of key, value pairs of properties to add to the System Properties.
     */
    @VisibleForTesting
    void injectToSystemProperties(final Map<String, String> properties) {

        // Get our current system properties:
        final Properties systemProperties = System.getProperties();

        for ( final Map.Entry<String, String> entry : properties.entrySet() ) {

            // If we have this key in our system already, we do NOT set it:
            if ( systemProperties.containsKey(entry.getKey()) ) {
                logger.debug("System property already exists.  Not overriding: " + entry.getKey());
                continue;
            }

            System.setProperty(entry.getKey(), entry.getValue());

            // Test for validation that it worked:
            final String propertyValueThatWasSet = System.getProperty(entry.getKey());
            if (propertyValueThatWasSet == null) {
                throw new GATKException("Unable to set System Property (" + entry.getKey() + "=" + entry.getValue() + ")!");
            }

            if (!propertyValueThatWasSet.equals(entry.getValue())) {
                throw new GATKException("System Property corrupted (" + entry.getKey() + "!=" + entry.getValue() + " -> " + propertyValueThatWasSet + ")!");
            }
        }
    }

    // =================================================================================================================

    /**
     * Quick way to get the GATK configuration.
     * @return The GATK Configuration.
     */
    public GATKConfig getGATKConfig() {
        return getOrCreate( GATKConfig.class );
    }

    /**
     * Dump the configuration to a file that can be easily read into {@link Properties}.
     * @param config Configuration instance to dump.
     * @param outFilePath {@link Path} to output location.
     * @param <T> Some configuration class that extends {@link Config}.
     */
    public static <T extends Config> void dumpConfigSettings(final T config, final Path outFilePath ) {
        final LinkedHashMap<String, Object> configMap = getConfigMap(config, false);

        final Properties properties = new Properties();
        properties.putAll(convertConfigMapToStringStringMap(configMap));

        final Date d = new Date();

        try ( final OutputStream outputStream = Files.newOutputStream(outFilePath, StandardOpenOption.CREATE_NEW) ) {
            properties.store(outputStream, "Created from " + config.getClass().getSimpleName() + " at " +
                    new SimpleDateFormat("HH.mm.ss").format(d) + " on " +
                    new SimpleDateFormat("yyyy.MM.dd").format(d));
        }
        catch (final Exception ex) {
            throw new GATKException("Could not write config (" + config.getClass().getTypeName() + ") to file: " + outFilePath, ex);
        }
    }

    /**
     * Wrapper around {@link org.aeonbits.owner.ConfigFactory#create(Class, Map[])} which will ensure that
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
    public <T extends Config> T create(final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(clazz);

        resolvePathVariables(clazz);

        return org.aeonbits.owner.ConfigFactory.create(clazz, imports);
    }

    /**
     * Wrapper around {@link ConfigCache#getOrCreate(Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * Gets from the cache or create, an instance of the given class using the given imports.
     * The factory used to create new instances is the static {@link org.aeonbits.owner.ConfigFactory#INSTANCE}.
     *
     * @param clazz     the interface extending from {@link Config} that you want to instantiate.
     * @param imports   additional variables to be used to resolve the properties.
     * @param <T>       type of the interface.
     * @return          an object implementing the given interface, that can be taken from the cache,
     *                  which maps methods to property values.
     */
    public <T extends Config> T getOrCreate(final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(clazz);

        resolvePathVariables(clazz);

        return ConfigCache.getOrCreate(clazz, imports);
    }

    /**
     * Wrapper around {@link ConfigCache#getOrCreate(Factory, Class, Map[])} which will ensure that
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
    public <T extends Config> T getOrCreate(final Factory factory, final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(factory);
        Utils.nonNull(clazz);

        resolvePathVariables(clazz);

        return ConfigCache.getOrCreate(factory, clazz, imports);
    }

    /**
     * Wrapper around {@link ConfigCache#getOrCreate(Object, Class, Map[])} which will ensure that
     * path variables specified in {@link org.aeonbits.owner.Config.Sources} annotations are resolved prior
     * to creation.
     *
     * Gets from the cache or create, an instance of the given class using the given imports.
     * The factory used to create new instances is the static {@link org.aeonbits.owner.ConfigFactory#INSTANCE}.
     *
     * @param key       the key object to be used to identify the instance in the cache.
     * @param clazz     the interface extending from {@link Config} that you want to instantiate.
     * @param imports   additional variables to be used to resolve the properties.
     * @param <T>       type of the interface.
     * @return          an object implementing the given interface, that can be taken from the cache,
     *                  which maps methods to property values.
     */
    public <T extends Config> T getOrCreate(final Object key, final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(key);
        Utils.nonNull(clazz);

        resolvePathVariables(clazz);

        return ConfigCache.getOrCreate(key, clazz, imports);
    }

    /**
     * Wrapper around {@link ConfigCache#getOrCreate(Factory, Object, Class, Map[])} which will ensure that
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
    public <T extends Config> T getOrCreate(final Factory factory, final Object key,
                                                   final Class<? extends T> clazz, final Map<?, ?>... imports) {

        Utils.nonNull(factory);
        Utils.nonNull(key);
        Utils.nonNull(clazz);

        resolvePathVariables(clazz);

        return ConfigCache.getOrCreate(key, clazz, imports);
    }

    private synchronized <T extends Config> void resolvePathVariables(final Class<? extends T> clazz) {
        if ( !alreadyResolvedPathVariables.contains(clazz) ) {
            checkFileNamePropertyExistenceAndSetConfigFactoryProperties(
                    getSourcesAnnotationPathVariables(clazz)
            );

            alreadyResolvedPathVariables.add(clazz);
        }
    }

    /**
     * Wrapper around {@link ConfigCache#get(Object)}.
     * This method is here to complete the interface for getting {@link Config} objects.
     *
     * Gets from the cache the {@link Config} instance identified by the given key.
     *
     * @param key       the key object to be used to identify the instance in the cache.
     * @param <T>       type of the interface.
     * @return          the {@link Config} object from the cache if exists, or <tt>null</tt> if it doesn't.
     */
    public <T extends Config> T get(final Object key) {
        Utils.nonNull(key);
        return ConfigCache.get(key);
    }

    /**
     * Get the configuration file name from the given arguments.
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

                // Get the config file name (ignoring other arguments):
                if ( ((i+1) < args.length) && (!args[i+1].startsWith("-")) ) {
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
     * Get the configuration filename from the command-line (if it exists) and create a configuration for it.
     * Configuration type defaults to {@link GATKConfig}
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     */
    public synchronized void initializeConfigurationsFromCommandLineArgs(final String[] argList,
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
    public synchronized <T extends Config> void initializeConfigurationsFromCommandLineArgs(final String[] argList,
                                                                                            final String configFileOption,
                                                                                            final Class<? extends T> configClass) {
        Utils.nonNull(argList);
        Utils.nonNull(configFileOption);
        Utils.nonNull(configClass);

        // Get main config from args:
        final String configFileName = getConfigFilenameFromArgs( argList, configFileOption );

        // Load the configuration from the given file:
        final T configuration = getOrCreateConfigFromFile(configFileName, configClass);

        // To start with we inject our system properties to ensure they are defined for downstream components:
        injectSystemPropertiesFromConfig( configuration );
    }

    @VisibleForTesting
    synchronized <T extends Config> T getOrCreateConfigFromFile(final String configFileName, final Class<? extends T> configClass) {

        // Set the config path if we've specified it:
        if ( configFileName != null ){
            org.aeonbits.owner.ConfigFactory.setProperty( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, configFileName );
        }

        // Set the config file to be the one we want to use from the command-line:
        return ConfigFactory.getInstance().getOrCreate(configClass);
    }

    @VisibleForTesting
    synchronized <T extends Config> T createConfigFromFile(final String configFileName, final Class<? extends T> configClass) {

        // Set the config path if we've specified it:
        if ( configFileName != null ){
            org.aeonbits.owner.ConfigFactory.setProperty( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, configFileName );
        }

        // Set the config file to be the one we want to use from the command-line:
        return ConfigFactory.getInstance().create(configClass);
    }

    /**
     * Injects system properties from the given configuration file.
     * System properties are specified by the presence of the {@link SystemProperty} annotation.
     * This will NOT override properties that already exist in the system.
     * @param config The {@link GATKConfig} object from which to inject system properties.
     */
    public synchronized <T extends Config> void injectSystemPropertiesFromConfig(final T config) {
        
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
        logConfigFields(config, Log.LogLevel.DEBUG);
    }

    /**
     * Gets all system properties from the given {@link Config}-derived object.
     * System properties are denoted via the presence of the {@link SystemProperty} annotation.
     * @param config A {@link Config}-derived object from which to read system properties.
     * @param <T> A {@link Config}-derived type from which to read System Properties.
     * @return A properties {@link Map} representing all System Properties in the given {@code config}.
     */
    @VisibleForTesting
    static <T extends Config> LinkedHashMap<String, String> getSystemPropertiesFromConfig(final T config) {

        Utils.nonNull(config);

        final LinkedHashMap<String, String> properties = new LinkedHashMap<>();

        for ( final Map.Entry<String, Object> entry : getConfigMap(config, true).entrySet() ) {
            properties.put(entry.getKey(), String.valueOf(entry.getValue()));
        }

        return properties;
    }

    /**
     * Logs all the parameters in the given {@link Config} object at the given {@link Log.LogLevel}
     * @param config A {@link Config} object from which to log all parameters and values.
     * @param logLevel The log {@link htsjdk.samtools.util.Log.LogLevel} at which to log the data in {@code config}
     * @param <T> any {@link Config} type to use to log all configuration information.
     */
    public static <T extends Config> void logConfigFields(final T config, final Log.LogLevel logLevel) {

        Utils.nonNull(config);
        Utils.nonNull(logLevel);

        final Level level = LoggingUtils.levelToLog4jLevel(logLevel);

        // Only continue in this method here if we would log the given level:
        if ( !logger.isEnabled(level) ) {
            return;
        }

        logger.log(level, "Configuration file values: ");
        for ( final Map.Entry<String, Object> entry : getConfigMap(config, false).entrySet() ) {
            logger.log(level, "\t" + entry.getKey() + " = " + entry.getValue());
        }
    }

    @VisibleForTesting
    static <T extends Config> LinkedHashMap<String, Object> getConfigMap( final T config, final boolean onlySystemProperties ) {
        final LinkedHashMap<String, Object> configMap = new LinkedHashMap<>();

        // This is gross and uses reflection to get all methods in the given config class
        // and then interrogates those methods for internal data on the config parameters.

        // We have to match our interfaces to the config interface that we're actually using.
        // It's not as simple as using getDeclaredMethods on the Class object because we'll get
        // a LOT of extraneous stuff that we don't care about.
        // So we make sure we have an interface that is a child of the OWNER Config interface.
        for ( final Class<?> classInterface : ClassUtils.getClassesOfType(Config.class, Arrays.asList(config.getClass().getInterfaces())) ) {

            // Now we cycle through our interface methods, resolve parameter names,
            // and log their values at the given level:
            for (final Method propertyMethod : classInterface.getDeclaredMethods()) {

                // Get the property name:
                String propertyName = propertyMethod.getName();

                // Get the real property name if we've overwritten it with a key:
                final Config.Key key = propertyMethod.getAnnotation(Config.Key.class);
                if (key != null) {
                    propertyName = key.value();
                }

                try {
                    if ( onlySystemProperties ) {
                        if ( propertyMethod.isAnnotationPresent(SystemProperty.class) ) {
                            configMap.put(propertyName, propertyMethod.invoke(config));
                        }
                    }
                    else {
                        configMap.put(propertyName, propertyMethod.invoke(config));
                    }
                } catch (final IllegalAccessException ex) {
                    throw new GATKException("Could not access the config getter: " +
                            config.getClass().getSimpleName() + "." +
                            propertyMethod.getName(), ex);

                } catch (final InvocationTargetException ex) {
                    throw new GATKException("Could not invoke the config getter: " +
                            config.getClass().getSimpleName() + "." +
                            propertyMethod.getName(), ex);
                }
            }
        }

        return configMap;
    }

    private static LinkedHashMap<String, String> convertConfigMapToStringStringMap(final LinkedHashMap<String, Object> configMap) {
        final LinkedHashMap<String, String> map = new LinkedHashMap<>();

        for ( final Map.Entry<String, Object> entry : configMap.entrySet() ){
            map.put( entry.getKey(), entry.getValue().toString() );
        }

        return map;
    }
}
