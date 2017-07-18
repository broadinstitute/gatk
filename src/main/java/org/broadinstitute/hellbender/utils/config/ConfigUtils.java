package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
 * A class to contain configuration utility methods.
 * Created by jonn on 7/19/17.
 */
public class ConfigUtils {

    static private final Logger logger = LogManager.getLogger(ConfigUtils.class);

    // This class just has static methods to help with configuration, so no need for this constructor.
    private ConfigUtils() {}

    /**
     * Whether we have already set the config factory variable defaults.
     */
    private static boolean hasSetConfigFactoryVariableDefaults = false;

    /**
     * Sets the {@link org.aeonbits.owner.ConfigFactory} variables so that it knows about
     * the variable paths for config files.
     */
    public static final void setConfigFactoryVariableDefaults() {

        // You only need to do this once.

        if ( !hasSetConfigFactoryVariableDefaults ) {

            ArrayList<String> propertyNames = new ArrayList<>();

            // Get the classes from which we need to look for sources:
            Class<?>[] configurationClasses = new Class<?>[] {
                    GATKConfig.class
            };

            // Create a regex to use to look for variables in the Sources annotation:
            Pattern p = Pattern.compile("\\$\\{(.*)}");

            // Loop through our classes and grab any sources with variables in there:
            for ( Class<?> clazz : configurationClasses ) {
                Set<Class<?>> interfaces = new HashSet<>(Arrays.asList(clazz.getInterfaces()));

                // Make sure that we get config classes here:
                if ( interfaces.contains(Config.class) ||
                         interfaces.contains(Accessible.class) ||
                         interfaces.contains(Mutable.class) ||
                         interfaces.contains(Reloadable.class) ) {
                    Config.Sources annotation = clazz.getAnnotation(Config.Sources.class);

                    String[] annotationValues = annotation.value();

                    for ( String val : annotationValues ) {

                        Matcher m = p.matcher(val);
                        if ( m.find() ) {
                            propertyNames.add(m.group(1));
                        }
                    }
                }
            }

            // Grab the system properties:
            Properties systemProperties = System.getProperties();

            // Grab the environment properties:
            Map<String, String> environmentProperties = System.getenv();

            // Make sure that if our property isn't in the system and environment
            // properties, that we set it to a neutral value that will not contain
            // anything (so that the property will fall back into the next value).
            for (String property : propertyNames) {

                if ((!environmentProperties.keySet().contains(property)) &&
                        (!systemProperties.containsKey(property))) {

                    ConfigFactory.setProperty(property, "/dev/null");
                }
            }

            hasSetConfigFactoryVariableDefaults = true;
        }
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
    public static final String getConfigFilenameFromArgs( final ArrayList<String> args, final String configFileOption ) {

        Utils.nonNull(args);
        Utils.nonNull(configFileOption);

        String configFileName = null;

        for ( int i = 0 ; i < args.size() ; ++i ) {
            if (args.get(i).compareTo(configFileOption) == 0) {

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

                    String message = "ERROR: Configuration file not given after config file option specified: " + configFileOption;
                    System.err.println(message);
                    throw new UserException.BadInput(message);
                }
            }
        }

        return configFileName;
    }

    /**
     * Get the given key from a {@link GATKConfig} class.
     * @param key The key in the {@link GATKConfig} to get.
     * @param <T> The type of the key to obtain.
     * @return The value for the specified {@code key}.
     */
    @SuppressWarnings("unchecked")
    public static <T> T getFromConfig(String key) {

        Utils.nonNull(key);

        // Get our configuration:
        final GATKConfig config = ConfigCache.getOrCreate( GATKConfig.class );

        try {
            for (Method m : GATKConfig.class.getDeclaredMethods()) {
                if (m.getName().equals(key)) {
                    return (T) m.invoke(config);
                }
            }
        }
        catch (IllegalAccessException e) {
            // We shouldn't need to worry about this one, but just in case:
            logger.fatal("Could not access the given key from GATKConfig class: " + key);
        }
        catch (InvocationTargetException e) {
            logger.fatal("Could not resolve the given key from GATKConfig class: " + key);
        }

        // This should never happen.
        // If it does, throw an exception:
        throw new RuntimeException("Unable to retrieve property, " + key + ",from configuration: " + GATKConfig.class.toString());
    }

    /**
     * Get the configuration filename the command-line (if it exists) and create a configuration for it.
     * Configuration type defaults to {@link GATKConfig}
     * Removes the configuration filenames and configuration file options from the given {@code argList}.
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     */
    public static final <T extends Config> void initializeConfigurationsFromCommandLineArgs(final ArrayList<String> argList,
                                                                                            String configFileOption) {
        initializeConfigurationsFromCommandLineArgs(
                argList,
                configFileOption,
                GATKConfig.class
        );
    }

    /**
     * Get the configuration filename the command-line (if it exists) and create a configuration for it of the given type.
     * Removes the configuration filenames and configuration file options from the given {@code argList}.
     * Also sets system-level properties from the system config file.
     * @param argList The list of arguments from which to read the config file.
     * @param configFileOption The command-line option specifying the main configuration file.
     * @param configClazz The class of the configuration file to instantiate.
     */
    public static final <T extends Config> void initializeConfigurationsFromCommandLineArgs(final ArrayList<String> argList,
                                                                                            String configFileOption,
                                                                                            Class<? extends T> configClazz) {
        Utils.nonNull(argList);
        Utils.nonNull(configFileOption);

        // Get main config from args:
        final String configFileName = getConfigFilenameFromArgs( argList, configFileOption );

        // Set the config file to be the one we want to use from the command-line:
        ConfigFactory.setProperty( GATKConfig.CONFIG_FILE_VARIABLE_NAME, configFileName );
        T gatkConfig = ConfigCache.getOrCreate(configClazz);

//        // Alternate way to load the config file:
//        T gatkConfig = ConfigUtils.initializeConfigurationAsProperties( configFileName, configClazz );

        // To start with we inject our system properties to ensure they are defined for downstream components:
        ConfigUtils.injectSystemPropertiesFromConfig( gatkConfig, GATKConfig.SYSTEM_PROPERTY_PREFIX );
    }

    /**
     * Initializes and returns the configuration as specified by {@code configFileName}
     * Also caches this configuration in the {@link ConfigCache} for use elsewhere.
     * @param configFileName The name of the file from which to initialize the configuration
     * @param configClass The type of configuration in which to interpret the given {@code configFileName}
     * @return The configuration instance implementing {@link GATKConfig} containing any overrides in the given file.
     */
    public static final <T extends Config> T initializeConfigurationAsProperties(final String configFileName, Class<? extends T> configClass) {

        Utils.nonNull(configClass);

        // Get a place to store our properties:
        final Properties userConfigFileProperties = new Properties();

        // Try to get the config from the specified file:
        if ( configFileName != null ) {

            try {
                final FileInputStream userConfigFileInputStream = new FileInputStream(configFileName);
                userConfigFileProperties.load(userConfigFileInputStream);

                if (configFileName != null) {
                    System.out.println("Found " + configClass.getSimpleName() + " Configuration File: " + configFileName);
                }

            } catch (final FileNotFoundException e) {
                System.err.println("WARNING: unable to find specified " + configClass.getSimpleName() + " configuration file: "
                        + configFileName + " - defaulting to built-in config settings.");
            }
            catch (final IOException e) {
                System.err.println("WARNING: unable to load specified " + configClass.getSimpleName() + " configuration file: "
                        + configFileName + " - defaulting to built-in config settings.");
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
     * @param config The {@link GATKConfig} object from which to inject system properties.
     */
    public static final <T extends Config> void injectSystemPropertiesFromConfig(T config) {
        injectSystemPropertiesFromConfig(config, "");
    }

    /**
     * Injects system properties from the given configuration file.
     * @param config The {@link GATKConfig} object from which to inject system properties.
     * @param propertyPrefix Prefix that is required for a property in {@code config} to be imported to System properties.
     */
    public static final <T extends Config> void injectSystemPropertiesFromConfig(T config, String propertyPrefix) {
        
        Utils.nonNull(config);
        Utils.nonNull(propertyPrefix);

        // This is gross and uses reflection to get all methods in the given config class
        // and then interrogates those methods for internal data on the config parameters.

        // We have to match our interfaces to the config interface that we're actually using.
        // It's not as simple as using getDeclaredMethods on the Class object because we'll get
        // a LOT of extraneous stuff that we don't care about.
        for ( Class<?> classInterface : config.getClass().getInterfaces() ){

            // If we have an interface that is a child of the OWNER Config interface, then
            // we must have the right object.
            if (Config.class.isAssignableFrom(classInterface)) {

                // Now we cycle through our interface methods, resolve parameter names,
                // and set the values in the system.
                for (Method propertyMethod : classInterface.getDeclaredMethods()) {

                    // Get the property name:
                    String propertyName = propertyMethod.getName();

                    Config.Key key = propertyMethod.getAnnotation(Config.Key.class);
                    if (key != null) {
                        propertyName = key.value();
                    }

                    // Only set properties beginning with the given prefix:
                    if (!propertyName.startsWith(propertyPrefix)) {
                        continue;
                    }

                    // Get the value of the property into a stringbuilder:
                    StringBuilder sb = new StringBuilder();
                    try {
                        sb.append(propertyMethod.invoke( config, new Object[]{} ));
                    } catch (IllegalAccessException ex) {
                        throw new RuntimeException("Could not access the config getter: " +
                                config.getClass().getSimpleName() + "." +
                                propertyMethod.getName(), ex);

                    } catch (InvocationTargetException ex) {
                        throw new RuntimeException("Could not invoke the config getter: " +
                                config.getClass().getSimpleName() + "." +
                                propertyMethod.getName(), ex);
                    }

                    // Set our property:
                    System.setProperty(propertyName, sb.toString());
                }
            }
        }
    }
}
