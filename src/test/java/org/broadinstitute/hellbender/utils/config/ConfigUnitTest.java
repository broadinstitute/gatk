package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Config;
import org.aeonbits.owner.ConfigCache;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit test for GATK configuration file handling.
 *
 * Created by jonn on 7/19/17.
 */
public class ConfigUnitTest extends GATKBaseTest {

    private final static String testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/config/";

    // ================================================================================
    // Helper Methods:
    // ================================================================================

    private static void listAndStoreConfigToStdOut(final Accessible config) throws IOException {

        config.list(System.out);
        System.out.println();

        config.store(System.out, "");
        System.out.println();

        config.storeToXML(System.out, "");
        System.out.println();
    }

    private static void validateAndClearProperty(final String key, final Object value) {

        Assert.assertEquals( String.valueOf( System.getProperty(key) ), String.valueOf(value) );

        System.clearProperty(key);
    }

    // ================================================================================
    // Tests:
    // ================================================================================

    @DataProvider
    Object[][] createArgsAndConfigFileOptionsBadInput() {
        return new Object[][] {
                // Main and Config specified, no file given
                {
                        new String[] {"main", "testArg", "--config"},
                        "--config",
                },
                // Main and Config specified, no file given
                {
                        new String[] {"main","--config"},
                        "--config",
                },
                // Config specified, no file given
                {
                        new String[] {"--config"},
                        "--config",
                },
                // Config specified other argument given
                {
                        new String[] {"main", "--config", "--otherArgument"},
                        "--config",
                },
                // Config specified other argument given
                {
                        new String[] {"main", "--config", "-X"},
                        "--config",
                }
        };
    }

    @Test(dataProvider = "createArgsAndConfigFileOptionsBadInput",
            expectedExceptions = UserException.BadInput.class)
    public void testGetConfigFilenameFromArgsBadInput( final String[] args,
                                        final String configFileOption) {

        ConfigFactory.getConfigFilenameFromArgs(args, configFileOption);
    }

    @Test
    public void testDumpConfigSettings() {
        // Test with our basic test class:
        final SystemTestConfig testConfig = ConfigFactory.getInstance().getOrCreate(SystemTestConfig.class);

        final File outFile = getSafeNonExistentFile("tmpConfigOutFile.config");

        // Dump the config:
        ConfigFactory.dumpConfigSettings(testConfig, outFile.toPath());

        // Read the config:
        final Properties properties = new Properties();
        try ( final InputStream inputStream = Files.newInputStream(outFile.toPath(), StandardOpenOption.READ) ) {
            properties.load(inputStream);
        }
        catch (final Exception ex) {
            throw new GATKException("ERROR OPENING TEST FILE: " + outFile, ex);
        }

        // Validate the properties
        // (as strings because the properties are Strings and the Config values are actual types)
        Assert.assertEquals(properties.getProperty("booleanDefTrue"),             String.valueOf(testConfig.booleanDefTrue()) );
        Assert.assertEquals(properties.getProperty("booleanDefFalse"),            String.valueOf(testConfig.booleanDefFalse()) );
        Assert.assertEquals(properties.getProperty("intDef207"),                  String.valueOf(testConfig.intDef207()) );
        Assert.assertEquals(properties.getProperty("listOfStringTest"),           String.valueOf(testConfig.listOfStringTest()) );
        Assert.assertEquals(properties.getProperty("system.Boolean.Def.True"),    String.valueOf(testConfig.systemBooleanDefTrue()) );
        Assert.assertEquals(properties.getProperty("system.Boolean.Def.False"),   String.valueOf(testConfig.systemBooleanDefFalse()) );
        Assert.assertEquals(properties.getProperty("system.Int.Def.207"),         String.valueOf(testConfig.systemIntDef207()) );
        Assert.assertEquals(properties.getProperty("system.List.Of.String.Test"), String.valueOf(testConfig.systemListOfStringTest()) );
    }

    @DataProvider
    Object[][] createArgsAndConfigFileOptions() {
        return new Object[][] {
                // incorrect config option given, no file name found
                {
                        new String[] {"main","--zonfigurati","DUMMY_FILE","END"},
                        "--config",
                        null,
                },
                // correct config option given, file given with trailing args
                {
                        new String[] {"main","--config","DUMMY_FILE","END"},
                        "--config",
                        "DUMMY_FILE",
                },
                // correct config option given, file given with no trailing args
                {
                        new String[] {"main","END","--config","DUMMY_FILE"},
                        "--config",
                        "DUMMY_FILE",
                },
        };
    }

    @Test(dataProvider= "createArgsAndConfigFileOptions")
    public void testGetConfigFilenameFromArgs( final String[] args,
                                        final String configFileOption,
                                        final String expectedFilename) {

        final String configFileName = ConfigFactory.getConfigFilenameFromArgs(args, configFileOption);

        Assert.assertEquals(configFileName, expectedFilename);
    }

    @Test
    public void testSystemConfigurationOwnerMethods() {
        // Test with our basic test class:
        final SystemTestConfig testConfig = org.aeonbits.owner.ConfigFactory.create(SystemTestConfig.class);

        // Make sure that none of these properties are set already:
        Assert.assertEquals( System.getProperty("booleanDefTrue"),      null );
        Assert.assertEquals( System.getProperty("booleanDefFalse"),     null );
        Assert.assertEquals( System.getProperty("intDef207"),           null );
        Assert.assertEquals( System.getProperty("listOfStringTest"),    null );
        Assert.assertEquals( System.getProperty("system.Boolean.Def.True"),   null );
        Assert.assertEquals( System.getProperty("system.Boolean.Def.False"),  null );
        Assert.assertEquals( System.getProperty("system.Int.Def.207"),        null );
        Assert.assertEquals( System.getProperty("system.List.Of.String.Test"),null );

        ConfigFactory.getInstance().injectSystemPropertiesFromConfig(testConfig);

        //Verify that the system contains the properties we expect:
        validateAndClearProperty("booleanDefTrue",       null);
        validateAndClearProperty("booleanDefFalse",      null);
        validateAndClearProperty("intDef207",            null);
        validateAndClearProperty("listOfStringTest",     null);

        validateAndClearProperty("system.Boolean.Def.True",    testConfig.systemBooleanDefTrue());
        validateAndClearProperty("system.Boolean.Def.False",   testConfig.systemBooleanDefFalse());
        validateAndClearProperty("system.Int.Def.207",         testConfig.systemIntDef207());
        validateAndClearProperty("system.List.Of.String.Test", testConfig.systemListOfStringTest());
    }

    @Test
    public void testInitializeConfigurationsFromCommandLineArgs() {

        // Make sure that none of these properties are set already:
        Assert.assertEquals( System.getProperty("booleanDefTrue"),      null );
        Assert.assertEquals( System.getProperty("booleanDefFalse"),     null );
        Assert.assertEquals( System.getProperty("intDef207"),           null );
        Assert.assertEquals( System.getProperty("listOfStringTest"),    null );
        Assert.assertEquals( System.getProperty("system.Boolean.Def.True"),   null );
        Assert.assertEquals( System.getProperty("system.Boolean.Def.False"),  null );
        Assert.assertEquals( System.getProperty("system.Int.Def.207"),        null );
        Assert.assertEquals( System.getProperty("system.List.Of.String.Test"),null );

        // Test with our basic test class:
        ConfigFactory.getInstance().initializeConfigurationsFromCommandLineArgs(
                new String[] {"--config", "dummyFile"},
                "--config",
                SystemTestConfig.class
        );

        final SystemTestConfig testConfig = ConfigFactory.getInstance().getOrCreate(SystemTestConfig.class);

        //Verify that the system contains the properties we expect:
        validateAndClearProperty("booleanDefTrue",       null);
        validateAndClearProperty("booleanDefFalse",      null);
        validateAndClearProperty("intDef207",            null);
        validateAndClearProperty("listOfStringTest",     null);

        validateAndClearProperty("system.Boolean.Def.True",    testConfig.systemBooleanDefTrue());
        validateAndClearProperty("system.Boolean.Def.False",   testConfig.systemBooleanDefFalse());
        validateAndClearProperty("system.Int.Def.207",         testConfig.systemIntDef207());
        validateAndClearProperty("system.List.Of.String.Test", testConfig.systemListOfStringTest());
    }

    @Test
    public void testOwnerConfiguration() throws IOException {

        // Test with our basic test class:
        final BasicTestConfig basicTestConfig = org.aeonbits.owner.ConfigFactory.create(BasicTestConfig.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfig);

        Assert.assertEquals(basicTestConfig.booleanDefFalse(), false);
        Assert.assertEquals(basicTestConfig.booleanDefTrue(), true);
        Assert.assertEquals(basicTestConfig.intDef207(), 207);
        Assert.assertEquals(basicTestConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string1", "string2", "string3", "string4")));

    }

    @Test
    public void testOwnerConfigurationWithClassPathOverrides() throws IOException {

        // Test with the class that overrides on the class path:
        final BasicTestConfigWithClassPathOverrides basicTestConfigWithClassPathOverrides =
                org.aeonbits.owner.ConfigFactory.create(BasicTestConfigWithClassPathOverrides.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfigWithClassPathOverrides);

        Assert.assertEquals(basicTestConfigWithClassPathOverrides.booleanDefFalse(), true);
        Assert.assertEquals(basicTestConfigWithClassPathOverrides.booleanDefTrue(), false);
        Assert.assertEquals(basicTestConfigWithClassPathOverrides.intDef207(), 702);
        Assert.assertEquals(basicTestConfigWithClassPathOverrides.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
    }

    @Test
    public void testOwnerConfigurationWithClassPathOverridesAndVariableFileInput() throws IOException {

        // Start with the name of the properties file to copy:
        final String overrideFilename = "AdditionalTestOverrides.properties";

        // Create a temporary folder in which to place the config file:
        final File outputDir = BaseTest.createTempDir("testOwnerConfigurationWithClassPathOverridesAndVariableFileInput");
        outputDir.deleteOnExit();

        // Put the known config file in the new directory:
        Files.copy(new File("src/test/resources/org/broadinstitute/hellbender/utils/config/" + overrideFilename).toPath(),
                new File(outputDir.getAbsolutePath() + File.separator + overrideFilename).toPath(),
                StandardCopyOption.REPLACE_EXISTING);

        // Assert that our config property is not set:
        Assert.assertEquals( org.aeonbits.owner.ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);

        // Set our file location here:
        org.aeonbits.owner.ConfigFactory.setProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME, outputDir.getAbsolutePath() + File.separator + overrideFilename);

        // Test with the class that overrides on the class path:
        final BasicTestConfigWithClassPathOverridesAndVariableFile basicTestConfigWithClassPathOverridesAndVariableFile =
                org.aeonbits.owner.ConfigFactory.create(BasicTestConfigWithClassPathOverridesAndVariableFile.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfigWithClassPathOverridesAndVariableFile);

        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefFalse(), true);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefTrue(), false);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.intDef207(), 999);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.customBoolean().booleanValue(), false);

        // Reset the config factory:
        org.aeonbits.owner.ConfigFactory.clearProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME);
        Assert.assertEquals(org.aeonbits.owner.ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);
    }

    @Test
    public void testOwnerConfigurationWithClassPathOverridesAndVariableFileInput_ParameterizedConfigFileDoesNotExist() throws IOException {

        // Start with the name of the properties file to copy:
        final String overrideFilename = "AdditionalTestOverrides.properties";

        // Create a temporary folder in which to place the config file:
        final File outputDir = BaseTest.createTempDir("testOwnerConfigurationWithClassPathOverridesAndVariableFileInput");
        outputDir.deleteOnExit();

        // Assert that our config property is not set:
        Assert.assertEquals( org.aeonbits.owner.ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);

        // Set our file location here:
        org.aeonbits.owner.ConfigFactory.setProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME, outputDir.getAbsolutePath() + File.separator + overrideFilename);

        // Test with the class that overrides on the class path:
        final BasicTestConfigWithClassPathOverridesAndVariableFile basicTestConfigWithClassPathOverridesAndVariableFile =
                org.aeonbits.owner.ConfigFactory.create(BasicTestConfigWithClassPathOverridesAndVariableFile.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfigWithClassPathOverridesAndVariableFile);

        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefFalse(), true);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefTrue(), false);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.intDef207(), 702);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.customBoolean().booleanValue(), true);

        // Reset the config factory:
        org.aeonbits.owner.ConfigFactory.clearProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME);
        Assert.assertEquals(org.aeonbits.owner.ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);
    }

    @DataProvider
    Object[][] createSourcesFilenamesAndSystemPropertiesFlags() {
        return new Object[][] {
                // Add no values to system properties
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(false, false, false)
                },
                // Add first value to system properties
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(true, false, false)
                },
                // Add first and second values to system properties
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(true, true, false)
                },
                // Add all values to system properties
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(true, true, true)
                },
        };
    }

    @Test(expectedExceptions = NullPointerException.class)
    public void testOwnerConfigFactorySetProperty() {
        org.aeonbits.owner.ConfigFactory.setProperty(UUID.randomUUID().toString(), null);
    }

    @Test(dataProvider = "createSourcesFilenamesAndSystemPropertiesFlags")
    public void testConfigUtilsCheckFileNamePropertyExistenceAndSetConfigFactoryProperties(final List<String> sourcesFilenames,
                                                                                    final List<Boolean> addToSystemProperties) {
        Assert.assertEquals(sourcesFilenames.size(), addToSystemProperties.size());

        // Ensure that we don't have any of these in our system properties or environment:
        for( final String fileNameVar : sourcesFilenames ) {
            Assert.assertEquals(System.getenv(fileNameVar), null);
            Assert.assertEquals(System.getProperty(fileNameVar), null);
        }

        // Add in the relevant values to our system properties:
        for ( int i = 0 ; i < sourcesFilenames.size() ; ++i ) {
            if ( addToSystemProperties.get(i) ) {
                System.setProperty(sourcesFilenames.get(i), sourcesFilenames.get(i));
            }
        }

        // Set up our config factory values:
        ConfigFactory.getInstance().checkFileNamePropertyExistenceAndSetConfigFactoryProperties(sourcesFilenames);

        // Now we check them:
        for ( int i = 0 ; i < sourcesFilenames.size() ; ++i ) {

            // If we added them to the properties, we don't need to set the names in
            // the ConfigFactory.  Therefore they should be null:
            if ( addToSystemProperties.get(i) ) {
                Assert.assertEquals(org.aeonbits.owner.ConfigFactory.getProperty(sourcesFilenames.get(i)), null);
            }
            // If the filename didn't exist in the system properties, then we did need to add them to the
            // ConfigFactory properties:
            else {
                Assert.assertEquals(org.aeonbits.owner.ConfigFactory.getProperty(sourcesFilenames.get(i)), ConfigFactory.NO_PATH_VARIABLE_VALUE);
            }
        }

        // Remove these from our System Properties and ConfigFactory properties:
        for ( int i = 0 ; i < sourcesFilenames.size() ; ++i ) {
            if ( addToSystemProperties.get(i) ) {
                System.clearProperty(sourcesFilenames.get(i));
            }
            else {
                org.aeonbits.owner.ConfigFactory.clearProperty(sourcesFilenames.get(i));
                Assert.assertEquals(org.aeonbits.owner.ConfigFactory.getProperty(sourcesFilenames.get(i)), null);
            }
        }
    }

    @DataProvider
    Object[][] createConfigListsAndResults() {
        return new Object[][] {
                // Empty List:
                {
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                // Zero Paths:
                {
                        Collections.singletonList( BasicTestConfig.class ),
                        Collections.emptyList()
                },
                // One Path:
                {
                        Arrays.asList( BasicTestConfig.class, BasicTestConfigWithClassPathOverridesAndVariableFile.class  ),
                        Collections.singletonList( BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME )
                },
                // Two Path:
                {
                        Collections.singletonList( GATKConfig.class ),
                        Arrays.asList( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, GATKConfig.CONFIG_FILE_VARIABLE_CLASS_PATH )
                },
                // Three Paths:
                {
                        Arrays.asList( BasicTestConfig.class, GATKConfig.class, BasicTestConfigWithClassPathOverridesAndVariableFile.class ),
                        Arrays.asList( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, GATKConfig.CONFIG_FILE_VARIABLE_CLASS_PATH, BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME )
                },
        };
    }

    @Test(dataProvider= "createConfigListsAndResults")
    public void testConfigUtilsGetConfigPathVariableNamesFromConfigClasses(final List<Class<?>> configurationClasses, final List<String> expectedFilePathVariables ) {
        Assert.assertEquals( ConfigFactory.getInstance().getConfigPathVariableNamesFromConfigClasses(configurationClasses), expectedFilePathVariables );
    }

    @DataProvider
    Object[][] createConfigFileAndVariableNames() {
        return new Object[][] {
                // Class with no @Sources annotation:
                {
                        BasicTestConfig.class,
                        Collections.emptyList()
                },
                // Class with one @Sources annotation:
                {
                        BasicTestConfigWithClassPathOverridesAndVariableFile.class,
                        Collections.singletonList( BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME )
                },
                // Class with @Sources annotation with a file and a class path:
                {
                        GATKConfig.class,
                        Arrays.asList( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, GATKConfig.CONFIG_FILE_VARIABLE_CLASS_PATH )
                },
        };
    }

    @Test(dataProvider = "createConfigFileAndVariableNames")
    public void testGetSourcesAnnotationPathVariables(final Class<? extends Config> configClass, final List<String> expectedFilePathVariables ) {
            Assert.assertEquals( ConfigFactory.getInstance().getSourcesAnnotationPathVariables(configClass), expectedFilePathVariables );
    }

    @DataProvider
    Object[][] createPropertiesToInjectIntoSystem() {

        final Map<String, String> propertyMap1 = new HashMap<>();

        final Map<String, String> propertyMap2 = new HashMap<>();
        propertyMap2.put("TESTKEY1", "TESTVAL1");

        final Map<String, String> propertyMap3 = new HashMap<>();
        propertyMap2.put("TESTKEY1", "TESTVAL1");
        propertyMap2.put("TESTKEY2", "TESTVAL2");

        final Map<String, String> propertyMap4 = new HashMap<>();
        propertyMap2.put("TESTKEY1", "TESTVAL1");
        propertyMap2.put("TESTKEY2", "TESTVAL2");
        propertyMap2.put("TESTKEY3", "TESTVAL3");

        // Four sets of properties to add in to the system to test with.
        // One set is empty, the others have key:value pairs:
        return new Object[][] {
                {
                        propertyMap1
                },
                {
                        propertyMap2
                },
                {
                        propertyMap3
                },
                {
                        propertyMap4
                },
        };
    }

    @Test(dataProvider = "createPropertiesToInjectIntoSystem")
    public void testConfigUtilsInjectToSystemProperties(final Map<String, String> propertiesToInject) {

        // Ensure our properties are not there to start with:
        for ( final String key : propertiesToInject.keySet() ) {
            Assert.assertEquals(System.getProperty(key), null);
        }

        // Inject the stuff:
        ConfigFactory.getInstance().injectToSystemProperties(propertiesToInject);

        // Verify our injection:
        for ( final String key : propertiesToInject.keySet() ) {
            Assert.assertEquals(System.getProperty(key), propertiesToInject.get(key));
        }

        // Reset our properties after:
        for ( final String key : propertiesToInject.keySet() ) {
            System.clearProperty(key);
            Assert.assertEquals(System.getProperty(key), null);
        }
    }

    @Test
    public void testDerivedClassProperties() throws IOException {

        // Test with our basic test class:
        final ChildClassConfig childClassConfig = org.aeonbits.owner.ConfigFactory.create(ChildClassConfig.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(childClassConfig);

        Assert.assertEquals(childClassConfig.booleanDefTrue(), false);
        Assert.assertEquals(childClassConfig.booleanDefFalse(), true);
        Assert.assertEquals(childClassConfig.intDef207(), 702);
        Assert.assertEquals(childClassConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(childClassConfig.customBoolean().booleanValue(), true);
        Assert.assertEquals(childClassConfig.newCustomBooleanThatDefaultsToTrue().booleanValue(), false);
    }

    @Test
    public void testDerivedClassAsParentClassProperties() throws IOException {

        // Test with our basic test class:
        final BasicTestConfigWithClassPathOverridesAndVariableFile parentClassConfig = org.aeonbits.owner.ConfigFactory.create(ChildClassConfig.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(parentClassConfig);

        Assert.assertEquals(parentClassConfig.booleanDefTrue(), false);
        Assert.assertEquals(parentClassConfig.booleanDefFalse(), true);
        Assert.assertEquals(parentClassConfig.intDef207(), 702);
        Assert.assertEquals(parentClassConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(parentClassConfig.customBoolean().booleanValue(), true);
    }

    @Test
    public void testDerivedClassMappingToParentAndSelf() throws IOException {

        Assert.assertEquals(
                org.aeonbits.owner.ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME),
                null
        );

        // This is `assertTrue` because of conflicting type comparisons:
        Assert.assertTrue(ConfigCache.get(ChildClassConfig.class) == null);
        Assert.assertTrue(ConfigCache.get(BasicTestConfigWithClassPathOverridesAndVariableFile.class) == null);

        final BasicTestConfigWithClassPathOverridesAndVariableFile parentClassConfig = ConfigFactory.getInstance().getOrCreate(ChildClassConfig.class);
        final ChildClassConfig childClassConfig = ConfigFactory.getInstance().getOrCreate(ChildClassConfig.class);

        ConfigCache.add(BasicTestConfigWithClassPathOverridesAndVariableFile.class, childClassConfig);

        // List parent properties for inspection:
        listAndStoreConfigToStdOut(parentClassConfig);

        Assert.assertEquals(parentClassConfig.booleanDefTrue(), false);
        Assert.assertEquals(parentClassConfig.booleanDefFalse(), true);
        Assert.assertEquals(parentClassConfig.intDef207(), 702);
        Assert.assertEquals(parentClassConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(parentClassConfig.customBoolean().booleanValue(), true);

        // List child properties for inspection:
        listAndStoreConfigToStdOut(childClassConfig);

        Assert.assertEquals(childClassConfig.booleanDefTrue(), false);
        Assert.assertEquals(childClassConfig.booleanDefFalse(), true);
        Assert.assertEquals(childClassConfig.intDef207(), 702);
        Assert.assertEquals(childClassConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(childClassConfig.customBoolean().booleanValue(), true);
        Assert.assertEquals(childClassConfig.newCustomBooleanThatDefaultsToTrue().booleanValue(), false);

        // Now check that getting them back produces the same object:
        // This is `assertTrue` because of conflicting type comparisons:
        Assert.assertTrue(
                ConfigCache.get(BasicTestConfigWithClassPathOverridesAndVariableFile.class) ==  ConfigCache.get(ChildClassConfig.class)
        );

        // Now we clean the cache and the path variables:
        org.aeonbits.owner.ConfigFactory.clearProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME);
        ConfigCache.remove(BasicTestConfigWithClassPathOverridesAndVariableFile.class);
        ConfigCache.remove(ChildClassConfig.class);

        // Verify that we've cleaned up:
        Assert.assertEquals(
                org.aeonbits.owner.ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME),
                null
        );

        // This is `assertTrue` because of conflicting type comparisons:
        Assert.assertTrue(ConfigCache.get(ChildClassConfig.class) == null);
        Assert.assertTrue(ConfigCache.get(BasicTestConfigWithClassPathOverridesAndVariableFile.class) == null);
    }

    @Test
    public void ensurePackagedGATKConfigDefaultsAreSameAsHardcodedDefaults() {

        final Map<String, String> defaultValueMap = getDefaultValuesFromConfig(GATKConfig.class);

        final GATKConfig gatkConfig = ConfigFactory.getInstance().createConfigFromFile(
                packageMainResourcesDir + "utils/config/" + "GATKConfig.properties",
                GATKConfig.class
        );

        final LinkedHashMap<String, Object> configFileDefaultsMap = ConfigFactory.getConfigMap(gatkConfig, false);

        for ( final Map.Entry<String, String> entry : defaultValueMap.entrySet() ) {
            Assert.assertTrue( configFileDefaultsMap.containsKey(entry.getKey()) );

            // Special case for lists:
            if ( Collection.class.isAssignableFrom( configFileDefaultsMap.get(entry.getKey()).getClass() ) ) {

                // Get our delimiter:
                final String delimiter = getFieldDelimiterFromGATKConfig( entry.getKey() );

                // Is this a regular collection or a sorted set?
                if ( SortedSet.class.isAssignableFrom( configFileDefaultsMap.get(entry.getKey()).getClass() ) ||
                     !Set.class.isAssignableFrom( configFileDefaultsMap.get(entry.getKey()).getClass() ) ) {

                    // Get our actual collection string without spaces between elements and delimiters:
                    @SuppressWarnings("unchecked")
                    final String actualCollectionString = ((Collection<Object>) configFileDefaultsMap.get(entry.getKey())).stream()
                            .map(Object::toString)
                            .collect(
                                Collectors.joining(delimiter)
                            );

                    // Get our expected collection string without spaces between elements and delimiters:
                    final String configCollectionDefaultsString = entry.getValue().replaceAll("\\s*" + delimiter + "\\s*", delimiter);

                    // Check for equality:
                    Assert.assertEquals( actualCollectionString, configCollectionDefaultsString );
                }
                else {
                    // REGULAR SET HERE:

                    // Check that each value exists in the set:
                    for ( final String value : entry.getValue().split(delimiter) ) {
                        final String strippedValue = value.replaceAll("^\\s+", "").replaceAll("\\s+$", "");

                        @SuppressWarnings("unchecked")
                        final Set<String> stringValueSet = ((Set<Object>)configFileDefaultsMap.get(entry.getKey())).stream()
                                .map(Object::toString)
                                .collect(Collectors.toSet());

                        Assert.assertTrue( stringValueSet.contains(strippedValue) );
                    }
                }
            }
            else {
                Assert.assertEquals( String.valueOf(configFileDefaultsMap.get( entry.getKey() )), entry.getValue() );
            }
        }
    }

    private Map<String, String> getDefaultValuesFromConfig(final Class<GATKConfig> gatkConfigClass) {

        final Map<String, String> defaultValueMap = new HashMap<>();

        // Now we cycle through our interface methods, resolve parameter names,
        // and log their values at the given level:
        for (final Method propertyMethod : gatkConfigClass.getDeclaredMethods()) {

            final String defaultValue;
            final Config.DefaultValue defVal = propertyMethod.getAnnotation( Config.DefaultValue.class );
            if (defVal != null) {

                defaultValue = defVal.value();

                // Get the property name:
                String propertyName = propertyMethod.getName();

                // Get the real property name if we've overwritten it with a key:
                final Config.Key key = propertyMethod.getAnnotation(Config.Key.class);
                if ( key != null ) {
                    propertyName = key.value();
                }

                defaultValueMap.put(propertyName, defaultValue);
            }
        }

        return defaultValueMap;
    }

    private String getFieldDelimiterFromGATKConfig(final String parameterName) {

        String delimiter = ",";

        for (final Method propertyMethod : GATKConfig.class.getDeclaredMethods()) {

            // Get the property name:
            String propertyName = propertyMethod.getName();

            // Get the real property name if we've overwritten it with a key:
            final Config.Key key = propertyMethod.getAnnotation(Config.Key.class);
            if ( key != null ) {
                propertyName = key.value();
            }

            // Check to see if this is the property we're interested in:
            if ( parameterName.equals(propertyName) ) {
                // Get the delimiter:
                final Config.Separator separator = propertyMethod.getAnnotation(Config.Separator.class);
                if ( separator != null ) {
                    delimiter = separator.value();
                }
            }
        }

        return delimiter;
    }
}
