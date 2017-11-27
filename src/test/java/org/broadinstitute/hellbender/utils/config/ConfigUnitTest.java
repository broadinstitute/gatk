package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Config;
import org.aeonbits.owner.ConfigCache;
import org.aeonbits.owner.ConfigFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;

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
    // Data Providers:
    // ================================================================================

    @DataProvider
    Object[][] createArgsAndConfigFileOptions() {
        return new Object[][] {
                {
                    new String[] {"main","--zonfigurati","DUMMY_FILE","END"},
                    "--config",
                    null,
                },
                {
                    new String[] {"main","--config","DUMMY_FILE","END"},
                    "--config",
                    "DUMMY_FILE",
                },
                {
                    new String[] {"main","END","--config","DUMMY_FILE"},
                    "--config",
                    "DUMMY_FILE",
                },
        };
    }

    @DataProvider
    Object[][] createArgsAndConfigFileOptionsBadInput() {
        return new Object[][] {
                {
                        new String[] {"main", "testArg", "--config"},
                        "--config",
                },
                {
                        new String[] {"main","--config"},
                        "--config",
                },
                {
                        new String[] {"--config"},
                        "--config",
                },
        };
    }

    @DataProvider
    Object[][] createArgsAndConfigFileOptionsWithExpectedArgs() {
        return new Object[][] {
                {
                        new String[] {"main","--zonfigurati","DUMMY_FILE","END"},
                        "--config",
                        new String[] {"main","--zonfigurati","DUMMY_FILE","END"},
                },
                {
                        new String[] {"main","--config","DUMMY_FILE","END"},
                        "--config",
                        new String[] {"main","END"},
                },
                {
                        new String[] {"main","END","--config","DUMMY_FILE"},
                        "--config",
                        new String[] {"main","END"},
                },
                {
                        new String[] {"--config","DUMMY_FILE","main","END"},
                        "--config",
                        new String[] {"main","END"},
                },
                {
                        new String[] {"main","--zonfigurati","DUMMY_FILE","END"},
                        null,
                        new String[] {"main","--zonfigurati","DUMMY_FILE","END"},
                },
                {
                        new String[] {"main", StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION,"DUMMY_FILE","END"},
                        null,
                        new String[] {"main","END"},
                },
                {
                        new String[] {"main","END",StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION,"DUMMY_FILE"},
                        null,
                        new String[] {"main","END"},
                },
                {
                        new String[] {StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION, "DUMMY_FILE","main","END"},
                        null,
                        new String[] {"main","END"},
                },
        };
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

    @DataProvider
    Object[][] createConfigFileAndVariableNames() {
        return new Object[][] {
                {
                        BasicTestConfig.class,
                        Collections.emptyList()
                },
                {
                        BasicTestConfigWithClassPathOverridesAndVariableFile.class,
                        Collections.singletonList( BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME )
                },
                {
                        GATKConfig.class,
                        Arrays.asList( GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME, GATKConfig.CONFIG_FILE_VARIABLE_CLASS_PATH )
                },
        };
    }

    @DataProvider
    Object[][] createSourcesFilenamesAndSystemPropertiesFlags() {
//        final List<String> sourcesFilenames,
//        final List<Boolean> addToSystemProperties
        return new Object[][] {
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(false, false, false)
                },
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(true, false, false)
                },
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(true, true, false)
                },
                {
                        Arrays.asList("testPropertyFileName1", "testPropertyFileName2", "testPropertyFileName3"),
                        Arrays.asList(true, true, true)
                },
        };
    }

//    <T extends Config> void testConfigUtilsGetSystemPropertiesFromConfig(T config, Map<String, String> expectedProperties) {
    @DataProvider
    Object[][] createConfigAndExpectedProperties() {

        final Map<String, String> propertyMap1 = new HashMap<>();

        propertyMap1.put("gatk_stacktrace_on_user_exception", "true");
        propertyMap1.put("samjdk.use_async_io_read_samtools", "false");
        propertyMap1.put("samjdk.use_async_io_write_samtools", "true");
        propertyMap1.put("samjdk.use_async_io_write_tribble", "false");
        propertyMap1.put("samjdk.compression_level", "1");
        propertyMap1.put("snappy.disable", "true");
        propertyMap1.put( "spark.kryoserializer.buffer.max", "512m" );
        propertyMap1.put( "spark.driver.maxResultSize", "0" );
        propertyMap1.put( "spark.driver.userClassPathFirst", "true" );
        propertyMap1.put( "spark.io.compression.codec", "lzf" );
        propertyMap1.put( "spark.yarn.executor.memoryOverhead", "600" );
        propertyMap1.put( "spark.driver.extraJavaOptions", "" );
        propertyMap1.put( "spark.executor.extraJavaOptions", "" );

        final Map<String, String> propertyMap2 = new HashMap<>();

        return new Object [][] {
                {
                        ConfigUtils.create(GATKConfig.class),
                        propertyMap1
                },
                {
                        ConfigUtils.create(BasicTestConfig.class),
                        propertyMap2
                },
        };
    }

//    void testConfigUtilsInjectToSystemProperties(Map<String, String> propertiesToInject) {
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

    // ================================================================================
    // Tests:
    // ================================================================================

    @Test(dataProvider = "createArgsAndConfigFileOptionsBadInput",
            expectedExceptions = UserException.BadInput.class)
    void testGetConfigFilenameFromArgsBadInput( final String[] args,
                                        final String configFileOption) {

        ConfigUtils.getConfigFilenameFromArgs(args, configFileOption);
    }

    @Test(dataProvider= "createArgsAndConfigFileOptions")
    void testGetConfigFilenameFromArgs( final String[] args,
                                        final String configFileOption,
                                        final String expectedFilename) {

        final String outFileName = ConfigUtils.getConfigFilenameFromArgs(args, configFileOption);

        Assert.assertEquals(expectedFilename, outFileName);
    }

    @Test
    void testInitializeConfiguration() throws IOException {
        final String inputPropertiesFile = testResourceDir + "AdditionalTestOverrides.properties";

        final BasicTestConfig basicTestConfig = ConfigUtils.initializeConfigurationAsProperties(inputPropertiesFile, BasicTestConfig.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfig);

        // Check for our values:
        Assert.assertEquals(basicTestConfig.booleanDefFalse(), false);
        Assert.assertEquals(basicTestConfig.booleanDefTrue(), true);
        Assert.assertEquals(basicTestConfig.intDef207(), 999);
        Assert.assertEquals(basicTestConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string1", "string2", "string3", "string4")));

    }

    @Test
    void testSystemConfiguration() {
        // Test with our basic test class:
        final SystemTestConfig testConfig = ConfigFactory.create(SystemTestConfig.class);

        // Make sure that none of these properties are set already:
        Assert.assertEquals( System.getProperty("booleanDefTrue"),      null );
        Assert.assertEquals( System.getProperty("booleanDefFalse"),     null );
        Assert.assertEquals( System.getProperty("intDef207"),           null );
        Assert.assertEquals( System.getProperty("listOfStringTest"),    null );
        Assert.assertEquals( System.getProperty("system.Boolean.Def.True"),   null );
        Assert.assertEquals( System.getProperty("system.Boolean.Def.False"),  null );
        Assert.assertEquals( System.getProperty("system.Int.Def.207"),        null );
        Assert.assertEquals( System.getProperty("system.List.Of.String.Test"),null );

        ConfigUtils.injectSystemPropertiesFromConfig(testConfig);

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
    void testOwnerConfiguration() throws IOException {

        // Test with our basic test class:
        final BasicTestConfig basicTestConfig = ConfigFactory.create(BasicTestConfig.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfig);

        Assert.assertEquals(basicTestConfig.booleanDefFalse(), false);
        Assert.assertEquals(basicTestConfig.booleanDefTrue(), true);
        Assert.assertEquals(basicTestConfig.intDef207(), 207);
        Assert.assertEquals(basicTestConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string1", "string2", "string3", "string4")));

    }

    @Test
    void testOwnerConfigurationWithClassPathOverrides() throws IOException {

        // Test with the class that overrides on the class path:
        final BasicTestConfigWithClassPathOverrides basicTestConfigWithClassPathOverrides =
                ConfigFactory.create(BasicTestConfigWithClassPathOverrides.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfigWithClassPathOverrides);

        Assert.assertEquals(basicTestConfigWithClassPathOverrides.booleanDefFalse(), true);
        Assert.assertEquals(basicTestConfigWithClassPathOverrides.booleanDefTrue(), false);
        Assert.assertEquals(basicTestConfigWithClassPathOverrides.intDef207(), 702);
        Assert.assertEquals(basicTestConfigWithClassPathOverrides.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
    }

    @Test
    void testOwnerConfigurationWithClassPathOverridesAndVariableFileInput() throws IOException {

        // Start with the name of the properties file to copy:
        final String overrideFilename = "AdditionalTestOverrides.properties";

        // Create a temporary folder in which to place the config file:
        final File outputDir = Files.createTempDirectory("testOwnerConfigurationWithClassPathOverridesAndVariableFileInput").toAbsolutePath().toFile();
        outputDir.deleteOnExit();

        // Put the known config file in the new directory:
        Files.copy(new File("src/test/resources/org/broadinstitute/hellbender/utils/config/" + overrideFilename).toPath(),
                new File(outputDir.getAbsolutePath() + File.separator + overrideFilename).toPath(),
                StandardCopyOption.REPLACE_EXISTING);

        // Assert that our config property is not set:
        Assert.assertEquals( ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);

        // Set our file location here:
        ConfigFactory.setProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME, outputDir.getAbsolutePath() + File.separator + overrideFilename);

        // Test with the class that overrides on the class path:
        final BasicTestConfigWithClassPathOverridesAndVariableFile basicTestConfigWithClassPathOverridesAndVariableFile =
                ConfigFactory.create(BasicTestConfigWithClassPathOverridesAndVariableFile.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfigWithClassPathOverridesAndVariableFile);

        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefFalse(), true);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefTrue(), false);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.intDef207(), 999);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.customBoolean().booleanValue(), false);

        // Reset the config factory:
        ConfigFactory.clearProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME);
        Assert.assertEquals(ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);
    }

    @Test
    void testOwnerConfigurationWithClassPathOverridesAndVariableFileInput_ParameterizedConfigFileDoesNotExist() throws IOException {

        // Start with the name of the properties file to copy:
        final String overrideFilename = "AdditionalTestOverrides.properties";

        // Create a temporary folder in which to place the config file:
        final File outputDir = Files.createTempDirectory("testOwnerConfigurationWithClassPathOverridesAndVariableFileInput").toAbsolutePath().toFile();
        outputDir.deleteOnExit();

        // Assert that our config property is not set:
        Assert.assertEquals( ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);

        // Set our file location here:
        ConfigFactory.setProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME, outputDir.getAbsolutePath() + File.separator + overrideFilename);

        // Test with the class that overrides on the class path:
        final BasicTestConfigWithClassPathOverridesAndVariableFile basicTestConfigWithClassPathOverridesAndVariableFile =
                ConfigFactory.create(BasicTestConfigWithClassPathOverridesAndVariableFile.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(basicTestConfigWithClassPathOverridesAndVariableFile);

        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefFalse(), true);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.booleanDefTrue(), false);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.intDef207(), 702);
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(basicTestConfigWithClassPathOverridesAndVariableFile.customBoolean().booleanValue(), true);

        // Reset the config factory:
        ConfigFactory.clearProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME);
        Assert.assertEquals(ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME), null);
    }

    @Test(dataProvider = "createSourcesFilenamesAndSystemPropertiesFlags")
    void testConfigUtilsCheckFileNamePropertyExistenceAndSetConfigFactoryProperties(final List<String> sourcesFilenames,
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
        ConfigUtils.checkFileNamePropertyExistenceAndSetConfigFactoryProperties(sourcesFilenames);

        // Now we check them:
        for ( int i = 0 ; i < sourcesFilenames.size() ; ++i ) {

            // If we added them to the properties, we don't need to set the names in
            // the ConfigFactory.  Therefore they should be null:
            if ( addToSystemProperties.get(i) ) {
                Assert.assertEquals(ConfigFactory.getProperty(sourcesFilenames.get(i)), null);
            }
            // If the filename didn't exist in the system properties, then we did need to add them to the
            // ConfigFactory properties:
            else {
                Assert.assertEquals(ConfigFactory.getProperty(sourcesFilenames.get(i)), ConfigUtils.NO_PATH_VARIABLE_VALUE);
            }
        }

        // Remove these from our System Properties and ConfigFactory properties:
        for ( int i = 0 ; i < sourcesFilenames.size() ; ++i ) {
            if ( addToSystemProperties.get(i) ) {
                System.clearProperty(sourcesFilenames.get(i));
            }
            else {
                ConfigFactory.clearProperty(sourcesFilenames.get(i));
                Assert.assertEquals(ConfigFactory.getProperty(sourcesFilenames.get(i)), null);
            }
        }
    }

    @Test(dataProvider= "createConfigListsAndResults")
    void testConfigUtilsGetConfigPathVariableNamesFromConfigClasses(final List<Class<?>> configurationClasses, final List<String> expectedFilePathVariables ) {
        Assert.assertEquals( ConfigUtils.getConfigPathVariableNamesFromConfigClasses(configurationClasses), expectedFilePathVariables );
    }

    @Test(dataProvider = "createConfigFileAndVariableNames")
    void testGetSourcesAnnotationPathVariables(final Class<? extends Config> configClass, final List<String> expectedFilePathVariables ) {
            Assert.assertEquals( ConfigUtils.getSourcesAnnotationPathVariables(configClass), expectedFilePathVariables );
    }

    @Test(dataProvider = "createConfigAndExpectedProperties")
    <T extends Config> void testConfigUtilsGetSystemPropertiesFromConfig(final T config, final Map<String, String> expectedProperties) {
        Assert.assertEquals( ConfigUtils.getSystemPropertiesFromConfig(config), expectedProperties );
    }

    @Test(dataProvider = "createPropertiesToInjectIntoSystem")
    void testConfigUtilsInjectToSystemProperties(final Map<String, String> propertiesToInject) {

        // Ensure our properties are not there to start with:
        for ( final String key : propertiesToInject.keySet() ) {
            Assert.assertEquals(System.getProperty(key), null);
        }

        // Inject the stuff:
        ConfigUtils.injectToSystemProperties(propertiesToInject);

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
    void testDerivedClassProperties() throws IOException {

        // Test with our basic test class:
        final ChildClassConfig childClassConfig = ConfigFactory.create(ChildClassConfig.class);

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
    void testDerivedClassAsParentClassProperties() throws IOException {

        // Test with our basic test class:
        final BasicTestConfigWithClassPathOverridesAndVariableFile parentClassConfig = ConfigFactory.create(ChildClassConfig.class);

        // List properties for inspection:
        listAndStoreConfigToStdOut(parentClassConfig);

        Assert.assertEquals(parentClassConfig.booleanDefTrue(), false);
        Assert.assertEquals(parentClassConfig.booleanDefFalse(), true);
        Assert.assertEquals(parentClassConfig.intDef207(), 702);
        Assert.assertEquals(parentClassConfig.listOfStringTest(), new ArrayList<>(Arrays.asList("string4", "string3", "string2", "string1")));
        Assert.assertEquals(parentClassConfig.customBoolean().booleanValue(), true);
    }

    @Test
    void testDerivedClassMappingToParentAndSelf() throws IOException {

        Assert.assertEquals(
                ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME),
                null
        );

        // This is `assertTrue` because of conflicting type comparisons:
        Assert.assertTrue(ConfigCache.get(ChildClassConfig.class) == null);
        Assert.assertTrue(ConfigCache.get(BasicTestConfigWithClassPathOverridesAndVariableFile.class) == null);

        final BasicTestConfigWithClassPathOverridesAndVariableFile parentClassConfig = ConfigUtils.getOrCreate(ChildClassConfig.class);
        final ChildClassConfig childClassConfig = ConfigUtils.getOrCreate(ChildClassConfig.class);

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
        ConfigFactory.clearProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME);
        ConfigCache.remove(BasicTestConfigWithClassPathOverridesAndVariableFile.class);
        ConfigCache.remove(ChildClassConfig.class);

        // Verify that we've cleaned up:
        Assert.assertEquals(
                ConfigFactory.getProperty(BasicTestConfigWithClassPathOverridesAndVariableFile.CONFIG_FILE_VARIABLE_FILE_NAME),
                null
        );

        // This is `assertTrue` because of conflicting type comparisons:
        Assert.assertTrue(ConfigCache.get(ChildClassConfig.class) == null);
        Assert.assertTrue(ConfigCache.get(BasicTestConfigWithClassPathOverridesAndVariableFile.class) == null);
    }

    @Test(dataProvider = "createArgsAndConfigFileOptionsWithExpectedArgs")
    void testRemoveConfigOptionAndFileFromArgs(final String[] args, final String configFileOption, final String[] expectedArgs) {

        final String[] newArgs;

        if ( configFileOption == null ) {
            newArgs = ConfigUtils.removeConfigOptionAndFileFromArgs(args);
        }
        else {
            newArgs = ConfigUtils.removeConfigOptionAndFileFromArgs(args, configFileOption);
        }

        Assert.assertEquals( newArgs, expectedArgs );
    }
    @Test(expectedExceptions = UserException.BadInput.class)
    void testRemoveConfigOptionAndFileFromArgsWithException() {
        ConfigUtils.removeConfigOptionAndFileFromArgs(new String[] {"main","END",StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION});
    }
}
