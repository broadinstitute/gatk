package org.broadinstitute.hellbender.utils.config;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Integration test
 * Created by jonn on 12/7/17.
 */
public class ConfigIntegrationTest extends CommandLineProgramTest {

    @CommandLineProgramProperties(
            summary = "Dummy empty command line that requires a reference .",
            oneLineSummary = "empty class",
            programGroup = TestProgramGroup.class
    )
    public static class DummyGatkTool extends GATKTool {

        public static final String SYSTEM_VAR_ARG_LONG_NAME = "SystemVar";
        public static final String SYSTEM_OUT_FILE_ARG_LONG_NAME = "SysOutFile";

        @Argument(
                doc = "Output file with only the configuration settings in it.",
                fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
        )
        private File outputFilePath;

        @Argument(
                doc = "System properties to save to system properties output file.",
                fullName = SYSTEM_VAR_ARG_LONG_NAME
        )
        private List<String> systemPropertiesToSave;

        @Argument(
                doc = "Output file with only the specified system properties in it.",
                fullName = SYSTEM_OUT_FILE_ARG_LONG_NAME
        )
        private File systemPropertiesOutputFile;

        @Override
        public void traverse() {

            // Create our config and dump the config settings:
            final GATKConfig gatkConfig = ConfigFactory.getInstance().getGATKConfig();
            ConfigFactory.dumpConfigSettings( gatkConfig, outputFilePath.toPath() );

            final Properties systemProperties = new Properties();
            for (final String propName : systemPropertiesToSave) {
                systemProperties.put( propName, System.getProperty(propName) );
            }

            try ( final PrintStream systemOutPrintStream = new PrintStream(systemPropertiesOutputFile)) {
                systemProperties.store(systemOutPrintStream, "");
            }
            catch ( final IOException ex ) {
                throw new GATKException("Could not open output file for system properties.", ex);
            }
        }
    }

    private void listProperties(final Properties props) {

        final List<String> propKeys = props.keySet().stream()
                .map(Object::toString)
                .sorted()
                .collect(Collectors.toCollection(ArrayList::new));

        System.out.println("-- listing properties --");
        for ( final String key : propKeys ) {
            System.out.print(key);
            System.out.print("=");
            System.out.println(props.get(key));
        }
    }

    @Override
    public String getTestedClassName() {
        return DummyGatkTool.class.getSimpleName();
    }

    @DataProvider
    private Object[][] provideForTestToolWithConfigValidation() {

        return new Object[][] {
                // Default config file (because nothing can be loaded):
                { getSafeNonExistentPath("nonExistentConfig.config") },

                // Config file with overrides:
                { Paths.get(publicTestDir + "org/broadinstitute/hellbender/utils/config/" + "TestGATKConfigOverrides.properties") },
        };
    }

    @Test(dataProvider = "provideForTestToolWithConfigValidation",
            singleThreaded = true)
    public void testToolWithConfigValidation(final Path configFilePath) {

        // Grab the properties so we can reset them later:
        // NOTE: CLONE HERE IS ABSOLUTELY NECESSARY OR THIS WILL FAIL VERY VERY HARD!
        final Properties systemProperties = (Properties)System.getProperties().clone();
//        System.out.println("================================================================================");
//        listProperties(systemProperties);
//        System.out.println("================================================================================");

        try {
            final File tmpConfigPropsFile = getSafeNonExistentFile("ConfigIntegrationTest.config.properties");
            final File tmpSystemPropsFile = getSafeNonExistentFile("ConfigIntegrationTest.system.properties");

            // Add in some system options.
            // None of these should be masked by the options of the same name in the config file.
            final Map<String, String> systemPropertiesToQuery = new HashMap<>();
            systemPropertiesToQuery.put("gatk_stacktrace_on_user_exception", "true");
            systemPropertiesToQuery.put("samjdk.compression_level", "7777");

            for (final Map.Entry<String, String> entry : systemPropertiesToQuery.entrySet()) {
                System.setProperty(entry.getKey(), entry.getValue());
            }

            // Create some arguments for our command:
            final ArgumentsBuilder args = new ArgumentsBuilder();

            args.addArgument(StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION, configFilePath.toUri().toString());
            args.addArgument(DummyGatkTool.SYSTEM_OUT_FILE_ARG_LONG_NAME, tmpSystemPropsFile.toString());
            args.addOutput(tmpConfigPropsFile);

            // Add in our system properties to check:
            for ( final String sysProp : systemPropertiesToQuery.keySet() ) {
                args.addArgument(DummyGatkTool.SYSTEM_VAR_ARG_LONG_NAME, sysProp);
            }

            // Run our command:
            runCommandLine(args);

            // =========================================================================================================
            // Now we get to read in the file's contents and compare them to our config settings:
            final Properties configProperties = new Properties();
            try ( final FileInputStream inputStream = new FileInputStream(tmpConfigPropsFile) ) {
                configProperties.load(inputStream);
            }
            catch ( final Exception ex ) {
                throw new GATKException("Test error!", ex);
            }

            final Properties systemPropertiesPostToolRun = new Properties();
            try ( final FileInputStream inputStream = new FileInputStream(tmpSystemPropsFile) ) {
                systemPropertiesPostToolRun.load(inputStream);
            }
            catch ( final Exception ex ) {
                throw new GATKException("Test error!", ex);
            }

            // Create a config object to compare with the properties:
            final GATKConfig config = ConfigFactory.getInstance().createConfigFromFile(configFilePath.toString(), GATKConfig.class);

            // Get the config properties:
            final LinkedHashMap<String, Object> configMap = ConfigFactory.getConfigMap(config, false);

            // Compare the configuration properties and assert their equality:
            for ( final Map.Entry<String, Object> entry : configMap.entrySet() ) {
                Assert.assertEquals(configProperties.getProperty(entry.getKey()), String.valueOf(entry.getValue()));
            }

            // Compare the system properties and assert their equality:
            for ( final String sysPropKey : systemPropertiesToQuery.keySet() ) {
                Assert.assertEquals(systemPropertiesPostToolRun.getProperty(sysPropKey), systemPropertiesToQuery.get(sysPropKey));
            }

            // Make sure they have the same size:
            Assert.assertEquals(configProperties.size(), configMap.size());

        }
        // =========================================================================================================
        // Now we have to reset our system options:
        finally {
            // Clear our system properties:
            final List<Object> keyList = new ArrayList<>(System.getProperties().keySet());
            for ( final Object key : keyList ) {
                System.clearProperty(key.toString());
            }

            // Set back the old System Properties:
            System.setProperties(systemProperties);
        }
    }
}
