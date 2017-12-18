package org.broadinstitute.hellbender.utils.config;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;

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

        @Argument(
                doc = "Output file with only the configuration settings in it.",
                fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
        )
        private File outputFilePath;

        @Override
        public void traverse() {
            final GATKConfig gatkConfig = ConfigFactory.getInstance().getGATKConfig();
            ConfigFactory.dumpConfigSettings( gatkConfig, outputFilePath.toPath() );
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

    @Test(dataProvider = "provideForTestToolWithConfigValidation")
    public void testToolWithConfigValidation(final Path configFilePath) {

        final File tmpFile = getSafeNonExistentFile("ConfigIntegrationTest.output.txt");
        tmpFile.deleteOnExit();

        // Create some arguments for our command:
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION);
        args.add(configFilePath.toUri().toString());
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(tmpFile);

        // Run our command:
        runCommandLine(args.getArgsArray());

        // Now we get to read in the file's contents and compare them to our config settings:
        final Properties properties = new Properties();
        try ( final FileInputStream inputStream = new FileInputStream(tmpFile) ) {
            properties.load(inputStream);
        }
        catch (final Exception ex) {
            throw new GATKException("Test error!", ex);
        }

        // Create a config object to compare with the properties:
        final GATKConfig config = ConfigFactory.getInstance().createConfigFromFile(configFilePath.toString(), GATKConfig.class);

        // Get the config properties:
        final LinkedHashMap<String, Object> configMap = ConfigFactory.getConfigMap(config, false);

        // Compare the properties and assert their equality:
        for ( final Map.Entry<String, Object> entry : configMap.entrySet() ) {
            Assert.assertEquals( properties.getProperty(entry.getKey()), String.valueOf(entry.getValue()) );
        }

        // Make sure they have the same size:
        Assert.assertEquals( properties.size(), configMap.size() );

        // =========================================================================================================
        // Now we have to reset our system options:

        // Run with normal options to reset any transient testing changes:
        final GATKConfig defaultConfig = ConfigFactory.getInstance().createConfigFromFile(null, GATKConfig.class);

        // Get our system properties and what they should be:
        final Map<String, Object> systemConfigOptions = ConfigFactory.getConfigMap(defaultConfig, true);

        // Reset our system properties:
        for ( final Map.Entry<String, Object> entry : systemConfigOptions.entrySet() ) {
            System.setProperty(entry.getKey(), entry.getValue().toString());
        }
    }
}
