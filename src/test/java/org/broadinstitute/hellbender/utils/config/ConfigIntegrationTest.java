package org.broadinstitute.hellbender.utils.config;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

/**
 * Integration test suite for Configuration tests.
 * Created by jonn on 10/4/17.
 */
public class ConfigIntegrationTest extends CommandLineProgramTest  {

    private static final String configFilePath = "src/main/resources/org/broadinstitute/hellbender/utils/config/GATKConfig.properties";

    @Override
    public String getTestedToolName() {
        return PrintReads.class.getSimpleName();
    }

    @Test
    public void testPrintReadsWithConfigFile() throws Exception {

        final File tmpDir = createTempDir("testPrintReadsWithConfigFile");

        final String inputFile = publicTestDir + "NA12878.chr17_69k_70k.dictFix.bam";
        final String outputFile = tmpDir + File.separator + "out.bam";

        // Create some arguments for our command:
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION);
        args.add(configFilePath);
        args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(inputFile);
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(outputFile);

        // Run our command:
        runCommandLine(args.getArgsArray());

        // Ensure the files are the same:
        IntegrationTestSpec.assertMatchingFiles(
                Collections.singletonList(new File(inputFile)),
                Collections.singletonList(outputFile),
                true,
                ValidationStringency.LENIENT
        );
    }
}
