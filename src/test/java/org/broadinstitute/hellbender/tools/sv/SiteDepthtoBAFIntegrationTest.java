package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class SiteDepthtoBAFIntegrationTest extends CommandLineProgramTest {
    public static final String sd2bafTestDir = toolsTestDir + "walkers/sv/SiteDepthToBAF/";

    @Override public String getTestedToolName() { return "SiteDepthtoBAF"; }

    @Test
    public void testStdDevTooBig() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(SiteDepthtoBAF.MIN_HET_PROBABILITY, "0.");
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, hg19_chr1_1M_dict);
        argsBuilder.add(SiteDepthtoBAF.BAF_SITES_VCF_LONG_NAME, sd2bafTestDir + "test1.vcf");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, sd2bafTestDir + "test1.sd.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(sd2bafTestDir + "test1.baf.txt"));
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test1: standard deviation too big", this);
    }

    @Test
    public void testAOK() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, publicTestDir + "hg19micro.dict");
        argsBuilder.add(SiteDepthtoBAF.BAF_SITES_VCF_LONG_NAME, sd2bafTestDir + "test2.vcf");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, sd2bafTestDir + "test2.sd.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(sd2bafTestDir + "test2.baf.txt"));
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test2: a-ok", this);
    }

    @Test
    public void testSitesMismatch() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, publicTestDir + "hg19micro.dict");
        argsBuilder.add(SiteDepthtoBAF.BAF_SITES_VCF_LONG_NAME, sd2bafTestDir + "test2a.vcf");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, sd2bafTestDir + "test2.sd.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(), 1, UserException.class);
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test2: a-ok", this);
    }

    @Test
    public void testAdjustMedian() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, hg19_chr1_1M_dict);
        argsBuilder.add(SiteDepthtoBAF.BAF_SITES_VCF_LONG_NAME, sd2bafTestDir + "test3.vcf");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, sd2bafTestDir + "test3.sd.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(sd2bafTestDir + "test3.baf.txt"));
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test3: adjust median", this);
    }
}
