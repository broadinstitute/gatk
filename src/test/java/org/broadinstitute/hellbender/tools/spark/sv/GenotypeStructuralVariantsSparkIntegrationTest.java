package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by valentin on 4/21/17.
 */
public class GenotypeStructuralVariantsSparkIntegrationTest extends CommandLineProgramTest {

    private static final File INPUT_VCF = new File(largeFileTestDir, "sv_evidence_for_variants_input.vcf.gz");


    @Test
    public void test() {
        final File outputVcf = createTempFile("output", ".vcf.gz");
        outputVcf.delete();
        final File outputVcfIndex = new File(outputVcf.getPath() + ".tbi");
        outputVcfIndex.delete();
        outputVcfIndex.deleteOnExit();
        final List<String> args = new ArrayList<String>() {{
            add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
            add(INPUT_VCF.getAbsolutePath());
            add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
            add(outputVcf.getAbsolutePath());
        }};
        Assert.assertTrue(outputVcf.isFile());
        Assert.assertTrue(outputVcfIndex.isFile());
        runCommandLine(args);
    }

}

