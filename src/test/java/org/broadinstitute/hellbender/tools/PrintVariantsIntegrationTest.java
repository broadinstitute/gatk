package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public final class PrintVariantsIntegrationTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Override
    public String getTestedClassName() {
        return PrintVariants.class.getSimpleName();
    }


    public String baseTestString(final String vcf) {
        return "--variant " + vcf + " --output %s";
    }

    @Test
    public void testFileToFile() throws Exception {
        final String f = getToolTestDataDir() + "print_variants.vcf";
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(f),
                Arrays.asList(f)
                );
        spec.executeTest("test good file", this);
    }

}