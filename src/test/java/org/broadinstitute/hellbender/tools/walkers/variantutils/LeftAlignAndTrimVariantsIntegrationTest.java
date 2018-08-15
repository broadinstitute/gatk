package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;

import static org.testng.Assert.*;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;


public class LeftAlignAndTrimVariantsIntegrationTest extends CommandLineProgramTest {
    final Path testDataDir=Paths.get(getToolTestDataDir());
    @DataProvider(name="LeftAlignDataProvider")
    public Object[][] LeftAlignTestData() {
        return new Object[][] {{testDataDir.resolve("test_left_align_hg38.vcf"),Paths.get(b38_reference_20_21),testDataDir.resolve("expected_left_align_hg38.vcf")}};
    }

    @Test(dataProvider = "LeftAlignDataProvider")
    public void testLeftAlignment(Path inputFile,Path ref,Path expectedOutputFile) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + ref.toString()
                + " -V " + inputFile
                + " -O %s"
                + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false"
                + " --suppress-reference-path ",
                Collections.singletonList(expectedOutputFile.toString())
        );
        spec.executeTest("testLeftAlignment--"+inputFile.toString(),this);
    }



}