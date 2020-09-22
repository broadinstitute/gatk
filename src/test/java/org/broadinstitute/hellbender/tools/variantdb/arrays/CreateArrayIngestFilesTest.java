package org.broadinstitute.hellbender.tools.variantdb.arrays;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import javax.validation.constraints.AssertTrue;
import java.io.File;
import java.io.IOException;

public class CreateArrayIngestFilesTest extends CommandLineProgramTest {

    final String resourceDir = getTestDataDir() + "/variantdb/arrays";

    @Test
    public void testArrayIngest() throws IOException {
        final File sampleMap = new File(resourceDir,"sampleMap.csv");
        final File arrayVCF = new File(resourceDir,"array.vcf");
        final File arrayMetrics = new File(resourceDir,"testSampleMetrics.tsv");
        final File probeInfo = new File(resourceDir,"probe_info.csv");
        final File expectedSampleTsv = new File(resourceDir, "sample_001_testSample.tsv");
        final File expectedDataTsv = new File(resourceDir, "raw_001_testSample.tsv");

        final String[] args = {
                "-V", arrayVCF.getAbsolutePath(),
                "--metrics-file", arrayMetrics.getAbsolutePath(),
                "-SNM", sampleMap.getAbsolutePath(),
                "--probe-info-file", probeInfo.getAbsolutePath(),
                "--ref-version", "37"
        };

        runCommandLine(args);

        File actualSampleTsv = new File("sample_001_testSample.tsv");
        File actualDataTsv = new File("raw_001_testSample.tsv");
        Assert.assertTrue(FileUtils.contentEquals(expectedDataTsv, actualDataTsv), "array data files differ");
        Assert.assertTrue(FileUtils.contentEquals(expectedSampleTsv, actualSampleTsv), "array sample files differ");

    }
}