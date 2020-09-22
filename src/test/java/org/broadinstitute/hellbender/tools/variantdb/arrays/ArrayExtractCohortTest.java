package org.broadinstitute.hellbender.tools.variantdb.arrays;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryLoadData;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class ArrayExtractCohortTest extends CommandLineProgramTest {

    final String resourceDir = getTestDataDir() + "/variantdb/arrays";

    private  void loadData() {
//        BigQueryLoadData.loadData("");
    }


    @Test
    public void testArrayIngest() throws IOException {
        final File hg19ref = new File(getTestDataDir(), "hg19mini.fasta");
        final File probeInfo = new File(resourceDir,"probe_info.csv");
        final File expectedVcf = new File(resourceDir, "extract.vcf");

        final String[] args = {
                "-O", "cohort.vcf",
                "--project-id", "broad-dsp-methods",
                "-R", hg19ref.getAbsolutePath(),
                "--cohort-extract-table", "",
                "--cohort-sample-table", "",
                "--probe-info-file", probeInfo.getAbsolutePath(),
        };

        runCommandLine(args);

        File cohortVcf = new File("cohort.vcf");
        Assert.assertTrue(FileUtils.contentEquals(expectedVcf, cohortVcf), "vcf extract files differ");

    }

}
