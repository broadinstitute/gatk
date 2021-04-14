package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;


class ExtractCohortTest extends CommandLineProgramTest {
  private final String prefix = getToolTestDataDir();
  private final String cohortAvroFileName = prefix +"000000000000.avro";
  private final String sampleFile = prefix +"sample_list";
  private final String outputFileName = prefix +"outputfile.vcf";
  private final String expectedFileName = prefix +"expected.vcf.gz";

  @Test
  public void testFinalVCFfromAvro() throws Exception {

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("query-mode", "LOCAL_SORT")
        .add("R", hg38Reference)
        .add("O", outputFileName)
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile)
        .add("project-id", "noProjectId");

    runCommandLine(args);

    final File expectedVCF = getTestFile(expectedFileName);
    final File outputVCF = getTestFile(outputFileName);

    IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
  }
}
