package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;


public class ExtractCohortEngineTest extends CommandLineProgramTest {

  private static final String cohortAvroFileName = // is it okay to just grab one?
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/AnVIL_WGS_100_exported_cohort/000000000000";
  private static final String sampleFile =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/AnVIL_WGS_100_exported_cohort/sample_list";
  private static final String outputFileName =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/outputfile.vcf";
  private static final String expectedFileName =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/AnVIL_WGS_100_exported_cohort/expected.vcf.gz";

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
        .add("project-id", "blah blah");


    runCommandLine(args);

    final File expectedVCF = createTempDir("expectedFileName");
    final File outputVCF = createTempDir("outputFileName");

    IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
  }
}
