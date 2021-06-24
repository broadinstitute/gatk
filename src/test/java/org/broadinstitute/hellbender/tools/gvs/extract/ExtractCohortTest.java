package org.broadinstitute.hellbender.tools.gvs.extract;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;


class ExtractCohortTest extends CommandLineProgramTest {
  private final String prefix = getToolTestDataDir();
  private final String cohortAvroFileName = prefix +"chr20_subset_3_samples.avro";
  private final String sampleFile = prefix +"sample_list";

  @Test
  public void testFinalVCFfromAvro() throws Exception {
    // To create the expectedVCF file (of expected output) --create a temp table in BQ with the following query
    // and then export it through the BQ GUI as an avro file into GCS.
    // SELECT * FROM `spec-ops-aou.anvil_100_for_testing.exported_cohort_all_samples`
    // where location < 20000000200000 and location >= 20000000100000
    // and (sample_name="HG00405" or sample_name="HG00418" or sample_name="HG00408")
    final File expectedVCF = getTestFile("expected_extract.vcf");

    // create a temporary file (that will get cleaned up after the test has run) to hold the output data in
    final File outputVCF = createTempFile("extract_output", "vcf");

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", outputVCF.getAbsolutePath())
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile)
        .add("emit-pls", false);

    runCommandLine(args);
    IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
  }

  @Test(expectedExceptions = UserException.class)
  public void testThrowFilterError() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile)
        .add("filter-set-info-table", "something")
        .add("filter-set-name", "something")
        .add("emit-pls", false);
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testNoFilteringThresholdsError() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile)
        .add("emit-pls", false)
        .add("filter-set-info-table", "foo")
        .add("vqslod-filter-by-site", true);
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testFakeFilteringError() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    // No filterSetInfoTableName included, so should throw a user error with the performSiteSpecificVQSLODFiltering flag
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile)
        .add("emit-pls", false)
        .add("filter-set-name", "foo")
        .add("vqslod-filter-by-site", true);
    runCommandLine(args);
  }
}
