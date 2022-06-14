package org.broadinstitute.hellbender.tools.gvs.extract;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;


class ExtractCohortTest extends CommandLineProgramTest {
  private final String prefix = getToolTestDataDir();
  private final String quickstart10mbRefRangesAvroFile = prefix + "quickstart_10mb_ref_ranges.avro";
  private final String quickstart10mbVetAvroFile = prefix + "quickstart_10mb_vet.avro";
  private final String quickstartSampleListFile = prefix + "quickstart.sample.list";

  @Test
  public void testFinalVCFfromRangesAvro() throws Exception {
    // To generate the Avro input files, create a table for export using the GVS QuickStart Data
    //
    // CREATE OR REPLACE TABLE `spec-ops-aou.terra_test_1.ref_ranges_for_testing` AS
    // SELECT * FROM `spec-ops-aou.terra_test_1.ref_ranges_001`
    // WHERE location >= (20 * 1000000000000) + 10000000 - 1001 AND location <= (20 * 1000000000000) + 20000000;
    //
    // Then export in GUI w/ Avro + Snappy
    //
    // And the same for the VET data:
    // CREATE OR REPLACE TABLE `spec-ops-aou.terra_test_1.vet_for_testing` AS
    // SELECT * FROM `spec-ops-aou.terra_test_1.vet_001`
    // WHERE location >= (20 * 1000000000000) + 10000000 - 1001 AND location <= (20 * 1000000000000) + 20000000
    //
    final File expectedVCF = getTestFile("ranges_extract.expected.vcf");

    // create a temporary file (that will get cleaned up after the test has run) to hold the output data in
    final File outputVCF = createTempFile("extract_output", "vcf");

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", outputVCF.getAbsolutePath())
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("L", "chr20:10000000-20000000");

    runCommandLine(args);
    IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
  }

  @Test(expectedExceptions = UserException.class)
  public void testThrowFilterError() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("filter-set-info-table", "something")
        .add("filter-set-name", "something")
        .add("emit-pls", false);
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testNoFilteringThresholdsError() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
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
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("emit-pls", false)
        .add("filter-set-name", "foo")
        .add("vqslod-filter-by-site", true);
    runCommandLine(args);
  }
}
