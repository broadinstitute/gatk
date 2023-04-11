package org.broadinstitute.hellbender.tools.gvs.extract;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.pgen.PgenEmptyPgenException;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;


public class ExtractCohortToPgenTest extends CommandLineProgramTest {
  private final String prefix = getToolTestDataDir();
  private final String quickstart10mbRefRangesAvroFile = prefix + "quickstart_10mb_ref_ranges.avro";
  private final String quickstart10mbVetAvroFile = prefix + "quickstart_10mb_vet.avro";
  private final String quickstartSampleListFile = prefix + "quickstart.sample.list";

  @AfterTest
  public void afterTest() {
    try {
      String [] filesToCleanUp = {"anything", "anything.idx"};
      for (String file : filesToCleanUp) {
        Files.deleteIfExists(Paths.get(file));
      }
    } catch (IOException e) {
      throw new RuntimeException("Failed cleaning up 'anything' outputs: ", e);
    }
  }

  @Test
  public void testFinalVQSRLitePgenfromRangesAvro() throws Exception {
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
    final File expectedPgen = getTestFile("ranges_extract.expected_vqsr_lite.pgen");
    final File expectedPsam = getTestFile("ranges_extract.expected_vqsr_lite.psam");
    final File expectedPvar = getTestFile("ranges_extract.expected_vqsr_lite.pvar.zst");

    // Create a temp dif for the output
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");
    final File outputPsam = new File(outputBasePath + ".psam");
    final File outputPvar = new File(outputBasePath + ".pvar.zst");


    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
            .add("ref-version", 38)
            .add("R", hg38Reference)
            .add("O", outputPgen)
            .add("local-sort-max-records-in-ram", 10000000)
            .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
            .add("vet-avro-file-name", quickstart10mbVetAvroFile)
            .add("sample-file", quickstartSampleListFile)
            .add("L", "chr20:10000000-20000000")
            .add("pgen-chromosome-code", "chrM");

    runCommandLine(args);
    Assert.assertEquals(Files.mismatch(outputPgen.toPath(), expectedPgen.toPath()), -1L);
    IntegrationTestSpec.assertEqualTextFiles(outputPsam, expectedPsam);
    Assert.assertEquals(Files.mismatch(outputPvar.toPath(), expectedPvar.toPath()), -1L);
  }

  @Test
  public void testFinalVQSRLitePgenfromRangesAvroLowAlleleCountMax() throws Exception {
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
    final File expectedPgen = getTestFile("ranges_extract.expected_low_AC.pgen");
    final File expectedPsam = getTestFile("ranges_extract.expected_low_AC.psam");
    final File expectedPvar = getTestFile("ranges_extract.expected_low_AC.pvar.zst");

    // Create a temp dif for the output
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");
    final File outputPsam = new File(outputBasePath + ".psam");
    final File outputPvar = new File(outputBasePath + ".pvar.zst");


    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
            .add("ref-version", 38)
            .add("R", hg38Reference)
            .add("O", outputPgen)
            .add("local-sort-max-records-in-ram", 10000000)
            .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
            .add("vet-avro-file-name", quickstart10mbVetAvroFile)
            .add("sample-file", quickstartSampleListFile)
            .add("L", "chr20:10000000-20000000")
            .add("pgen-chromosome-code", "chrM")
            .add("max-alt-alleles", 9);

    runCommandLine(args);
    Assert.assertEquals(Files.mismatch(outputPgen.toPath(), expectedPgen.toPath()), -1L);
    IntegrationTestSpec.assertEqualTextFiles(outputPsam, expectedPsam);
    Assert.assertEquals(Files.mismatch(outputPvar.toPath(), expectedPvar.toPath()), -1L);
  }

  @Test
  public void testFinalVQSRLitePgenfromRangesAvroSeparateIndex() throws Exception {
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
    final File expectedPgen = getTestFile("ranges_extract.expected_separate_index.pgen");
    final File expectedPgi = getTestFile("ranges_extract.expected_separate_index.pgen.pgi");
    final File expectedPsam = getTestFile("ranges_extract.expected_separate_index.psam");
    final File expectedPvar = getTestFile("ranges_extract.expected_separate_index.pvar.zst");

    // Create a temp dif for the output
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");
    final File outputPgi = new File(outputBasePath + ".pgen.pgi");
    final File outputPsam = new File(outputBasePath + ".psam");
    final File outputPvar = new File(outputBasePath + ".pvar.zst");


    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
            .add("ref-version", 38)
            .add("R", hg38Reference)
            .add("O", outputPgen)
            .add("local-sort-max-records-in-ram", 10000000)
            .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
            .add("vet-avro-file-name", quickstart10mbVetAvroFile)
            .add("sample-file", quickstartSampleListFile)
            .add("L", "chr20:10000000-20000000")
            .add("pgen-chromosome-code", "chrM")
            .add("write-mode", "WRITE_SEPARATE_INDEX");

    runCommandLine(args);
    Assert.assertEquals(Files.mismatch(outputPgen.toPath(), expectedPgen.toPath()), -1L);
    Assert.assertEquals(Files.mismatch(outputPgi.toPath(), expectedPgi.toPath()), -1L);
    IntegrationTestSpec.assertEqualTextFiles(outputPsam, expectedPsam);
    Assert.assertEquals(Files.mismatch(outputPvar.toPath(), expectedPvar.toPath()), -1L);
  }

  @Test
  public void testFinalVQSRLitePgenfromRangesAvroEmptyPgenAllow() throws Exception {
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
    // Create a temp dif for the output
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");
    final File outputPsam = new File(outputBasePath + ".psam");
    final File outputPvar = new File(outputBasePath + ".pvar.zst");


    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
            .add("ref-version", 38)
            .add("R", hg38Reference)
            .add("O", outputPgen)
            .add("local-sort-max-records-in-ram", 10000000)
            .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
            .add("vet-avro-file-name", quickstart10mbVetAvroFile)
            .add("sample-file", quickstartSampleListFile)
            .add("L", "chr20:10000000-20000000")
            .add("XL", "chr20:10000005-20000000")
            .add("pgen-chromosome-code", "chrM")
            .add("allow-empty-pgen", true);

    runCommandLine(args);
    Assert.assertTrue(outputPgen.exists());
    Assert.assertEquals(outputPgen.length(), 0L);
    Assert.assertTrue(outputPsam.exists());
    Assert.assertEquals(outputPsam.length(), 0L);
    Assert.assertTrue(outputPvar.exists());
    Assert.assertEquals(outputPvar.length(), 0L);
  }

  @Test
  public void testFinalVQSRClassicPgenfromRangesAvro() throws Exception {
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
    final File expectedPgen = getTestFile("ranges_extract.expected_vqsr_classic.pgen");
    final File expectedPsam = getTestFile("ranges_extract.expected_vqsr_classic.psam");
    final File expectedPvar = getTestFile("ranges_extract.expected_vqsr_classic.pvar.zst");

    // create a temporary directory for the output
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");
    final File outputPsam = new File(outputBasePath + ".psam");
    final File outputPvar = new File(outputBasePath + ".pvar.zst");

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("use-vqsr-classic-scoring", true)
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", outputPgen)
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("L", "chr20:10000000-20000000")
        .add("pgen-chromosome-code", "chrM");

    runCommandLine(args);
    Assert.assertEquals(Files.mismatch(outputPgen.toPath(), expectedPgen.toPath()), -1L);
    IntegrationTestSpec.assertEqualTextFiles(outputPsam, expectedPsam);
    Assert.assertEquals(Files.mismatch(outputPvar.toPath(), expectedPvar.toPath()), -1L);
  }

  @Test(expectedExceptions = PgenEmptyPgenException.class)
  public void testEmptyPgenExceptionVQSRLite() throws Exception {
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");

    final ArgumentsBuilder args = new ArgumentsBuilder();
    // No filterSetInfoTableName included, so should throw a user error with the performSiteSpecificVQSLODFiltering flag
    args
            .add("ref-version", 38)
            .add("R", hg38Reference)
            .add("O", outputPgen)
            .add("local-sort-max-records-in-ram", 10000000)
            .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
            .add("vet-avro-file-name", quickstart10mbVetAvroFile)
            .add("sample-file", quickstartSampleListFile)
            .add("L", "chr20:10000000-20000000")
            .add("XL", "chr20:10000005-20000000")
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testThrowFilterErrorVQSRLite() throws Exception {
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
            .add("emit-pls", false)
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }
  @Test(expectedExceptions = UserException.class)
  public void testThrowFilterErrorVQSRClassic() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("use-vqsr-classic-scoring", true)
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("filter-set-info-table", "something")
        .add("filter-set-name", "something")
        .add("emit-pls", false)
        .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testNoFilteringThresholdsErrorVQSRLite() throws Exception {
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
            .add("vqsr-score-filter-by-site", true)
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }
  @Test(expectedExceptions = UserException.class)
  public void testNoFilteringThresholdsErrorVQSRClassic() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("use-vqsr-classic-scoring", true)
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("emit-pls", false)
        .add("filter-set-info-table", "foo")
        .add("vqsr-score-filter-by-site", true)
        .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testFakeFilteringErrorVQSRLite() throws Exception {
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
            .add("vqsr-score-filter-by-site", true)
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testFakeFilteringErrorVQSRClassic() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    // No filterSetInfoTableName included, so should throw a user error with the performSiteSpecificVQSLODFiltering flag
    args
        .add("use-vqsr-classic-scoring", true)
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("emit-pls", false)
        .add("filter-set-name", "foo")
        .add("vqsr-score-filter-by-site", true)
        .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }
}
