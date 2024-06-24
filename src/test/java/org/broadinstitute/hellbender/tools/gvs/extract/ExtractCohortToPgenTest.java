package org.broadinstitute.hellbender.tools.gvs.extract;

import com.github.luben.zstd.ZstdInputStream;
import htsjdk.samtools.util.BufferedLineReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.pgen.PgenEmptyPgenException;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;


public class ExtractCohortToPgenTest extends CommandLineProgramTest {
  private final String prefix = getToolTestDataDir();
  private final String quickstart10mbRefRangesAvroFile = prefix + "quickstart_10mb_ref_ranges.avro";
  private final String quickstart10mbVetAvroFile = prefix + "quickstart_10mb_vet.avro";
  private final String quickstartSampleListFile = prefix + "quickstart.sample.list";


  private File decompressPvar(final File pvarFile) {
    try {
      // Open an input stream for reading from the compressed pvar
      final ZstdInputStream pvarZstdStream = new ZstdInputStream(new FileInputStream(pvarFile));
      final BufferedLineReader pvarZstdReader = new BufferedLineReader(pvarZstdStream);
      // Open an output stream to write the decompressed pvar to
      final File decompressedPvar = createTempFile(pvarFile.getName(), ".pvar");
      final BufferedWriter decompressedPvarWriter = new BufferedWriter(new FileWriter(decompressedPvar));

      // Read from the Zstd reader and write to the decompressed Pvar
      while (pvarZstdReader.ready()) {
        final String nextLine = pvarZstdReader.readLine();
        // Skip the source header line because some weirdness in the way the gatk build works makes that not match
        if(nextLine.startsWith("##source"))
          continue;

        decompressedPvarWriter.write(nextLine);
        decompressedPvarWriter.newLine();
      }

      pvarZstdReader.close();
      decompressedPvarWriter.close();

      return decompressedPvar;
    }
    catch(IOException e) {
      throw new RuntimeException("Failed to decompress pvar.zst file: ", e);
    }
  }

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
  public void testFinalVETSPgenfromRangesAvro() throws Exception {
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
    final File expectedPgen = getTestFile("ranges_extract.expected_vets.pgen");
    final File expectedPsam = getTestFile("ranges_extract.expected_vets.psam");
    final File expectedPvar = getTestFile("ranges_extract.expected_vets.pvar");

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
    // Decompress the pvar for validation
    final File decompressedPvar = decompressPvar(outputPvar);
    IntegrationTestSpec.assertEqualTextFiles(decompressedPvar, expectedPvar);
  }

  @Test
  public void testFinalVETSPgenfromRangesAvroLowAlleleCountMax() throws Exception {
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
    final File expectedPvar = getTestFile("ranges_extract.expected_low_AC.pvar");

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
    // Decompress the pvar for validation
    final File decompressedPvar = decompressPvar(outputPvar);
    IntegrationTestSpec.assertEqualTextFiles(decompressedPvar, expectedPvar);
  }

  @Test
  public void testFinalVETSPgenfromRangesAvroSeparateIndex() throws Exception {
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
    final File expectedPvar = getTestFile("ranges_extract.expected_separate_index.pvar");

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
    // Decompress the pvar for validation
    final File decompressedPvar = decompressPvar(outputPvar);
    IntegrationTestSpec.assertEqualTextFiles(decompressedPvar, expectedPvar);
  }

  @Test
  public void testFinalVETSPgenfromRangesAvroEmptyPgenAllow() throws Exception {
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
  public void testFinalVQSRPgenfromRangesAvro() throws Exception {
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
    final File expectedPgen = getTestFile("ranges_extract.expected_vqsr.pgen");
    final File expectedPsam = getTestFile("ranges_extract.expected_vqsr.psam");
    final File expectedPvar = getTestFile("ranges_extract.expected_vqsr.pvar");

    // create a temporary directory for the output
    final File outputDir = createTempDir("extract_output");
    final String outputBasePath = outputDir.getAbsolutePath() + "/extract_output";
    final File outputPgen = new File(outputBasePath + ".pgen");
    final File outputPsam = new File(outputBasePath + ".psam");
    final File outputPvar = new File(outputBasePath + ".pvar.zst");

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("use-vqsr-scoring", true)
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
    // Decompress the pvar for validation
    final File decompressedPvar = decompressPvar(outputPvar);
    IntegrationTestSpec.assertEqualTextFiles(decompressedPvar, expectedPvar);
  }

  @Test(expectedExceptions = PgenEmptyPgenException.class)
  public void testEmptyPgenExceptionVETS() throws Exception {
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
  public void testThrowErrorIfNoCallingFilteredGtsAndFilteringBySite() {
    // Verifies that an exception is thrown if you try to --convert-filtered-genotypes-to-no-calls, but are using site filtering
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
            .add("filter-set-site-table", "something")
            .add("filter-set-name", "something")
            .add("vqs-score-filter-by-site", true)
            .add("convert-filtered-genotypes-to-no-calls", true)
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testThrowFilterErrorVETS() throws Exception {
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
  public void testThrowFilterErrorVQSR() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("use-vqsr-scoring", true)
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
  public void testNoFilteringThresholdsErrorVETS() throws Exception {
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
            .add("vqs-score-filter-by-site", true)
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }
  @Test(expectedExceptions = UserException.class)
  public void testNoFilteringThresholdsErrorVQSR() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("use-vqsr-scoring", true)
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("emit-pls", false)
        .add("filter-set-info-table", "foo")
        .add("vqs-score-filter-by-site", true)
        .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testFakeFilteringErrorVETS() throws Exception {
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
            .add("vqs-score-filter-by-site", true)
            .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }

  @Test(expectedExceptions = UserException.class)
  public void testFakeFilteringErrorVQSR() throws Exception {
    final ArgumentsBuilder args = new ArgumentsBuilder();
    // No filterSetInfoTableName included, so should throw a user error with the performSiteSpecificVQSLODFiltering flag
    args
        .add("use-vqsr-scoring", true)
        .add("ref-version", 38)
        .add("R", hg38Reference)
        .add("O", "anything")
        .add("local-sort-max-records-in-ram", 10000000)
        .add("ref-ranges-avro-file-name", quickstart10mbRefRangesAvroFile)
        .add("vet-avro-file-name", quickstart10mbVetAvroFile)
        .add("sample-file", quickstartSampleListFile)
        .add("emit-pls", false)
        .add("filter-set-name", "foo")
        .add("vqs-score-filter-by-site", true)
        .add("pgen-chromosome-code", "chrM");
    runCommandLine(args);
  }
}
