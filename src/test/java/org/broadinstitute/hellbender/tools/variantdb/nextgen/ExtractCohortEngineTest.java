package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

public class ExtractCohortEngineTest extends CommandLineProgramTest {

  private static final String cohortAvroFileName = // is it okay to just grab one?
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/AnVIL_WGS_100_exported_cohort/000000000000";
  private static final String outputFile =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/outputfile.vcf";
  private static final String sampleFile =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/AnVIL_WGS_100_exported_cohort/sample_list";

  @Test
  public void testFinalVCFfromAvro() throws Exception {

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("query-mode", "LOCAL_SORT")
        .add("R", hg38Reference)
        .add("O", outputFile) // TODO questionable?
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile)
        .add("project-id", "blah blah");


    runCommandLine(args);

    final File newOutputVCF = createTempDir("");
    final File expectedOutputVCF = createTempDir("");

    IntegrationTestSpec.assertEqualTextFiles(newOutputVCF, expectedOutputVCF);


    // *manual_nocalls.vcf.gz modified manually to introduce a no-call (./.) GT at chr20 61417
    // no-call to test marking ./. as GQ 0 in PET and abstaining from putting variant into VET
    // original GT at chr20 61417 = 0/1
    final String input_vcf_file = "NA12878.haplotypeCalls.reblocked.chr20.100k.manual_nocall.vcf.gz";
    final String interval_list_file = "wgs_calling_regions.hg38.chr20.100k.interval_list";
    final String sample_map_file = "test_sample_map.tsv";
    final File outputDir = createTempDir("output_dir");
    final List<String> expectedOutputFiles = new ArrayList<>(Arrays.asList(
        getToolTestDataDir() + "expected.metadata_001_NA12878.tsv",
        getToolTestDataDir() + "expected.pet_001_NA12878.tsv",
        getToolTestDataDir() + "expected.vet_001_NA12878.tsv"));
;

  }
}
