package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;


public class ExtractCohortEngineTest extends CommandLineProgramTest {
  private static final Logger logger = LogManager.getLogger(ExtractCohortEngineTest.class);

  private final String dirPath = getToolTestDataDir(); // TODO this goes to src/test/resources....
  // which makes me wonder if I should be putting them someplace else
  private static final String prefix =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/";
  private static final String cohortAvroFileName = prefix +"AnVIL_WGS_100_exported_cohort/000000000000";
  private static final String sampleFile = prefix +"AnVIL_WGS_100_exported_cohort/sample_list";
  private static final String outputFileName = prefix +"outputfile.vcf";
  private static final String expectedFileName = prefix +"AnVIL_WGS_100_exported_cohort/expected.vcf.gz";

  @Test
  public void testFinalVCFfromAvro() throws Exception {
    logger.info(dirPath); // TODO just to check this. but it looks like it goes to resources

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args
        .add("mode", "GENOMES")
        .add("ref-version", 38)
        .add("query-mode", "LOCAL_SORT")
        .add("R", hg38Reference)
        .add("O", outputFileName)
        .add("local-sort-max-records-in-ram", 10000000)
        .add("cohort-avro-file-name", cohortAvroFileName)
        .add("sample-file", sampleFile);


    runCommandLine(args);

    final File expectedVCF = createTempDir("expectedFileName");
    final File outputVCF = createTempDir("outputFileName");

    IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
  }
}
