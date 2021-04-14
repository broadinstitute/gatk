package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;


class ExtractCohortTest extends CommandLineProgramTest {
  private static final Logger logger = LogManager.getLogger(ExtractCohortTest.class);

  private final String dirPath = getToolTestDataDir(); // TODO this goes to src/test/resources....
  // which makes me wonder if I should be putting them someplace else
  private static final String prefix =
      "src/test/java/org/broadinstitute/hellbender/tools/variantdb/nextgen/AnVIL_WGS_100_exported_cohort/";
  private static final String cohortAvroFileName = prefix +"000000000000";
  private static final String sampleFile = prefix +"sample_list";
  private static final String outputFileName = prefix +"outputfile.vcf";
  private static final String expectedFileName = prefix +"expected.vcf.gz";

  @Test
  public void testFinalVCFfromAvro() throws Exception {
    logger.info("yooooo" + dirPath); // TODO just to check this. but it looks like it goes to resources
    logger.info("yooooo"); // TODO just to check this. but it looks like it goes to resources

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

    final File expectedVCF = createTempFile("expectedFileName", "vcf"); // TODO wait this ne3ds to be a file not a directoryh
    final File outputVCF = createTempFile("outputFileName", "vcf");

    IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
  }
}
