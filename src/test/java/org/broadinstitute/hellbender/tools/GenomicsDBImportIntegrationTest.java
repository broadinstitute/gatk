package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class GenomicsDBImportIntegrationTest extends CommandLineProgramTest {

  private static final File TEST_DATA_DIR = getTestDataDir();
  private static final String TEST_OUTPUT_DIRECTORY =
    publicTestDir + "org/broadinstitute/hellbender/tools/genomicsdb";
  private static final File GENOMICSDB_WORKSPACE =
    new File(TEST_OUTPUT_DIRECTORY + "/tiledb-ws");

  @Override
  public String getTestedClassName() {
    return GenomicsDBImport.class.getSimpleName();
  }

  @Test
  public void testGenomicsDBImporter() throws IOException {
  }
}