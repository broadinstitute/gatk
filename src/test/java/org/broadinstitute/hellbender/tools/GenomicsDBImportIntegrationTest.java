package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.SimpleIntervalUnitTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
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

  private static final String hg00096 = "src/test/resources/large/gvcfs/HG00096.g.vcf.gz";
  private static final String hg00268 = "src/test/resources/large/gvcfs/HG00268.g.vcf.gz";
  private static final String na19625 = "src/test/resources/large/gvcfs/NA19625.g.vcf.gz";

  @Override
  public String getTestedClassName() {
    return GenomicsDBImport.class.getSimpleName();
  }

  @Test
  public void testGenomicsDBImporter() throws IOException {

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args.add("-GW"); args.add(GENOMICSDB_WORKSPACE.getAbsoluteFile());

    SimpleInterval simpleInterval = new SimpleInterval("20", 17959479, 82403646);
    args.add("-L"); args.add(simpleInterval);

    args.add("-V"); args.add(hg00096); args.add(hg00268); args.add(na19625);

    runCommandLine(args);
  }
}