package org.broadinstitute.hellbender.tools;

import com.googlecode.protobuf.format.JsonFormat;
import com.intel.genomicsdb.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.GenomicsDBFeatureReader;
import com.intel.genomicsdb.GenomicsDBImportConfiguration;
import com.intel.genomicsdb.GenomicsDBImporter;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.SimpleIntervalUnitTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.util.parsing.json.JSONFormat;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.BiConsumer;

import static com.googlecode.protobuf.format.JsonFormat.printToString;
import static org.broadinstitute.hellbender.engine.GenomicsDBIntegrationTest.assertCondition;

public final class GenomicsDBImportIntegrationTest extends CommandLineProgramTest {

  private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/genomicsdb";
  private static final File GENOMICSDB_WORKSPACE =new File(TEST_OUTPUT_DIRECTORY + "/tiledb-ws");
  private static final String GENOMICSDB_ARRAYNAME = "gatk4-genomicsdb-test-0";
  private static final File TEMP_CALLSETMAP_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/callset.json");
  private static final File TEMP_VIDMAP_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/vidmap.json");
  private static final File TEST_LOADER_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/loader.json");
  private static final File TEST_QUERY_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/query.json");

  private static final String hg00096 = "src/test/resources/large/gvcfs/HG00096.g.vcf.gz";
  private static final String hg00268 = "src/test/resources/large/gvcfs/HG00268.g.vcf.gz";
  private static final String na19625 = "src/test/resources/large/gvcfs/NA19625.g.vcf.gz";
  private static final File TEST_REFERENCE_GENOME = getTestDataDir() + "/large/";

  @Override
  public String getTestedClassName() {
    return GenomicsDBImport.class.getSimpleName();
  }

  @Test
  public void testGenomicsDBImporter() throws IOException {

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args.add("-GW"); args.add(GENOMICSDB_WORKSPACE.getAbsoluteFile());

    SimpleInterval simpleInterval = new SimpleInterval("20", 69491, 69521);
    args.add("-L"); args.add(simpleInterval);
    args.add("-V"); args.add(hg00096); args.add(hg00268); args.add(na19625);

    runCommandLine(args);

    File loaderJSON = newImportConfigurationFile();
    File queryJSON = newExportConfiguration();

    GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> featureReader =
      new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(loaderJSON, queryJSON);

    final Iterable<VariantContext> actualVcs = new FeatureDataSource<>(output);
    final Iterable<VariantContext> expectedVcs = new FeatureDataSource<>(expected);

    assertCondition(actualVcs, expectedVcs, (a,e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a,e, Collections.emptyList()));

    IOUtils.deleteRecursivelyOnExit(GENOMICSDB_WORKSPACE);
  }

  File newImportConfigurationFile() {
    GenomicsDBImportConfiguration.ImportConfiguration.Builder builder =
      GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();

    GenomicsDBImportConfiguration.Partition.Builder pBuilder =
      GenomicsDBImportConfiguration.Partition.newBuilder();

    GenomicsDBImportConfiguration.Partition partition =
      pBuilder
        .setArray(GENOMICSDB_ARRAYNAME)
        .setWorkspace(GENOMICSDB_WORKSPACE.getAbsolutePath())
        .setBegin(0)
        .build();

    List<GenomicsDBImportConfiguration.Partition> partitionList = new ArrayList<>();
    partitionList.add(partition);

    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      builder
        .setProduceTiledbArray(true)
        .setCallsetMappingFile(TEMP_CALLSETMAP_JSON_FILE.getAbsolutePath())
        .setVidMappingFile(TEMP_VIDMAP_JSON_FILE.getAbsolutePath())
        .setSizePerColumnPartition(1000L)
        .setSegmentSize(1048576)
        .addAllColumnPartitions(partitionList)
        .build();

    return GenomicsDBImporter.printLoaderJSONFile(
      importConfiguration, TEST_LOADER_JSON_FILE.getAbsolutePath());
  }

  File newExportConfiguration() {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder eBuilder =
      GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();

    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration =
      eBuilder
        .setTiledbWorkspace(GENOMICSDB_WORKSPACE.getAbsolutePath())
        .setTiledbArrayName(GENOMICSDB_ARRAYNAME)
        .setReferenceGenome(TEST_REFERENCE_GENOME)
        .build();

    return printQueryJSONFile(exportConfiguration, TEST_QUERY_JSON_FILE.getAbsolutePath());
  }

  private static <T> void assertCondition(Iterable<T> actual, Iterable<T> expected, BiConsumer<T,T> assertion){
    final Iterator<T> iterActual = actual.iterator();
    final Iterator<T> iterExpected = expected.iterator();
    while(iterActual.hasNext() && iterExpected.hasNext()){
      assertion.accept(iterActual.next(), iterExpected.next());
    }
    if (iterActual.hasNext()){
      Assert.fail("actual is longer than expected with at least one additional element: " + iterActual.next());
    }
    if (iterExpected.hasNext()){
      Assert.fail("actual is shorter than expected, missing at least one element: " + iterExpected.next());
    }

  }

  File printQueryJSONFile(
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration,
    String filename) {
    String queryJSONString = JsonFormat.printToString(exportConfiguration);

    File tempQueryJSONFile = new File(filename);

    try( PrintWriter out = new PrintWriter(tempQueryJSONFile)  ){
      out.println(queryJSONString);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    return tempQueryJSONFile;
  }
}