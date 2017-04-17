package org.broadinstitute.hellbender.tools;

import com.intel.genomicsdb.GenomicsDBFeatureReader;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

public final class GenomicsDBImportIntegrationTest extends CommandLineProgramTest {

  private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools";
  private static final File GENOMICSDB_WORKSPACE = new File(TEST_OUTPUT_DIRECTORY + "/tiledb-ws");
  private static final String GENOMICSDB_ARRAYNAME = "gatk4-genomicsdb-test-0";
  private static final File TEST_CALLSETMAP_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/callset.json");
  private static final File TEST_VIDMAP_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/vidmap.json");

  private static final String hg00096 = publicTestDir + "large/gvcfs/HG00096.g.vcf.gz";
  private static final String hg00268 = publicTestDir + "large/gvcfs/HG00268.g.vcf.gz";
  private static final String na19625 = publicTestDir + "large/gvcfs/NA19625.g.vcf.gz";
  private static final String combined = publicTestDir + "large/gvcfs/combined.gatk3.7.g.vcf.gz";

  private static final File TEST_REFERENCE_GENOME = new File(publicTestDir + "/large/Homo_sapiens_assembly38.20.21.fasta");

  @Override
  public String getTestedClassName() {
    return GenomicsDBImport.class.getSimpleName();
  }

  @Test
  public void testGenomicsDBImporter() throws IOException {

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args.add("-GW"); args.add(GENOMICSDB_WORKSPACE.getAbsolutePath());
    args.add("-GCWS"); args.add(true);
    args.add("-GA"); args.add(GENOMICSDB_ARRAYNAME);
    args.add("-R"); args.add(TEST_REFERENCE_GENOME);

    SimpleInterval simpleInterval = new SimpleInterval("chr20", 17960187, 17981446);
    args.add("-L"); args.add(simpleInterval);
    args.add("-V"); args.add(hg00096);
    args.add("-V"); args.add(hg00268);
    args.add("-V"); args.add(na19625);                  
    args.add("-GVID"); args.add(TEST_VIDMAP_JSON_FILE.getAbsolutePath());
    args.add("-GCS"); args.add(TEST_CALLSETMAP_JSON_FILE.getAbsolutePath());

    runCommandLine(args);

    GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
      new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(
        TEST_VIDMAP_JSON_FILE.getAbsolutePath(),
        TEST_CALLSETMAP_JSON_FILE.getAbsolutePath(),
        GENOMICSDB_WORKSPACE.getAbsolutePath(),
        GENOMICSDB_ARRAYNAME,
        TEST_REFERENCE_GENOME.getAbsolutePath(), null, new BCF2Codec());

    CloseableTribbleIterator<VariantContext> actualVcs = genomicsDBFeatureReader.query(simpleInterval.getContig(), simpleInterval.getStart(), simpleInterval.getEnd());

    AbstractFeatureReader<VariantContext, LineIterator> combinedVCFReader = AbstractFeatureReader.getFeatureReader(combined, new VCFCodec(), true);
    CloseableTribbleIterator<VariantContext> expectedVcs = combinedVCFReader.query(simpleInterval.getContig(), simpleInterval.getStart(), simpleInterval.getEnd());

    GenomicsDBTestUtils.assertCondition(actualVcs, expectedVcs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a,e, Collections.emptyList()));

    actualVcs.close();
    expectedVcs.close();
    IOUtils.deleteRecursivelyOnExit(GENOMICSDB_WORKSPACE);
    FileUtils.deleteQuietly(TEST_CALLSETMAP_JSON_FILE);
    FileUtils.deleteQuietly(TEST_VIDMAP_JSON_FILE);
  }
}
