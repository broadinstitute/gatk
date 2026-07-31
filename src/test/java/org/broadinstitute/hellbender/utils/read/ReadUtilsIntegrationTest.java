package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ReadUtilsIntegrationTest extends GATKBaseTest {


  @DataProvider(name="createSAMWriter")
  public Object[][] createSAMWriterData() {
    return new Object[][] {
        // Note: We expect to silently fail to create an index if createIndex is true but sort order is not coord.
        {getTestFile("print_reads.sam"),     getTestFile("print_reads.fasta"), ".cram", false,  true, true, false},
        {getTestFile("query_sorted.bam"),     null, ".bam", false,  true, true, false},
        {getTestFile("coordinate_sorted.bam"),null, ".bam", false,  true, true, true},
        {getTestFile("query_sorted.bam"),     null, ".bam", true,   true, true, false},
        {getTestFile("coordinate_sorted.bam"),null, ".bam", true,   true, true, true},
        {getTestFile("query_sorted.bam"),     null, ".bam", true,   true, false, false},
        {getTestFile("coordinate_sorted.bam"),null, ".bam", true,   true, false, true},
        {getTestFile("coordinate_sorted.bam"),null, ".bam", true,   false, false, false}        };
  }

  /**
   * Test writing a SAM file to GCS.
   *
   */
  @Test(dataProvider="createSAMWriter", groups = {"bucket"})
  public void testCreatePathSAMWriter(
      final File bamFile,
      final File referenceFile,
      final String outputExtension,
      final boolean preSorted,
      final boolean createIndex,
      final boolean createMD5,
      final boolean expectIndex) throws Exception {

    final String outputPathName = BucketUtils.getTempFilePath(getGCPTestStaging() + "samWriterTest", outputExtension);
    final Path outputPath = BucketUtils.getPathOnGcs(outputPathName);
    final Path md5Path = BucketUtils.getPathOnGcs(outputPathName + ".md5");
    final Path referencePath = referenceFile == null ? null : referenceFile.toPath();

    try (final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(bamFile)) {

      final SAMFileHeader header = samReader.getFileHeader();
      if (expectIndex) { // ensure test condition
        Assert.assertEquals(expectIndex, header.getSortOrder() == SAMFileHeader.SortOrder.coordinate);
      }

      try (final SAMFileWriter samWriter = ReadUtils.createCommonSAMWriter
          (outputPath, referencePath, samReader.getFileHeader(), preSorted, createIndex, createMD5)) {
        final Iterator<SAMRecord> samRecIt = samReader.iterator();
        while (samRecIt.hasNext()) {
          samWriter.addAlignment(samRecIt.next());
        }
      }

      Assert.assertEquals(expectIndex, null != SamFiles.findIndex(outputPath));
      Assert.assertEquals(createMD5, Files.exists(md5Path));
    }

    // Write to a local temp file to use as the comparison baseline. For CRAM output, TLEN is
    // elided during encode and recomputed from alignment positions on decode (htsjdk 5.0.0+
    // matches htslib behavior), so the round-tripped TLEN may differ from the original SAM.
    // Comparing against a locally-written file in the same format ensures both sides undergo
    // identical encoding/decoding, producing consistent TLEN values.
    final File localTempFile = File.createTempFile("samWriterTest", outputExtension);
    localTempFile.deleteOnExit();
    try (final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(bamFile)) {
      try (final SAMFileWriter localWriter = ReadUtils.createCommonSAMWriter(
          localTempFile.toPath(), referencePath, samReader.getFileHeader(), preSorted, false, false)) {
        final Iterator<SAMRecord> samRecIt = samReader.iterator();
        while (samRecIt.hasNext()) {
          localWriter.addAlignment(samRecIt.next());
        }
      }
    }

    // now check the GCS output matches the locally-written file
    try (final SamReader localReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(localTempFile);
        final SamReader outputReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(outputPath)) {
      Assert.assertEquals(localReader.iterator(), outputReader.iterator());
    }
  }

}

