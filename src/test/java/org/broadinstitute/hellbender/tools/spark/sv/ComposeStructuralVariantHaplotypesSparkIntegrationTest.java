package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.BAMFileReader;
import htsjdk.samtools.BamFileIoUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.picard.sam.CompareSAMs;
import org.broadinstitute.hellbender.utils.gcs.BamBucketIoUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.seqdoop.hadoop_bam.SAMRecordReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;

/**
 * Created by valentin on 8/21/17.
 */
public class ComposeStructuralVariantHaplotypesSparkIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_FILES_DIR = new File(CommandLineProgramTest.publicTestDir,
            ComposeStructuralVariantHaplotypesSpark.class.getPackage().getName().replace(".", File.separator) + File.separator + "integration").getAbsolutePath();
    private static final String INPUT_VCF = new File(TEST_FILES_DIR, "compose_genotype_contigs_input.vcf.gz").getAbsolutePath();
    private static final String INPUT_BAM = new File(TEST_FILES_DIR, "compose_genotype_contigs_input.bam").getAbsolutePath();
    private static final String EXPECTED_OUTPUT_BAM = new File(TEST_FILES_DIR, "compose_genotype_contigs_output.bam").getAbsolutePath();

    @Test
    public void testComposeSVContigs() throws IOException {
        final File outputBamFile = File.createTempFile("test-output", ".bam");
        final File alignedOutputBamFile = File.createTempFile("test-output-aln", ".bam");
        //final File outputBamFile = createTempFile("test-output", ".bam");
        final IntegrationTestSpec spec = new IntegrationTestSpec(String.format("-%s %s -%s %s -%s %s -%s %s -%s %s",
                StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b38_reference_20_21,
                StandardArgumentDefinitions.VARIANT_SHORT_NAME, INPUT_VCF,
                ComposeStructuralVariantHaplotypesSpark.CONTIGS_FILE_SHORT_NAME, INPUT_BAM,
                StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputBamFile,
                ComposeStructuralVariantHaplotypesSpark.ALIGNED_OUTPUT_SHORT_NAME, alignedOutputBamFile), Collections.emptyList());
        spec.executeTest("testComposeSVContigs", this);
        final ReadsDataSource expected = new ReadsDataSource(new File(EXPECTED_OUTPUT_BAM).toPath());
        final ReadsDataSource actual = new ReadsDataSource(outputBamFile.toPath());
        final Iterator<GATKRead> expectedReads = expected.iterator();
        final Iterator<GATKRead> actualReads = actual.iterator();
        while (expectedReads.hasNext() && actualReads.hasNext()) {
            final GATKRead expectedRead = expectedReads.next();
            final GATKRead actualRead = actualReads.next();
            assertReadsAreEqual(expectedRead, actualRead);
        }
        Assert.assertEquals(actualReads.hasNext(), expectedReads.hasNext());
        expected.close();
        actual.close();
        //Assert.assertTrue(alignedOutputBamFile.delete());
        //Assert.assertTrue(outputBamFile.delete());
    }

    private void assertReadsAreEqual(final GATKRead expectedRead, final GATKRead actualRead) {
        Assert.assertEquals(expectedRead.getName(), actualRead.getName());
        Assert.assertEquals(expectedRead.getCigar(), actualRead.getCigar());
        Assert.assertEquals(expectedRead.getContig(), actualRead.getContig());
        Assert.assertEquals(expectedRead.getStart(), actualRead.getStart());
        Assert.assertEquals(expectedRead.getBases(), actualRead.getBases());
    }
}
