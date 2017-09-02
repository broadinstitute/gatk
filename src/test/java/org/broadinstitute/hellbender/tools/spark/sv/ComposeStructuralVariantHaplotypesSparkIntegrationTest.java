package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by valentin on 8/21/17.
 */
public class ComposeStructuralVariantHaplotypesSparkIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_FILES_DIR = new File(CommandLineProgramTest.publicTestDir,
            ComposeStructuralVariantHaplotypesSpark.class.getPackage().getName().replace(".", File.separator) + File.separator + "integration").getAbsolutePath();
    private static final String INPUT_VCF = new File(TEST_FILES_DIR, "compose_genotype_contigs_input.vcf.gz").getAbsolutePath();
    private static final String INPUT_BAM = new File(TEST_FILES_DIR, "compose_genotype_contigs_input.bam").getAbsolutePath();
    private static final String EXPECTED_OUTPUT_BAM = new File(TEST_FILES_DIR, "compose_genotype_contigs_output.bam").getAbsolutePath();
    private static final String EXPECTED_ALIGNED_OUTPUT_BAM = new File(TEST_FILES_DIR, "compose_genotype_contigs_aligned_output.bam").getAbsolutePath();

    @Test
    public void testComposeSVContigs() throws IOException {
        File outputBamFile = null, alignedOutputBamFile = null;
        try {
            outputBamFile = createTempFile("test-output", ".bam");
            alignedOutputBamFile = createTempFile("test-output-aln", ".bam");
            final IntegrationTestSpec spec = new IntegrationTestSpec(String.format("-%s %s -%s %s -%s %s -%s %s -%s %s",
                    StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b38_reference_20_21,
                    StandardArgumentDefinitions.VARIANT_SHORT_NAME, INPUT_VCF,
                    ComposeStructuralVariantHaplotypesSpark.CONTIGS_FILE_SHORT_NAME, INPUT_BAM,
                    StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputBamFile,
                    ComposeStructuralVariantHaplotypesSpark.ALIGNED_OUTPUT_SHORT_NAME, alignedOutputBamFile), Collections.emptyList());
            spec.executeTest("testComposeSVContigs", this);
            assertBamsAreEqual(outputBamFile, new File(EXPECTED_OUTPUT_BAM));
            assertBamsAreEqual(alignedOutputBamFile, new File(EXPECTED_ALIGNED_OUTPUT_BAM));
        } finally {
            try { if (alignedOutputBamFile != null) alignedOutputBamFile.delete(); } catch (final Throwable th) {};
            try { if (outputBamFile != null) outputBamFile.delete(); } catch (final Throwable th) {};
        }
    }

    private static void assertBamsAreEqual(final File actualFile, final File expectedFile) {
        final ReadsDataSource expected = new ReadsDataSource(expectedFile.toPath());
        final ReadsDataSource actual = new ReadsDataSource(actualFile.toPath());
        final SAMFileHeader expectedHeader = expected.getHeader();
        final SAMFileHeader actualHeader = actual.getHeader();
        final Iterator<GATKRead> expectedReads = expected.iterator();
        final Iterator<GATKRead> actualReads = actual.iterator();
        while (expectedReads.hasNext() && actualReads.hasNext()) {
            final SAMRecord expectedRead = expectedReads.next().convertToSAMRecord(expectedHeader);
            final SAMRecord actualRead = actualReads.next().convertToSAMRecord(actualHeader);
            assertReadsAreEqual(expectedRead, actualRead);
        }
        Assert.assertEquals(actualReads.hasNext(), expectedReads.hasNext());
        expected.close();
        actual.close();
    }

    private static void assertReadsAreEqual(final SAMRecord expectedRead, final SAMRecord actualRead) {
        Assert.assertEquals(expectedRead.getReadName(), actualRead.getReadName());
        Assert.assertEquals(expectedRead.getCigar(), actualRead.getCigar());
        Assert.assertEquals(expectedRead.getContig(), actualRead.getContig());
        Assert.assertEquals(expectedRead.getStart(), actualRead.getStart());
        Assert.assertEquals(expectedRead.getReadBases(), actualRead.getReadBases());
        Assert.assertEquals(expectedRead.getFlags(), actualRead.getFlags());
        Assert.assertEquals(expectedRead.getReferenceName(), actualRead.getReferenceName());
        Assert.assertEquals(expectedRead.getStart(), actualRead.getStart());
        Assert.assertEquals(
                composeReadAttributeSortedList(expectedRead),
                composeReadAttributeSortedList(actualRead));
    }

    private static List<Pair<String, Object>> composeReadAttributeSortedList(SAMRecord expectedRead) {
        return expectedRead.getAttributes().stream()
                .map(x -> new ImmutablePair<>(x.tag, x.value))
                .sorted(Comparator.comparing(Pair::getLeft))
                .collect(Collectors.toList());
    }
}
