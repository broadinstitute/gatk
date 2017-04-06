package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by valentin on 4/7/17.
 */
public class EvidenceSamForVariantsSparkIntegrationTest extends CommandLineProgramTest {

    private static final File SMALL_TEST_FILE_ROOT = new File(publicTestDir, "org/broadinstitute/hellbender/tools/spark/sv/EvidenceSamForVariantsSpark");
    private static final File INTERVALS_FILE = new File(SMALL_TEST_FILE_ROOT, "test_intervals.list");
    private static final File INPUT_VARIANT_FILE = new File(largeFileTestDir, "sv_evidence_sam_for_variants_input.vcf.gz");
    private static final File INPUT_BAM_FILE = new File(largeFileTestDir, "sv_evidence_for_variants_input.bam");
    private static final File FASTQ_DIR = new File(largeFileTestDir, "sv_evidence_for_variants_fastqs");

    @Test
    public void testExampleVariantWalker() throws IOException {
        final File readsIntervalsFile = createTempFile("read-intervals-", ".list");
        final File outputBamFile = createTempFile("output", ".bam");
        readsIntervalsFile.delete();
        outputBamFile.delete();
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                String.format(" -%s %s -%s %s -%s %s -%s %s -%s %s -%s %s -%s %s",
                        StandardArgumentDefinitions.VARIANT_SHORT_NAME, INPUT_VARIANT_FILE.getAbsolutePath(),
                        EvidenceSamForVariantsSpark.ASSEMBLY_DIRECTORY_SHORT_NAME, FASTQ_DIR.getAbsolutePath(),
                        StandardArgumentDefinitions.REFERENCE_SHORT_NAME, BaseTest.b37_reference_20_21,
                        IntervalArgumentCollection.INCLUDED_INTERVALS_SHORT_NAME, INTERVALS_FILE.getAbsolutePath(),
                        StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_BAM_FILE.getAbsolutePath(),
                        EvidenceSamForVariantsSpark.READ_INTERVALS_SHORT_NAME, readsIntervalsFile,
                        StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputBamFile), Collections.emptyList());
        testSpec.executeTest("testExampleVariantWalker", this);
        Assert.assertTrue(readsIntervalsFile.isFile());
        Assert.assertTrue(outputBamFile.isFile());
        assertEquals(outputBamFile, INPUT_BAM_FILE);
    }

    private void assertEquals(final File actual, final File expected) throws IOException {

        final SamReader actualReader = SamReaderFactory.make().open(actual);
        final SamReader expectedReader = SamReaderFactory.make().open(expected);
        final List<SAMRecord> actualRecords = new ArrayList<>();
        final List<SAMRecord> expectedRecords = new ArrayList<>();
        actualReader.forEach(actualRecords::add);
        expectedReader.forEach(expectedRecords::add);
        Assert.assertEquals(expectedRecords.size(), actualRecords.size());
        Assert.assertEquals(actualRecords, expectedRecords);
    }

    private File indexFile(final File bamFile) {
        return new File(bamFile.getParent(), bamFile.getName().replace("\\.bam$", ".bai"));
    }
}
