package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqFilterSparkIntegrationTest extends CommandLineProgramTest {

    static final String imagePath = publicTestDir + "hg19mini.fasta.img";
    static final String libraryPath = publicTestDir + PathSeqFilterSpark.class.getPackage().getName().replace(".", "/") + "/hg19mini.hss";

    @DataProvider(name = "pathseqFilterTestData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"basic_input.bam",
                    "basic_output.paired.bam",
                    "basic_output.unpaired.bam",
                    "basic_filter.metrics",
                    false, false, false, false, false},
                {"basic_input.bam",
                    "dup_output.paired.bam",
                    "dup_output.unpaired.bam",
                    "dup_filter.metrics",
                    false, true, true, false, false},
                {"basic_input.bam",
                    "kmer_output.paired.bam",
                    "kmer_output.unpaired.bam",
                    "kmer_filter.metrics",
                    false, true, false, true, false},
                {"basic_input.bam",
                    "bwa_output.paired.bam",
                    "bwa_output.unpaired.bam",
                    "bwa_filter.metrics",
                    false, true, false, false, true},
                {"aligned_input.bam",
                    "aligned_output.paired.bam",
                    null,
                    "aligned_filter.metrics",
                    true, true, false, false, false},
                {"basic_input.bam",
                    "all_output.paired.bam",
                    "all_output.unpaired.bam",
                    "all_filter.metrics",
                    false, false, true, true, true},
        };
    }

    @Override
    public String getTestedClassName() {
        return PathSeqFilterSpark.class.getSimpleName();
    }

    @Test(dataProvider = "pathseqFilterTestData")
    private void testFilterTool(final String inputBamFilename, final String expectedPairedBamFilename,
                                final String expectedUnpairedBamFilename, final String expectedMetricsFilename,
                                final boolean isHostAligned, final boolean skipFilters,
                                final boolean filterDuplicates, final boolean useKmerFilter,
                                final boolean useBwaFilter) throws Exception {
        final File inputBamFile = getTestFile(inputBamFilename);
        final File expectedMetricsFile = getTestFile(expectedMetricsFilename);

        final File outputPairedBamFile = createTempFile("output_paired", ".bam");
        final File outputUnpairedBamFile = createTempFile("output_unpaired", ".bam");
        final File outputMetricsFile = createTempFile("metrics", ".txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addBooleanArgument(PSFilterArgumentCollection.IS_HOST_ALIGNED_LONG_NAME, isHostAligned);
        args.addBooleanArgument(PSFilterArgumentCollection.SKIP_FILTERS_LONG_NAME, skipFilters);
        args.addBooleanArgument(PSFilterArgumentCollection.FILTER_DUPLICATES_LONG_NAME, filterDuplicates);
        args.addInput(inputBamFile);
        args.addFileArgument(PathSeqFilterSpark.PAIRED_OUTPUT_LONG_NAME, outputPairedBamFile);
        args.addFileArgument(PathSeqFilterSpark.UNPAIRED_OUTPUT_LONG_NAME, outputUnpairedBamFile);
        args.addFileArgument(PSFilterArgumentCollection.FILTER_METRICS_FILE_LONG_NAME, outputMetricsFile);
        if (useKmerFilter) {
            args.addFileArgument(PSFilterArgumentCollection.KMER_FILE_PATH_LONG_NAME, new File(libraryPath));
        }
        if (useBwaFilter) {
            args.addFileArgument(PSFilterArgumentCollection.FILTER_BWA_IMAGE_LONG_NAME, new File(imagePath));
        }

        this.runCommandLine(args.getArgsArray());

        if (expectedPairedBamFilename != null) {
            SamAssertionUtils.assertSamsEqual(outputPairedBamFile, getTestFile(expectedPairedBamFilename), ValidationStringency.LENIENT, null);
        }
        if (expectedUnpairedBamFilename != null) {
            SamAssertionUtils.assertSamsEqual(outputUnpairedBamFile, getTestFile(expectedUnpairedBamFilename), ValidationStringency.LENIENT, null);
        }
        Assert.assertTrue(MetricsFile.areMetricsEqual(outputMetricsFile, expectedMetricsFile));
    }

}
