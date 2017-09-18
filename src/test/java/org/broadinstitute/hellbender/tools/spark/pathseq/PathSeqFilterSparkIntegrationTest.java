package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqFilterSparkIntegrationTest extends CommandLineProgramTest {

    static final String imagePath = "src/test/resources/" + PathSeqFilterSpark.class.getPackage().getName().replace(".", "/") + "/hg19mini.fasta.img";
    static final String libraryPath = "src/test/resources/" + PathSeqFilterSpark.class.getPackage().getName().replace(".", "/") + "/hg19mini.hss";

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

        final File tmpDir = createTempDir("tmp_pathseqFilterTest");
        final File outputBamFilebase = new File(tmpDir.getAbsolutePath(), "output");
        final String outputBamBasePath = outputBamFilebase.getAbsolutePath();
        final File outputMetricsFile = createTempFile("metrics", ".txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addBooleanArgument("isHostAligned", isHostAligned);
        args.addBooleanArgument("skipFilters", skipFilters);
        args.addBooleanArgument("filterDuplicates", filterDuplicates);
        args.addInput(inputBamFile);
        args.addArgument("output", outputBamBasePath);
        args.addFileArgument("filterMetricsFile", outputMetricsFile);
        if (useKmerFilter) {
            args.addFileArgument("kmerLibraryPath", new File(libraryPath));
        }
        if (useBwaFilter) {
            args.addFileArgument("filterBwaImage", new File(imagePath));
        }

        this.runCommandLine(args.getArgsArray());

        if (expectedPairedBamFilename != null) {
            SamAssertionUtils.assertSamsEqual(new File(outputBamBasePath + ".paired.bam"), getTestFile(expectedPairedBamFilename), ValidationStringency.LENIENT, null);
        }
        if (expectedUnpairedBamFilename != null) {
            SamAssertionUtils.assertSamsEqual(new File(outputBamBasePath + ".unpaired.bam"), getTestFile(expectedUnpairedBamFilename), ValidationStringency.LENIENT, null);
        }
        Assert.assertTrue(MetricsFile.areMetricsEqual(outputMetricsFile, expectedMetricsFile));
    }

}
