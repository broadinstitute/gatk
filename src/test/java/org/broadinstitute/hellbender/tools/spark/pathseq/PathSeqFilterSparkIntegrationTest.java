package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqFilterSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PathSeqFilterSpark.class.getSimpleName();
    }

    private void testFilterTool(final ArgumentsBuilder args, final String inputFilename, final String expectedFilePaired, final String expectedFileUnpaired) throws Exception {
        final File inputFile = getTestFile(inputFilename);
        final File tmpDir = createTempDir("tmp");
        if (!tmpDir.delete()) {
            Assert.fail();
        }
        final File outputFile = new File(tmpDir.getAbsolutePath(), "output");
        final String outputBasePath = outputFile.getAbsolutePath();
        final File metricsFile = createTempFile("metrics", ".txt");
        if (!metricsFile.delete()) {
            Assert.fail();
        }
        args.addInput(inputFile);
        args.addArgument("output", outputBasePath);
        args.addFileArgument("metricsFile", metricsFile);

        this.runCommandLine(args.getArgsArray());

        if (expectedFilePaired != null) {
            SamAssertionUtils.assertSamsEqual(new File(outputBasePath + ".paired.bam"), getTestFile(expectedFilePaired), ValidationStringency.LENIENT, null);
        }
        if (expectedFileUnpaired != null) {
            SamAssertionUtils.assertSamsEqual(new File(outputBasePath + ".unpaired.bam"), getTestFile(expectedFileUnpaired), ValidationStringency.LENIENT, null);
        }
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testFilterBasic() throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addBooleanArgument("filterDuplicates", false);
        testFilterTool(args, "basic_input.bam", "basic_output.paired.bam", "basic_output.unpaired.bam");
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testFilterDuplicates() throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addBooleanArgument("filterDuplicates", true);
        args.addBooleanArgument("skipFilters", true);
        testFilterTool(args, "basic_input.bam", "dup_output.paired.bam", "dup_output.unpaired.bam");
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testFilterKmer() throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final String libraryPath = "src/test/resources/" + PathSeqBuildKmers.class.getPackage().getName().replace(".", "/") + "/hg19mini.hss";
        args.addArgument("kmerLibraryPath", libraryPath);
        args.addBooleanArgument("skipFilters", true);
        args.addBooleanArgument("filterDuplicates", false);
        testFilterTool(args, "basic_input.bam", "kmer_output.paired.bam", "kmer_output.unpaired.bam");
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testFilterBwa() throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final String imagePath = "src/test/resources/" + PathSeqBuildKmers.class.getPackage().getName().replace(".", "/") + "/hg19mini.fasta.bwa_image";
        args.addArgument("bwamemIndexImage", imagePath);
        args.addBooleanArgument("skipFilters", true);
        args.addBooleanArgument("filterDuplicates", false);
        testFilterTool(args, "basic_input.bam", "bwa_output.paired.bam", "bwa_output.unpaired.bam");
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testFilterPrealignedOnly() throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addBooleanArgument("isHostAligned", true);
        args.addBooleanArgument("skipFilters", true);
        args.addBooleanArgument("filterDuplicates", false);
        testFilterTool(args, "aligned_input.bam", "aligned_output.paired.bam", null);
    }

}
