package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.spark.pipelines.PrintReadsSpark;
import org.broadinstitute.hellbender.tools.spark.pipelines.PrintVariantsSpark;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.DataprocTestUtils;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Iterator;

/**
 * tests that actually run spark tools on dataproc
 */
public class DataprocIntegrationTest extends CommandLineProgramTest{
    private String clusterName;

    @BeforeClass(groups = "cloud")
    public void startCluster(){
        clusterName = DataprocTestUtils.getTestCluster();
    }

    @DataProvider
    public Object[][] getCloudPaths(){
        return new Object[][]{
                {"org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam"},
                {"large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam"},
        };
    }

    @Test(dataProvider = "getCloudPaths", groups = "cloud")
    public void printReadSparkOnDataproc(final String input) throws IOException {
        final String gcsInputPath = getGCPTestInputPath() + input;
        final String outputPath = BucketUtils.getTempFilePath(getGCPTestStaging(), ".bam");

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addArgument("input", gcsInputPath)
                .addArgument("output", outputPath)
                .addArgument("bamPartitionSize", String.valueOf(10*1024*1024));
        DataprocTestUtils.launchGatkTool(PrintReadsSpark.class.getSimpleName(), argBuilder.getArgsList(), clusterName);
        final File expected = copyLocally(gcsInputPath, "expected");
        final File actual = copyLocally(outputPath, "actual");
        IntegrationTestSpec.assertMatchingFiles(Collections.singletonList(actual), Collections.singletonList(expected.toString()), true, ValidationStringency.LENIENT);
        assertReadsAreInCoordinatishOrder(actual);
    }

    private static void assertReadsAreInCoordinatishOrder(final File bam) {
        try(final ReadsDataSource reads = new ReadsDataSource(bam.toPath())){
            final Iterator<GATKRead> iter = reads.iterator();
            Locatable previous = null;
            while(iter.hasNext()){
                final GATKRead current = iter.next();
                if ( previous != null && previous.contigsMatch(current)){
                    Assert.assertTrue(previous.getStart() <= current.getStart());
                }
                previous = current;
            }
        }
    }

    //disabled due to https://github.com/broadinstitute/gatk/issues/3840
    @Test(groups = "cloud", enabled=false)
    public void printVariantsOnDataproc() throws IOException {
        final String gcsInputPath = getGCPTestInputPath() + "large/gvcfs/gatk3.7_30_ga4f720357.24_sample.21.expected.vcf";
        final String outputPath = BucketUtils.getTempFilePath(getGCPTestStaging(), ".vcf");

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addArgument("V", gcsInputPath)
                .addArgument("output", outputPath);
        DataprocTestUtils.launchGatkTool(PrintVariantsSpark.class.getSimpleName(), argBuilder.getArgsList(), clusterName);
        final File expected = copyLocally(gcsInputPath, "expected");
        final File actual = copyLocally(outputPath, "actual");
        IntegrationTestSpec.assertMatchingFiles(Collections.singletonList(actual), Collections.singletonList(expected.toString()), false, ValidationStringency.LENIENT);
    }

    private static File copyLocally(final String gcsInputPath, final String name) throws IOException {
        final Path input = IOUtils.getPath(gcsInputPath);
        final Path output = Files.createTempFile(name, input.getFileName().toString());
        Files.delete(output);
        Files.copy(input, output);
        return output.toFile();
    }
}
