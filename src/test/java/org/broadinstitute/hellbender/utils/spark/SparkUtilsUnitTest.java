package org.broadinstitute.hellbender.utils.spark;

import htsjdk.samtools.*;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SparkUtilsUnitTest extends BaseTest {

    @Test
    public void testConvertHeaderlessHadoopBamShardToBam() throws Exception {
        final File bamShard = new File(publicTestDir + "org/broadinstitute/hellbender/utils/spark/reads_data_source_test1.bam.headerless.part-r-00000");
        final File output = createTempFile("testConvertHadoopBamShardToBam", ".bam");
        final File headerSource = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final int expectedReadCount = 11;

        boolean shardIsNotValidBam = false;
        try ( final ReadsDataSource readsSource = new ReadsDataSource(bamShard.toPath()) ) {
            for ( final GATKRead read : readsSource ) {}
        }
        catch ( SAMFormatException e ) {
            shardIsNotValidBam = true;
        }

        Assert.assertTrue(shardIsNotValidBam, "Input shard should not be a valid BAM");

        SAMFileHeader header = null;
        try ( final SamReader headerReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(headerSource) ) {
            header = headerReader.getFileHeader();
        }
        catch ( IOException e ) {
            throw new UserException("Error reading header from " + headerSource.getAbsolutePath(), e);
        }

        SparkUtils.convertHeaderlessHadoopBamShardToBam(bamShard, header, output);

        int actualCount = 0;
        try ( final ReadsDataSource readsSource = new ReadsDataSource(output.toPath()) ) {
            for ( final GATKRead read : readsSource ) { ++actualCount; }
        }

        Assert.assertEquals(actualCount, expectedReadCount, "Wrong number of reads in final BAM file");
    }

    @Test
    public void testPathExists() throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            //use the HDFS on the mini cluster
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path tempPath = new Path(workingDirectory, "testFileExists.txt");
            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

            Assert.assertFalse(SparkUtils.pathExists(ctx, tempPath));
            final FileSystem fs = tempPath.getFileSystem(ctx.hadoopConfiguration());
            final FSDataOutputStream fsOutStream = fs.create(tempPath);
            fsOutStream.close();
            fs.deleteOnExit(tempPath);
            Assert.assertTrue(SparkUtils.pathExists(ctx, tempPath));
        });

    }
}
