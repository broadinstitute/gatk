package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ConvertHeaderlessHadoopBamShardToBamIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testConvertHeaderlessHadoopBamShardToBam() {
        final File bamShard = new File(publicTestDir + "org/broadinstitute/hellbender/utils/spark/reads_data_source_test1.bam.headerless.part-r-00000");
        final File output = createTempFile("testConvertHeaderlessHadoopBamShardToBam", ".bam");
        final File headerSource = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final int expectedReadCount = 11;

        List<String> args = Arrays.asList(
            "--" + ConvertHeaderlessHadoopBamShardToBam.BAM_SHARD_LONG_NAME, bamShard.getAbsolutePath(),
            "--" + ConvertHeaderlessHadoopBamShardToBam.BAM_WITH_HEADER_LONG_NAME, headerSource.getAbsolutePath(),
            "-" + ConvertHeaderlessHadoopBamShardToBam.OUTPUT_SHORT_NAME, output.getAbsolutePath()
        );
        runCommandLine(args);

        int actualCount = 0;
        try ( final ReadsDataSource readsSource = new ReadsDataSource(output.toPath()) ) {
            for ( final GATKRead read : readsSource ) { ++actualCount; }
        }

        Assert.assertEquals(actualCount, expectedReadCount, "Wrong number of reads in final BAM file");
    }
}
