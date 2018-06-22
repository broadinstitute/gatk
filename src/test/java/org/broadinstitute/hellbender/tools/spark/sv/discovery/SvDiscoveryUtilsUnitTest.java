package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.ValidationStringency;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class SvDiscoveryUtilsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void test() {

        ReadsSparkSource readSource = new ReadsSparkSource(SparkContextFactory.getTestSparkContext(), ValidationStringency.LENIENT);
        final String referenceFile = b38_reference_20_21;

        final String inputBam = toolsTestDir + "/spark/sv/utils/testSAMWriter_chr20.bam";

        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, referenceFile);

        final List<GATKRead> temp = rddParallelReads.collect();
        final byte[] bases = temp.get(0).getBases();
        Assert.assertNotNull(bases);
        SAMFileHeader header = readSource.getHeader(inputBam, referenceFile);

        final File tempDirNew = BaseTest.createTempDir("new");
        tempDirNew.deleteOnExit();

        // test output empty set
        String outputPath = tempDirNew.getAbsolutePath() + "/out.bam";
        SvDiscoveryUtils.writeSAMRecords(rddParallelReads, Collections.emptySet(), outputPath, header);
        Assert.assertTrue(readSource.getParallelReads(outputPath, referenceFile).collect().isEmpty());

        // test output BAM
        final HashSet<String> names = new HashSet<>(Arrays.asList("asm000004:tig00016",
                "asm028919:tig00012",
                "asm016781:tig00011",
                "asm028306:tig00006",
                "asm028311:tig00040",
                "asm000004:tig00007",
                "asm000004:tig00032",
                "asm016781:tig00011",
                "asm028317:tig00003",
                "asm000004:tig00022"));
        SvDiscoveryUtils.writeSAMRecords(rddParallelReads, names, outputPath, header);
        Assert.assertEquals(readSource.getParallelReads(outputPath, referenceFile).collect().size(), 16);

        // test output SAM
        outputPath = tempDirNew.getAbsolutePath() + "/out.sam";
        SvDiscoveryUtils.writeSAMRecords(rddParallelReads, names, outputPath, header);
        Assert.assertEquals(readSource.getParallelReads(outputPath, referenceFile).collect().size(), 16);
    }
}
