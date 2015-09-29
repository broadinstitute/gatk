package org.broadinstitute.hellbender.tools.spark.pipelines;


import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReadsPipelineSparkIntegrationTest extends CommandLineProgramTest {
    private String getResourceDir(){
        return getTestDataDir() + "/" + "BQSR" + "/";
    }

    @DataProvider(name = "EndToEndTestData")
    public Object[][] getEndToEndTestData() {
        final String GRCh37Ref = RefAPISource.GRCH37_REF_ID;
        final String HiSeqBam =  getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        return new Object[][] {
                // This test case checks that the input and output bams are equal (which should be the case with stub tool implementations)
                { HiSeqBam, GRCh37Ref, dbSNPb37, new File(HiSeqBam) }
        };
    }
    
    @Test(dataProvider = "EndToEndTestData")
    void dummyTest( final String inputBam, final String reference, final String knownSites, final File expectedOutput ) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        File outFile = new File("output.bam");
        IOUtils.deleteRecursivelyOnExit(outFile);
        String output = "file:///" + outFile.getCanonicalPath();

        ReadsPipelineSpark readsPipelineSpark = new ReadsPipelineSpark();
        readsPipelineSpark.run(ctx, inputBam, new OptionalIntervalArgumentCollection(), Arrays.asList(knownSites),
                RefAPISource.URL_PREFIX + reference, output, false, getDataflowTestApiKey());
    }
}