package org.broadinstitute.hellbender.tools.spark.transforms;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.JoinReadsWithVariants;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.dataflow.pipelines.BaseRecalibratorDataflow;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.spark.JoinStrategy;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

public class BaseRecalibratorSparkFnUnitTest extends BaseTest {

    private String getResourceDir() {
        return "src/test/resources/org/broadinstitute/hellbender/tools/BQSR/";
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        final String localResources = getResourceDir();
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        return new Object[][]{
                // local computation and files (except for the reference)
                {hg18Reference, HiSeqBam, dbSNPb37, getResourceDir() + "expected.NA12878.chr17_69k_70k.txt"},
//                {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "-knownSites " + moreSites, localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.2inputs.recal.txt")},
//                {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--indels_context_size 4",  localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.indels_context_size_4.recal.txt")},
//                {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--low_quality_tail 5",     localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.low_quality_tail_5.recal.txt")},
//                {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--quantizing_levels 6",    localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.quantizing_levels_6.recal.txt")},
//                {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--mismatches_context_size 4", localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.mismatches_context_size_4.recal.txt")},
//                //// //{new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},
        };
    }

    @Test(dataProvider = "BQSRTest")
    public void testApply(String reference, String bam, String vcf, String expected) throws Exception{
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        JavaRDD<GATKRead> initialReads = readSource.getParallelReads(bam, intervals);

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<Variant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(vcf);

        final ReferenceDataflowSource referenceDataflowSource = new ReferenceDataflowSource(null, reference, BaseRecalibratorDataflow.BQSR_REFERENCE_WINDOW_FUNCTION);

        final JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(initialReads.filter(r -> !r.isUnmapped()), referenceDataflowSource, bqsrKnownVariants, JoinStrategy.SHUFFLE);
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, readsHeader, referenceDataflowSource.getReferenceSequenceDictionary(null), new RecalibrationArgumentCollection());
        final File report = BaseTest.createTempFile("report", "txt");
        try( PrintStream out = new PrintStream(report)){
            bqsrReport.createGATKReport().print(out);
        }
        IntegrationTestSpec.assertEqualTextFiles(report, new File(expected));
    }

}