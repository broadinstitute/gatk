package org.broadinstitute.hellbender.tools.spark.transforms;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

public class ApplyBQSRSparkFnUnitTest extends BaseTest {

    final String resourceDir =  "src/test/resources/org/broadinstitute/hellbender/tools/BQSR/";
    final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    final String naBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        ApplyBQSRArgumentCollection qq1 = new ApplyBQSRArgumentCollection();
        qq1.quantizationLevels = 1;

        ApplyBQSRArgumentCollection qq6 = new ApplyBQSRArgumentCollection();
        qq6.quantizationLevels = 6;

        ApplyBQSRArgumentCollection diq = new ApplyBQSRArgumentCollection();
        diq.disableIndelQuals = true;

        tests.add(new Object[]{hiSeqBam, new ApplyBQSRArgumentCollection(), resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.bqsr.alternate_allaligned.bam"});
        tests.add(new Object[]{hiSeqBam, qq1, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.bqsr.qq-1.alternate_allaligned.bam"});
        tests.add(new Object[]{hiSeqBam, qq6, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.bqsr.qq6.alternate_allaligned.bam"});
        tests.add(new Object[]{hiSeqBam, diq, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.bqsr.DIQ.alternate_allaligned.bam"});

        // TODO: add test inputs with some unaligned reads

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest")
    public void testPR(String bam, ApplyBQSRArgumentCollection args, String expected) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        final JavaRDD<GATKRead> reads = readSource.getParallelReads(bam);
        final Broadcast<RecalibrationReport> recalReport = ctx.broadcast(new RecalibrationReport(new File(resourceDir + "HiSeq.20mb.1RG.table.gz")));
        final JavaRDD<GATKRead> applied = ApplyBQSRSparkFn.apply(reads, recalReport, readsHeader, args);
        final SortedSet<GATKRead> recalibratedReads = new TreeSet<>(new ReadCoordinateComparator(readsHeader));
        recalibratedReads.addAll(applied.collect());

        final ReadsDataSource readsSource = new ReadsDataSource(new File(expected));
        final SortedSet<GATKRead> expectedReads = new TreeSet<>(new ReadCoordinateComparator(readsHeader));
        Iterators.addAll(expectedReads, readsSource.iterator());

        Assert.assertEquals(recalibratedReads, expectedReads);
    }
}