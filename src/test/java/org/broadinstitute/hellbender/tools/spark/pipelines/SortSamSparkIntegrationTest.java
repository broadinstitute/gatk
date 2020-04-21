package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

public final class SortSamSparkIntegrationTest extends CommandLineProgramTest {

    public static final String COUNT_READS_SAM = "count_reads.sam";
    public static final String COORDINATE_SAM = "count_reads_sorted.sam";
    public static final String QUERY_NAME_BAM = "count_reads.bam";
    public static final String COORDINATE_BAM = "count_reads_sorted.bam";
    public static final String COORDINATE_CRAM = "count_reads_sorted.cram";
    public static final String QUERY_NAME_CRAM = "count_reads.cram";
    public static final String REF = "count_reads.fasta";
    public static final String CRAM = ".cram";
    public static final String BAM = ".bam";
    public static final String SAM = ".sam";

    @DataProvider(name="sortbams")
    public Object[][] sortBAMData() {
        return new Object[][] {
                {COUNT_READS_SAM, COORDINATE_SAM, null, SAM, SAMFileHeader.SortOrder.coordinate},
                {QUERY_NAME_BAM, COORDINATE_BAM, null, BAM, SAMFileHeader.SortOrder.coordinate},
                {QUERY_NAME_CRAM, COORDINATE_CRAM, REF, BAM, SAMFileHeader.SortOrder.coordinate},
                {QUERY_NAME_CRAM, COORDINATE_CRAM, REF, CRAM, SAMFileHeader.SortOrder.coordinate},
                {QUERY_NAME_BAM, COORDINATE_BAM, REF, CRAM, SAMFileHeader.SortOrder.coordinate},

                {COORDINATE_SAM, COUNT_READS_SAM, null, SAM, SAMFileHeader.SortOrder.queryname},
                {COORDINATE_BAM, QUERY_NAME_BAM, null, BAM, SAMFileHeader.SortOrder.queryname},
                {COORDINATE_CRAM, QUERY_NAME_CRAM, REF, BAM, SAMFileHeader.SortOrder.queryname},
                {COORDINATE_CRAM, QUERY_NAME_CRAM, REF, CRAM, SAMFileHeader.SortOrder.queryname},
                {COORDINATE_BAM, QUERY_NAME_BAM, REF, CRAM, SAMFileHeader.SortOrder.queryname},
        };
    }

    @Test(dataProvider="sortbams", groups="spark")
    public void testSortBAMs(
            final String inputFileName,
            final String expectedOutputFileName,
            final String referenceFileName,
            final String outputExtension,
            final SAMFileHeader.SortOrder sortOrder) throws Exception {
        final File inputFile =  getTestFile(inputFileName);
        final File expectedOutputFile =  getTestFile(expectedOutputFileName);
        final File actualOutputFile = createTempFile("sort_sam", outputExtension);
        File referenceFile = null == referenceFileName ? null : getTestFile(referenceFileName);

        final SamReaderFactory factory = SamReaderFactory.makeDefault();

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputFile);
        args.addOutput(actualOutputFile);
        if (null != referenceFile) {
            args.addReference(referenceFile);
            factory.referenceSequence(referenceFile);
        }
        args.add(StandardArgumentDefinitions.SORT_ORDER_LONG_NAME, sortOrder.name());

        this.runCommandLine(args);

        //test files are exactly equal
        SamAssertionUtils.assertSamsEqual(actualOutputFile, expectedOutputFile, ValidationStringency.DEFAULT_STRINGENCY, referenceFile);

        //test sorting matches htsjdk
        try(ReadsDataSource in = new ReadsPathDataSource(IOUtils.toGATKPath(actualOutputFile), factory )) {
            BaseTest.assertSorted(Utils.stream(in).map(read -> read.convertToSAMRecord(in.getHeader())).iterator(), sortOrder.getComparatorInstance());
        }
    }

    // This test is disabled until https://github.com/broadinstitute/gatk/issues/5881 is fixed
    @Test(enabled = false, dataProvider="sortbams", groups="spark")
    public void testSortBAMsSharded(
            final String inputFileName,
            final String unused,
            final String referenceFileName,
            final String outputExtension,
            final SAMFileHeader.SortOrder sortOrder) {
        final File inputFile = getTestFile(inputFileName);
        final File actualOutputFile = createTempFile("sort_sam", outputExtension);
        File referenceFile = null == referenceFileName ? null : getTestFile(referenceFileName);
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputFile);
        args.addOutput(actualOutputFile);
        if (null != referenceFile) {
            args.addReference(referenceFile);
        }
        args.add(StandardArgumentDefinitions.SORT_ORDER_LONG_NAME, sortOrder.name());
        args.add(GATKSparkTool.SHARDED_OUTPUT_LONG_NAME,true);
        args.add(GATKSparkTool.NUM_REDUCERS_LONG_NAME, "2");

        this.runCommandLine(args);

        final ReadsSparkSource source = new ReadsSparkSource(SparkContextFactory.getTestSparkContext());
        final JavaRDD<GATKRead> reads = source.getParallelReads(new GATKPath(actualOutputFile.getAbsolutePath()), referenceFile == null ? null : new GATKPath(referenceFile.getAbsolutePath()));

        final SAMFileHeader header = source.getHeader(new GATKPath(actualOutputFile.getAbsolutePath()),
                referenceFileName == null ? null : new GATKPath(referenceFile.getAbsolutePath()));

        final List<SAMRecord> reloadedReads = reads.collect().stream().map(read -> read.convertToSAMRecord(header)).collect(Collectors.toList());
        BaseTest.assertSorted(reloadedReads.iterator(), sortOrder.getComparatorInstance(),   reloadedReads.stream().map(SAMRecord::getSAMString).collect(Collectors.joining("\n")));
    }

    @DataProvider
    public Object[][] getInvalidSortOrders(){
        return new Object[][]{
                {SAMFileHeader.SortOrder.unknown},
                {SAMFileHeader.SortOrder.unsorted},
                {SAMFileHeader.SortOrder.duplicate}
        };
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class, dataProvider = "getInvalidSortOrders")
    public void testBadSortOrders(SAMFileHeader.SortOrder badOrder){
        final File unsortedBam = new File(getTestDataDir(), QUERY_NAME_BAM);
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(unsortedBam);
        args.addOutput(createTempFile("sort_bam_spark", BAM));
        args.add(StandardArgumentDefinitions.SORT_ORDER_LONG_NAME, badOrder.toString());

        this.runCommandLine(args);
    }
}
