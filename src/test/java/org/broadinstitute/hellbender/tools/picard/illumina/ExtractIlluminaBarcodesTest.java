package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ClusterData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProvider;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProviderFactory;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.broadinstitute.hellbender.utils.text.parsers.BasicInputParser;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import static htsjdk.samtools.util.IOUtil.*;
import static java.util.Arrays.asList;
import static java.util.Arrays.sort;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.*;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

/**
 * @author alecw@broadinstitute.org
 */
public class ExtractIlluminaBarcodesTest extends CommandLineProgramTest {
    private static final File SINGLE_DATA_DIR = new File(getTestDataDir(), "picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_DATA_DIR = new File(getTestDataDir(), "picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final String[] BARCODES = {
            "CAACTCTC",
            "CAACTCTG", // This one is artificial -- one edit away from the first one
            "ACAGGTAT",
            "GACCGTTG",
            "ATTATCAA",
            "TGCTGCTG",
            "AACAATGG",
            "TGTAATCA",
            "GCCGTCGA",
            "GTCCACAG",
            "TTGTCTAT",
            "GTGGAGAC",
            "TTGCAAAT"
    };

    private File basecallsDir;
    private File dual;
    private File qual;

    public String getTestedClassName() {
        return ExtractIlluminaBarcodes.class.getSimpleName();
    }

    @BeforeTest
    private void setUp() throws Exception {
        basecallsDir = File.createTempFile("eib.", ".tmp");
        assertTrue(basecallsDir.delete());
        assertTrue(basecallsDir.mkdir());
        copyDirectoryTree(SINGLE_DATA_DIR, basecallsDir);
        dual = File.createTempFile("eib_dual", ".tmp");
        assertTrue(dual.delete());
        assertTrue(dual.mkdir());
        copyDirectoryTree(DUAL_DATA_DIR, dual);
        qual = File.createTempFile("eib_qual", ".tmp");
        assertTrue(qual.delete());
        assertTrue(qual.mkdir());
        copyDirectoryTree(DUAL_DATA_DIR, qual);
    }

    @AfterTest
    private void tearDown() {
        deleteDirectoryTree(basecallsDir);
        deleteDirectoryTree(dual);
        deleteDirectoryTree(qual);
    }

    @Test
    public void testSingleEndWithBarcodeAtStart() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(1, "8B25T");
        assertEquals(metricsFile.getMetrics().get(11).PERFECT_MATCHES, 1);
    }

    @Test
    public void testSingleEndWithBarcodeAtEnd() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(1, "25T8B");
        assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 5);
    }

    @Test
    public void testPairedEndWithBarcodeOnFirstEnd() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(1, "25T8B25T");
        assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 5);
    }

    @Test
    public void testPairedEndWithBarcodeOnSecondEnd() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(1, "25T25T8B");
        assertEquals(metricsFile.getMetrics().get(12).PERFECT_MATCHES, 1);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonWritableOutputFile() throws Exception {
        final File existingFile = new File(basecallsDir, "s_1_1101_barcode.txt.gz");
        try {
            existingFile.setReadOnly();
            final String readStructure = "25T8B25T";
            final int lane = 1;

            final File metricsFile = File.createTempFile("eib.", ".metrics");
            metricsFile.deleteOnExit();

            final List<String> args = new ArrayList<>(asList(
                    "--BASECALLS_DIR", basecallsDir.getPath(),
                    "--LANE", Integer.toString(lane),
                    "--READ_STRUCTURE", readStructure,
                    "--METRICS_FILE", metricsFile.getPath(),
                    "--COMPRESS_OUTPUTS", "true"
            ));
            for (final String barcode : BARCODES) {
                args.add("--BARCODE");
                args.add(barcode);
            }
            assertEquals(runCommandLine(args), Main.ANY_OTHER_EXCEPTION_EXIT_VALUE);
        } finally {
            existingFile.setWritable(true);
        }

    }

    /**
     * 4 cases tested:
     * * exact match to ACAGTG
     * * inexact match within threshold to TGACCA
     * * inexact match not within threshold to TGACCA
     * * inexact match where the next match is too close to ACAGTG
     *
     * @throws Exception
     */
    @Test
    public void testBarcodeMatching() throws Exception {
        final int lane = 1;
        final int barcodePosition = 26;
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(lane, "25T8B25T");

        ExtractIlluminaBarcodes.BarcodeMetric metricOne = null;
        ExtractIlluminaBarcodes.BarcodeMetric metricTwo = null;
        ExtractIlluminaBarcodes.BarcodeMetric metricNoMatch = null;
        for (final ExtractIlluminaBarcodes.BarcodeMetric metric : metricsFile.getMetrics()) {
            if (metric.BARCODE.equals(BARCODES[0])) {
                metricOne = metric;
            } else if (metric.BARCODE.equals(BARCODES[2])) {
                metricTwo = metric;
            } else if (metric.BARCODE.equals("NNNNNNNN")) {
                metricNoMatch = metric;
            }
        }
        assertEquals(metricOne.PERFECT_MATCHES, 5);
        assertEquals(metricOne.ONE_MISMATCH_MATCHES, 0);
        assertEquals(metricOne.PF_READS, 3);
        assertEquals(metricOne.READS, 5);

        // one inexact match
        assertEquals(metricTwo.READS, 4);
        assertEquals(metricTwo.ONE_MISMATCH_MATCHES, 0);

        assertEquals(metricNoMatch.READS, 140);
        assertEquals(metricNoMatch.PF_READS, 112);

        // Check the barcode files themselves
        final File[] barcodeFiles = getFilesMatchingRegexp(basecallsDir, "s_" + lane + "_\\d{4}_barcode.txt");
        sort(barcodeFiles);

        final BasicInputParser barcodeParser = new BasicInputParser(true, barcodeFiles);

        // Exact match
        String[] illuminaFields = barcodeParser.next();
        assertEquals(illuminaFields[1], "Y");
        assertEquals(illuminaFields[2], "CAACTCTC");

        // Inexact match
        illuminaFields = barcodeParser.next();
        assertEquals(illuminaFields[1], "Y");
        assertEquals(illuminaFields[2], "ACAGGTAT");

        // Too many mismatches
        illuminaFields = barcodeParser.next();
        assertEquals(illuminaFields[1], "N");

        barcodeParser.close();

        // Tack on test of barcode-informed Illumina Basecall parsing
        final ReadStructure rs = new ReadStructure("25T8B25T");
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDir, lane, rs,
                new BclQualityEvaluationStrategy(ILLUMINA_ALLEGED_MINIMUM_QUALITY),
                BaseCalls, QualityScores, Barcodes);
        testParsing(factory, rs, metricOne, barcodePosition);
    }

    @Test
    public void testDualBarcodes() throws Exception {
        final File metricsFile = File.createTempFile("dual.", ".metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "--BASECALLS_DIR", dual.getAbsolutePath(),
                "--LANE", "1",
                "--METRICS_FILE", metricsFile.getPath(),
                "--READ_STRUCTURE", "25T8B8B25T",
                "--BARCODE", "CAATAGTCCGACTCTC"
        };

        assertEquals(runCommandLine(args), 0);
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> result = new MetricsFile<>();
        result.read(new FileReader(metricsFile));
        assertEquals(result.getMetrics().get(0).PERFECT_MATCHES, 1, "Got wrong number of perfect matches");
        assertEquals(result.getMetrics().get(0).ONE_MISMATCH_MATCHES, 0, "Got wrong number of one-mismatch matches");
    }

    /**
     * Testing the quality thresholding. Looking at a single barcode (ACAGTG) with a min quality of 25 and no mismatches
     */
    @Test(dataProvider = "qualityBarcodeData")
    public void testQualityBarcodes(final int quality,
                                    final int maxMismatches, final int perfectMatches, final int oneMismatch,
                                    final String testName) throws Exception {
        final File metricsFile = File.createTempFile("qual.", ".metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "--BASECALLS_DIR", qual.getPath(),
                "--LANE", "1",
                "--READ_STRUCTURE", "25T8B25T",
                "--METRICS_FILE", metricsFile.getPath(),
                "--MINIMUM_BASE_QUALITY", Integer.toString(quality),
                "--MAX_MISMATCHES", Integer.toString(maxMismatches),
                "--BARCODE", "CAATAGTC"
        };

        assertEquals(runCommandLine(args), 0);
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> result = new MetricsFile<>();
        result.read(new FileReader(metricsFile));
        assertEquals(result.getMetrics().get(0).PERFECT_MATCHES, perfectMatches, "Got wrong number of perfect matches for test: '" + testName + "'");
        assertEquals(result.getMetrics().get(0).ONE_MISMATCH_MATCHES, oneMismatch, "Got wrong number of one-mismatch matches for test: '" + testName + "'");
    }

    @DataProvider(name = "qualityBarcodeData")
    public Object[][] getQualityTestData() {
        return new Object[][]{
                {16, 0, 1, 0, "Barcode has good quality, 1 match"},
                {25, 0, 0, 0, "Barcode has quality failures, no matches"}
        };
    }

    private void testParsing(final IlluminaDataProviderFactory factory, final ReadStructure readStructure, final ExtractIlluminaBarcodes.BarcodeMetric metricACAGTG, final int barcodePosition) {

        int numReads = 0;

        final IlluminaDataProvider dataProvider = factory.makeDataProvider();
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();

            if (metricACAGTG.BARCODE.equals(cluster.getMatchedBarcode())) {
                ++numReads;
            }

            assertEquals(cluster.getRead(readStructure.templates.getIndices()[0]).getQualities().length, barcodePosition - 1);
            assertEquals(cluster.getRead(readStructure.templates.getIndices()[0]).getBases().length, barcodePosition - 1);
        }
        assertEquals(numReads, metricACAGTG.READS);
        dataProvider.close();
    }

    private MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> runIt(final int lane, final String readStructure)
            throws Exception {
        final File metricsFile = File.createTempFile("eib.", ".metrics");
        metricsFile.deleteOnExit();

        final List<String> args = new ArrayList<>(asList(
                "--BASECALLS_DIR", basecallsDir.getPath(),
                "--LANE", Integer.toString(lane),
                "--READ_STRUCTURE", readStructure,
                "--METRICS_FILE", metricsFile.getPath()
        ));
        for (final String barcode : BARCODES) {
            args.add("--BARCODE");
            args.add(barcode);
        }
        return runIt(args, metricsFile);
    }

    private MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> runIt(final List<String> args, final File metricsFile) throws Exception {
        // Generate _barcode.txt files and metrics file.
        assertEquals(runCommandLine(args), 0);

        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> retval = new MetricsFile<>();
        retval.read(new FileReader(metricsFile));
        return retval;
    }
}
