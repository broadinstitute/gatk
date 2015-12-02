package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;

/**
 * Tests CollectAlignmentSummaryStatistics
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public final class CollectAlignmentSummaryMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectAlignmentSummaryMetrics");

    public String getTestedClassName() {
        return CollectAlignmentSummaryMetrics.class.getSimpleName();
    }

    @DataProvider(name="metricsfiles")
    public Object[][] alignmentMetricsFiles() {
        return new Object[][] {
                {"summary_alignment_stats_test.sam", "summary_alignment_stats_test.fasta"},
                {"summary_alignment_stats_test.bam", "summary_alignment_stats_test.fasta"},
                {"summary_alignment_stats_test.cram", "summary_alignment_stats_test.fasta"}
        };
    }

    @Test(dataProvider="metricsfiles")
    public void test(final String inputFile, final String referenceFile) throws IOException {
        final File input = new File(TEST_DATA_DIR, inputFile);
        final File reference = new File(TEST_DATA_DIR, referenceFile);
        final File outfile = BaseTest.createTempFile("alignmentMetrics", ".txt");
        final String[] args = new String[] {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--reference", reference.getAbsolutePath()
        };
        runCommandLine(args);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));
        
        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 7);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 3);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 59);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 19.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 303);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*58D/303D*/0.191419);
                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                break;
            case SECOND_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 9);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 7);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 239);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 707);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*19D/707D*/0.026874);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 18);
                Assert.assertEquals(metrics.PF_READS, 16);
                Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 10);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 298);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 1010);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*77D/1010D*/0.076238);
                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testBisulfite() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_bisulfite_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile = BaseTest.createTempFile("alignmentMetrics", ".txt");
        final String[] args = new String[] {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--reference", reference.getAbsolutePath(),
                "--IS_BISULFITE_SEQUENCED", "true"
        };
        runCommandLine(args);

        final NumberFormat format =  NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                // 19 no-calls, one potentially methylated base, one mismatch at a potentially methylated base
                Assert.assertEquals(metrics.TOTAL_READS, 1);
                Assert.assertEquals(metrics.PF_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 20.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 20D / 100D);
                Assert.assertEquals(metrics.BAD_CYCLES, 20);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(20 / (double) 100));
                break;
            case SECOND_OF_PAIR:
                // Three no-calls, two potentially methylated bases
                Assert.assertEquals(metrics.TOTAL_READS, 1);
                Assert.assertEquals(metrics.PF_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*3D/99D*/0.030303);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(3 / (double) 99));
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 2);
                Assert.assertEquals(metrics.PF_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 202);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 11.5);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 202);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*23D/199D*/0.115578);
                Assert.assertEquals(metrics.BAD_CYCLES, 23);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(23 / (double) 199));
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }


    @Test
    public void testNoReference() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File outfile = BaseTest.createTempFile("alignmentMetrics", ".txt");
        final String[] args = new String[] {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath()
        };
        runCommandLine(args);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 7);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                break;
            case SECOND_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 9);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 18);
                Assert.assertEquals(metrics.PF_READS, 16);
                Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testZeroLengthReads() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test2.sam");
        final File outfile = BaseTest.createTempFile("alignmentMetrics", ".txt");
        final String[] args = new String[] {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath()
        };
        runCommandLine(args);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));
        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            // test that it doesn't blow up
        }
    }

    @Test
    public void testMultipleLevelsOfMetrics() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test_multiple.sam");
        final File outfile = BaseTest.createTempFile("alignmentMetrics", ".txt");
        final String[] args = new String[] {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--METRIC_ACCUMULATION_LEVEL", "ALL_READS",
                "--METRIC_ACCUMULATION_LEVEL", "SAMPLE",
                "--METRIC_ACCUMULATION_LEVEL", "LIBRARY",
                "--METRIC_ACCUMULATION_LEVEL", "READ_GROUP"
        };
        runCommandLine(args);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            if (metrics.SAMPLE == null) {
                switch (metrics.CATEGORY) {
                case FIRST_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 7);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                    Assert.assertEquals(metrics.BAD_CYCLES, 19);
                    break;
                case SECOND_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 9);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                    Assert.assertEquals(metrics.BAD_CYCLES, 3);
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 18);
                    Assert.assertEquals(metrics.PF_READS, 16);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                    Assert.assertEquals(metrics.BAD_CYCLES, 22);
                    break;
                case UNPAIRED:
                default:
                    Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                }
            }
            else if (metrics.SAMPLE.equals("Ma")) {
                // There's only one library and one read group for this sample so the metrics for
                // every level should be identical
                switch (metrics.CATEGORY) {
                    case FIRST_OF_PAIR:
                        Assert.assertEquals(metrics.TOTAL_READS, 5);
                        Assert.assertEquals(metrics.PF_READS, 3);
                        Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                        Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                        Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                        Assert.assertEquals(metrics.BAD_CYCLES, 24);
                        break;
                    case SECOND_OF_PAIR:
                        Assert.assertEquals(metrics.TOTAL_READS, 5);
                        Assert.assertEquals(metrics.PF_READS, 5);
                        Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                        Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                        Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                        Assert.assertEquals(metrics.BAD_CYCLES, 3);
                        break;
                    case PAIR:
                        Assert.assertEquals(metrics.TOTAL_READS, 10);
                        Assert.assertEquals(metrics.PF_READS, 8);
                        Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                        Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                        Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                        Assert.assertEquals(metrics.BAD_CYCLES, 27);
                        break;
                    case UNPAIRED:
                    default:
                        Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                }
            }
            else if (metrics.SAMPLE.equals("Pa")) {
                // Two libraries and three read groups for this sample
                if (metrics.LIBRARY == null) {
                    switch (metrics.CATEGORY) {
                        case FIRST_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 4);
                            Assert.assertEquals(metrics.PF_READS, 4);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                            Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                            Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                            Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 19);
                            break;
                        case SECOND_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 4);
                            Assert.assertEquals(metrics.PF_READS, 4);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                            Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                            Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                            Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 3);
                            break;
                        case PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 8);
                            Assert.assertEquals(metrics.PF_READS, 8);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                            Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                            Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                            Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 22);
                            break;
                        case UNPAIRED:
                        default:
                            Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                    }
                }
                else if (metrics.LIBRARY.equals("lib1")) {
                    // Only one read group in this library so library and RG metrics should be identical
                    switch (metrics.CATEGORY) {
                        case FIRST_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 2);
                            Assert.assertEquals(metrics.PF_READS, 2);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 19);
                            break;
                        case SECOND_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 2);
                            Assert.assertEquals(metrics.PF_READS, 2);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.BAD_CYCLES, 3);
                            break;
                        case PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 4);
                            Assert.assertEquals(metrics.PF_READS, 4);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.BAD_CYCLES, 22);
                            break;
                        case UNPAIRED:
                        default:
                            Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                    }

                }
                else if (metrics.LIBRARY.equals("lib2")) {
                    if (metrics.READ_GROUP == null) {
                        switch (metrics.CATEGORY) {
                            case FIRST_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                                break;
                            case SECOND_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                                break;
                            case PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 4);
                                Assert.assertEquals(metrics.PF_READS, 4);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                                break;
                            case UNPAIRED:
                            default:
                                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                        }
                    }
                    else if (metrics.READ_GROUP.equals("i")) {
                        switch (metrics.CATEGORY) {
                            case FIRST_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                                break;
                            case SECOND_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                                break;
                            case PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                                break;
                            case UNPAIRED:
                            default:
                                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                        }
                    }
                    else if (metrics.READ_GROUP.equals("i2")) {
                        switch (metrics.CATEGORY) {
                            case FIRST_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 27);
                                break;
                            case SECOND_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                                break;
                            case PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 30);
                                break;
                            case UNPAIRED:
                            default:
                                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                        }
                    }
                    else {
                        Assert.fail("Data does not contain this read group: " + metrics.READ_GROUP);
                    }

                }
                else {
                    Assert.fail("Data does not contain this library: " + metrics.LIBRARY);
                }
            }
            else {
                Assert.fail("Data does not contain this sample: " + metrics.SAMPLE);
            }
        }

    }

}
