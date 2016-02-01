package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link FilterSamples}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class FilterSamplesIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return FilterSamples.class.getSimpleName();
    }

    @Test(dataProvider = "sampleExtremeMeanCoverageFilterData")
    public void testSampleExtremeMeanCoverage(final List<SampleCoverageStats> sampleStats,
                                              final double minMean, final double maxMean)
        throws IOException
    {
        final File samplesFile = createSampleStatsFile(sampleStats);
        final List<SampleCoverageStats> leftInSamples = sampleStats.stream()
                .filter(stats -> stats.mean >= minMean && stats.mean <= maxMean)
                .collect(Collectors.toList());
        final List<SampleCoverageStats> leftOutSamples = sampleStats.stream()
                .filter(stats -> stats.mean < minMean || stats.mean > maxMean)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, samplesFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterSamples.MINIMUM_MEAN_COVERAGE_SHORT_NAME, String.valueOf(minMean),
                "-" + FilterSamples.MAXIMUM_MEAN_COVERAGE_SHORT_NAME, String.valueOf(maxMean),
                "-" + FilterSamples.REJECTION_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInSamples(leftInSamples, outputFile);
        checkRejectedSamples(leftOutSamples, rejectionFile, FilterSamples.SampleFilter.ExtremeCoverageMean);
    }

    @Test(dataProvider = "sampleExtremeMeanCoverageFilterData")
    public void testSampleExtremeMeanCoverageWithoutRejectionFile(final List<SampleCoverageStats> sampleStats,
                                                                  final double minMean, final double maxMean)
            throws IOException
    {
        final File samplesFile = createSampleStatsFile(sampleStats);
        final List<SampleCoverageStats> leftInSamples = sampleStats.stream()
                .filter(stats -> stats.mean >= minMean && stats.mean <= maxMean)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        runCommandLine(new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, samplesFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterSamples.MINIMUM_MEAN_COVERAGE_SHORT_NAME, String.valueOf(minMean),
                "-" + FilterSamples.MAXIMUM_MEAN_COVERAGE_SHORT_NAME, String.valueOf(maxMean),
        });

        checkLeftInSamples(leftInSamples, outputFile);
    }

    @Test(dataProvider = "sampleExtremeCoverageVarianceFilterData")
    public void testSampleExtremeCoverageVariance(final List<SampleCoverageStats> sampleStats,
                                                  final double min, final double max)
            throws IOException
    {
        final File samplesFile = createSampleStatsFile(sampleStats);
        final List<SampleCoverageStats> leftInSamples = sampleStats.stream()
                .filter(stats -> stats.variance >= min && stats.variance <= max)
                .collect(Collectors.toList());
        final List<SampleCoverageStats> leftOutSamples = sampleStats.stream()
                .filter(stats -> stats.variance < min || stats.variance > max)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, samplesFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterSamples.MINIMUM_COVERAGE_VARIANCE_SHORT_NAME, String.valueOf(min),
                "-" + FilterSamples.MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME, String.valueOf(max),
                "-" + FilterSamples.REJECTION_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInSamples(leftInSamples, outputFile);
        checkRejectedSamples(leftOutSamples, rejectionFile, FilterSamples.SampleFilter.ExtremeCoverageVariance);
    }

    @Test(dataProvider = "sampleExtremeCoverageVarianceFilterData")
    public void testSampleExtremeCoverageVarianceWithoutRejectionFile(final List<SampleCoverageStats> sampleStats,
                                                                      final double min, final double max)
            throws IOException
    {
        final File samplesFile = createSampleStatsFile(sampleStats);
        final List<SampleCoverageStats> leftInSamples = sampleStats.stream()
                .filter(stats -> stats.variance >= min && stats.variance <= max)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        runCommandLine(new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, samplesFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterSamples.MINIMUM_COVERAGE_VARIANCE_SHORT_NAME, String.valueOf(min),
                "-" + FilterSamples.MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME, String.valueOf(max),
        });

        checkLeftInSamples(leftInSamples, outputFile);
    }

    private void checkLeftInSamples(final List<SampleCoverageStats> expectedSampleStats, final File outputFile) throws IOException {
        final TableReader<String> outputReader = new TableReader<String>(outputFile) {

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                Assert.assertTrue(columns.containsExactly(FilterSamples.OUTPUT_NAME_COLUMN));
            }

            @Override
            protected String createRecord(final DataLine dataLine) {
                return dataLine.get(0);
            }
        };
        final List<String> outputSamples = outputReader.stream().collect(Collectors.toList());
        final List<String> expectedSamples = expectedSampleStats.stream().map(stats -> stats.sample).collect(Collectors.toList());
        Assert.assertEquals(outputSamples, expectedSamples);
    }

    private void checkRejectedSamples(final List<SampleCoverageStats> expectedSampleStats, final File rejectionFile,
                                      final FilterSamples.SampleFilter filter) throws IOException {
        try (final TableReader<String[]> rejectionReader = new TableReader<String[]>(rejectionFile) {
            @Override
            protected String[] createRecord(DataLine dataLine) {
                return dataLine.toArray();
            }
        }) {
            Assert.assertTrue(rejectionReader.columns().containsAll(
                    FilterSamples.REJECTION_SAMPLE_NAME_COLUMN,
                    FilterSamples.REJECTION_FILTER_NAME_COLUMN,
                    FilterSamples.REJECTION_REASON_COLUMN
            ));
            for (final SampleCoverageStats stats : expectedSampleStats) {
                final String[] values = rejectionReader.readRecord();
                Assert.assertEquals(values[rejectionReader.columns().indexOf(FilterSamples.REJECTION_SAMPLE_NAME_COLUMN)],
                        stats.sample);
                Assert.assertEquals(values[rejectionReader.columns().indexOf(FilterSamples.REJECTION_FILTER_NAME_COLUMN)],
                        filter.toString());
            }
        }
    }

    private File createSampleStatsFile(final List<SampleCoverageStats> sampleStats) {
        final File result = createTempFile("filter-test-samples", ".tab");
        try (final SampleCoverageStatsWriter writer = new SampleCoverageStatsWriter(result)) {
            for (final SampleCoverageStats stats : sampleStats) {
                writer.writeRecord(stats);
            }
            return result;
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    @DataProvider(name = "sampleExtremeMeanCoverageFilterData")
    public Object[][] sampleExtremeMeanCoverageFilterData() {
        final Random rdn = new Random(1313);
        final List<SampleCoverageStats> stats = IntStream.range(1, 1001)
                .mapToObj(i -> new SampleCoverageStats("sample_" + (i + 1), rdn.nextDouble() * 200 , 100))
                .collect(Collectors.toList());

        return new Object[][] {
                { stats,  0, 10},
                { stats,  10, 50},
                { stats,  90, 110},
                { stats,  -1.0, Double.POSITIVE_INFINITY},
                { stats, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY},
        };
    }

    @DataProvider(name = "sampleExtremeCoverageVarianceFilterData")
    public Object[][] sampleExtremeCoverageVarianceFilterData() {
        final Random rdn = new Random(1313);
        final List<SampleCoverageStats> stats = IntStream.range(1, 1001)
                .mapToObj(i -> new SampleCoverageStats("sample_" + (i + 1), 100, rdn.nextDouble() * 200))
                .collect(Collectors.toList());

        return new Object[][] {
                { stats,  0, 10},
                { stats,  10, 50},
                { stats,  90, 110},
                { stats,  -1.0, Double.POSITIVE_INFINITY},
                { stats, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY},
        };
    }
}
