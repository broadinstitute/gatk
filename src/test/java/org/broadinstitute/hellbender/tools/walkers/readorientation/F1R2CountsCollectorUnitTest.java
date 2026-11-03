package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.sun.management.UnixOperatingSystemMXBean;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.metrics.VersionHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class F1R2CountsCollectorUnitTest extends GATKBaseTest {
    // write a histogram file with the given headers and read it as LearnReadOrientationModel does
    private static MetricsFile<?, Integer> writeAndRead(final List<Header> headers) {
        final MetricsFile<?, Integer> metricsFile = new MetricsFile<>();
        headers.forEach(metricsFile::addHeader);
        final File file = createTempFile("histograms", F1R2CountsCollector.REF_HIST_EXTENSION);
        metricsFile.write(file);
        return LearnReadOrientationModel.readMetricsFile(file);
    }

    @Test
    public void testSampleHeader() {
        // LearnReadOrientationModel used to take a histogram file's first header as its sample
        final VersionHeader otherHeader = new VersionHeader();
        otherHeader.setVersionedItem("GATK");
        otherHeader.setVersionString("4.0");
        final MetricsFile<?, Integer> metricsFile = new MetricsFile<>();
        metricsFile.addHeader(otherHeader);
        F1R2CountsCollector.addSampleHeader(metricsFile, "tumor");

        final MetricsFile<?, Integer> read = writeAndRead(metricsFile.getHeaders());
        Assert.assertEquals(F1R2CountsCollector.getSample(read, new File("histograms")), "tumor");
    }

    @DataProvider
    public Object[][] headersWithoutOneSample() {
        return new Object[][] { { List.of() }, { List.of(new StringHeader("tumor"), new StringHeader("normal")) } };
    }

    @Test(dataProvider = "headersWithoutOneSample", expectedExceptions = UserException.BadInput.class)
    public void testHistogramFileWithoutOneSampleHeader(final List<Header> headers) {
        F1R2CountsCollector.getSample(writeAndRead(headers), new File("histograms"));
    }

    // A writer for each sample and 3-base context used to stay open throughout collection, so that 16 samples took
    // 1024 file descriptors for the alt tables alone.
    @Test
    public void testOpenFilesDuringCollectionDoNotScaleWithContexts() {
        final OperatingSystemMXBean os = ManagementFactory.getOperatingSystemMXBean();
        if (!(os instanceof UnixOperatingSystemMXBean)) {
            throw new SkipException("counting open files requires a Unix operating system");
        }
        final UnixOperatingSystemMXBean unixOs = (UnixOperatingSystemMXBean) os;

        final int numSamples = 20;
        final List<String> samples = IntStream.range(0, numSamples).mapToObj(n -> "sample" + n).collect(Collectors.toList());
        final File outputTarGz = createTempFile("f1r2", ".tar.gz");

        final long openFilesBefore = unixOs.getOpenFileDescriptorCount();
        final F1R2CountsCollector collector = new F1R2CountsCollector(new CollectF1R2CountsArgumentCollection(),
                new SAMFileHeader(), outputTarGz, samples);
        final long openFilesDuringCollection = unixOs.getOpenFileDescriptorCount() - openFilesBefore;
        collector.writeHistograms();
        collector.closeAndArchiveFiles();

        // one alt table per sample, not one per sample and context
        Assert.assertTrue(openFilesDuringCollection <= 2 * numSamples,
                String.format("%d files were open during collection for %d samples", openFilesDuringCollection, numSamples));

        // the archive holds an alt table for each sample and context, and two histograms for each sample
        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(outputTarGz.toPath(), extractedDir.toPath());
        final Set<String> expectedAltTables = samples.stream()
                .flatMap(sample -> F1R2FilterConstants.ALL_KMERS.stream().map(context -> F1R2CountsCollector.getAltTableName(sample, context)))
                .collect(Collectors.toSet());
        final Set<String> altTables = F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir).stream()
                .map(File::getName).collect(Collectors.toSet());
        Assert.assertEquals(altTables, expectedAltTables);
        Assert.assertEquals(extractedDir.listFiles().length, expectedAltTables.size() + 2 * numSamples,
                Arrays.toString(extractedDir.list()));
    }
}
