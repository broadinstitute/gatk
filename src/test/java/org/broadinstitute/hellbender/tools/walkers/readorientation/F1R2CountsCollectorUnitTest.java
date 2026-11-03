package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.sun.management.UnixOperatingSystemMXBean;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.SkipException;
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
