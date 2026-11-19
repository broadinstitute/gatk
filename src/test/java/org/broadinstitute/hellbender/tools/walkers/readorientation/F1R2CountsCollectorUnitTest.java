package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.sun.management.UnixOperatingSystemMXBean;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class F1R2CountsCollectorUnitTest extends GATKBaseTest {
    @Test
    public void testMedianIsAtLeast() {
        // agrees with MathUtils.median for odd and even numbers of values, including the mean of middle values on either
        // side of the threshold, and ties at the threshold
        final Random random = new Random(1);
        for (int trial = 0; trial < 10000; trial++) {
            final List<Integer> values = IntStream.range(0, 1 + random.nextInt(8)).mapToObj(n -> 45 + random.nextInt(10))
                    .collect(Collectors.toList());
            for (int threshold = 44; threshold <= 55; threshold++) {
                final int t = threshold;
                final int numAtLeast = (int) values.stream().filter(v -> v >= t).count();
                final int maxBelow = values.stream().filter(v -> v < t).max(Integer::compare).orElse(Integer.MIN_VALUE);
                final int minAtLeast = values.stream().filter(v -> v >= t).min(Integer::compare).orElse(Integer.MAX_VALUE);
                Assert.assertEquals(F1R2CountsCollector.medianIsAtLeast(values.size(), numAtLeast, maxBelow, minAtLeast, t),
                        MathUtils.median(values) >= t, values + " with threshold " + t);
            }
        }
    }

    @Test
    public void testDepthBins() {
        // exact depths up to 24, then rounded to the nearest multiple of 10, at most the maximum depth
        final int maxDepth = 200;
        IntStream.rangeClosed(1, 24).forEach(depth -> Assert.assertEquals(F1R2FilterUtils.binDepth(depth, maxDepth), depth));
        Assert.assertEquals(F1R2FilterUtils.binDepth(25, maxDepth), 30);
        Assert.assertEquals(F1R2FilterUtils.binDepth(34, maxDepth), 30);
        Assert.assertEquals(F1R2FilterUtils.binDepth(35, maxDepth), 40);
        Assert.assertEquals(F1R2FilterUtils.binDepth(194, maxDepth), 190);
        Assert.assertEquals(F1R2FilterUtils.binDepth(195, maxDepth), 200);
        Assert.assertEquals(F1R2FilterUtils.binDepth(1000, maxDepth), 200);
        Assert.assertEquals(F1R2FilterUtils.binDepth(204, 205), 200);
        Assert.assertEquals(F1R2FilterUtils.binDepth(205, 205), 205);
        Assert.assertEquals(F1R2FilterUtils.binDepth(15, 12), 12);

        final int[] bins = F1R2FilterUtils.getDepthBins(maxDepth);
        Assert.assertEquals(bins.length, 24 + 18);
        Assert.assertEquals(bins[23], 24);
        Assert.assertEquals(bins[24], 30);
        Assert.assertEquals(bins[bins.length - 1], maxDepth);
    }

    // Each block's histogram tables have a row for each context, and alt base and orientation, with sites, and the counts
    // are reset for the next block
    @Test
    public void testBlockSiteCounts() throws IOException {
        final int maxDepth = 200;
        final BlockSiteCounts counts = new BlockSiteCounts(maxDepth);
        final Path refTable = createTempPath("counts", F1R2CountsCollector.REF_HIST_EXTENSION);
        final Path altTable = createTempPath("counts", F1R2CountsCollector.ALT_HIST_EXTENSION);
        try (TableWriter<F1R2HistogramTable.Row> refWriter = F1R2HistogramTable.writer(refTable, "tumor", false, counts.getBinDepths());
             TableWriter<F1R2HistogramTable.Row> altWriter = F1R2HistogramTable.writer(altTable, "tumor", true, counts.getBinDepths())) {
            counts.addRefSite("ACG", 10);
            counts.addRefSite("ACG", 33);   // in the bin of depth 30
            counts.addRefSite("ACG", 500);  // in the bin of the maximum depth
            counts.addOneAltReadSite("ACG", Nucleotide.T, ReadOrientation.F1R2, 27);
            counts.writeRowsAndReset(refWriter, altWriter, "chr1", 1);
            counts.addRefSite("TTT", 10);
            counts.writeRowsAndReset(refWriter, altWriter, "chr1", 1_000_001);
        }

        final List<F1R2HistogramTable.Row> refRows = readRows(refTable);
        Assert.assertEquals(refRows.stream().map(row -> row.contig() + ":" + row.blockStart() + " " + row.histogramLabel()).collect(Collectors.toList()),
                List.of("chr1:1 ACG", "chr1:1000001 TTT"));
        Assert.assertEquals(refRows.get(0).depths(), F1R2FilterUtils.getDepthBins(maxDepth));
        Assert.assertEquals(countsByDepth(refRows.get(0)), Map.of(10, 1, 30, 1, 200, 1));
        Assert.assertEquals(countsByDepth(refRows.get(1)), Map.of(10, 1));

        final List<F1R2HistogramTable.Row> altRows = readRows(altTable);
        Assert.assertEquals(altRows.size(), 1);
        Assert.assertEquals(altRows.get(0).histogramLabel(), "ACG_T_F1R2");
        Assert.assertEquals(countsByDepth(altRows.get(0)), Map.of(30, 1));
    }

    private static List<F1R2HistogramTable.Row> readRows(final Path table) throws IOException {
        try (TableReader<F1R2HistogramTable.Row> reader = F1R2HistogramTable.reader(table)) {
            Assert.assertEquals(reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG), "tumor");
            return reader.toList();
        }
    }

    private static Map<Integer, Integer> countsByDepth(final F1R2HistogramTable.Row row) {
        return IntStream.range(0, row.depths().length).filter(n -> row.counts()[n] > 0).boxed()
                .collect(Collectors.toMap(n -> row.depths()[n], n -> row.counts()[n]));
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

        // a ref and an alt histogram table and an alt table per sample, not an alt table per sample and context
        Assert.assertTrue(openFilesDuringCollection <= 3 * numSamples,
                String.format("%d files were open during collection for %d samples", openFilesDuringCollection, numSamples));

        // the archive holds an alt table for each sample and context, and two histogram tables for each sample
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
