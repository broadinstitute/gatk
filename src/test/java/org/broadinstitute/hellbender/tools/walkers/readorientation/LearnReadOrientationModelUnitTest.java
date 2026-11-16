package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class LearnReadOrientationModelUnitTest extends CommandLineProgramTest {
    @Test
    public void testCombineHistograms(){
        final File unscatteredCountsTarGz = createTempFile("unscattered", ".tar.gz");
        final File scatteredDir = createTempDir("scattered");

        // Step 0: Run CollectF1R2Counts without scatter-gather
        runCommandLine(Arrays.asList(
                "-R", b37_reference_20_21,
                "-I", LearnReadOrientationModelIntegrationTest.hapmapBamSnippet,
                "-L", LearnReadOrientationModelIntegrationTest.intervalList,
                "-O", unscatteredCountsTarGz.getAbsolutePath()),
                CollectF1R2Counts.class.getSimpleName());

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(unscatteredCountsTarGz.toPath(), extractedDir.toPath());

        // Step 1: SplitIntervals
        final File intervalDir = createTempDir("intervals");
        final int scatterCount = 5;
        runCommandLine(Arrays.asList(
            "-R", b37_reference_20_21,
            "-L", LearnReadOrientationModelIntegrationTest.intervalList,
            "-O", intervalDir.getAbsolutePath(),
            "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount)),
            SplitIntervals.class.getSimpleName());

        // Step 2: CollectF1R2Counts
        final File[] intervals = intervalDir.listFiles();
        final List<File> extractedDirs = IntStream.range(0, intervals.length).mapToObj(i -> createTempDir("extracted_" + i)).collect(Collectors.toList());
        for (int i = 0; i < intervals.length; i++){
            final File scatteredTarGz = new File(scatteredDir, "scatter_" + i + ".tar.gz");
            runCommandLine(Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-I", LearnReadOrientationModelIntegrationTest.hapmapBamSnippet,
                    "-L", intervals[i].getAbsolutePath(),
                    "-O", scatteredTarGz.getAbsolutePath()),
                    CollectF1R2Counts.class.getSimpleName());

            IOUtils.extractTarGz(scatteredTarGz.toPath(), extractedDirs.get(i).toPath());
        }

        final List<File> scatteredRefHistogramTables = extractedDirs.stream()
                .flatMap(dir -> F1R2CountsCollector.getRefHistogramsFromExtractedTar(dir).stream()).collect(Collectors.toList());

        final List<File> scatteredAltHistogramTables = extractedDirs.stream()
                .flatMap(dir -> F1R2CountsCollector.getAltHistogramsFromExtractedTar(dir).stream()).collect(Collectors.toList());

        final List<File> scatteredAltTableFiles = extractedDirs.stream()
                .flatMap(dir -> F1R2CountsCollector.getAltTablesFromExtractedTar(dir).stream()).collect(Collectors.toList());

        final List<Histogram<Integer>> ref = sumHistogramTables(scatteredRefHistogramTables, false);
        final List<Histogram<Integer>> alt = sumHistogramTables(scatteredAltHistogramTables, true);
        final List<AltSiteRecord> altSites = LearnReadOrientationModel.gatherAltSiteRecords(scatteredAltTableFiles);

        final List<File> altTablesUnscattered = F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir);
        final List<AltSiteRecord> altSitesTruth = LearnReadOrientationModel.gatherAltSiteRecords(altTablesUnscattered);

        final List<Histogram<Integer>> refTruth = sumHistogramTables(F1R2CountsCollector.getRefHistogramsFromExtractedTar(extractedDir), false);
        final List<Histogram<Integer>> altTruth = sumHistogramTables(F1R2CountsCollector.getAltHistogramsFromExtractedTar(extractedDir), true);

        for (Histogram<Integer> truth : refTruth){
            final Histogram<Integer> eval = ref.stream().filter(h -> h.getValueLabel().equals(truth.getValueLabel())).findAny().get();
            Assert.assertEquals(eval.getSumOfValues(), truth.getSumOfValues());
        }

        for (Histogram<Integer> truth : altTruth){
            final Histogram<Integer> eval = alt.stream().filter(h -> h.getValueLabel().equals(truth.getValueLabel())).findAny().get();
            Assert.assertEquals(eval.getSum(), truth.getSum());
            Assert.assertEquals(eval.getSumOfValues(), truth.getSumOfValues());
        }

        Assert.assertEquals(altSites.size(), altSitesTruth.size());

    }

    // the histograms of the one sample of the tables
    private static List<Histogram<Integer>> sumHistogramTables(final List<File> tables, final boolean alt) {
        final Map<String, List<Histogram<Integer>>> histogramsBySample =
                LearnReadOrientationModel.sumHistogramTables(tables, alt, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
        Assert.assertEquals(histogramsBySample.size(), 1);
        return histogramsBySample.values().iterator().next();
    }
}
