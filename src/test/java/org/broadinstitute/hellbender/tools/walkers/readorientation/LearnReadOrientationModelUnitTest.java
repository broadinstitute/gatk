package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class LearnReadOrientationModelUnitTest extends CommandLineProgramTest {
    @Test
    public void testCombineHistograms(){
        final File refHistTruthFile = createTempFile("ref", "metrics");
        final File altHistTruthFile = createTempFile("alt", "metrics");
        final File altTableTruthFile = createTempFile("alt","tsv");

        final File refMetricsDir = createTempDir("rh");
        final File altMetricsDir = createTempDir("ah");
        final File altTableDir = createTempDir("at");

        // Step 0: Run CollectF1R2Counts without scatter-gather
        runCommandLine(Arrays.asList(
                "-R", b37_reference_20_21,
                "-I", ReadOrientationModelIntegrationTest.hapmapBamSnippet,
                "-L", ReadOrientationModelIntegrationTest.intervalList,
                "--" + CollectF1R2Counts.ALT_DATA_TABLE_LONG_NAME, altTableTruthFile.getAbsolutePath(),
                "--" + CollectF1R2Counts.REF_SITE_METRICS_LONG_NAME, refHistTruthFile.getAbsolutePath(),
                "--" + CollectF1R2Counts.ALT_DEPTH1_HISTOGRAM_LONG_NAME, altHistTruthFile.getAbsolutePath()),
                CollectF1R2Counts.class.getSimpleName());


        // Step 1: SplitIntervals
        final File intervalDir = createTempDir("intervals");
        final int scatterCount = 5;
        runCommandLine(Arrays.asList(
            "-R", b37_reference_20_21,
            "-L", ReadOrientationModelIntegrationTest.intervalList,
            "-O", intervalDir.getAbsolutePath(),
            "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount)),
            SplitIntervals.class.getSimpleName());

        // Step 2: CollectF1R2Counts
        final File[] intervals = intervalDir.listFiles();
        for (int i = 0; i < intervals.length; i++){
            runCommandLine(Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-I", ReadOrientationModelIntegrationTest.hapmapBamSnippet,
                    "-L", intervals[i].getAbsolutePath(),
                    "--" + CollectF1R2Counts.ALT_DATA_TABLE_LONG_NAME, altTableDir.getAbsolutePath() + "/" + i + ".tsv",
                    "--" + CollectF1R2Counts.REF_SITE_METRICS_LONG_NAME, refMetricsDir.getAbsolutePath() + "/" + i + ".metrics",
                    "--" + CollectF1R2Counts.ALT_DEPTH1_HISTOGRAM_LONG_NAME, altMetricsDir.getAbsolutePath() + "/" + i + ".metrics"),
                    CollectF1R2Counts.class.getSimpleName());
        }

        final List<Histogram<Integer>> ref = LearnReadOrientationModel.sumHistogramsFromFiles(Arrays.asList(refMetricsDir.listFiles()), true);
        final List<Histogram<Integer>> alt = LearnReadOrientationModel.sumHistogramsFromFiles(Arrays.asList(altMetricsDir.listFiles()), false);
        final List<AltSiteRecord> altSites = LearnReadOrientationModel.gatherAltSiteRecords(Arrays.asList(altTableDir.listFiles())).getRight();

        final List<Histogram<Integer>> refTruth = LearnReadOrientationModel.readMetricsFile(refHistTruthFile).getAllHistograms();
        final List<Histogram<Integer>> altTruth = LearnReadOrientationModel.readMetricsFile(altHistTruthFile).getAllHistograms();
        final List<AltSiteRecord> altSitesTruth = AltSiteRecord.readAltSiteRecords(altTableTruthFile.toPath()).getRight();


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
}
