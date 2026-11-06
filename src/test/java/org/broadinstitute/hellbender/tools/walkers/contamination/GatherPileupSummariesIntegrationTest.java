package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public class GatherPileupSummariesIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testGather() throws IOException {
        final Random rng = new Random();
        final String sampleName = "sample1";

        final List<Pair<String, Integer>> pairs = Arrays.asList(
                new ImmutablePair<>("21", 10_000_000),
                new ImmutablePair<>("1", 10_000_000),
                new ImmutablePair<>("14", 10_000_000),
                new ImmutablePair<>("MT", 10_000),
                new ImmutablePair<>("14", 20_000_000),
                new ImmutablePair<>("1", 20_000_000),
                new ImmutablePair<>("MT", 12_000),
                new ImmutablePair<>("14", 30_000_000),
                new ImmutablePair<>("X", 20_000_000));
        final List<String> contigsOrdered = Arrays.asList("1", "14", "21", "X", "MT");
        final int numEntriesPerFile = 100;

        final Path directory = Files.createTempDirectory("GPS-test");
        for (final Pair<String, Integer> pair : pairs) {
            final int depth = rng.nextInt(100) + 1;
            final int refCount = rng.nextInt(depth);
            final int altCount = depth - refCount;
            final double af = rng.nextDouble();

            final String contig = pair.getLeft();
            final int firstPosition = pair.getRight();

            int position = firstPosition;
            final List<PileupSummary> records = new ArrayList<>();
            for (int i = 0; i < numEntriesPerFile; i++){
                records.add(new PileupSummary(contig, position, refCount, altCount, 0, af));
                position += 1;
            }

            final Path file = Files.createTempFile(directory, "pileupSummary", ".tsv");
            PileupSummary.writeToFile(sampleName, records, file.toFile());
        }

        final File combinedPileupSummary = createTempFile("combined", "tsv");
        final File[] files = directory.toFile().listFiles();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG19_DICT);
        args.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, combinedPileupSummary.getAbsolutePath());
        Arrays.stream(files).forEach(f -> {
            args.add(StandardArgumentDefinitions.INPUT_SHORT_NAME, f.getAbsolutePath());
        });

        runCommandLine(args.getArgsList(), GatherPileupSummaries.class.getSimpleName());

        // Use PileupSummaryComparator to test
        final List<PileupSummary> scanned = PileupSummary.readFromFile(combinedPileupSummary).getRight();
        Assert.assertTrue(IntStream.range(0, scanned.size()-1).allMatch(i ->
                contigsOrdered.indexOf(scanned.get(i).getContig()) < contigsOrdered.indexOf(scanned.get(i+1).getContig()) ||
                        (contigsOrdered.indexOf(scanned.get(i).getContig()) == contigsOrdered.indexOf(scanned.get(i+1).getContig()) && scanned.get(i).getStart() <= scanned.get(i+1).getStart())));
    }

    @Test
    public void testAllelicCountsOutput() {
        // two scatters, given out of order, the first with a site that no reads cover, and an empty scatter
        final File emptyScatter = createTempFile("pileupSummary", ".tsv");
        PileupSummary.writeToFile("sample1", Collections.emptyList(), emptyScatter);
        final File firstScatter = createTempFile("pileupSummary", ".tsv");
        PileupSummary.writeToFile("sample1", Arrays.asList(
                new PileupSummary("1", 100, 10, 5, 1, 0.3, Nucleotide.A, Nucleotide.G),
                new PileupSummary("1", 200, 0, 0, 0, 0.1, Nucleotide.C, Nucleotide.T)), firstScatter);
        final File secondScatter = createTempFile("pileupSummary", ".tsv");
        PileupSummary.writeToFile("sample1", Arrays.asList(
                new PileupSummary("14", 50, 7, 8, 0, 0.5, Nucleotide.G, Nucleotide.A)), secondScatter);

        final File combinedPileupSummary = createTempFile("combined", ".tsv");
        final File allelicCounts = createTempFile("allelic-counts", ".tsv");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG19_DICT)
                .add(StandardArgumentDefinitions.INPUT_SHORT_NAME, secondScatter.getAbsolutePath())
                .add(StandardArgumentDefinitions.INPUT_SHORT_NAME, emptyScatter.getAbsolutePath())
                .add(StandardArgumentDefinitions.INPUT_SHORT_NAME, firstScatter.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, combinedPileupSummary.getAbsolutePath())
                .add(GatherPileupSummaries.ALLELIC_COUNTS_OUTPUT_LONG_NAME, allelicCounts.getAbsolutePath());
        runCommandLine(args.getArgsList(), GatherPileupSummaries.class.getSimpleName());

        // alt counts include other alts, as in CollectAllelicCounts
        final AllelicCountCollection result = new AllelicCountCollection(allelicCounts);
        Assert.assertEquals(result.getMetadata().getSampleName(), "sample1");
        Assert.assertTrue(CopyNumberArgumentValidationUtils.isSameDictionary(result.getMetadata().getSequenceDictionary(),
                ReferenceUtils.loadFastaDictionary(new File(FULL_HG19_DICT))));
        Assert.assertEquals(result.getRecords(), Arrays.asList(
                new AllelicCount(new SimpleInterval("1", 100, 100), 10, 6, Nucleotide.A, Nucleotide.G),
                new AllelicCount(new SimpleInterval("1", 200, 200), 0, 0, Nucleotide.C, Nucleotide.T),
                new AllelicCount(new SimpleInterval("14", 50, 50), 7, 8, Nucleotide.G, Nucleotide.A)));
    }

    @Test
    public void testScatterGather() throws IOException {
        final File intervalDir = createTempDir("intervals");
        final File pileupSummaryDir = createTempDir("pileupSummary");
        final String scatterCount = "20";

        // Step 1: SplitIntervals
        runCommandLine(Arrays.asList(
                "-R", b37_reference_20_21,
                "-L", NA12878_20_21_covered_regions,
                "-O", intervalDir.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount),
                SplitIntervals.class.getSimpleName());

        // Step 2: GetPileupSummaries
        final File[] intervals = intervalDir.listFiles();
        IntStream.range(0, intervals.length).forEach(i ->
            runCommandLine(Arrays.asList(
                    "-I", NA12878_20_21_WGS_bam,
                    "-V", thousandGenomes,
                    "-L", intervals[i].getAbsolutePath(),
                    "-O", pileupSummaryDir.getAbsolutePath() + "/" + i + ".tsv",
                    "-" + GetPileupSummaries.MAX_SITE_AF_SHORT_NAME, "0.9"),
                    GetPileupSummaries.class.getSimpleName())
        );

        // Step 3: GatherPileupSummaries
        final File combinedPileupSummary = createTempFile("combined", "tsv");
        final File[] pileupSummaries = pileupSummaryDir.listFiles();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG19_DICT);
        args.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, combinedPileupSummary.getAbsolutePath());
        Arrays.stream(pileupSummaries).forEach(f ->
                args.add(StandardArgumentDefinitions.INPUT_SHORT_NAME, f.getAbsolutePath()));

        runCommandLine(args.getArgsList(), GatherPileupSummaries.class.getSimpleName());

        final List<PileupSummary> combined = PileupSummary.readFromFile(combinedPileupSummary).getRight();

        Assert.assertTrue(IntStream.range(0, combined.size()-1).allMatch(i ->
                combined.get(i).getStart() <= combined.get(i+1).getStart() ||
                        Integer.valueOf(combined.get(i).getContig()) < Integer.valueOf(combined.get(i+1).getContig())));
    }
}