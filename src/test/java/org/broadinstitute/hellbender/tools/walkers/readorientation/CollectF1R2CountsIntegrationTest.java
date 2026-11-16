package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2TestingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by tsato on 3/14/18.
 */
public class CollectF1R2CountsIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testOnSyntheticBam() throws IOException {
        final File outputTarGz = createTempFile("f1r2", ".tar.gz");
        final String sample = "SAMPLE";

        final int numAltReads = 30;
        final int numRefReads = 70;
        final int depth = numAltReads + numRefReads;
        final File samFile = createSyntheticSam(numRefReads, numAltReads, sample);

        final String[] args = {
                "-R", hg19_chr1_1M_Reference,
                "-I", samFile.getAbsolutePath(),
                "-O", outputTarGz.getAbsolutePath()
        };

        runCommandLine(args);

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(outputTarGz.toPath(), extractedDir.toPath());

        final File refHistogramTable = F1R2CountsCollector.getRefHistogramsFromExtractedTar(extractedDir).stream().findFirst().get();

        // the reads are all in the first block of contig 1
        Assert.assertEquals(getBlocks(refHistogramTable), Set.of("1:1"));

        final List<Histogram<Integer>> histograms = LearnReadOrientationModel.sumHistogramTables(List.of(refHistogramTable),
                false, F1R2FilterConstants.DEFAULT_MAX_DEPTH).get(sample);

        /** Expected result
         *
         * (position chr1: {@code alignmentStart})
         *
         * coord.:    99,990             100,000             100,010             100,020
         *              |                   |                   |                   |
         * reference:   T C A T C A C A C T C A C T A A G C A C A C A G A G A A T A A T
         * alt read:    T C A T C A C A C T C T C T G A C C A A A G A G A A T A A T A A
         *                                    *     *   *     *
         *
         * At each alt site, we have 75 ref and 25 alt reads
         *
         * context, alt counts, alt f1r2 counts, depth
         * -------------------------------------------
         * CAC, {70, 0, 0, 30}, {35, 0, 0, 15}, 100
         * TAA, {70, 0, 30, 0}, {35, 0, 15, 0}, 100
         * AGC, {0, 30, 70, 0}, {0, 15, 35, 0}, 100
         * ACA, {30, 70, 0, 0}, {15, 35, 0, 0}, 100
         *
         *
         * Ref: 6 sites at depth 100
         * TCA, CAT, ATC, TCA, CAC, ACA, CAC, ACT, CTC, TCA, ACT, CTA, AAG, GCA, CAC, CAC, ACA, CAG, AGA, GAG, AGA, GAA, AAT, ATA, TAA, AAT
         * TCA*3, CAT, ATC, CAC*4, ACA*2, ACT*2, CTC, CTA, AAG, GCA, CAG, AGA*2, GAG, GAA, AAT*2, ATA, TAA
         *
         **/

        // alt site tests
        final int totalAltRecords = LearnReadOrientationModel.gatherAltSiteRecords(F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir)).size();
        Assert.assertEquals(totalAltRecords, 4);

        // the reads start at position 99,991 of contig 1
        final Map<String, Integer> positions = Map.of("CAC", 100_001, "TAA", 100_004, "AGC", 100_006, "ACA", 100_009);
        for (final String context : List.of("CAC", "TAA", "AGC", "ACA")) {
            final File altTableFile = new File(extractedDir, F1R2CountsCollector.getAltTableName(sample, context));
            final List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(altTableFile.toPath(), context).getRight();
            Assert.assertEquals(altDesignMatrix.size(), 1);
            final AltSiteRecord record = altDesignMatrix.get(0);
            Assert.assertEquals(record.getReferenceContext(), context);
            Assert.assertEquals(record.getRefCount(), 70);
            Assert.assertEquals(record.getAltCount(), 30);
            Assert.assertEquals(record.getRefF1R2(), 35);
            Assert.assertEquals(record.getAltF1R2(), 15);
            Assert.assertEquals(record.getContig(), "1");
            Assert.assertEquals(record.getPosition(), (int) positions.get(context));
        }

        // check ref histograms
        // ImmutableMap.of can take up to 5 elements, and ImmutableMap.Builder is ugly
        final ImmutableMap<String, Integer> expectedReferenceCounts1 = ImmutableMap.of(
                "TCA", 3,
                "CAT", 1,
                "ATC", 1,
                "CAC", 4,
                "ACA", 2);
        final ImmutableMap<String, Integer>  expectedReferenceCounts2 = ImmutableMap.of(
                "ACT", 2,
                "CTC", 1,
                "CTA", 1,
                "AAG", 1,
                "GCA", 1);
        final ImmutableMap<String, Integer>  expectedReferenceCounts3 = ImmutableMap.of(
                "CAG", 1,
                "AGA", 2,
                "GAG", 1,
                "GAA", 1,
                "AAT", 2);
        final ImmutableMap<String, Integer>  expectedReferenceCounts4 = ImmutableMap.of("ATA", 1, "TAA", 1);

        for (ImmutableMap<String, Integer> expectedContextCountPairs : Arrays.asList(
                expectedReferenceCounts1, expectedReferenceCounts2, expectedReferenceCounts3, expectedReferenceCounts4)){
            expectedContextCountPairs.entrySet().stream().forEach(entrySet -> {
                final String context = entrySet.getKey();
                final double expectedCount = entrySet.getValue();
                Histogram<Integer> histogram = histograms.stream()
                        .filter(hist -> hist.getValueLabel().equals(context))
                        .findFirst().get();
                Assert.assertEquals(histogram.get(depth).getValue(), expectedCount);
                Assert.assertEquals(histogram.getSumOfValues(), expectedCount );
            });
        }
    }

    @Test
    public void testHistograms() throws IOException {
        final File outputTarGz = createTempFile("f1r2", ".tar.gz");
        final String sample = "SAMPLE";
        final File sam = createSyntheticSam(30, 1, sample);

        final String[] args = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam.getAbsolutePath(),
                "-O", outputTarGz.getAbsolutePath()
        };

        runCommandLine(args);

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(outputTarGz.toPath(), extractedDir.toPath());

        final File altHistogramTable = F1R2CountsCollector.getAltHistogramsFromExtractedTar(extractedDir).stream().findFirst().get();

        final List<Histogram<Integer>> histograms = LearnReadOrientationModel.sumHistogramTables(List.of(altHistogramTable),
                true, F1R2FilterConstants.DEFAULT_MAX_DEPTH).get(sample);

        // TODO: should there be 64*3*2 = 384 histograms or just the non-zero ones?
        final String[] expectedTransitions = new String[]{"CAC_T_F2R1", "TAA_G_F2R1", "AGC_C_F2R1", "ACA_A_F2R1"};
        // Assert.assertEquals(histograms.size(), expectedTransitions.length, "alt histogram must contain the expected number of contexts");

        for (String transition : expectedTransitions ){
            Optional<Histogram<Integer>> histogram = histograms.stream()
                    .filter(h -> h.getValueLabel().equals(transition))
                    .findFirst();
            Assert.assertTrue(histogram.isPresent(), "histogram must exist");
            Assert.assertEquals((int) histogram.get().getSumOfValues(), 1, "histogram must only contain one read");
        }
    }

    @Test
    public void testTwoSamples() throws IOException {
        final File outputTarGZ1 = createTempFile("f1r21", ".tar.gz");
        final File outputTarGZ2 = createTempFile("f1r22",".tar.gz");
        final File outputTarGZ12 = createTempFile("f1r212",".tar.gz");
        final String sample1 = "SAMPLE1";
        final String sample2 = "SAMPLE2";
        final File sam1 = createSyntheticSam(10, 1, sample1);
        final File sam2 = createSyntheticSam(20, 2, sample2);

        final String[] args1 = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam1.getAbsolutePath(),
                "-O", outputTarGZ1.getAbsolutePath()
        };

        final String[] args2 = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam2.getAbsolutePath(),
                "-O", outputTarGZ2.getAbsolutePath()
        };

        final String[] args12 = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam1.getAbsolutePath(),
                "-I", sam2.getAbsolutePath(),
                "-O", outputTarGZ12.getAbsolutePath()
        };

        runCommandLine(args1);
        runCommandLine(args2);
        runCommandLine(args12);
        final File extractedDir1 = createTempDir("extracted1");
        IOUtils.extractTarGz(outputTarGZ1.toPath(), extractedDir1.toPath());

        final File extractedDir2 = createTempDir("extracted2");
        IOUtils.extractTarGz(outputTarGZ2.toPath(), extractedDir2.toPath());

        final File extractedDir12 = createTempDir("extracted12");
        IOUtils.extractTarGz(outputTarGZ12.toPath(), extractedDir12.toPath());

        // each sample's counts are the same whether it is collected alone or together with the other sample
        for (final Map.Entry<String, File> sampleAndDir : Map.of(sample1, extractedDir1, sample2, extractedDir2).entrySet()) {
            final String sample = sampleAndDir.getKey();
            for (final String context : F1R2FilterConstants.ALL_KMERS) {
                final String altTable = F1R2CountsCollector.getAltTableName(sample, context);
                IntegrationTestSpec.assertEqualTextFiles(new File(sampleAndDir.getValue(), altTable), new File(extractedDir12, altTable));
            }

            for (final String extension : List.of(F1R2CountsCollector.REF_HIST_EXTENSION, F1R2CountsCollector.ALT_HIST_EXTENSION)) {
                IntegrationTestSpec.assertEqualTextFiles(new File(sampleAndDir.getValue(), sample + extension),
                        new File(extractedDir12, sample + extension));
            }
        }
    }

    // The histogram tables have rows for each block of the genome with data, and the alt tables have the sites' positions
    @Test
    public void testBlocks() throws IOException {
        final File outputTarGz = createTempFile("f1r2", ".tar.gz");
        runCommandLine(new ArgumentsBuilder().addReference(b37_reference_20_21)
                .addInput(new File(LearnReadOrientationModelIntegrationTest.hapmapBamSnippet))
                .addIntervals(new File(LearnReadOrientationModelIntegrationTest.intervalList))
                .addOutput(outputTarGz));

        final File extractedDir = createTempDir("extracted");
        IOUtils.extractTarGz(outputTarGz.toPath(), extractedDir.toPath());

        // the intervals of chromosome 20, each in its own block
        final List<SimpleInterval> intervals = List.of(new SimpleInterval("20", 5296760, 5297096), new SimpleInterval("20", 23420905, 23421244),
                new SimpleInterval("20", 25755641, 25755987), new SimpleInterval("20", 34144590, 34144918), new SimpleInterval("20", 62165386, 62165688));
        final Set<String> blocks = Set.of("20:5000001", "20:23000001", "20:25000001", "20:34000001", "20:62000001");

        Assert.assertEquals(getBlocks(F1R2CountsCollector.getRefHistogramsFromExtractedTar(extractedDir).get(0)), blocks);
        Assert.assertEquals(getBlocks(F1R2CountsCollector.getAltHistogramsFromExtractedTar(extractedDir).get(0)), blocks);

        final List<AltSiteRecord> altSites = LearnReadOrientationModel.gatherAltSiteRecords(F1R2CountsCollector.getAltTablesFromExtractedTar(extractedDir));
        Assert.assertFalse(altSites.isEmpty());
        for (final AltSiteRecord altSite : altSites) {
            final SimpleInterval locus = new SimpleInterval(altSite.getContig(), altSite.getPosition(), altSite.getPosition());
            Assert.assertTrue(intervals.stream().anyMatch(interval -> interval.contains(locus)), locus.toString());
        }
    }

    // the blocks of the genome of the rows of a histogram table, as contig:block start
    private static Set<String> getBlocks(final File histogramTable) throws IOException {
        try (TableReader<F1R2HistogramTable.Row> reader = F1R2HistogramTable.reader(histogramTable.toPath())) {
            return reader.stream().map(row -> row.contig() + ":" + row.blockStart()).collect(Collectors.toSet());
        }
    }

    static File createSyntheticSam(final int refDepth, final int altDepth, final String sampleName) throws IOException {
        final File samFile = IOUtils.createTempFile("synthetic", ".bam");
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader(sampleName);
        final SAMFileGATKReadWriter writer = M2TestingUtils.getBareBonesSamWriter(samFile, samHeader);
        //            Ref Sequence: "CATCACACTCACTAAGCACACAGAGAATAAT".getBytes();
        //          SNPs positions:            *  * *  *
        final byte[] altReadBases = "CATCACACTCTCTGACCAAACAGAGAATAAT".getBytes();
        final List<GATKRead> refReads = M2TestingUtils.createReads(refDepth, M2TestingUtils.DEFAULT_REF_BASES, samHeader, (byte)30, "ref");
        final List<GATKRead> alt1Reads = M2TestingUtils.createReads(altDepth, altReadBases, samHeader, (byte)30, "alt");
        refReads.forEach(writer::addRead);
        alt1Reads.forEach(writer::addRead);
        writer.close(); // closing the writer writes to the file

        return samFile;
    }

}
