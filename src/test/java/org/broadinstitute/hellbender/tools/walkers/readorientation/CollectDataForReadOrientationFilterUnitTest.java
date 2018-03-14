package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static org.testng.Assert.*;

/**
 * Created by tsato on 3/14/18.
 */
public class CollectDataForReadOrientationFilterUnitTest extends CommandLineProgramTest {
    @Test
    public void testOnSyntheticBam() throws IOException {
        final File refMetrics = createTempFile("ref", ".metrics");
        final File altMetrics = createTempFile("alt", ".metrics");
        final File altTable = createTempFile("alt", ".table");

        final int numAltReads = 30;
        final int numRefReads = 70;
        final int depth = numAltReads + numRefReads;
        final File samFile = createSyntheticSam(numRefReads, numAltReads);

        final String[] args = {
                "-R", hg19_chr1_1M_Reference,
                "-I", samFile.getAbsolutePath(),
                "-" + CollectDataForReadOrientationFilter.ALT_DATA_TABLE_SHORT_NAME, altTable.getAbsolutePath(),
                "-" + CollectDataForReadOrientationFilter.REF_SITE_METRICS_SHORT_NAME, refMetrics.getAbsolutePath(),
                "-" + "alt-hist", altMetrics.getAbsolutePath()
        };

        runCommandLine(args);

        final MetricsFile<?, Integer> referenceSiteMetrics = new MetricsFile<>();
        final Reader in = IOUtil.openFileForBufferedReading(refMetrics);
        referenceSiteMetrics.read(in);
        CloserUtil.close(in);

        List<Histogram<Integer>> histograms = referenceSiteMetrics.getAllHistograms();
        List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(altTable);

        /** Expected result
         *
         * (position chr1: {@code alignmentStart})
         *
         * coord.:    99,990             100,000             100,010             100,020
         *              |                   |                   |                   |
         * reference: | T C A T C A C A C T C A C T A A G C A C A C A G A G A A T A A T
         * alt read:  | T C A T C A C A C T C T C T G A C C A A A G A G A A T A A T A A
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
        Assert.assertEquals(altDesignMatrix.size(), 4);
        AltSiteRecord recordCAC = altDesignMatrix.stream().
                filter(rec -> rec.getReferenceContext().equals("CAC")).findFirst().get();
        ArrayAsserts.assertArrayEquals(recordCAC.getBaseCounts(), new int[]{70, 0, 0, 30});
        ArrayAsserts.assertArrayEquals(recordCAC.getF1R2Counts(), new int[]{35, 0, 0, 15});

        AltSiteRecord recordTAA = altDesignMatrix.stream().
                filter(rec -> rec.getReferenceContext().equals("TAA")).findFirst().get();
        ArrayAsserts.assertArrayEquals(recordTAA.getBaseCounts(), new int[]{70, 0, 30, 0});
        ArrayAsserts.assertArrayEquals(recordTAA.getF1R2Counts(), new int[]{35, 0, 15, 0});

        AltSiteRecord recordAGC = altDesignMatrix.stream().
                filter(rec -> rec.getReferenceContext().equals("AGC")).findFirst().get();
        ArrayAsserts.assertArrayEquals(recordAGC.getBaseCounts(), new int[]{0, 30, 70, 0});
        ArrayAsserts.assertArrayEquals(recordAGC.getF1R2Counts(), new int[]{0, 15, 35, 0});

        AltSiteRecord recordACA = altDesignMatrix.stream().
                filter(rec -> rec.getReferenceContext().equals("ACA")).findFirst().get();
        ArrayAsserts.assertArrayEquals(recordACA.getBaseCounts(), new int[]{30, 70, 0, 0});
        ArrayAsserts.assertArrayEquals(recordACA.getF1R2Counts(), new int[]{15, 35, 0, 0});

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
        final File refMetrics = createTempFile("ref", ".metrics");
        final File altMetrics = createTempFile("alt", ".metrics");
        final File altTable = createTempFile("alt", ".table");
        final File sam = createSyntheticSam(30, 1);


        final String[] args = {
                "-R", hg19_chr1_1M_Reference,
                "-I", sam.getAbsolutePath(),
                "-" + CollectDataForReadOrientationFilter.ALT_DATA_TABLE_SHORT_NAME, altTable.getAbsolutePath(),
                "-" + CollectDataForReadOrientationFilter.REF_SITE_METRICS_SHORT_NAME, refMetrics.getAbsolutePath(),
                "-" + "alt-hist", altMetrics.getAbsolutePath()
        };

        runCommandLine(args);

        final MetricsFile<?, Integer> altSiteMetrics = new MetricsFile<>();
        final Reader altMetricsReader = IOUtil.openFileForBufferedReading(altMetrics);
        altSiteMetrics.read(altMetricsReader);
        CloserUtil.close(altMetricsReader);

        List<Histogram<Integer>> histograms = altSiteMetrics.getAllHistograms();

        // TODO: should there be 64 histograms in alt?
        final String[] expectedTransitions = new String[]{"CAC_T_F2R1", "TAA_G_F2R1", "AGC_C_F2R1", "ACA_A_F2R1"};
        Assert.assertEquals(histograms.size(), expectedTransitions.length, "alt histogram must contain the expected number of contexts");

        for (String transition : expectedTransitions ){
            Optional<Histogram<Integer>> histogram = histograms.stream()
                    .filter(h -> h.getValueLabel().equals(transition))
                    .findFirst();
            Assert.assertTrue(histogram.isPresent(), "histogram must exist");
            Assert.assertEquals(histogram.get().getSumOfValues() == 1, "histogram must only contain one read");
        }
    }

    private File createSyntheticSam(final int refDepth, final int altDepth) throws IOException {
        // create a sam header
        final int alignmentStart = 99_991;
        final int numChromosomes = 1;
        final int startingChromosome = 1;
        final int chromosomeSize = 1_000_000;
        final String readGroupName = "HELLO.1";
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(
                numChromosomes, startingChromosome, chromosomeSize);
        samHeader.addReadGroup(new SAMReadGroupRecord(readGroupName));

        // create a sample list
        final int chromosomeIndex = 0;
        final String sampleName = "samthreetree";
        final SampleList sampleList = new IndexedSampleList(sampleName);

        // specify characteristics of reads
        final int depth = altDepth + refDepth;

        final byte[] refReadBases = "CATCACACTCACTAAGCACACAGAGAATAAT".getBytes();
        final byte[] altReadBases = "CATCACACTCTCTGACCAAACAGAGAATAAT".getBytes();
        final int readLength = refReadBases.length;
        final byte baseq = 26;
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);
        final int mapq = 60;

        // create reads
        final List<GATKRead> reads = new ArrayList<>(depth);
        for (int i = 0; i < depth; i++) {
            final byte[] bases = i < altDepth ? altReadBases : refReadBases;
            final String cigar = refReadBases.length + "M";
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "Read" + i, chromosomeIndex,
                    alignmentStart, bases, quals, cigar);
            read.setReadGroup(readGroupName);
            read.setMappingQuality(mapq);
            read.setIsFirstOfPair();
            read.setIsReverseStrand(i % 2 == 0);
            reads.add(read);
        }

        final File samFile = File.createTempFile("synthetic", ".sam");
        final SAMFileGATKReadWriter writer = new SAMFileGATKReadWriter(
                ReadUtils.createCommonSAMWriter(samFile, null, samHeader, true, false, false));
        reads.forEach(writer::addRead);
        writer.close(); // closing the writer writes to the file

        return samFile;
    }

}