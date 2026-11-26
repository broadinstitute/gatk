package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.readorientation.CollectF1R2CountsArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2CountsCollector;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Mutect2 accumulates the evidence of each locus one read at a time; check it against the pileups of
 * {@link LocusIteratorByState}, from which Mutect2 used to read the same evidence one locus at a time.
 */
public class LocusEvidenceUnitTest extends GATKBaseTest {
    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final int PCR_SNV_QUAL = 40;
    private static final int SUBSTITUTION_QUAL_CORRECTION = (int) Math.round(10 * Math.log10(3));

    @DataProvider
    public Object[][] realReads() {
        return new Object[][] {
                // tumor, normal (or null), reference, interval
                {DREAM_BAMS_DIR + "tumor_4.bam", DREAM_BAMS_DIR + "normal_4.bam", b37_reference_20_21, new SimpleInterval("20", 240_000, 480_000)},
                {DREAM_BAMS_DIR + "tumor_3.bam", DREAM_BAMS_DIR + "normal_3.bam", b37_reference_20_21, new SimpleInterval("20", 1_100_000, 1_400_000)},
                {DREAM_BAMS_DIR + "tumor.bam", DREAM_BAMS_DIR + "normal.bam", b37_reference_20_21, new SimpleInterval("20", 10_000_000, 10_100_000)},
                {largeFileTestDir + "mutect/highDPMTsnippet.bam", null, toolsTestDir + "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta", new SimpleInterval("chrM", 1, 2_000)}
        };
    }

    @Test(dataProvider = "realReads")
    public void testRealReads(final String tumorBam, final String normalBam, final String referencePath, final SimpleInterval interval) {
        final List<java.nio.file.Path> bams = new ArrayList<>(List.of(IOUtils.getPath(tumorBam)));
        if (normalBam != null) {
            bams.add(IOUtils.getPath(normalBam));
        }
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(bams);
             final ReferenceDataSource reference = ReferenceDataSource.of(IOUtils.getPath(referencePath))) {
            final SAMFileHeader header = readsSource.getHeader();
            final ReadFilter readFilter = ReadFilter.fromList(Mutect2Engine.makeStandardMutect2ReadFilters(), header);
            final List<GATKRead> reads = Utils.stream(readsSource.query(interval)).filter(readFilter).collect(Collectors.toList());
            final Set<String> normalSamples = normalBam == null ? Set.of() : Set.of(onlySample(normalBam));
            Assert.assertTrue(reads.size() > 1000);
            compareToPileups(reads, header, reference, normalSamples, interval, true);
            compareToPileups(reads, header, reference, normalSamples, interval, false);
        }
    }

    private static String onlySample(final String bam) {
        try (final ReadsDataSource source = new ReadsPathDataSource(IOUtils.getPath(bam))) {
            return ReadUtils.getSamplesFromHeader(source.getHeader()).iterator().next();
        }
    }

    /**
     * Random alignments that exercise the rules: soft clips with low-quality neighbors, insertions and deletions next to
     * each other and to clips, reference skips, hard clips, short fragments with adaptor bases, mates that overlap, base
     * qualities of 128 or more, and N bases.
     */
    @Test
    public void testRandomReads() {
        try (final ReferenceDataSource reference = ReferenceDataSource.of(IOUtils.getPath(b37_reference_20_21))) {
            final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(reference.getSequenceDictionary());
            for (final String sample : List.of("tumor", "normal")) {
                final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(sample + "_rg");
                readGroup.setSample(sample);
                header.addReadGroup(readGroup);
            }

            final Random rng = new Random(13);
            final int intervalStart = 1_000_000;
            final int intervalEnd = 1_003_000;
            final List<GATKRead> reads = new ArrayList<>();
            for (int n = 0; n < 4000; n++) {
                reads.add(randomRead(rng, header, reference, "read" + n, intervalStart - 200 + rng.nextInt(intervalEnd - intervalStart + 200)));
            }
            reads.sort(Comparator.comparingInt(GATKRead::getStart));
            for (final boolean collectF1R2 : List.of(true, false)) {
                compareToPileups(reads, header, reference, Set.of("normal"), new SimpleInterval("20", intervalStart, intervalEnd), collectF1R2);
            }
        }
    }

    private static GATKRead randomRead(final Random rng, final SAMFileHeader header, final ReferenceDataSource reference,
                                       final String name, final int start) {
        final StringBuilder cigar = new StringBuilder();
        final List<Character> operators = new ArrayList<>();
        final List<Integer> lengths = new ArrayList<>();
        if (rng.nextDouble() < 0.1) { operators.add('H'); lengths.add(1 + rng.nextInt(5)); }
        if (rng.nextDouble() < 0.3) { operators.add('S'); lengths.add(1 + rng.nextInt(10)); }
        if (rng.nextDouble() < 0.05) { operators.add('I'); lengths.add(1 + rng.nextInt(3)); }
        operators.add('M'); lengths.add(1 + rng.nextInt(60));
        final int numIndels = rng.nextInt(4);
        for (int k = 0; k < numIndels; k++) {
            final double r = rng.nextDouble();
            operators.add(r < 0.4 ? 'I' : (r < 0.9 ? 'D' : 'N'));
            lengths.add(1 + rng.nextInt(r < 0.9 ? 6 : 50));
            if (rng.nextDouble() < 0.1) {   // occasionally two adjacent indels
                operators.add(rng.nextBoolean() ? 'I' : 'D');
                lengths.add(1 + rng.nextInt(3));
            }
            operators.add('M'); lengths.add(1 + rng.nextInt(60));
        }
        if (rng.nextDouble() < 0.3) { operators.add('S'); lengths.add(1 + rng.nextInt(10)); }
        if (rng.nextDouble() < 0.1) { operators.add('H'); lengths.add(1 + rng.nextInt(5)); }

        int readLength = 0;
        int referenceLength = 0;
        for (int k = 0; k < operators.size(); k++) {
            cigar.append(lengths.get(k)).append(operators.get(k));
            final char operator = operators.get(k);
            readLength += (operator == 'M' || operator == 'I' || operator == 'S') ? lengths.get(k) : 0;
            referenceLength += (operator == 'M' || operator == 'D' || operator == 'N') ? lengths.get(k) : 0;
        }

        // mostly reference bases, so that substitutions are a minority, and qualities that straddle the thresholds
        final byte[] referenceBases = reference.queryAndPrefetch("20", start, start + referenceLength).getBases();
        final byte[] bases = new byte[readLength];
        final byte[] quals = new byte[readLength];
        final byte[] alphabet = {'A', 'C', 'G', 'T', 'N'};
        for (int i = 0; i < readLength; i++) {
            bases[i] = rng.nextDouble() < 0.8 ? referenceBases[Math.min(i, referenceBases.length - 1)] : alphabet[rng.nextInt(alphabet.length)];
            final double r = rng.nextDouble();
            quals[i] = (byte) (r < 0.1 ? rng.nextInt(8) : (r < 0.12 ? 128 + rng.nextInt(100) : 10 + rng.nextInt(40)));
        }

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, start, bases, quals, cigar.toString());
        read.setReadGroup(rng.nextBoolean() ? "tumor_rg" : "normal_rg");
        read.setMappingQuality(20 + rng.nextInt(41));
        read.setIsReverseStrand(rng.nextBoolean());
        if (rng.nextDouble() < 0.85) {
            read.setIsPaired(true);
            read.setIsProperlyPaired(rng.nextDouble() < 0.8);
            if (rng.nextBoolean()) {
                read.setIsFirstOfPair();
            } else {
                read.setIsSecondOfPair();
            }
            if (rng.nextDouble() < 0.05) {
                read.setMateIsUnmapped();
            } else {
                read.setMatePosition("20", Math.max(1, start - 150 + rng.nextInt(300)));
                read.setMateIsReverseStrand(rng.nextDouble() < 0.8 != read.isReverseStrand());
                // many short fragments, whose adaptor bases pileups leave out
                final int fragmentLength = rng.nextDouble() < 0.4 ? 20 + rng.nextInt(90) : 100 + rng.nextInt(300);
                read.setFragmentLength(read.isReverseStrand() ? -fragmentLength : fragmentLength);
            }
        }
        return read;
    }

    /**
     * Accumulate the evidence of reads as Mutect2's traversal does, adding the reads that start at or before each locus
     * before reading its evidence, and check every locus of an interval against its pileup
     *
     * @param collectF1R2 whether to collect, and check, the F1R2 statistics of the tumor samples, which changes how
     *                    the evidence of their reads is accumulated
     */
    private static void compareToPileups(final List<GATKRead> reads, final SAMFileHeader header, final ReferenceDataSource reference,
                                         final Set<String> normalSamples, final SimpleInterval interval, final boolean collectF1R2) {
        final List<String> samples = new ArrayList<>(ReadUtils.getSamplesFromHeader(header));
        final List<String> tumorSamples = samples.stream().filter(s -> !normalSamples.contains(s)).collect(Collectors.toList());
        final CollectF1R2CountsArgumentCollection f1r2Args = new CollectF1R2CountsArgumentCollection();
        final F1R2CountsCollector f1r2Collector = collectF1R2 ? new F1R2CountsCollector(f1r2Args, header, createTempFile("f1r2", ".tar.gz"), tumorSamples) : null;
        final LocusEvidence evidence = new LocusEvidence(header, interval.getContig(), interval.getStart(), reference, null,
                normalSamples, PCR_SNV_QUAL, SUBSTITUTION_QUAL_CORRECTION, f1r2Collector, null);

        final LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), DownsamplingMethod.NONE, samples, header, true);
        final IntervalAlignmentContextIterator pileups = new IntervalAlignmentContextIterator(libs,
                new IntervalLocusIterator(List.of(interval).iterator()), header.getSequenceDictionary());

        int nextRead = 0;
        int lociWithEvents = 0;
        for (int locus = interval.getStart(); locus <= interval.getEnd(); locus++) {
            evidence.discardLociBefore(locus);
            while (nextRead < reads.size() && reads.get(nextRead).getStart() <= locus) {
                evidence.add(reads.get(nextRead++));
            }

            final AlignmentContext context = pileups.next();
            Assert.assertEquals(context.getStart(), locus);
            final ReadPileup pileup = context.getBasePileup();
            final byte refBase = reference.queryAndPrefetch(interval.getContig(), locus, locus).getBases()[0];
            final String where = interval.getContig() + ":" + locus;
            for (final boolean normal : List.of(false, true)) {
                final ReadPileup samplesPileup = pileup.getPileupForSamples(s -> normalSamples.contains(s) == normal, header);
                Assert.assertEquals(evidence.depth(locus, normal), samplesPileup.size(), where);
                final List<ByteArrayList> expected = isActiveBuffers(samplesPileup, refBase);
                final Mutect2Engine.PileupQualBuffer altQuals = evidence.altQuals(locus, normal);
                for (int altAllele = 0; altAllele < expected.size(); altAllele++) {
                    Assert.assertEquals(altQuals.quals(altAllele), expected.get(altAllele), where + " allele " + altAllele);
                    lociWithEvents += expected.get(altAllele).isEmpty() ? 0 : 1;
                }
            }

            for (final String sample : collectF1R2 ? tumorSamples : List.<String>of()) {
                final F1R2CountsCollector.PileupStatistics expected = f1r2Collector.makePileupStatistics();
                for (final PileupElement pe : pileup.getPileupForSample(sample, header)) {
                    if (pe.getQual() > f1r2Args.minBaseQuality) {
                        expected.add(pe.isDeletion(), pe.getBase(), ReadUtils.isF1R2(pe.getRead()),
                                pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart(), pe.getMappingQual());
                    }
                }
                Assert.assertEquals(evidence.f1r2Statistics(sample, locus).toString(), expected.toString(), where);
            }
        }
        Assert.assertFalse(pileups.hasNext());
        // so that the test doesn't pass trivially
        Assert.assertTrue(lociWithEvents > 10);
    }

    /**
     * The qualities that Mutect2's isActive accumulated from a pileup for each alt allele -- A, C, G, T, another
     * substitution and an indel -- copied from its PileupQualBuffer
     */
    private static List<ByteArrayList> isActiveBuffers(final ReadPileup pileup, final byte refBase) {
        final int otherSubstitution = 4;
        final int indel = 5;
        final List<ByteArrayList> buffers = new ArrayList<>();
        for (int n = 0; n < 6; n++) {
            buffers.add(new ByteArrayList());
        }
        final int position = pileup.getLocation().getStart();
        for (final PileupElement pe : pileup) {
            final int indelLength = pe.isDeletion() ? pe.getCurrentCigarElement().getLength() : pe.getLengthOfImmediatelyFollowingIndel();
            if (indelLength > 0) {
                buffers.get(indel).add(Mutect2Engine.indelQual(indelLength));
            } else if (isNextToUsefulSoftClip(pe)) {
                buffers.get(indel).add(Mutect2Engine.indelQual(1));
            } else if (pe.getBase() != refBase && pe.getQual() > Mutect2Engine.MINIMUM_BASE_QUALITY) {
                final GATKRead read = pe.getRead();
                final int mateStart = (!read.isProperlyPaired() || read.mateIsUnmapped()) ? Integer.MAX_VALUE : read.getMateStart();
                final boolean overlapsMate = mateStart <= position && position < mateStart + read.getLength();
                final byte qual = overlapsMate ? (byte) FastMath.min(pe.getQual(), PCR_SNV_QUAL / 2) : pe.getQual();
                final int index = BaseUtils.simpleBaseToBaseIndex(pe.getBase());
                if (index == -1) {
                    buffers.get(otherSubstitution).add(qual);
                } else {
                    buffers.get(index).add((byte) FastMath.min(qual + (double) SUBSTITUTION_QUAL_CORRECTION, QualityUtils.MAX_QUAL));
                }
            }
        }
        return buffers;
    }

    private static boolean isNextToUsefulSoftClip(final PileupElement pe) {
        final int offset = pe.getOffset();
        return pe.getQual() > Mutect2Engine.MINIMUM_BASE_QUALITY &&
                ((pe.isBeforeSoftClip() && pe.getRead().getBaseQuality(offset + 1) > Mutect2Engine.MINIMUM_BASE_QUALITY)
                        || (pe.isAfterSoftClip() && pe.getRead().getBaseQuality(offset - 1) > Mutect2Engine.MINIMUM_BASE_QUALITY));
    }
}
