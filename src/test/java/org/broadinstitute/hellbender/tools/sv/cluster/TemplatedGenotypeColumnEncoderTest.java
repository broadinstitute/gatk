package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * Byte-identity tests for {@link TemplatedGenotypeColumnEncoder}: for randomized sparse records the VCF line
 * produced from the spliced, lazily-attached column string must equal, character for character, the line the
 * encoder produces from a record whose every sample was filled with a real {@link Genotype}.
 */
public class TemplatedGenotypeColumnEncoderTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICT = SVTestUtils.hg38Dict;
    private static final String FLOAT_FORMAT = "TSTF";
    private static final String STRING_LIST_FORMAT = "TSTS";
    private static final int NUM_SAMPLES = 48;

    private static VCFHeader header;
    private static TreeSet<String> samples;
    private static PloidyTable ploidyTable;

    private static synchronized void setUp() {
        if (header != null) {
            return;
        }
        samples = new TreeSet<>();
        final Map<String, Map<String, Integer>> ploidy = new HashMap<>();
        for (int i = 0; i < NUM_SAMPLES; i++) {
            final String s = String.format("S%03d", i);
            samples.add(s);
            final Map<String, Integer> contigPloidy = new HashMap<>();
            contigPloidy.put("chr1", 2);
            contigPloidy.put("chrX", i % 3 == 0 ? 1 : 2);   // two classes
            contigPloidy.put("chrY", i % 3 == 0 ? 1 : 0);   // includes ploidy 0
            ploidy.put(s, contigPloidy);
        }
        ploidyTable = new PloidyTable(ploidy);

        final Set<VCFHeaderLine> lines = new HashSet<>();
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
        lines.add(new VCFFormatHeaderLine(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "ECN"));
        lines.add(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "CN"));
        lines.add(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "RD_CN"));
        lines.add(new VCFFormatHeaderLine(FLOAT_FORMAT, 1, VCFHeaderLineType.Float, "test float"));
        lines.add(new VCFFormatHeaderLine(STRING_LIST_FORMAT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "test strings"));
        header = new VCFHeader(lines, samples);
        header.setSequenceDictionary(DICT);
    }

    private static VCFEncoder encoder() {
        // Missing INFO header lines are allowed so the test can focus on genotype columns; the same flag is
        // passed to the encoder under test, matching how the walker mirrors the writer's --lenient setting.
        return new VCFEncoder(header, true, false);
    }

    private static SVCallRecord makeRecord(final String id, final String contig, final int start,
                                           final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                           final List<Genotype> carriers) {
        final Allele alt;
        final int end;
        switch (type) {
            case DEL: alt = Allele.SV_SIMPLE_DEL; end = start + 5000; break;
            case DUP: alt = Allele.SV_SIMPLE_DUP; end = start + 5000; break;
            default: alt = Allele.SV_SIMPLE_INS; end = start; break;
        }
        final boolean ins = type == GATKSVVCFConstants.StructuralVariantAnnotationType.INS;
        return new SVCallRecord(id, contig, start, ins ? Boolean.TRUE : null, contig, end, ins ? Boolean.FALSE : null, type, null,
                Collections.emptyList(), type == GATKSVVCFConstants.StructuralVariantAnnotationType.INS ? 300 : null,
                Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, alt), carriers,
                Collections.singletonMap("TSTI", type.name()), Collections.emptySet(), null, DICT);
    }

    /** Random carrier genotype with a random subset of FORMAT fields. */
    private static Genotype randomCarrier(final String sample, final Allele alt, final Random rng) {
        final GenotypeBuilder gb = new GenotypeBuilder(sample);
        switch (rng.nextInt(4)) {
            case 0: gb.alleles(Arrays.asList(Allele.REF_N, alt)); break;
            case 1: gb.alleles(Arrays.asList(alt, alt)); break;
            case 2: gb.alleles(Collections.singletonList(alt)); break;
            default: gb.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)); break;
        }
        if (rng.nextBoolean()) gb.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1 + rng.nextInt(2));
        if (rng.nextInt(3) == 0) gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, rng.nextInt(4));
        if (rng.nextInt(3) == 0) gb.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, rng.nextInt(4));
        if (rng.nextInt(3) == 0) gb.GQ(rng.nextInt(100));
        if (rng.nextInt(4) == 0) gb.DP(rng.nextInt(60));
        if (rng.nextInt(5) == 0) gb.AD(new int[]{rng.nextInt(20), rng.nextInt(20)});
        if (rng.nextInt(5) == 0) gb.PL(new int[]{0, rng.nextInt(99), rng.nextInt(999)});
        if (rng.nextInt(5) == 0) gb.filter("LOW");
        if (rng.nextInt(4) == 0) gb.attribute(FLOAT_FORMAT, rng.nextDouble() * 10);
        if (rng.nextInt(4) == 0) gb.attribute(STRING_LIST_FORMAT, Arrays.asList("a", "b" + rng.nextInt(9)));
        if (rng.nextInt(6) == 0) gb.phased(true);
        return gb.make();
    }

    private static VariantContext fill(final SVCallRecord record, final Set<String> fillSamples, final boolean refDefault) {
        final GenotypesContext filled = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                record, fillSamples, refDefault, ploidyTable, header);
        return SVCallRecordUtils.getVariantBuilder(SVCallRecordUtils.copyCallWithNewGenotypes(record, filled)).make();
    }

    @DataProvider(name = "scenarios")
    public Object[][] scenarios() {
        return new Object[][]{
                {"chr1", true, 0.3, 1L},
                {"chr1", false, 0.3, 2L},
                {"chrX", true, 0.3, 3L},
                {"chrX", false, 0.5, 4L},
                {"chrY", true, 0.4, 5L},
                {"chr1", true, 0.0, 6L},   // no carriers at all
                {"chr1", true, 1.0, 7L},   // every sample a carrier: no representatives possible
        };
    }

    @Test(dataProvider = "scenarios")
    public void testExpandedLineIsByteIdentical(final String contig, final boolean refDefault,
                                                final double carrierFraction, final long seed) {
        setUp();
        final Random rng = new Random(seed);
        final TemplatedGenotypeColumnEncoder templated = new TemplatedGenotypeColumnEncoder(header, ploidyTable, true);
        final VCFEncoder reference = encoder();
        final GATKSVVCFConstants.StructuralVariantAnnotationType[] types = {
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS};
        for (int r = 0; r < 60; r++) {
            final GATKSVVCFConstants.StructuralVariantAnnotationType type = types[r % types.length];
            final Allele alt = type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL ? Allele.SV_SIMPLE_DEL
                    : type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP ? Allele.SV_SIMPLE_DUP : Allele.SV_SIMPLE_INS;
            final List<Genotype> carriers = new ArrayList<>();
            for (final String s : samples) {
                if (rng.nextDouble() < carrierFraction) {
                    carriers.add(randomCarrier(s, alt, rng));
                }
            }
            final SVCallRecord record = makeRecord("rec" + r, contig, 10_000 + 100 * r, type, carriers);

            // Reference: every sample gets a real genotype object, then the full encoder runs over all of them.
            final String expected = reference.encode(fill(record, samples, refDefault));

            // Under test: representatives only, then template splice.
            final Set<String> representatives = templated.representativeSamples(contig, record.getGenotypes().getSampleNames());
            for (final String rep : representatives) {
                Assert.assertFalse(record.getGenotypes().containsSample(rep), "representative must not be a carrier");
            }
            final VariantContext sparse = fill(record, representatives, refDefault);
            Assert.assertEquals(sparse.getNSamples(), carriers.size() + representatives.size());
            final VariantContext expanded = templated.expand(sparse, representatives);
            final String actual = reference.encode(expanded);

            Assert.assertEquals(actual, expected, "record " + r + " on " + contig);
            // Site fields untouched
            Assert.assertEquals(expanded.getContig(), record.getContigA());
            Assert.assertEquals(expanded.getStart(), record.getPositionA());
            Assert.assertEquals(expanded.getID(), record.getId());
        }
    }

    @Test
    public void testRepresentativesCoverEachPloidyClassOnce() {
        setUp();
        final TemplatedGenotypeColumnEncoder templated = new TemplatedGenotypeColumnEncoder(header, ploidyTable, true);
        // chrX has ploidy 1 (i % 3 == 0) and ploidy 2 samples
        final Set<String> reps = templated.representativeSamples("chrX", Collections.emptySet());
        Assert.assertEquals(reps.size(), 2);
        final Set<Integer> ploidies = new HashSet<>();
        for (final String s : reps) {
            ploidies.add(ploidyTable.get(s, "chrX"));
        }
        Assert.assertEquals(ploidies, new HashSet<>(Arrays.asList(1, 2)));
        // First non-carrier in header order is chosen: exclude S000 (ploidy 1) and expect S003
        final Set<String> reps2 = templated.representativeSamples("chrX", Collections.singleton("S000"));
        Assert.assertTrue(reps2.contains("S003"), reps2.toString());
        Assert.assertFalse(reps2.contains("S000"));
        // A class made entirely of carriers gets no representative
        final Set<String> allPloidy1 = new HashSet<>();
        for (final String s : samples) {
            if (ploidyTable.get(s, "chrX") == 1) {
                allPloidy1.add(s);
            }
        }
        final Set<String> reps3 = templated.representativeSamples("chrX", allPloidy1);
        Assert.assertEquals(reps3.size(), 1);
        Assert.assertEquals((int) ploidyTable.get(reps3.iterator().next(), "chrX"), 2);
    }

    @Test
    public void testFallbackDecodeMatchesFullGenotypes() {
        setUp();
        final Random rng = new Random(11L);
        final TemplatedGenotypeColumnEncoder templated = new TemplatedGenotypeColumnEncoder(header, ploidyTable, true);
        final List<Genotype> carriers = new ArrayList<>();
        for (final String s : samples) {
            if (rng.nextDouble() < 0.3) {
                carriers.add(randomCarrier(s, Allele.SV_SIMPLE_DEL, rng));
            }
        }
        final SVCallRecord record = makeRecord("fb", "chrX", 50_000,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, carriers);
        final VariantContext full = fill(record, samples, true);
        final Set<String> reps = templated.representativeSamples("chrX", record.getGenotypes().getSampleNames());
        final VariantContext expanded = templated.expand(fill(record, reps, true), reps);

        // Forcing a decode goes through the VCFCodec-based parser; compare against the full record after a
        // VCF round trip of its own so both sides have codec-normalized attribute types.
        final VariantContext fullRoundTrip = VariantContextTestUtils.readEntireVCFIntoMemory(
                writeSingleRecordVcf(full)).getValue().get(0);
        Assert.assertEquals(expanded.getNSamples(), NUM_SAMPLES);
        for (final String s : samples) {
            final Genotype a = expanded.getGenotype(s);
            final Genotype b = fullRoundTrip.getGenotype(s);
            Assert.assertNotNull(a, s);
            Assert.assertEquals(a.getAlleles(), b.getAlleles(), s);
            Assert.assertEquals(a.getExtendedAttributes(), b.getExtendedAttributes(), s);
            Assert.assertEquals(a.getGQ(), b.getGQ(), s);
            Assert.assertEquals(a.getDP(), b.getDP(), s);
        }
    }

    private static String writeSingleRecordVcf(final VariantContext vc) {
        final java.io.File out = createTempFile("templated_full", ".vcf");
        try (final htsjdk.variant.variantcontext.writer.VariantContextWriter writer =
                     new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
                             .setOutputFile(out).clearOptions()
                             .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                             .build()) {
            writer.writeHeader(header);
            writer.add(vc);
        }
        return out.getAbsolutePath();
    }

    /**
     * A sample absent from the ploidy table is legal in the full fill path as long as it always carries a genotype.
     * The templated path must accept exactly the same records and reject exactly the same ones.
     */
    @Test
    public void testSampleMissingFromPloidyTableMatchesFillPath() {
        setUp();
        final Map<String, Map<String, Integer>> partial = new HashMap<>();
        for (final String s : samples) {
            if (!s.equals("S005")) {
                partial.put(s, Collections.singletonMap("chr1", 2));
            }
        }
        final PloidyTable partialTable = new PloidyTable(partial);
        final TemplatedGenotypeColumnEncoder templated = new TemplatedGenotypeColumnEncoder(header, partialTable, true);
        final VCFEncoder reference = encoder();

        // S005 is a carrier: both paths succeed and agree
        final List<Genotype> carriers = new ArrayList<>();
        carriers.add(new GenotypeBuilder("S005", Arrays.asList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make());
        carriers.add(new GenotypeBuilder("S010", Arrays.asList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).GQ(50).make());
        final SVCallRecord ok = makeRecord("ok", "chr1", 20_000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, carriers);
        final GenotypesContext fullFilled = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                ok, samples, true, partialTable, header);
        final String expected = reference.encode(SVCallRecordUtils.getVariantBuilder(
                SVCallRecordUtils.copyCallWithNewGenotypes(ok, fullFilled)).make());
        final Set<String> reps = templated.representativeSamples("chr1", ok.getGenotypes().getSampleNames());
        Assert.assertFalse(reps.contains("S005"));
        final GenotypesContext sparseFilled = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                ok, reps, true, partialTable, header);
        final VariantContext sparse = SVCallRecordUtils.getVariantBuilder(
                SVCallRecordUtils.copyCallWithNewGenotypes(ok, sparseFilled)).make();
        Assert.assertEquals(reference.encode(templated.expand(sparse, reps)), expected);

        // S005 is not a carrier: the fill path throws, and so must the templated path
        final SVCallRecord bad = makeRecord("bad", "chr1", 30_000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                Collections.singletonList(carriers.get(1)));
        Assert.assertThrows(IllegalArgumentException.class, () ->
                SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(bad, samples, true, partialTable, header));
        final Set<String> badReps = templated.representativeSamples("chr1", bad.getGenotypes().getSampleNames());
        final GenotypesContext badSparseFilled = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                bad, badReps, true, partialTable, header);
        final VariantContext badSparse = SVCallRecordUtils.getVariantBuilder(
                SVCallRecordUtils.copyCallWithNewGenotypes(bad, badSparseFilled)).make();
        Assert.assertThrows(IllegalArgumentException.class, () -> templated.expand(badSparse, badReps));
    }
}
