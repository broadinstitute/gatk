package org.broadinstitute.hellbender.tools.exome;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration test for {@link GenotypeCopyNumberTriStateSegments}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeCopyNumberTriStateSegmentsIntegrationTest extends CopyNumberTriStateSegmentsCallerIntegrationTest {

    @Override
    public String getTestedClassName() {
        return GenotypeCopyNumberTriStateSegments.class.getSimpleName();
    }

    @Test(dataProvider = "simulateChainData")
    public void testRun(final HiddenMarkovModelChain chain)
      throws IOException {
        final DiscoverCopyNumberTriStateSegmentsIntegrationTest discovery = new DiscoverCopyNumberTriStateSegmentsIntegrationTest();
        final File inputFile = CopyNumberTriStateSegmentsCallerIntegrationTest.writeChainInTempFile(chain);
        final File segmentsFile = createTempFile("segments", ".tab");
        final File outputFile = createTempFile("output", ".vcf");
        discovery.runCommandLine(chain, inputFile, segmentsFile);
        Assert.assertTrue(segmentsFile.exists());
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.getAbsolutePath());
        arguments.add("-" + GenotypeCopyNumberTriStateSegments.DISCOVERY_FILE_SHORT_NAME);
        arguments.add(segmentsFile.getAbsolutePath());
        arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
        arguments.add(CopyNumberTriStateSegmentsCallerIntegrationTest.REALISTIC_TARGETS_FILE.getAbsolutePath());
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        loadModelArguments(chain, arguments);
        arguments.add(String.valueOf(CopyNumberTriStateSegmentCaller.ZScoreDimension.NONE.toString()));
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        assertSegmentsAndOutputConcordance(segmentsFile, outputFile);
    }

    private void assertSegmentsAndOutputConcordance(final File segmentsFile, final File outputFile) throws IOException {
        final VCFFileReader outputReader = new VCFFileReader(outputFile, false);
        final List<CopyNumberTriStateSegmentRecord> segments = readSegmentRecords(segmentsFile);
        assertHeader(outputReader, segments);
        final List<CopyNumberTriStateSegmentRecord> variantSegments = segments.stream()
                .filter(s -> s.getSegment().getCall() != CopyNumberTriState.NEUTRAL)
                .collect(Collectors.toList());
        final List<VariantContext> variants = new ArrayList<>();
        outputReader.forEach(variants::add);
        assertVariantsOrder(variants);
        assertVariantsPlGtAndGQAreConsistent(variants);
        assertVariantsInfoFieldsAreConsistent(variants);
        assertVariantSegmentsAreCovered(variants, variantSegments);
        assertVariantsAreCoveredBySegments(variants, variantSegments);
    }

    private void assertVariantsAreCoveredBySegments(final List<VariantContext> variants, final List<CopyNumberTriStateSegmentRecord> variantSegments) {
        for (final VariantContext variant : variants) {
            final List<CopyNumberTriStateSegmentRecord> matches =
                    variantSegments.stream()
                    .filter(s -> new SimpleInterval(variant).equals(s.getSegment().getInterval()))
                    .collect(Collectors.toList());
            Assert.assertFalse(matches.isEmpty());
            for (final Genotype genotype : variant.getGenotypes()) {
                final boolean discovery = genotype.getExtendedAttribute(GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY).toString().equals(GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE);
                if (discovery) {
                    Assert.assertTrue(matches.stream().anyMatch(s -> s.getSampleName().equals(genotype.getSampleName())));
                } else {
                    Assert.assertTrue(matches.stream().noneMatch(s -> s.getSampleName().equals(genotype.getSampleName())));
                }
            }
        }
    }

    private void assertVariantsInfoFieldsAreConsistent(final List<VariantContext> variants) {
        for (final VariantContext variant : variants) {
            final int expectedAN = variant.getGenotypes().size();
            final int[] expectedAC = new int[CopyNumberTriState.values().length];
            final List<Allele> alleles = variant.getAlleles();
            for (final Genotype genotype : variant.getGenotypes()) {
                Assert.assertEquals(genotype.getAlleles().size(), 1);
                final int alleleIndex = alleles.indexOf(genotype.getAllele(0));
                Assert.assertTrue(alleleIndex >= 0);
                expectedAC[alleleIndex]++;
            }
            Assert.assertEquals(variant.getAlleles(), GenotypeCopyNumberTriStateSegments.CNVAllele.ALL_ALLELES);
            Assert.assertTrue(variant.hasAttribute(GenotypeCopyNumberTriStateSegments.NUMBER_OF_TARGETS_KEY));
            Assert.assertTrue(variant.hasAttribute(VCFConstants.ALLELE_COUNT_KEY));
            Assert.assertTrue(variant.hasAttribute(VCFConstants.END_KEY));
            Assert.assertTrue(variant.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
            Assert.assertTrue(variant.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY));
            final int expectedANAnotherWay = IntStream.of(expectedAC).sum();
            Assert.assertEquals(expectedANAnotherWay, expectedAN);
            Assert.assertEquals(variant.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, -1), expectedAN);
            final double[] expectedAF = IntStream.of(expectedAC).mapToDouble(c -> c / (double) expectedAN).toArray();
            final double[] observedAFWithoutRef = variant.getAttributeAsList(VCFConstants.ALLELE_FREQUENCY_KEY).stream()
                    .mapToDouble(o -> Double.parseDouble(String.valueOf(o))).toArray();
            Assert.assertEquals(observedAFWithoutRef.length, expectedAF.length - 1);
            for (int i = 0; i < observedAFWithoutRef.length; i++) {
                Assert.assertEquals(observedAFWithoutRef[i], expectedAF[i+1], 0.001);
            }
            final int[] observedACWithoutRef = variant.getAttributeAsList(VCFConstants.ALLELE_COUNT_KEY).stream()
                    .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
            Assert.assertEquals(observedACWithoutRef.length, expectedAC.length - 1);
            for (int i = 0; i < observedACWithoutRef.length; i++) {
                Assert.assertEquals(observedACWithoutRef[i], expectedAC[i+1]);
            }
            Assert.assertEquals(variant.getAttributeAsInt(GenotypeCopyNumberTriStateSegments.NUMBER_OF_TARGETS_KEY, -1),
                    CopyNumberTriStateSegmentsCallerIntegrationTest.REALISTIC_TARGETS.targetCount(variant));
        }
    }

    private void assertVariantsPlGtAndGQAreConsistent(final List<VariantContext> variants) {
        for (final VariantContext vc :variants) {
            for (final Genotype gt : vc.getGenotypes()) {
                final int[] PL = gt.getPL();
                Assert.assertNotNull(PL);
                final int minPLIndex = MathUtils.minElementIndex(PL);
                Assert.assertEquals(vc.getAlleles().indexOf(gt.getAlleles().get(0)), minPLIndex);
                final int secondPLIndex = IntStream.range(0, PL.length)
                        .filter(i -> i != minPLIndex)
                        .boxed()
                        .sorted((a, b) -> Integer.compare(PL[a], PL[b]))
                        .findFirst().get();
                final int expectedGQ = Math.min(GenotypeCopyNumberTriStateSegments.MAX_GQ, PL[secondPLIndex] - PL[minPLIndex]);
                Assert.assertEquals(gt.getGQ(), expectedGQ);
            }
        }
    }

    private void assertVariantSegmentsAreCovered(final List<VariantContext> variants, final List<CopyNumberTriStateSegmentRecord> variantSegments) {
        for (final CopyNumberTriStateSegmentRecord variantSegment : variantSegments) {
            final Optional<VariantContext> match = variants.stream().filter(vc -> new SimpleInterval(vc).equals(variantSegment.getSegment().getInterval()))
                    .findFirst();
            Assert.assertTrue(match.isPresent());
            final VariantContext matchedVariant = match.get();
            final Genotype genotype = matchedVariant.getGenotype(variantSegment.getSampleName());
            final String discovery = genotype.getAnyAttribute(GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY).toString();
            Assert.assertTrue(discovery.equals(GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE));
            final CopyNumberTriState call = variantSegment.getSegment().getCall();

            final List<Allele> gt = genotype.getAlleles();
            Assert.assertEquals(gt.size(), 1);
            // The call may not be the same for case where the event-quality is relatively low.
            if (variantSegment.getSegment().getEventQuality() > 10) {
                Assert.assertEquals(GenotypeCopyNumberTriStateSegments.CNVAllele.valueOf(gt.get(0)).state, call, genotype.toString());
            }
            final String[] SQ = genotype.getAnyAttribute(GenotypeCopyNumberTriStateSegments.SOME_QUALITY_KEY).toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
            final double someQual = variantSegment.getSegment().getSomeQuality();
            Assert.assertEquals(Double.parseDouble(SQ[call == CopyNumberTriState.DELETION ? 0 : 1]), someQual, GenotypeCopyNumberTriStateSegments.PHRED_SCORE_PRECISION, variantSegment.getSampleName() + " => " + genotype.toString());

            final String[] LQ = genotype.getAnyAttribute(GenotypeCopyNumberTriStateSegments.START_QUALITY_KEY).toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
            final double startQuality = variantSegment.getSegment().getStartQuality();
            Assert.assertEquals(Double.parseDouble(LQ[call == CopyNumberTriState.DELETION ? 0 : 1]), startQuality, GenotypeCopyNumberTriStateSegments.PHRED_SCORE_PRECISION, variantSegment.getSampleName() + " => " + genotype.toString());

            final String[] RQ = genotype.getAnyAttribute(GenotypeCopyNumberTriStateSegments.END_QUALITY_KEY).toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
            final double endQuality = variantSegment.getSegment().getEndQuality();
            Assert.assertEquals(Double.parseDouble(RQ[call == CopyNumberTriState.DELETION ? 0 : 1]), endQuality, GenotypeCopyNumberTriStateSegments.PHRED_SCORE_PRECISION, variantSegment.getSampleName() + " => " + genotype.toString());

            // Check the PL.
            final int[] PL = genotype.getPL();
            final int observedGQFromPL = Math.min(GenotypeCopyNumberTriStateSegments.MAX_GQ, PL[GenotypeCopyNumberTriStateSegments.CNVAllele.REF.ordinal()] - PL[GenotypeCopyNumberTriStateSegments.CNVAllele.valueOf(call).ordinal()]);
            final double expectedCallPL = CopyNumberTriStateSegmentCaller.roundPhred(QualityUtils.phredScaleErrorRate(QualityUtils.qualToProb(variantSegment.getSegment().getExactQuality())));
            final double expectedRefPL = CopyNumberTriStateSegmentCaller.roundPhred(QualityUtils.phredScaleCorrectRate(QualityUtils.qualToProb(variantSegment.getSegment().getEventQuality())));
            final int expectedGQFromPL = Math.min(GenotypeCopyNumberTriStateSegments.MAX_GQ, (int) Math.round(expectedRefPL - expectedCallPL));
            Assert.assertTrue(Math.abs(observedGQFromPL - expectedGQFromPL) <= 1, genotype.toString() + " " + variantSegment.getSegment().getEventQuality());
        }
    }

    private void assertVariantsOrder(List<VariantContext> variants) {
        // check variant order.
        final Map<String, Integer> contigs = new HashMap<>();
        for (int i = 1; i < variants.size(); i++) {
            final VariantContext prev = variants.get(i - 1);
            final VariantContext next = variants.get(i);
            // add the contigs to the map; the value will be the contigs index in order of appearance.
            contigs.putIfAbsent(prev.getContig(), contigs.size());
            contigs.putIfAbsent(next.getContig(), contigs.size());
            if (prev.getContig().equals(next.getContig())) {
                Assert.assertFalse(prev.getStart() > next.getStart(), "variant context in same contig out of POS order");
                Assert.assertTrue(prev.getStart() < next.getStart() ||
                 prev.getEnd() != next.getEnd(), "consecutive variant context have exactly the same coordinates");
            } else {
                Assert.assertTrue(contigs.get(prev.getContig()) < contigs.get(next.getContig()), "contigs are mixed in output");
            }
        }
    }

    private void assertHeader(final VCFFileReader outputReader, final List<CopyNumberTriStateSegmentRecord> segments) {
        final VCFHeader header = outputReader.getFileHeader();
        // Check the sample names
        Assert.assertEquals(new HashSet<>(header.getSampleNamesInOrder()),
                new HashSet<>(segments.stream().map(CopyNumberTriStateSegmentRecord::getSampleName).collect(Collectors.toSet())));
        assertFormatHeaderLine(header, GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY, VCFHeaderLineCount.INTEGER, 1, VCFHeaderLineType.Character);
        assertFormatHeaderLine(header, VCFConstants.GENOTYPE_PL_KEY, VCFHeaderLineCount.G, -1, VCFHeaderLineType.Integer);
        assertFormatHeaderLine(header, VCFConstants.GENOTYPE_KEY, VCFHeaderLineCount.INTEGER, 1, VCFHeaderLineType.String);
        assertFormatHeaderLine(header, VCFConstants.GENOTYPE_QUALITY_KEY, VCFHeaderLineCount.INTEGER, 1, VCFHeaderLineType.Integer);
        assertFormatHeaderLine(header, GenotypeCopyNumberTriStateSegments.SOME_QUALITY_KEY, VCFHeaderLineCount.A, -1, VCFHeaderLineType.Float);
        assertFormatHeaderLine(header, GenotypeCopyNumberTriStateSegments.START_QUALITY_KEY, VCFHeaderLineCount.A, -1, VCFHeaderLineType.Float);
        assertFormatHeaderLine(header, GenotypeCopyNumberTriStateSegments.END_QUALITY_KEY, VCFHeaderLineCount.A, -1, VCFHeaderLineType.Float);
        assertInfoHeaderLine(header, GenotypeCopyNumberTriStateSegments.NUMBER_OF_TARGETS_KEY, VCFHeaderLineCount.INTEGER, 1, VCFHeaderLineType.Integer);
        assertInfoHeaderLine(header, VCFConstants.END_KEY, VCFHeaderLineCount.INTEGER, 1, VCFHeaderLineType.Integer);
        assertInfoHeaderLine(header, VCFConstants.ALLELE_COUNT_KEY, VCFHeaderLineCount.A, -1, VCFHeaderLineType.Integer);
        assertInfoHeaderLine(header, VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, -1, VCFHeaderLineType.Float);
        assertInfoHeaderLine(header, VCFConstants.ALLELE_NUMBER_KEY, VCFHeaderLineCount.INTEGER, 1, VCFHeaderLineType.Integer);
    }

    private void assertInfoHeaderLine(final VCFHeader header, final String id, final VCFHeaderLineCount countType, final int fixedCount, final VCFHeaderLineType lineType) {
        final VCFInfoHeaderLine line = header.getInfoHeaderLine(id);
        Assert.assertNotNull(line);
        Assert.assertEquals(line.getCountType(), countType);
        if (countType == VCFHeaderLineCount.INTEGER) {
            Assert.assertEquals(line.getCount(), fixedCount);
        }
        Assert.assertEquals(line.getType(), lineType);
    }

    private void assertFormatHeaderLine(final VCFHeader header, final String id, final VCFHeaderLineCount count, final int number, final VCFHeaderLineType lineType) {
        final VCFFormatHeaderLine line = header.getFormatHeaderLine(id);
        Assert.assertNotNull(line);
        Assert.assertEquals(line.getCountType(), count);
        if (count == VCFHeaderLineCount.INTEGER) {
            Assert.assertEquals(line.getCount(), number);
        }
        Assert.assertEquals(line.getType(), lineType);
    }

    private List<CopyNumberTriStateSegmentRecord> readSegmentRecords(final File segmentsFile) throws IOException {
        try (final CopyNumberTriStateSegmentRecordReader reader = new CopyNumberTriStateSegmentRecordReader(segmentsFile)) {
            return reader.toList();
        }
    }

    @DataProvider(name = "simulateChainData")
    public static Object[][] simulateChainDataProvider() {
        return CopyNumberTriStateSegmentsCallerIntegrationTest.simulateChainDataProvider();
    }
}
