package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class BreakpointRefinerTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final double ERROR_TOL = 1e-6;

    @DataProvider(name = "refineSplitReadSiteTestData")
    public Object[][] refineSplitReadSiteTestData() {
        return new Object[][]{
                // Empty site
                {
                    Collections.emptyList(),
                        Collections.emptyList(),
                        Collections.emptyList(),
                        Collections.emptyMap(),
                        30,
                        1,
                        new SplitReadSite(1,
                                Collections.emptyMap(),
                                null)
                },
                // 0 carrier / 1 background sample
                {
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 10, 1, true)),
                        Collections.emptyList(),
                        Collections.singletonList("sample1"),
                        Collections.singletonMap("sample1", 22.),
                        30,
                        1,
                        new SplitReadSite(1,
                                Collections.emptyMap(),
                                null)
                },
                // 1 carrier / 0 background sample
                {
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.emptyList(),
                        Collections.singletonList(new HashMap.SimpleEntry<>("sample1", 22.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Collections.singletonMap("sample1", 1),
                                new EvidenceStatUtils.PoissonTestResult(0.2557291599131005, 1.3636363636363638, 0.0))
                },
                // 1 carrier / 1 background sample
                {
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Collections.singletonMap("sample1", 1),
                                new EvidenceStatUtils.PoissonTestResult(0.2557291599131005, 1.3636363636363638, 0.0))
                },
                // 1 carrier / 1 background sample with evidence
                {
                        Lists.newArrayList(new SplitReadEvidence("sample1", "chr21", 10, 1, true),
                                new SplitReadEvidence("sample2", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 1), new HashMap.SimpleEntry<>("sample2", 1)).stream()
                                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                                new EvidenceStatUtils.PoissonTestResult(0.8422154564080215, 1.3636363636363638, 1.5))
                },
                // 1 carrier / 1 background; multiple loci with evidence
                {
                        Lists.newArrayList(new SplitReadEvidence("sample1", "chr21", 5, 1, true),
                                new SplitReadEvidence("sample2", "chr21", 5, 1, true),
                                new SplitReadEvidence("sample1", "chr21", 10, 6, true),
                                new SplitReadEvidence("sample2", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 6), new HashMap.SimpleEntry<>("sample2", 1)).stream()
                                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                                new EvidenceStatUtils.PoissonTestResult(0.011929713129970071, 8.181818181818182, 1.5))
                },
        };
    }

    @Test(dataProvider= "refineSplitReadSiteTestData")
    public void refineSplitReadSiteTest(final List<SplitReadEvidence> sortedEvidence,
                                        final Collection<String> carrierSamples,
                                        final Collection<String> backgroundSamples,
                                        final Map<String, Double> sampleCoverageMap,
                                        final int representativeDepth,
                                        final int defaultPosition,
                                        final SplitReadSite expected) {
        final SplitReadSite test = BreakpointRefiner.refineSplitReadSite(sortedEvidence, carrierSamples,
                backgroundSamples, sampleCoverageMap, representativeDepth, defaultPosition);
        Assert.assertEquals(test.getPosition(), expected.getPosition());
        final Set<String> samples = new HashSet<>(carrierSamples);
        samples.addAll(backgroundSamples);
        for (final String s : samples) {
            Assert.assertEquals(test.getCount(s), expected.getCount(s));
        }
        Assert.assertTrue((test.getP() == null && expected.getP() == null)
                || Math.abs(test.getP() - expected.getP()) <= ERROR_TOL);
        Assert.assertTrue((test.getCarrierSignal() == null && expected.getCarrierSignal() == null)
                || Math.abs(test.getCarrierSignal() - expected.getCarrierSignal()) <= ERROR_TOL);
        Assert.assertTrue((test.getBackgroundSignal() == null && expected.getBackgroundSignal() == null)
                || Math.abs(test.getBackgroundSignal() - expected.getBackgroundSignal()) <= ERROR_TOL);
    }

    @DataProvider(name = "refineCallTestData")
    public Object[][] refineCallTestData() {
        return new Object[][]{
                // No evidence
                {
                        Collections.emptySet(),
                        Collections.emptySet(),
                        Collections.emptySet(),
                        Collections.emptyList(),
                        Collections.emptyList(),
                        1000,
                        2000,
                        Collections.emptyMap(),
                        Collections.emptyMap()
                },
                // Single carrier with evidence at new coordinates
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 1005, 5, true)),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 2005, 4, false)),
                        1005,
                        2005,
                        Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 5), new HashMap.SimpleEntry<>("sample2", 0)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 4), new HashMap.SimpleEntry<>("sample2", 0)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                },
                // Single carrier with evidence at new coordinates, but excluded
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.singleton("sample1"),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 1005, 5, true)),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 2005, 5, false)),
                        1000,
                        2000,
                        Collections.emptyMap(),
                        Collections.emptyMap()
                },
                // Single carrier and background sample with evidence at new coordinates
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1005, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1005, 1, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 2005, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 2005, 1, false)
                        ),
                        1005,
                        2005,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                // Single carrier and background sample with invalid end evidence upstream of refined start
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1500, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1500, 1, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1400, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 1400, 1, false)
                        ),
                        1500,
                        2000,
                        Collections.singletonMap("sample1", 5),
                        Collections.emptyMap()
                },
                // Single carrier and background sample with valid end evidence upstream of original start but not refined start
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 900, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 900, 1, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 950, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 950, 1, false)
                        ),
                        900,
                        950,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                // Single carrier and background sample with multiple evidence sites
                //   --Refined locus 3'
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1005, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1005, 1, true),
                                new SplitReadEvidence("sample1", "chr21", 1010, 6, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 2005, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 2005, 1, false),
                                new SplitReadEvidence("sample1", "chr21", 2010, 6, false)
                        ),
                        1010,
                        2010,
                        Collections.singletonMap("sample1", 6),
                        Collections.singletonMap("sample1", 6)
                },
                //   --Refined locus 5'
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 995, 6, true),
                                new SplitReadEvidence("sample1", "chr21", 1005, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1005, 1, true)
                        ),
                        Lists.newArrayList(

                                new SplitReadEvidence("sample1", "chr21", 1995, 6, false),
                                new SplitReadEvidence("sample1", "chr21", 2005, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 2005, 1, false)
                        ),
                        995,
                        1995,
                        Collections.singletonMap("sample1", 6),
                        Collections.singletonMap("sample1", 6)
                },
                //   --3' evidence weaker
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1005, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1005, 1, true),
                                new SplitReadEvidence("sample1", "chr21", 1010, 1, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 2005, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 2005, 1, false),
                                new SplitReadEvidence("sample1", "chr21", 2010, 1, false)
                        ),
                        1005,
                        2005,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                //   --Equal evidence, refine closest to original coordinates
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 990, 5, true),
                                new SplitReadEvidence("sample1", "chr21", 992, 5, true),
                                new SplitReadEvidence("sample1", "chr21", 995, 5, true),
                                new SplitReadEvidence("sample1", "chr21", 1001, 5, true),
                                new SplitReadEvidence("sample1", "chr21", 1010, 5, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 990, 5, false),
                                new SplitReadEvidence("sample1", "chr21", 992, 5, false),
                                new SplitReadEvidence("sample1", "chr21", 995, 5, false),
                                new SplitReadEvidence("sample1", "chr21", 2001, 5, false),
                                new SplitReadEvidence("sample1", "chr21", 2010, 5, false)
                        ),
                        1001,
                        2001,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                //   --Ignore excluded sample
                {
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.singleton("sample2"),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1005, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1005, 5, true),
                                new SplitReadEvidence("sample1", "chr21", 1010, 4, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 2005, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 2005, 5, false),
                                new SplitReadEvidence("sample1", "chr21", 2010, 4, false)
                        ),
                        1005,
                        2005,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
        };
    }

    @Test(dataProvider= "refineCallTestData")
    public void refineCallTest(final Set<String> carrierSamples,
                               final Set<String> backgroundSamples,
                               final Set<String> excludedSamples,
                               final List<SplitReadEvidence> startEvidence,
                               final List<SplitReadEvidence> endEvidence,
                               final int expectedPositionA,
                               final int expectedPositionB,
                               final Map<String, Integer> expectedStartSampleCounts,
                               final Map<String, Integer> expectedEndSampleCounts) {
        final GenotypesContext genotypes = GenotypesContext.create();
        carrierSamples.stream()
                .forEach(s -> genotypes.add(new GenotypeBuilder(s).alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1).make()));
        backgroundSamples.stream()
                .forEach(s -> genotypes.add(new GenotypeBuilder(s).alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make()));
        final Map<String, Object> attributes = Collections.singletonMap("TEST_KEY", "TEST_VALUE");
        final SVCallRecord record = new SVCallRecord("call1", "chr21", 1000, true, "chr21", 2000, false, StructuralVariantType.DEL,
                null, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), genotypes, attributes, DICTIONARY);
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        sampleCoverageMap.put("sample1", 35.);
        sampleCoverageMap.put("sample2", 25.);

        final BreakpointRefiner refiner = new BreakpointRefiner(sampleCoverageMap, 20, DICTIONARY);
        final BreakpointRefiner.RefineResult result = refiner.testRecord(record, startEvidence, endEvidence,
                Sets.difference(carrierSamples, excludedSamples), Sets.difference(backgroundSamples, excludedSamples), null);
        final SVCallRecord test = refiner.applyToRecord(record, result);
        Assert.assertEquals(test.getId(), "call1");
        Assert.assertEquals(test.getPositionA(), 1000);
        Assert.assertEquals(test.getPositionB(), 2000);
        Assert.assertEquals(test.getAttributes().get(GATKSVVCFConstants.START_SPLIT_POSITION_ATTRIBUTE), expectedPositionA);
        Assert.assertEquals(test.getAttributes().get(GATKSVVCFConstants.END_SPLIT_POSITION_ATTRIBUTE), expectedPositionB);
        Assert.assertEquals(test.getType(), StructuralVariantType.DEL);
        Assert.assertEquals(test.getAlgorithms(), Collections.singletonList("pesr"));
        Assert.assertEquals(test.getAttributes().get("TEST_KEY"), "TEST_VALUE");
        Assert.assertTrue(test.getGenotypes().containsSamples(carrierSamples));
        Assert.assertTrue(test.getGenotypes().containsSamples(backgroundSamples));
        for (final String s : expectedStartSampleCounts.keySet()) {
            final Integer count = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get(s),
                    GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, -1);
            Assert.assertEquals(count, expectedStartSampleCounts.get(s));
        }
        for (final String s : expectedEndSampleCounts.keySet()) {
            final Integer count = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get(s),
                    GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, -1);
            Assert.assertEquals(count, expectedEndSampleCounts.get(s));
        }
    }

    @Test
    public void refineInsertionCallTest() {
        final int maxInsertionSplitReadCrossDistance = 20;
        final GenotypesContext genotypes = GenotypesContext.create();
        genotypes.add(new GenotypeBuilder("sample1").alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make());
        genotypes.add(new GenotypeBuilder("sample2").alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make());
        final SVCallRecord record = new SVCallRecord("call1", "chr21", 1000, true, "chr21", 1001, false, StructuralVariantType.INS,
                500, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS), genotypes, Collections.emptyMap(), DICTIONARY);
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        sampleCoverageMap.put("sample1", 35.);
        sampleCoverageMap.put("sample2", 25.);

        final List<SplitReadEvidence> startEvidence = Lists.newArrayList(
                new SplitReadEvidence("sample1", "chr21", 995, 5, true)
        );
        // First case has most evidence, but is beyond the max crossover distance
        // Second case has second-most evidence and is in valid range
        // Third case has least evidence
        final List<SplitReadEvidence> endEvidence = Lists.newArrayList(
                new SplitReadEvidence("sample1", "chr21", 995 - maxInsertionSplitReadCrossDistance - 1, 10, false),
                new SplitReadEvidence("sample1", "chr21", 995 - maxInsertionSplitReadCrossDistance, 9, false),
                new SplitReadEvidence("sample1", "chr21", 2010, 5, false)
        );

        final BreakpointRefiner refiner = new BreakpointRefiner(sampleCoverageMap, maxInsertionSplitReadCrossDistance, DICTIONARY);
        final BreakpointRefiner.RefineResult result = refiner.testRecord(record, startEvidence, endEvidence,
                Collections.singleton("sample1"), Collections.emptySet(), null);
        final SVCallRecord test = refiner.applyToRecord(record, result);
        Assert.assertEquals(test.getId(), "call1");
        Assert.assertEquals(test.getPositionA(), 1000);
        Assert.assertEquals(test.getPositionB(), 1001);
        Assert.assertEquals(test.getAlgorithms(), Collections.singletonList("pesr"));
        Assert.assertEquals(test.getLength(), Integer.valueOf(500));
        Assert.assertEquals((int) test.getAttributes().get(GATKSVVCFConstants.START_SPLIT_POSITION_ATTRIBUTE), 995);
        Assert.assertEquals((int) test.getAttributes().get(GATKSVVCFConstants.END_SPLIT_POSITION_ATTRIBUTE), 995 - maxInsertionSplitReadCrossDistance);
        final Integer startCount = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get("sample1"),
                GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, -1);
        Assert.assertEquals(startCount, Integer.valueOf(5));
        final Integer endCount = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get("sample1"),
                GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, -1);
        Assert.assertEquals(endCount, Integer.valueOf(9));
    }

    @Test
    public void refineInterchromosomalCallTest() {
        final Allele bndAllele = Allele.create("<BND>", false);
        final GenotypesContext genotypes = GenotypesContext.create();
        genotypes.add(new GenotypeBuilder("sample1").alleles(Lists.newArrayList(Allele.REF_N, bndAllele))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make());
        genotypes.add(new GenotypeBuilder("sample2").alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make());
        final SVCallRecord record = new SVCallRecord("call1", "chr21", 1000, true, "chr22", 8000, false, StructuralVariantType.BND,
                null, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, bndAllele), genotypes, Collections.emptyMap(), DICTIONARY);
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        sampleCoverageMap.put("sample1", 35.);
        sampleCoverageMap.put("sample2", 25.);

        final List<SplitReadEvidence> startEvidence = Collections.emptyList();
        // No limits on end position since it's interchromosomal
        final List<SplitReadEvidence> endEvidence = Lists.newArrayList(
                new SplitReadEvidence("sample1", "chr22", 1, 10, false),
                new SplitReadEvidence("sample1", "chr22", 8000, 9, false)
        );

        final BreakpointRefiner refiner = new BreakpointRefiner(sampleCoverageMap, 20, DICTIONARY);
        final BreakpointRefiner.RefineResult result = refiner.testRecord(record, startEvidence, endEvidence,
                Collections.singleton("sample1"), Collections.emptySet(), null);
        final SVCallRecord test = refiner.applyToRecord(record, result);
        Assert.assertEquals(test.getId(), "call1");
        Assert.assertEquals(test.getContigA(), "chr21");
        Assert.assertEquals(test.getPositionA(), 1000);
        Assert.assertEquals(test.getAttributes().get(GATKSVVCFConstants.START_SPLIT_POSITION_ATTRIBUTE), 1000);
        Assert.assertEquals(test.getContigB(), "chr22");
        Assert.assertEquals(test.getPositionB(), 8000);
        Assert.assertEquals(test.getAttributes().get(GATKSVVCFConstants.END_SPLIT_POSITION_ATTRIBUTE), 1);
        Assert.assertEquals(test.getAlgorithms(), Collections.singletonList("pesr"));
        final Integer startCount = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get("sample1"),
                GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, -1);
        Assert.assertEquals(startCount, Integer.valueOf(0));
        final Integer endCount = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get("sample1"),
                GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, -1);
        Assert.assertEquals(endCount, Integer.valueOf(10));
    }

}