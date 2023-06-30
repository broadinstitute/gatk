package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class SplitReadEvidenceTesterTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final double ERROR_TOL = 1e-6;
    private static final int MAX_SR_CROSS_DISTANCE = 20;

    private final SVCallRecord TEST_DEL_RECORD = new SVCallRecord("call1", "chr21", 1000, true, "chr21", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
            null, null, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);

    private final SVCallRecord TEST_INS_RECORD = new SVCallRecord("call1", "chr21", 1000, true, "chr21", 1000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
            null, 500, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);

    private final SVCallRecord TEST_BND_RECORD = new SVCallRecord("call1", "chr21", 1000, true, "chr22", 8000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
            null, null, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.create("<BND>", false)), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);

    private final SVCallRecord TEST_INV_RECORD = new SVCallRecord("call1", "chr21", 1000, true, "chr21", 8000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
            null, null, Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INV), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);

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
                        true,
                        new SplitReadSite("chr21",
                                1,
                                true,
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
                        true,
                        new SplitReadSite("chr21",
                                1,
                                true,
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
                        true,
                        new SplitReadSite("chr21",
                                10,
                                true,
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
                        true,
                        new SplitReadSite("chr21",
                                10,
                                true,
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
                        true,
                        new SplitReadSite("chr21",
                                10,
                                true,
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
                        true,
                        new SplitReadSite("chr21",
                                10,
                                true,
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
                                        final boolean strand,
                                        final SplitReadSite expected) {
        final SplitReadSite test = SplitReadEvidenceTester.testSplitReadSite(sortedEvidence, strand, carrierSamples,
                backgroundSamples, sampleCoverageMap, representativeDepth, "chr21", defaultPosition);
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
                        TEST_DEL_RECORD,
                        Collections.emptySet(),
                        Collections.emptySet(),
                        Collections.emptySet(),
                        Collections.emptyList(),
                        Collections.emptyList(),
                        1000,
                        2000,
                        1000,
                        2000,
                        Collections.emptyMap(),
                        Collections.emptyMap()
                },
                // Single carrier with evidence at new coordinates
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 1005, 5, true)),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 2005, 4, false)),
                        1000,
                        2000,
                        1005,
                        2005,
                        Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 5), new HashMap.SimpleEntry<>("sample2", 0)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 4), new HashMap.SimpleEntry<>("sample2", 0)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                },
                // Single carrier with evidence at new coordinates, but excluded
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.singleton("sample1"),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 1005, 5, true)),
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 2005, 5, false)),
                        1000,
                        2000,
                        1000,
                        2000,
                        Collections.emptyMap(),
                        Collections.emptyMap()
                },
                // Single carrier and background sample with evidence at new coordinates
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        1005,
                        2005,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                // Single carrier and background sample with invalid end evidence upstream of refined start
                // Start evidence better
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1500, 6, true),
                                new SplitReadEvidence("sample2", "chr21", 1500, 1, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1400, 5, false),
                                new SplitReadEvidence("sample2", "chr21", 1400, 1, false)
                        ),
                        1000,
                        2000,
                        1500,
                        2000,
                        Collections.singletonMap("sample1", 6),
                        Collections.emptyMap()
                },
                // Single carrier and background sample with invalid end evidence upstream of refined start
                // End evidence better
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1500, 5, true),
                                new SplitReadEvidence("sample2", "chr21", 1500, 1, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1400, 6, false),
                                new SplitReadEvidence("sample2", "chr21", 1400, 1, false)
                        ),
                        1000,
                        2000,
                        1000,
                        1400,
                        Collections.emptyMap(),
                        Collections.singletonMap("sample1", 6)
                },
                // Single carrier and background sample with valid end evidence upstream of original start but not refined start
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        900,
                        950,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                // Single carrier and background sample with multiple evidence sites
                //   --Refined locus 3'
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        1010,
                        2010,
                        Collections.singletonMap("sample1", 6),
                        Collections.singletonMap("sample1", 6)
                },
                //   --Refined locus 5'
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        995,
                        1995,
                        Collections.singletonMap("sample1", 6),
                        Collections.singletonMap("sample1", 6)
                },
                //   --3' evidence weaker
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        1005,
                        2005,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                //   --Equal evidence, refine closest to original coordinates
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        1001,
                        2001,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                //   --Ignore excluded sample
                {
                        TEST_DEL_RECORD,
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
                        1000,
                        2000,
                        1005,
                        2005,
                        Collections.singletonMap("sample1", 5),
                        Collections.singletonMap("sample1", 5)
                },
                // INS simple case
                {
                        TEST_INS_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 995, 12, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 1001, 10, false)
                        ),
                        1000,
                        1000,
                        995,
                        1001,
                        Collections.singletonMap("sample1", 12),
                        Collections.singletonMap("sample1", 10)
                },
                // INS
                // First end has most evidence, but is beyond the max crossover distance
                // Second end has second-most evidence and is in valid range, so should be chosen
                // Third end has least evidence
                {
                        TEST_INS_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 995, 12, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 995 - MAX_SR_CROSS_DISTANCE - 1, 10, false),
                                new SplitReadEvidence("sample1", "chr21", 995 - MAX_SR_CROSS_DISTANCE, 9, false),
                                new SplitReadEvidence("sample1", "chr21", 2010, 5, false)
                        ),
                        1000,
                        1000,
                        995,
                        995 - MAX_SR_CROSS_DISTANCE,
                        Collections.singletonMap("sample1", 12),
                        Collections.singletonMap("sample1", 9)
                },
                // INS, reverse crossover case
                {
                        TEST_INS_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 995 + MAX_SR_CROSS_DISTANCE + 1, 10, false),
                                new SplitReadEvidence("sample1", "chr21", 995 + MAX_SR_CROSS_DISTANCE, 9, false),
                                new SplitReadEvidence("sample1", "chr21", 2010, 5, false)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 995, 12, true)
                        ),
                        1000,
                        1000,
                        995 + MAX_SR_CROSS_DISTANCE,
                        995,
                        Collections.singletonMap("sample1", 9),
                        Collections.singletonMap("sample1", 12)
                },
                // BND, interchromosomal
                {
                        TEST_BND_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 900, 4, false)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr22", 1, 10, false),
                                new SplitReadEvidence("sample1", "chr22", 8000, 9, false)
                        ),
                        1000,
                        8000,
                        900,
                        1,
                        Collections.singletonMap("sample1", 4),
                        Collections.singletonMap("sample1", 10)
                },
                // INV, with same evidence at POS=4000
                {
                        TEST_INV_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 4000, 4, true)
                        ),
                        Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr21", 4000, 5, true),
                                new SplitReadEvidence("sample1", "chr21", 2000, 3, true)
                        ),
                        1000,
                        8000,
                        4000,
                        2000,
                        Collections.singletonMap("sample1", 4),
                        Collections.singletonMap("sample1", 3)
                },
        };
    }

    @Test(dataProvider= "refineCallTestData")
    public void refineCallTest(final SVCallRecord baseRecord,
                               final Set<String> carrierSamples,
                               final Set<String> backgroundSamples,
                               final Set<String> excludedSamples,
                               final List<SplitReadEvidence> startEvidence,
                               final List<SplitReadEvidence> endEvidence,
                               final int expectedPositionA,
                               final int expectedPositionB,
                               final int expectedFirstSRPosition,
                               final int expectedSecondSRPosition,
                               final Map<String, Integer> expectedStartSampleCounts,
                               final Map<String, Integer> expectedEndSampleCounts) {
        final GenotypesContext genotypes = GenotypesContext.create();
        carrierSamples.stream()
                .forEach(s -> genotypes.add(new GenotypeBuilder(s).alleles(Lists.newArrayList(baseRecord.getRefAllele(), baseRecord.getAltAlleles().get(0)))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make()));
        backgroundSamples.stream()
                .forEach(s -> genotypes.add(new GenotypeBuilder(s).alleles(Lists.newArrayList(baseRecord.getRefAllele(), baseRecord.getRefAllele()))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make()));

        final Map<String, Object> attributes = Collections.singletonMap("TEST_KEY", "TEST_VALUE");

        final SVCallRecord record = SVCallRecordUtils.copyCallWithNewGenotypes(SVCallRecordUtils.copyCallWithNewAttributes(baseRecord, attributes), genotypes);

        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        sampleCoverageMap.put("sample1", 35.);
        sampleCoverageMap.put("sample2", 25.);

        final SplitReadEvidenceTester refiner = new SplitReadEvidenceTester(sampleCoverageMap, MAX_SR_CROSS_DISTANCE, DICTIONARY);
        final SplitReadEvidenceTester.SplitReadTestResult result = refiner.test(record, startEvidence, endEvidence,
                Sets.difference(carrierSamples, excludedSamples), Sets.difference(backgroundSamples, excludedSamples));
        final SVCallRecord test = refiner.applyToRecord(record, result);
        Assert.assertEquals(test.getId(), record.getId());
        Assert.assertEquals(test.getContigA(), record.getContigA());
        Assert.assertEquals(test.getContigB(), record.getContigB());
        Assert.assertEquals(test.getPositionA(), expectedPositionA);
        Assert.assertEquals(test.getPositionB(), expectedPositionB);
        Assert.assertEquals(test.getAttributes().get(GATKSVVCFConstants.FIRST_SPLIT_POSITION_ATTRIBUTE), expectedFirstSRPosition);
        Assert.assertEquals(test.getAttributes().get(GATKSVVCFConstants.SECOND_SPLIT_POSITION_ATTRIBUTE), expectedSecondSRPosition);
        Assert.assertEquals(test.getStrandA(), record.getStrandA());
        Assert.assertEquals(test.getStrandB(), record.getStrandB());
        Assert.assertEquals(test.getType(), record.getType());
        Assert.assertEquals(test.getAlgorithms(), record.getAlgorithms());
        Assert.assertEquals(test.getAttributes().get("TEST_KEY"), "TEST_VALUE");
        Assert.assertTrue(test.getGenotypes().containsSamples(carrierSamples));
        Assert.assertTrue(test.getGenotypes().containsSamples(backgroundSamples));
        for (final String s : expectedStartSampleCounts.keySet()) {
            final Integer count = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get(s),
                    GATKSVVCFConstants.FIRST_SPLIT_READ_COUNT_ATTRIBUTE, -1);
            Assert.assertEquals(count, expectedStartSampleCounts.get(s));
        }
        for (final String s : expectedEndSampleCounts.keySet()) {
            final Integer count = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get(s),
                    GATKSVVCFConstants.SECOND_SPLIT_READ_COUNT_ATTRIBUTE, -1);
            Assert.assertEquals(count, expectedEndSampleCounts.get(s));
        }
    }
}