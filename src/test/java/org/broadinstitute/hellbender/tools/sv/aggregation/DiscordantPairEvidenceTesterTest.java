package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class DiscordantPairEvidenceTesterTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;

    private final SVCallRecord TEST_DEL_RECORD = new SVCallRecord("call1", "chr21", 1000, true, "chr21", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
            null, Collections.emptyList(), null,  Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);

    @DataProvider(name = "poissonTestTestData")
    public Object[][] poissonTestTestData() {
        return new Object[][]{
                // Empty site
                {
                        Collections.emptyList(),
                        Collections.emptyList(),
                        Collections.emptyList(),
                        Collections.emptyMap(),
                        new DiscordantPairEvidenceTester.DiscordantPairTestResult(
                                new EvidenceStatUtils.PoissonTestResult(1., 0, 0.),
                                Collections.emptyMap(),
                                Collections.emptyList()
                        )
                },
                // 0 carrier / 1 background sample
                {
                        Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false)),
                        Collections.emptyList(),
                        Collections.singletonList("sample1"),
                        Collections.singletonMap("sample1", 22.),
                        new DiscordantPairEvidenceTester.DiscordantPairTestResult(
                                new EvidenceStatUtils.PoissonTestResult(1, 0, 3.),
                                Collections.singletonMap("sample1", 1),
                                Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false))
                        )
                },
                // 1 carrier / 0 background sample
                {
                        Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false)),
                        Collections.singletonList("sample1"),
                        Collections.emptyList(),
                        Collections.singletonMap("sample1", 22.),
                        new DiscordantPairEvidenceTester.DiscordantPairTestResult(
                                new EvidenceStatUtils.PoissonTestResult(0.04978706836786395, 3.0, 0),
                                Collections.singletonMap("sample1", 1),
                                Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false))
                        )
                },
                // 1 carrier / 1 background sample
                {
                        Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        new DiscordantPairEvidenceTester.DiscordantPairTestResult(
                                new EvidenceStatUtils.PoissonTestResult(0.04978706836786395, 3., 0),
                                Collections.singletonMap("sample1", 1),
                                Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false))
                        )
                },

                // 1 carrier / 2 background samples with evidence
                {
                        Lists.newArrayList(
                                new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false),
                                new DiscordantPairEvidence("sample1", "chr21", 1005, true, "chr21", 2005, false),
                                new DiscordantPairEvidence("sample2", "chr21", 1000, true, "chr21", 2000, false)
                        ),
                        Collections.singletonList("sample1"),
                        Lists.newArrayList("sample2", "sample3"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.), new AbstractMap.SimpleEntry<>("sample3", 21.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        new DiscordantPairEvidenceTester.DiscordantPairTestResult(
                                new EvidenceStatUtils.PoissonTestResult(0.1246520194830811, 5.0, 1.5),
                                Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 2), new AbstractMap.SimpleEntry<>("sample2", 1)).stream()
                                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                                Lists.newArrayList(
                                        new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false),
                                        new DiscordantPairEvidence("sample1", "chr21", 1005, true, "chr21", 2005, false),
                                        new DiscordantPairEvidence("sample2", "chr21", 1000, true, "chr21", 2000, false)
                                )
                        )
                }
        };
    }

    @Test(dataProvider= "poissonTestTestData")
    public void poissonTestTest(final List<DiscordantPairEvidence> sortedEvidence,
                                        final Collection<String> carrierSamples,
                                        final Collection<String> backgroundSamples,
                                        final Map<String, Double> sampleCoverageMap,
                                        final DiscordantPairEvidenceTester.DiscordantPairTestResult expected) {
        final DiscordantPairEvidenceTester tester = new DiscordantPairEvidenceTester(sampleCoverageMap, DICTIONARY);
        final DiscordantPairEvidenceTester.DiscordantPairTestResult test =
                tester.poissonTest(sortedEvidence, carrierSamples, backgroundSamples);
        Assert.assertNotNull(test);
        Assert.assertEquals(test.getSampleCounts(), expected.getSampleCounts());
        Assert.assertEquals(test.getDiscordantPairEvidence(), expected.getDiscordantPairEvidence());

        final EvidenceStatUtils.PoissonTestResult testPoisson = test.getTest();
        final EvidenceStatUtils.PoissonTestResult expectedPoisson = expected.getTest();
        SVTestUtils.assertFloatWithinTolerance(testPoisson.getP(), expectedPoisson.getP());
        SVTestUtils.assertFloatWithinTolerance(testPoisson.getCarrierSignal(), expectedPoisson.getCarrierSignal());
        SVTestUtils.assertFloatWithinTolerance(testPoisson.getBackgroundSignal(), expectedPoisson.getBackgroundSignal());
    }


    @DataProvider(name = "applyToRecordTest")
    public Object[][] applyToRecordTest() {
        return new Object[][]{
                // No evidence
                {
                        TEST_DEL_RECORD,
                        Collections.emptySet(),
                        Collections.emptySet(),
                        Collections.emptySet(),
                        Collections.emptyList(),
                        Collections.emptyMap(),
                        0.,
                        0.
                },
                // Single carrier with evidence
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false)),
                        Collections.singletonMap("sample1", 1),
                        100.,
                        8.685889638065037
                },
                // Single carrier but excluded
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Sets.newHashSet("sample1"),
                        Collections.singletonList(new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false)),
                        Collections.singletonMap("sample1", 1),
                        0.,
                        0.
                },
                // Background evidence
                {
                        TEST_DEL_RECORD,
                        Sets.newHashSet("sample1"),
                        Sets.newHashSet("sample2"),
                        Collections.emptySet(),
                        Lists.newArrayList(
                                new DiscordantPairEvidence("sample1", "chr21", 1000, true, "chr21", 2000, false),
                                new DiscordantPairEvidence("sample1", "chr21", 1001, true, "chr21", 2001, false),
                                new DiscordantPairEvidence("sample2", "chr21", 1000, true, "chr21", 2000, false)
                        ),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 2), new AbstractMap.SimpleEntry<>("sample2", 1)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        60.,
                        3.734645199954624
                },
        };
    }

    @Test(dataProvider= "applyToRecordTest")
    public void applyToRecordTest(final SVCallRecord baseRecord,
                                  final Set<String> carrierSamples,
                                  final Set<String> backgroundSamples,
                                  final Set<String> excludedSamples,
                                  final List<DiscordantPairEvidence> evidence,
                                  final Map<String, Integer> expectedSampleCounts,
                                  final Double expectedCarrierSignal,
                                  final Double expectedQuality) {
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

        final DiscordantPairEvidenceTester tester = new DiscordantPairEvidenceTester(sampleCoverageMap, DICTIONARY);
        final DiscordantPairEvidenceTester.DiscordantPairTestResult result = tester.test(record, evidence,
                Sets.difference(carrierSamples, excludedSamples), Sets.difference(backgroundSamples, excludedSamples));
        final SVCallRecord test = tester.applyToRecord(record, result);
        Assert.assertEquals(test.getId(), record.getId());
        Assert.assertEquals(test.getContigA(), record.getContigA());
        Assert.assertEquals(test.getContigB(), record.getContigB());
        Assert.assertEquals(test.getPositionA(), record.getPositionA());
        Assert.assertEquals(test.getPositionB(), record.getPositionB());
        Assert.assertEquals(test.getStrandA(), record.getStrandA());
        Assert.assertEquals(test.getStrandB(), record.getStrandB());
        Assert.assertEquals(test.getType(), record.getType());
        Assert.assertEquals(test.getAlgorithms(), record.getAlgorithms());
        Assert.assertEquals(test.getAttributes().get("TEST_KEY"), "TEST_VALUE");
        Assert.assertTrue(test.getGenotypes().containsSamples(carrierSamples));
        Assert.assertTrue(test.getGenotypes().containsSamples(backgroundSamples));
        SVTestUtils.assertFloatWithinTolerance((Double) test.getAttributes().get(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE), expectedCarrierSignal);
        SVTestUtils.assertFloatWithinTolerance((Double) test.getAttributes().get(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE), expectedQuality);
        for (final String s : expectedSampleCounts.keySet()) {
            final Integer count = VariantContextGetters.getAttributeAsInt(test.getGenotypes().get(s),
                    GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, -1);
            Assert.assertEquals(count, expectedSampleCounts.get(s));
        }
    }

}