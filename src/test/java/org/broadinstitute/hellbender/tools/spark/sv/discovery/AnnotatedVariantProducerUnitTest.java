package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple4;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;


public class AnnotatedVariantProducerUnitTest extends BaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVDiscoveryTestDataProvider.testDataInitialized) {
            new SVDiscoveryTestDataProvider();
        }
    }


    // -----------------------------------------------------------------------------------------------
    // Evidence summary annotation
    // -----------------------------------------------------------------------------------------------
    /**
     * Not an exhaustive test on all attributes, only tests:
     * MAPPING_QUALITIES, ALIGNMENT_LENGTH
     */
    @Test(groups = "sv", dataProvider = "forEvidenceAnnotation")
    public void testEvidenceAnnotation(final Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData,
                                       final String[] expectedMappingQualitiesAsStrings,
                                       final String[] expectedAlignmentLengthsAsStrings) {

        final AlignmentInterval region1 = testData._1();
        final AlignmentInterval region2 = testData._2();
        final byte[] contigSeq = null; // hack, as the contig sequence is really not necessary for this test purpose

        final Map<String, Object> attributeMap =
                AnnotatedVariantProducer.getEvidenceRelatedAnnotations(Collections.singletonList(new ChimericAlignment(region1, region2, Collections.emptyList(), testData._4())));

        Assert.assertEquals(((String)attributeMap.get(MAPPING_QUALITIES)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedMappingQualitiesAsStrings);
        Assert.assertEquals(((String)attributeMap.get(ALIGN_LENGTHS)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedAlignmentLengthsAsStrings);
    }

    @DataProvider(name = "forEvidenceAnnotation")
    private Object[][] dataForEvidenceAnnotation() {
        final List<Object[]> data = new ArrayList<>(20);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint, new String[]{"60"}, new String[]{String.valueOf(1984)}});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_minus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_plus, new String[]{"60"}, new String[]{String.valueOf(100)}});

        // long range substitution
        data.add(new Object[]{forLongRangeSubstitution_minus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_plus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_minus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_plus, new String[]{"60"}, new String[]{String.valueOf(50)}});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_minus, new String[]{"60"}, new String[]{String.valueOf(137)}});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_plus, new String[]{"60"}, new String[]{String.valueOf(127)}});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_minus, new String[]{"60"}, new String[]{String.valueOf(31)}});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_plus, new String[]{"60"}, new String[]{String.valueOf(31)}});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_minus, new String[]{"60"}, new String[]{String.valueOf(127)}});

        return data.toArray(new Object[data.size()][]);
    }

    // -----------------------------------------------------------------------------------------------
    // Integrative test
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv", dataProvider = "forIntegrativeTest")
    public void testIntegrative(final Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData,
                                final List<String> expectedAttributeKeys,
                                final Broadcast<ReferenceMultiSource> referenceBroadcast) throws IOException {

        final AlignmentInterval region1 = testData._1();
        final AlignmentInterval region2 = testData._2();

        final Iterable<ChimericAlignment> evidence = Collections.singletonList(new ChimericAlignment(region1, region2, Collections.emptyList(), testData._4()));

        final NovelAdjacencyReferenceLocations breakpoints = testData._3();

        final VariantContext variantContext =
                AnnotatedVariantProducer.produceAnnotatedVcFromInferredTypeAndRefLocations(breakpoints.leftJustifiedLeftRefLoc,
                        breakpoints.leftJustifiedRightRefLoc.getStart(), breakpoints.complication, SvTypeInference.inferFromNovelAdjacency(breakpoints),
                        null, evidence, SparkContextFactory.getTestSparkContext().broadcast(SVDiscoveryTestDataProvider.reference));

        final List<String> attributeKeys = variantContext.getAttributes().keySet().stream().sorted().collect(Collectors.toList());

        Assert.assertEquals(attributeKeys, expectedAttributeKeys);
    }

    @DataProvider(name = "forIntegrativeTest")
    private Object[][] dataForIntegrativeTest() {
        final List<Object[]> data = new ArrayList<>(20);

        final Broadcast<ReferenceMultiSource> referenceBroadcast = SparkContextFactory.getTestSparkContext().broadcast(reference);

        final Set<String> commonAttributes = Sets.newHashSet(VCFConstants.END_KEY, SVLEN, SVTYPE, TOTAL_MAPPINGS,
                HQ_MAPPINGS, MAPPING_QUALITIES, ALIGN_LENGTHS, MAX_ALIGN_LENGTH, CONTIG_NAMES);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(INV33, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()), referenceBroadcast});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_minus, commonAttributes.stream().sorted().collect(Collectors.toList()), referenceBroadcast});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_plus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(INSERTED_SEQUENCE).stream()).sorted().collect(Collectors.toList()), referenceBroadcast});

        // long range substitution
        data.add(new Object[]{forLongRangeSubstitution_minus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(INSERTED_SEQUENCE).stream()).sorted().collect(Collectors.toList()), referenceBroadcast});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_plus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()), referenceBroadcast});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_minus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_CONTRACTION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_plus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUP_SEQ_CIGARS, DUPLICATION_NUMBERS, DUP_ORIENTATIONS).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_minus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUP_SEQ_CIGARS, DUPLICATION_NUMBERS, DUP_ORIENTATIONS, INSERTED_SEQUENCE).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_plus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_minus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_CONTRACTION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_plus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_CONTRACTION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_minus, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                referenceBroadcast});

        return data.toArray(new Object[data.size()][]);
    }


    @DataProvider(name = "CIIntervals")
    public Object[][] getCIIntervalTests() {
        return new Object[][] {
                new Object[] { 200, new SVInterval(1, 190, 225), "-10,25", null},
                new Object[] { 200, new SVInterval(1, 200, 225), "0,25", null},
                new Object[] { 200, new SVInterval(1, 201, 225), null, new IllegalStateException("Interval must contain point")}
        };
    }

    @Test(dataProvider = "CIIntervals")
    public void testProduceCIInterval(final int point, final SVInterval interval, final String expected, final Exception expectedException) throws Exception {
        if (expectedException == null) {
            Assert.assertEquals(AnnotatedVariantProducer.produceCIInterval(point, interval), expected);
        } else {
            try {
                AnnotatedVariantProducer.produceCIInterval(point, interval);
                Assert.fail("did not throw expected exception " + expectedException);
            } catch (Throwable e) {
                Assert.assertEquals(e.getClass(), expectedException.getClass());
                Assert.assertEquals(e.getMessage(), expectedException.getMessage());
            }
        }
    }
}