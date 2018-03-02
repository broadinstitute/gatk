package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.IntegerCopyNumberSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.IntegerCopyNumberSegmentCollectionUnitTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Unit test for {@link GermlineCNVSegmentVariantComposer}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVSegmentVariantComposerUnitTest extends GATKBaseTest {
    @Test(dataProvider = "variantCompositionSettings")
    public void testVariantComposition(final int refAutosomalCopyNumber,
                                       final Set<String> allosomalContigs) {
        /* read a segments collection */
        final IntegerCopyNumberSegmentCollection collection = new IntegerCopyNumberSegmentCollection(
                IntegerCopyNumberSegmentCollectionUnitTest.TEST_INTEGER_COPY_NUMBER_SEGMENTS_FILE);
        final File segmentsOutputVCF = createTempFile("test-write-segments", ".vcf");
        final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(segmentsOutputVCF,
                collection.getMetadata().getSequenceDictionary(), false);
        final GermlineCNVSegmentVariantComposer variantComposer = new GermlineCNVSegmentVariantComposer(
                writer, IntegerCopyNumberSegmentCollectionUnitTest.EXPECTED_SAMPLE_NAME,
                new IntegerCopyNumberState(refAutosomalCopyNumber), allosomalContigs);

        /* compose segments and assert correctness */
        for (final IntegerCopyNumberSegment segment: collection.getRecords()) {
            final VariantContext var = variantComposer.composeVariantContext(segment);
            Assert.assertEquals(var.getContig(), segment.getContig());
            Assert.assertEquals(var.getStart(), segment.getStart());
            Assert.assertEquals(var.getEnd(), segment.getEnd());
            Assert.assertEquals(var.getAlleles(), GermlineCNVSegmentVariantComposer.ALL_ALLELES);

            final Genotype gt = var.getGenotype(IntegerCopyNumberSegmentCollectionUnitTest.EXPECTED_SAMPLE_NAME);

            /* assert allele correctness */
            final Allele actualAllele = gt.getAlleles().get(0);
            final int refCopyNumber = allosomalContigs.contains(segment.getContig())
                    ? segment.getBaselineIntegerCopyNumberState().getCopyNumber()
                    : refAutosomalCopyNumber;
            final Allele expectedAllele;
            if (segment.getCallIntegerCopyNumberState().getCopyNumber() > refCopyNumber) {
                expectedAllele = GermlineCNVSegmentVariantComposer.DUP_ALLELE;
            } else if (segment.getCallIntegerCopyNumberState().getCopyNumber() < refCopyNumber) {
                expectedAllele = GermlineCNVSegmentVariantComposer.DEL_ALLELE;
            } else {
                expectedAllele = GermlineCNVSegmentVariantComposer.REF_ALLELE;
            }
            Assert.assertEquals(actualAllele, expectedAllele);

            /* assert correctness of quality metrics */
            Assert.assertEquals(
                    (long) gt.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS),
                    FastMath.round(segment.getQualitySomeCalled()));
            Assert.assertEquals(
                    (long) gt.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QA),
                    FastMath.round(segment.getQualityAllCalled()));
            Assert.assertEquals(
                    (long) gt.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSS),
                    FastMath.round(segment.getQualityStart()));
            Assert.assertEquals(
                    (long) gt.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSE),
                    FastMath.round(segment.getQualityEnd()));
        }
    }

    @DataProvider(name = "variantCompositionSettings")
    public Object[][] getVariantCompositionSettings() {
        return new Object[][] {
                {2, new HashSet<>(Arrays.asList("X", "Y"))},
                {2, new HashSet<>()},
                {3, new HashSet<>()},
                {1, new HashSet<>(Arrays.asList("1", "X"))},
        };
    }
}
