package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.IntegerCopyNumberSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.IntegerCopyNumberSegmentCollectionUnitTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
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
                                       final Set<String> allosomalContigs,
                                       final ReferenceSequenceFile reference) {
        /* read a segments collection */
        final IntegerCopyNumberSegmentCollection collection = new IntegerCopyNumberSegmentCollection(
                IntegerCopyNumberSegmentCollectionUnitTest.TEST_INTEGER_COPY_NUMBER_SEGMENTS_FILE);
        final File segmentsOutputVCF = createTempFile("test-write-segments", ".vcf");
        final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(segmentsOutputVCF.toPath(),
                collection.getMetadata().getSequenceDictionary(), false);
        final GermlineCNVSegmentVariantComposer variantComposer = new GermlineCNVSegmentVariantComposer(
                writer, IntegerCopyNumberSegmentCollectionUnitTest.EXPECTED_SAMPLE_NAME,
                new IntegerCopyNumberState(refAutosomalCopyNumber), allosomalContigs, reference);

        /* compose segments and assert correctness */
        for (final IntegerCopyNumberSegment segment: collection.getRecords()) {
            final VariantContext var = variantComposer.composeVariantContext(segment);
            Assert.assertEquals(var.getContig(), segment.getContig());
            Assert.assertEquals(var.getStart(), segment.getStart());
            Assert.assertEquals(var.getEnd(), segment.getEnd());

            final Genotype gt = var.getGenotype(IntegerCopyNumberSegmentCollectionUnitTest.EXPECTED_SAMPLE_NAME);

            /* assert allele correctness */
            final int refCopyNumber = allosomalContigs.contains(segment.getContig())
                    ? segment.getBaselineIntegerCopyNumberState().getCopyNumber()
                    : refAutosomalCopyNumber;
            if (segment.getBaselineIntegerCopyNumberState().getCopyNumber() != 0) {
                Assert.assertEquals(gt.getPloidy(), refCopyNumber);  //ploidy should match the autosomal reference copy number in the event that there's an aneuploidy
            }
            final Allele expectedAllele;
            if (segment.getCallIntegerCopyNumberState().getCopyNumber() == refCopyNumber) {
                Assert.assertEquals(var.getAlternateAlleles().size(), 0);
                Assert.assertTrue(gt.isHomRef());
            } else if (segment.getCallIntegerCopyNumberState().getCopyNumber() > refCopyNumber) {
                final Allele actualAllele = var.getAlternateAllele(0);
                Assert.assertEquals(var.getAlternateAlleles().size(), 1);
                expectedAllele = GATKSVVCFConstants.DUP_ALLELE;
                Assert.assertEquals(actualAllele, expectedAllele);
                if (refCopyNumber > 1) {
                    Assert.assertTrue(gt.getAlleles().stream().allMatch(a -> a.equals(Allele.NO_CALL)));
                } else {
                    Assert.assertEquals(gt.getAllele(0), GATKSVVCFConstants.DUP_ALLELE);
                }
            } else { //if (segment.getCallIntegerCopyNumberState().getCopyNumber() < refCopyNumber) {
                final Allele actualAllele = var.getAlternateAllele(0);
                Assert.assertEquals(var.getAlternateAlleles().size(), 1);
                expectedAllele = GATKSVVCFConstants.DEL_ALLELE;
                Assert.assertEquals(actualAllele, expectedAllele);
            }

            Assert.assertEquals(var.getAlleles().size(), gt.isHomRef() ? 1 : 2);

            //if a reference is supplied, use that to get the reference allele
            if (reference == null) {
                Assert.assertEquals(var.getReference(), Allele.REF_N);
            } else {
                Assert.assertFalse(var.getReference().isSymbolic());
                Assert.assertEquals(var.getReference().getBaseString().length(), 1);
                Assert.assertEquals(var.getReference().getBases(), ReferenceUtils.getRefBaseAtPosition(reference, var.getContig(), var.getStart()));
            }

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
                {2, new HashSet<>(Arrays.asList("X", "Y")), null},
                {2, new HashSet<>(), null},
                {3, new HashSet<>(), null},
                {1, new HashSet<>(Arrays.asList("1", "X")), null},
                {2, new HashSet<>(Arrays.asList("X","Y")), ReferenceUtils.createReferenceReader(new GATKPath(GATKBaseTest.b37Reference))}
        };
    }
}
