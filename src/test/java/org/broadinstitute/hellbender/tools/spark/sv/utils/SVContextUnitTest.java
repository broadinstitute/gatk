package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.integration.SVIntegrationTestDataProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


public class SVContextUnitTest extends GATKBaseTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "spark/sv/utils/";

    private static final File VALID_VARIANTS_FILE = new File(TEST_SUB_DIR, "SVContext.vcf.gz");

    private static final File REFERENCE_FILE = new File(b38_reference_20_21);

    /**
     * Tests {@link SVContext#of}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", groups = "sv")
    public void testCreate(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) {
        final SVContext svc = SVContext.of(vc);
        Assert.assertNotNull(svc);
    }

    /**
     * Tests {@link SVContext#getStructuralVariantLength()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"}, groups = "sv")
    public void testLength(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) {
        final SVContext svc = SVContext.of(vc);
        final int length = svc.getStructuralVariantLength();
        if (!vc.hasAttribute(GATKSVVCFConstants.SVLEN)) {
            Assert.assertEquals(length, -1);
        } else {
            Assert.assertEquals(length, Math.abs(vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)));
        }
    }

    /**
     * Tests {@link SVContext#getStructuralVariantType()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"}, groups = "sv")
    public void testType(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) {
        final SVContext svc = SVContext.of(vc);
        final StructuralVariantType type = svc.getStructuralVariantType();
        Assert.assertNotNull(type);
        final Allele alternativeAllele = svc.getAlternateAllele(0);
        final String alternativeAlleleName = alternativeAllele.getDisplayString();
        Assert.assertEquals("<" + type.name() + ">" , alternativeAlleleName);
    }

    /**
     * Tests {@link SVContext#getInsertedSequence()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"}, groups = "sv")
    public void testInsertedSequence(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) {
        final SVContext svc = SVContext.of(vc);
        final byte[] actual = svc.getInsertedSequence();
        Assert.assertEquals(actual == null ? "<null>" : new String(actual), vc.getAttributeAsString(GATKSVVCFConstants.INSERTED_SEQUENCE, "<null>"));
    }

    /**
     * Tests {@link SVContext#getEnd()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"}, groups = "sv")
    public void testEnd(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) {
        final SVContext svc = SVContext.of(vc);
        Assert.assertEquals(svc.getEnd(), vc.getEnd());
    }

    /**
     * Tests {@link SVContext#composeHaplotypeBasedOnReference(int, int, ReferenceMultiSparkSource)}} when used
     * to obtain the reference haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"}, groups = "sv")
    public void testComposeReferenceHaplotype(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) throws IOException {
        final SVContext svc = SVContext.of(vc);
        final int paddingSize = 10;
        final Haplotype refHaplotype = svc.composeHaplotypeBasedOnReference(0, paddingSize, reference);
        Assert.assertNotNull(refHaplotype);
        Assert.assertTrue(refHaplotype.isReference());
        final SimpleInterval expectedInterval = new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getEnd() + paddingSize);
        Assert.assertEquals(new SimpleInterval(refHaplotype.getGenomeLocation()), expectedInterval, svc.getContig() + ":" + svc.getStart() + " " + svc.getStructuralVariantType() + " " + svc.getStructuralVariantLength());
        Assert.assertEquals(refHaplotype.getBases(), reference.getReferenceBases(expectedInterval).getBases());
        Assert.assertEquals(refHaplotype.getCigar(), new Cigar(Collections.singletonList(new CigarElement(expectedInterval.size(), CigarOperator.M))));
    }

    /**
     * Tests {@link SVContext#composeHaplotypeBasedOnReference(int, int, ReferenceMultiSparkSource)}} when used
     * to obtain the reference alternative haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"}, groups = "sv")
    public void testComposeAlternativeHaplotype(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) throws IOException {
        final SVContext svc = SVContext.of(vc);
        if (svc.getStructuralVariantType() != StructuralVariantType.INS && svc.getStructuralVariantType() != StructuralVariantType.DEL) {
            throw new SkipException("unsupported type; skipped for now");
        }
        final int paddingSize = 10;
        final Haplotype altHaplotype = svc.composeHaplotypeBasedOnReference(1, paddingSize, reference);
        Assert.assertNotNull(altHaplotype);
        Assert.assertFalse(altHaplotype.isReference());
        final SimpleInterval expectedInterval = new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getEnd() + paddingSize);
        Assert.assertEquals(new SimpleInterval(altHaplotype.getGenomeLocation()), expectedInterval, svc.getContig() + ":" + svc.getStart() + " " + svc.getStructuralVariantType() + " " + svc.getStructuralVariantLength());
        final byte[] expectedBases;
        final Cigar expectedCigar;
        if (svc.getStructuralVariantType() == StructuralVariantType.INS) {
            expectedBases = Utils.concat(reference.getReferenceBases(
                    new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getStart())).getBases(),
                    svc.getInsertedSequence(),
                    reference.getReferenceBases(new SimpleInterval(vc.getContig(), vc.getStart() + 1, vc.getStart() + paddingSize)).getBases());
            expectedCigar = new Cigar(Arrays.asList(new CigarElement(paddingSize, CigarOperator.M), new CigarElement(svc.getStructuralVariantLength(), CigarOperator.I), new CigarElement(paddingSize, CigarOperator.M)));
        } else { // must be DEL.
            expectedBases = Utils.concat(reference.getReferenceBases(
                    new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getStart())).getBases(),
                    reference.getReferenceBases(
                            new SimpleInterval(vc.getContig(), vc.getStart() + svc.getStructuralVariantLength() + 1, vc.getStart() + svc.getStructuralVariantLength() + paddingSize)).getBases());
            expectedCigar = new Cigar(Arrays.asList(new CigarElement(paddingSize, CigarOperator.M), new CigarElement(svc.getStructuralVariantLength(), CigarOperator.D), new CigarElement(paddingSize, CigarOperator.M)));
        }
        Assert.assertEquals(altHaplotype.getCigar(), expectedCigar);
        Assert.assertEquals(altHaplotype.getBases(), expectedBases, svc.getStructuralVariantType().name() + " " + new String(altHaplotype.getBases()) + " vs " + new String(expectedBases));
    }


    /**
     * Tests {@link SVContext#getBreakPointIntervals(int, SAMSequenceDictionary, boolean)}.
     * to obtain the reference haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"}, groups = "sv")
    public void testGetBreakPoints(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) throws IOException {
        testGetBreakPoints(vc, reference, 0);
        testGetBreakPoints(vc, reference, 10);
    }

    private void testGetBreakPoints( final VariantContext vc, final ReferenceMultiSparkSource reference, final int paddingSize) throws IOException {
        final SVContext svc = SVContext.of(vc);
        if (svc.getStructuralVariantType() != StructuralVariantType.INS && svc.getStructuralVariantType() != StructuralVariantType.DEL) {
            throw new SkipException("unsupported type; skipped for now");
        }
        final List<SimpleInterval> breakPoints = svc.getBreakPointIntervals(paddingSize, reference.getReferenceSequenceDictionary(null), false);
        final int contigLength = reference.getReferenceSequenceDictionary(null).getSequence(vc.getContig()).getSequenceLength();
        final List<Integer> expectedOffsets = new ArrayList<>();
        if (svc.getStructuralVariantType() == StructuralVariantType.INS) {
            expectedOffsets.add(vc.getStart());
        } else if (svc.getStructuralVariantType() == StructuralVariantType.DEL) {
            expectedOffsets.add(vc.getStart() + 1);
            expectedOffsets.add(vc.getEnd());
        }
        final List<SimpleInterval> expectedBreakPoints = expectedOffsets.stream()
                .map(i -> new SimpleInterval(vc.getContig(), Math.max(1, paddingSize > 0 ? (i - paddingSize + 1) : i),
                        Math.min(contigLength, i + paddingSize)))
                .collect(Collectors.toList());
        Assert.assertEquals(breakPoints, expectedBreakPoints);
    }

    /**
     * Tests {@link SVContext#getSupportingContigIds()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"}, groups = "sv")
    public void testContigNames(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSparkSource reference) {
        final SVContext svc = SVContext.of(vc);
        final List<String> actual = svc.getSupportingContigIds();
        final List<String> expected = vc.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, null);
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "outputVariantTestFilesData", dependsOnMethods = {"testCreate"}, groups = "sv")
    public void testTestOutputFileContent( final VariantContext vc, final ReferenceMultiSparkSource reference, final String file) throws IOException {
        testCreate(vc, reference);
        testLength(vc, reference);
        testEnd(vc, reference);
        testType(vc, reference);
        testContigNames(vc, reference);
        testInsertedSequence(vc, reference);
        // these test might skip for records with yet unsupported SV types.
        // We just silence the skips.
        try { testComposeAlternativeHaplotype(vc, reference); } catch (final SkipException ex) {}
        try { testComposeReferenceHaplotype(vc, reference); } catch (final SkipException ex) {}
    }

    @DataProvider(name="validVariantContexts")
    public Object[][] validVariantContexts(){
        final ReferenceMultiSparkSource reference = referenceMultiSource(REFERENCE_FILE.getAbsolutePath());
        try (final VCFFileReader reader = new VCFFileReader(VALID_VARIANTS_FILE, false)) {
            return Utils.stream(reader)
                    .map(vc -> new Object[]{vc, reference})
                    .toArray(Object[][]::new);
        } catch (final Throwable ex) {
            throw new TestException("could not load the valid context file: " + VALID_VARIANTS_FILE.getAbsolutePath(), ex);
        }
    }

    @DataProvider(name = "outputVariantTestFilesData")
    public Object[][] outputVariantTestFilesData() {

        final List<Tuple2<String, String>> outputFilesAndReference = new ArrayList<>();
        outputFilesAndReference.add(
                new Tuple2<>(SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF, SVIntegrationTestDataProvider.reference.getAbsolutePath()));
        outputFilesAndReference.add(
                new Tuple2<>(SVIntegrationTestDataProvider.EXPECTED_SIMPLE_INV_VCF, SVIntegrationTestDataProvider.reference.getAbsolutePath()));

        final List<Object[]> result = new ArrayList<>();
        for (final Tuple2<String, String> outputAndReference : outputFilesAndReference) {
            final String file = outputAndReference._1();
            final String referenceName = outputAndReference._2();
            final ReferenceMultiSparkSource reference = referenceMultiSource(referenceName);
            try (final VCFFileReader reader = new VCFFileReader(new File(file), false)) {
                reader.forEach(vc -> result.add(new Object[] {vc, reference, file}));
            } catch (final Throwable ex) {
                throw new TestException("could not load the valid context file: " + VALID_VARIANTS_FILE.getAbsolutePath(), ex);
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    private static ReferenceMultiSparkSource referenceMultiSource( final String fastaFileName) {
        return new ReferenceMultiSparkSource(fastaFileName,
                                        (r) -> new SimpleInterval(r.getContig(), r.getAssignedStart(), r.getEnd()));
    }
}
