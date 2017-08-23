package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.integration.SVIntegrationTestDataProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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

/**
 * Unit tests for {@link StructuralVariantContext}.
 */
public class StructuralVariantContextUnitTest {

    private static final File VALID_VARIANTS_FILE = new File(BaseTest.largeFileTestDir, "svc_test.vcf.gz");

    private static final File REFERENCE_FILE = new File(CommandLineProgramTest.b38_reference_20_21);

    /**
     * Tests {@link StructuralVariantContext#create}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts")
    public void testCreate(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        Assert.assertNotNull(svc);
    }

    /**
     * Tests {@link StructuralVariantContext#getStructuralVariantLength()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testLength(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final int length = svc.getStructuralVariantLength();
        if (!vc.hasAttribute(GATKSVVCFConstants.SVLEN)) {
            Assert.assertEquals(length, -1);
        } else {
            Assert.assertEquals(length, Math.abs(vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)));
        }
    }

    /**
     * Tests {@link StructuralVariantContext#getStructuralVariantType()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testType(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final StructuralVariantType type = svc.getStructuralVariantType();
        Assert.assertNotNull(type);
        final Allele alternativeAllele = svc.getAlternateAllele(0);
        final String alternativeAlleleName = alternativeAllele.getDisplayString();
        Assert.assertEquals("<" + type.name() + ">" , alternativeAlleleName);
    }

    /**
     * Tests {@link StructuralVariantContext#getInsertedSequence()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testInsertedSequence(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final byte[] actual = svc.getInsertedSequence();
        Assert.assertEquals(actual == null ? "<null>" : new String(actual), vc.getAttributeAsString(GATKSVVCFConstants.INSERTED_SEQUENCE, "<null>"));
    }

    /**
     * Tests {@link StructuralVariantContext#getEnd()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testEnd(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        Assert.assertEquals(svc.getEnd(), vc.getEnd());
    }

    /**
     * Tests {@link StructuralVariantContext#composeHaplotypeBasedOnReference(int, int, ReferenceMultiSource, PipelineOptions)}} when used
     * to obtain the reference haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"})
    public void testComposeReferenceHaplotype(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) throws IOException {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final int paddingSize = 10;
        final Haplotype refHaplotype = svc.composeHaplotypeBasedOnReference(0, paddingSize, reference, null);
        Assert.assertNotNull(refHaplotype);
        Assert.assertTrue(refHaplotype.isReference());
        final SimpleInterval expectedInterval = new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getEnd() + paddingSize);
        Assert.assertEquals(new SimpleInterval(refHaplotype.getGenomeLocation()), expectedInterval, svc.getContig() + ":" + svc.getStart() + " " + svc.getStructuralVariantType() + " " + svc.getStructuralVariantLength());
        Assert.assertEquals(refHaplotype.getBases(), reference.getReferenceBases(null, expectedInterval).getBases());
        Assert.assertEquals(refHaplotype.getCigar(), new Cigar(Collections.singletonList(new CigarElement(expectedInterval.size(), CigarOperator.M))));
    }

    /**
     * Tests {@link StructuralVariantContext#composeHaplotypeBasedOnReference(int, int, ReferenceMultiSource, PipelineOptions)}} when used
     * to obtain the reference alternative haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"})
    public void testComposeAlternativeHaplotype(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) throws IOException {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        if (svc.getStructuralVariantType() != StructuralVariantType.INS && svc.getStructuralVariantType() != StructuralVariantType.DEL) {
            throw new SkipException("unsupported type; skipped for now");
        }
        final int paddingSize = 10;
        final Haplotype altHaplotype = svc.composeHaplotypeBasedOnReference(1, paddingSize, reference, null);
        Assert.assertNotNull(altHaplotype);
        Assert.assertFalse(altHaplotype.isReference());
        final SimpleInterval expectedInterval = new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getEnd() + paddingSize);
        Assert.assertEquals(new SimpleInterval(altHaplotype.getGenomeLocation()), expectedInterval, svc.getContig() + ":" + svc.getStart() + " " + svc.getStructuralVariantType() + " " + svc.getStructuralVariantLength());
        final byte[] expectedBases;
        final Cigar expectedCigar;
        if (svc.getStructuralVariantType() == StructuralVariantType.INS) {
            expectedBases = Utils.concat(reference.getReferenceBases((null), new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getStart())).getBases(),
                    svc.getInsertedSequence(),
                    reference.getReferenceBases(null, new SimpleInterval(vc.getContig(), vc.getStart() + 1, vc.getStart() + paddingSize)).getBases());
            expectedCigar = new Cigar(Arrays.asList(new CigarElement(paddingSize, CigarOperator.M), new CigarElement(svc.getStructuralVariantLength(), CigarOperator.I), new CigarElement(paddingSize, CigarOperator.M)));
        } else { // must be DEL.
            expectedBases = Utils.concat(reference.getReferenceBases(null, new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getStart())).getBases(),
                    reference.getReferenceBases(null, new SimpleInterval(vc.getContig(), vc.getStart() + svc.getStructuralVariantLength() + 1, vc.getStart() + svc.getStructuralVariantLength() + paddingSize)).getBases());
            expectedCigar = new Cigar(Arrays.asList(new CigarElement(paddingSize, CigarOperator.M), new CigarElement(svc.getStructuralVariantLength(), CigarOperator.D), new CigarElement(paddingSize, CigarOperator.M)));
        }
        Assert.assertEquals(altHaplotype.getCigar(), expectedCigar);
        Assert.assertEquals(altHaplotype.getBases(), expectedBases, svc.getStructuralVariantType().name() + " " + new String(altHaplotype.getBases()) + " vs " + new String(expectedBases));
    }


    /**
     * Tests {@link StructuralVariantContext#getBreakPointIntervals(int, SAMSequenceDictionary)}.
     * to obtain the reference haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"})
    public void testGetBreakPoints(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) throws IOException {
        testGetBreakPoints(vc, reference, 0);
        testGetBreakPoints(vc, reference, 10);
    }

    private void testGetBreakPoints(final VariantContext vc, final ReferenceMultiSource reference, final int paddingSize) throws IOException {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        if (svc.getStructuralVariantType() != StructuralVariantType.INS && svc.getStructuralVariantType() != StructuralVariantType.DEL) {
            throw new SkipException("unsupported type; skipped for now");
        }
        final List<SimpleInterval> breakPoints = svc.getBreakPointIntervals(paddingSize, reference.getReferenceSequenceDictionary(null));
        final int contigLength = reference.getReferenceSequenceDictionary(null).getSequence(vc.getContig()).getSequenceLength();
        final List<Integer> expectedOffsets = new ArrayList<>();
        if (svc.getStructuralVariantType() == StructuralVariantType.INS) {
            expectedOffsets.add(vc.getStart());
        } else if (svc.getStructuralVariantType() == StructuralVariantType.DEL) {
            expectedOffsets.add(vc.getStart());
            expectedOffsets.add(vc.getEnd());
        }
        final List<SimpleInterval> expectedBreakPoints = expectedOffsets.stream()
                .map(i -> new SimpleInterval(vc.getContig(), Math.max(1, paddingSize > 0 ? (i - paddingSize + 1) : i),
                        Math.min(contigLength, i + paddingSize)))
                .collect(Collectors.toList());
        Assert.assertEquals(breakPoints, expectedBreakPoints);
    }

    /**
     * Tests {@link StructuralVariantContext#getSupportingContigIds()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testContigNames(final VariantContext vc, @SuppressWarnings("unused") final ReferenceMultiSource reference) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final List<String> actual = svc.getSupportingContigIds();
        final List<String> expected = vc.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, null);
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "outputVariantTestFilesData", dependsOnMethods = {"testCreate"})
    public void testTestOutputFileContent(final VariantContext vc, final ReferenceMultiSource reference, final String file) throws IOException {
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
        final ReferenceMultiSource reference = referenceMultiSource(REFERENCE_FILE.getAbsolutePath());
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
            final ReferenceMultiSource reference = referenceMultiSource(referenceName);
            try (final VCFFileReader reader = new VCFFileReader(new File(file), false)) {
                reader.forEach(vc -> result.add(new Object[] {vc, reference, file}));
            } catch (final Throwable ex) {
                throw new TestException("could not load the valid context file: " + VALID_VARIANTS_FILE.getAbsolutePath(), ex);
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    private static ReferenceMultiSource referenceMultiSource(final String fastaFileName) {
        return new ReferenceMultiSource((PipelineOptions) null, fastaFileName,
                (r) -> new SimpleInterval(r.getContig(), r.getAssignedStart(), r.getEnd()));
    }
}
