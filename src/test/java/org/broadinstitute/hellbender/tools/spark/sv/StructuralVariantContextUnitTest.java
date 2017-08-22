package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.SimpleIntervalUnitTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

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
public class StructuralVariantContextUnitTest extends BaseTest {

    private static final File VALID_VARIANTS_FILE = new File(largeFileTestDir, "svc_test.vcf.gz");

    private static final File REFERENCE_FILE = new File(CommandLineProgramTest.b38_reference_20_21);

    /**
     * Tests {@link StructuralVariantContext#create}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts")
    public void testCreate(final VariantContext vc) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        Assert.assertNotNull(svc);
    }

    /**
     * Tests {@link StructuralVariantContext#getStructuralVariantLength()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testLength(final VariantContext vc) {
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
    public void testType(final VariantContext vc) {
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
    public void testInsertedSequence(final VariantContext vc) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final byte[] actual = svc.getInsertedSequence();
        Assert.assertEquals(actual == null ? "<null>" : new String(actual), vc.getAttributeAsString(GATKSVVCFConstants.INSERTED_SEQUENCE, "<null>"));
    }

    /**
     * Tests {@link StructuralVariantContext#getEnd()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testEnd(final VariantContext vc) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        Assert.assertEquals(svc.getEnd(), vc.getEnd());
    }

    /**
     * Tests {@link StructuralVariantContext#composeHaplotypeBasedOnReference(int, int, ReferenceMultiSource)}} when used
     * to obtain the reference haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"})
    public void testComposeReferenceHaplotype(final VariantContext vc) throws IOException {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final ReferenceMultiSource reference = referenceMultiSource();
        final int paddingSize = 10;
        final Haplotype refHaplotype = svc.composeHaplotypeBasedOnReference(0, paddingSize, reference);
        Assert.assertNotNull(refHaplotype);
        Assert.assertTrue(refHaplotype.isReference());
        final SimpleInterval expectedInterval = new SimpleInterval(vc.getContig(), vc.getStart() + 1 - paddingSize, vc.getEnd() + paddingSize);
        Assert.assertEquals(new SimpleInterval(refHaplotype.getGenomeLocation()), expectedInterval, svc.getContig() + ":" + svc.getStart() + " " + svc.getStructuralVariantType() + " " + svc.getStructuralVariantLength());
        Assert.assertEquals(refHaplotype.getBases(), reference.getReferenceBases(null, expectedInterval).getBases());
        Assert.assertEquals(refHaplotype.getCigar(), new Cigar(Collections.singletonList(new CigarElement(expectedInterval.size(), CigarOperator.M))));
    }

    /**
     * Tests {@link StructuralVariantContext#composeHaplotypeBasedOnReference(int, int, ReferenceMultiSource)}} when used
     * to obtain the reference alternative haplotype.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate", "testType", "testLength"})
    public void testComposeAlternativeHaplotype(final VariantContext vc) throws IOException {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        if (svc.getStructuralVariantType() != StructuralVariantType.INS && svc.getStructuralVariantType() != StructuralVariantType.DEL) {
            throw new SkipException("unsupported type; skipped for now");
        }
        final ReferenceMultiSource reference = referenceMultiSource();
        final int paddingSize = 10;
        final Haplotype altHaplotype = svc.composeHaplotypeBasedOnReference(1, paddingSize, reference);
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
    public void testGetBreakPoints(final VariantContext vc) throws IOException {
        testGetBreakPoints(vc, 0);
        testGetBreakPoints(vc, 10);
    }

    private void testGetBreakPoints(final VariantContext vc, final int paddingSize) throws IOException {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        if (svc.getStructuralVariantType() != StructuralVariantType.INS && svc.getStructuralVariantType() != StructuralVariantType.DEL) {
            throw new SkipException("unsupported type; skipped for now");
        }
        final ReferenceMultiSource reference = referenceMultiSource();
        final List<SimpleInterval> breakPoints = svc.getBreakPointIntervals(paddingSize, reference.getReferenceSequenceDictionary(null));
        final int contigLength = reference.getReferenceSequenceDictionary(null).getSequence(vc.getContig()).getSequenceLength();
        final List<Integer> expectedOffsets = new ArrayList<>();
        if (svc.getStructuralVariantType() == StructuralVariantType.INS) {
            expectedOffsets.add(vc.getStart());
        } else if (svc.getStructuralVariantType() == StructuralVariantType.DEL) {
            expectedOffsets.add(vc.getStart());
            expectedOffsets.add(vc.getStart() + svc.getStructuralVariantLength());
        }
        final List<SimpleInterval> expectedBreakPoints = expectedOffsets.stream()
                .map(i -> new SimpleInterval(vc.getContig(), Math.max(1, paddingSize > 0 ? (i - paddingSize + 1) : i),
                        Math.min(contigLength, i + paddingSize)))
                .collect(Collectors.toList());
        Assert.assertEquals(breakPoints, expectedBreakPoints);
    }

    /**
     * Tests {@link StructuralVariantContext#getContigNames()}.
     * @param vc input variant context.
     */
    @Test(dataProvider="validVariantContexts", dependsOnMethods = {"testCreate"})
    public void testContigNames(final VariantContext vc) {
        final StructuralVariantContext svc = StructuralVariantContext.create(vc);
        final List<String> actual = svc.getContigNames();
        final List<String> expected = vc.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, null);
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name="validVariantContexts")
    public Object[][] validVariantContexts(){
        try (final VCFFileReader reader = new VCFFileReader(VALID_VARIANTS_FILE, false)) {
            return Utils.stream(reader)
                    .map(vc -> new Object[]{vc})
                    .toArray(Object[][]::new);
        } catch (final Throwable ex) {
            throw new TestException("could not load the valid context file: " + VALID_VARIANTS_FILE.getAbsolutePath(), ex);
        }
    }

    public static ReferenceMultiSource referenceMultiSource() {
        return new ReferenceMultiSource((PipelineOptions) null, REFERENCE_FILE.getAbsolutePath(), (r) -> new SimpleInterval(r.getContig(), r.getAssignedStart(), r.getEnd()));
    }
}
