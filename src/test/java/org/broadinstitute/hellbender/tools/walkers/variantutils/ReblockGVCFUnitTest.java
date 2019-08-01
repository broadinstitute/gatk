package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ReblockGVCFUnitTest {
    private final static Allele LONG_REF = Allele.create("ACTG", true);
    private final static Allele DELETION = Allele.create("A", false);
    private final static Allele SHORT_REF = Allele.create("A", true);
    private final static Allele LONG_SNP = Allele.create("TCTA", false);

    @Test
    public void testCleanUpHighQualityVariant() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        //We need an annotation engine for cleanUpHighQualityVariant()
        reblocker.createAnnotationEngine();
        reblocker.dropLowQuals = true;
        reblocker.doQualApprox = true;

        final Genotype g0 = makeG("sample1", LONG_REF, DELETION, 41, 0, 37, 200, 100, 200, 400, 600, 800);
        final Genotype g = addAD(g0,0,13,17,0);
        final VariantContext extraAlt0 = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, LONG_SNP, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final Map<String, Object> attr = new HashMap<>();
        attr.put(VCFConstants.DEPTH_KEY, 32);
        final VariantContext extraAlt = addAttributes(extraAlt0, attr);
        //we'll call this with the same VC again under the assumption that STAND_CALL_CONF is zero so no alleles/GTs change
        final VariantContext cleaned1 = reblocker.cleanUpHighQualityVariant(extraAlt, extraAlt);
        Assert.assertTrue(cleaned1.getAlleles().size() == 3);
        Assert.assertTrue(cleaned1.getAlleles().contains(LONG_REF));
        Assert.assertTrue(cleaned1.getAlleles().contains(DELETION));
        Assert.assertTrue(cleaned1.getAlleles().contains(Allele.NON_REF_ALLELE));
        Assert.assertTrue(cleaned1.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        Assert.assertTrue(cleaned1.getAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY).equals(41));
        Assert.assertTrue(cleaned1.hasAttribute(GATKVCFConstants.VARIANT_DEPTH_KEY));
        Assert.assertTrue(cleaned1.getAttribute(GATKVCFConstants.VARIANT_DEPTH_KEY).equals(30));
        Assert.assertTrue(cleaned1.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        Assert.assertTrue(cleaned1.getAttributeAsString(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,"").split(",")[1].equals("32"));

        final Genotype hetNonRef = makeG("sample2", DELETION, LONG_SNP, 891,879,1128,84,0,30,891,879,84,891);
        final VariantContext keepAlts = makeDeletionVC("keepAllAlts", Arrays.asList(LONG_REF, DELETION, LONG_SNP, Allele.NON_REF_ALLELE), LONG_REF.length(), hetNonRef);
        Assert.assertTrue(keepAlts.getAlleles().size() == 4);
        Assert.assertTrue(keepAlts.getAlleles().contains(LONG_REF));
        Assert.assertTrue(keepAlts.getAlleles().contains(DELETION));
        Assert.assertTrue(keepAlts.getAlleles().contains(LONG_SNP));
        Assert.assertTrue(keepAlts.getAlleles().contains(Allele.NON_REF_ALLELE));
    }

    @Test
    public void testLowQualVariantToGQ0HomRef() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        reblocker.dropLowQuals = true;
        final Genotype g = makeG("sample1", LONG_REF, Allele.NON_REF_ALLELE, 200, 100, 200, 11, 0, 37);
        final VariantContext toBeNoCalled = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final VariantContext dropped = reblocker.lowQualVariantToGQ0HomRef(toBeNoCalled, toBeNoCalled);
        Assert.assertEquals(dropped, null);

        reblocker.dropLowQuals = false;
        final VariantContext modified = reblocker.lowQualVariantToGQ0HomRef(toBeNoCalled, toBeNoCalled);
        Assert.assertTrue(modified.getAttributes().containsKey(VCFConstants.END_KEY));
        Assert.assertTrue(modified.getAttributes().get(VCFConstants.END_KEY).equals(13));
        Assert.assertTrue(modified.getReference().equals(SHORT_REF));
        Assert.assertTrue(modified.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE));
        Assert.assertTrue(!modified.filtersWereApplied());
        Assert.assertTrue(modified.getLog10PError() == VariantContext.NO_LOG10_PERROR);

        //No-calls were throwing NPEs.  Now they're not.
        final Genotype g2 = makeG("sample1", Allele.NO_CALL,Allele.NO_CALL);
        final VariantContext noData = makeDeletionVC("noData", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g2);
        final VariantContext notCrashing = reblocker.lowQualVariantToGQ0HomRef(noData, noData);
        Assert.assertTrue(notCrashing.getGenotype(0).isNoCall());
    }

    @Test
    public void testChangeCallToGQ0HomRef() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final Genotype g = makeG("sample1", LONG_REF, Allele.NON_REF_ALLELE, 200, 100, 200, 11, 0, 37);
        final VariantContext toBeNoCalled = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final Map<String, Object> noAttributesMap = new HashMap<>();
        final GenotypeBuilder noCalled = reblocker.changeCallToGQ0HomRef(toBeNoCalled, noAttributesMap);
        final Genotype newG = noCalled.make();
        Assert.assertTrue(noAttributesMap.containsKey(VCFConstants.END_KEY));
        Assert.assertTrue(noAttributesMap.get(VCFConstants.END_KEY).equals(13));
        Assert.assertTrue(newG.getAllele(0).equals(SHORT_REF));
        Assert.assertTrue(newG.getAllele(1).equals(SHORT_REF));
        Assert.assertTrue(!newG.hasAD());
    }

    @Test  //no-calls can be dropped or reblocked just like hom-refs, i.e. we don't have to preserve them like variants
    public void testNoCalls() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final Genotype g2 = makeG("sample1", Allele.NO_CALL,Allele.NO_CALL);
        final VariantContext noData = makeDeletionVC("noData", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g2);
        Assert.assertTrue(reblocker.shouldBeReblocked(noData));
    }

    //TODO: these are duplicated from PosteriorProbabilitiesUtilsUnitTest but PR #4947 modifies VariantContextTestUtils, so I'll do some refactoring before the second of the two is merged
    private Genotype makeG(final String sample, final Allele a1, final Allele a2, final int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private VariantContext makeDeletionVC(final String source, final List<Allele> alleles, final int refLength, final Genotype... genotypes) {
        final int start = 10;
        final int stop = start+refLength-1;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(Arrays.asList(genotypes)).unfiltered().make();
    }

    private Genotype addAD(final Genotype g, final int... ads) {
        return new GenotypeBuilder(g).AD(ads).make();
    }

    private VariantContext addAttributes(final VariantContext vc, final Map<String, Object> attributes) {
        return new VariantContextBuilder(vc).attributes(attributes).make();
    }

}