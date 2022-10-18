package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY;

/**
 * Tests {@link ReferenceConfidenceVariantContextMerger}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReferenceConfidenceVariantContextMergerUnitTest extends GATKBaseTest {
    private final Allele Aref = Allele.create("A", true);
    private final Allele C = Allele.create("C");
    private final Allele G = Allele.create("G");
    private final Allele ATC = Allele.create("ATC");
    private final Allele del = Allele.SPAN_DEL;
    private final Allele ATCref = Allele.create("ATC", true);


    private static VariantAnnotatorEngine getAnnotationEngine() {
        return new VariantAnnotatorEngine(VariantContextTestUtils.getAllAnnotations(), null, Collections.emptyList(), false, false);
    }

    @Test(dataProvider = "referenceConfidenceMergeData")
    public void testReferenceConfidenceMerge(final String testID, final List<VariantContext> toMerge, final Locatable loc,
                                             final boolean returnSiteEvenIfMonomorphic, final boolean uniquifySamples, final VariantContext expectedResult) {
        ReferenceConfidenceVariantContextMerger merger = new ReferenceConfidenceVariantContextMerger(getAnnotationEngine(), new VCFHeader());
        final VariantContext result = merger.merge(toMerge, loc, returnSiteEvenIfMonomorphic ? (byte) 'A' : null, true, uniquifySamples);
        if ( result == null ) {
            Assert.assertTrue(expectedResult == null);
            return;
        }
        Assert.assertEquals(result.getAlleles(), expectedResult.getAlleles(),testID);
        Assert.assertEquals(result.getNSamples(), expectedResult.getNSamples(),testID);
        for ( final Genotype expectedGenotype : expectedResult.getGenotypes() ) {
            Assert.assertTrue(result.hasGenotype(expectedGenotype.getSampleName()), "Missing " + expectedGenotype.getSampleName());
            VariantContextTestUtils.assertGenotypesAreEqual(result.getGenotype(expectedGenotype.getSampleName()), expectedGenotype);
        }
    }

    @DataProvider
    public Object[][] getVariousDepths() {
        Genotype baseGenotype = new GenotypeBuilder("sample", Arrays.asList(C, G)).make();
        return new Object[][]{
                {baseGenotype, 0},
                {new GenotypeBuilder(baseGenotype).DP(10).attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, 5).make(), 5},
                {new GenotypeBuilder(baseGenotype).DP(10).attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, "5").make(), 5},
                {new GenotypeBuilder(baseGenotype).DP(10).make(), 10}
        };
    }

    @Test(dataProvider = "getVariousDepths")
    public void testGetBestDepthValue(final Genotype genotype, final int expectedDepth){
        Assert.assertEquals(ReferenceConfidenceVariantContextMerger.getBestDepthValue(genotype), expectedDepth);
    }

    @Test
    public void testGenerateADWithNewAlleles() {

        final int[] originalAD = new int[] {1,2,0};
        final int[] indexesOfRelevantAlleles = new int[] {0,1,2,2};

        final int[] newAD = AlleleSubsettingUtils.generateAD(originalAD, indexesOfRelevantAlleles);
        Assert.assertEquals(newAD, new int[]{1,2,0,0});
    }

    @DataProvider(name = "referenceConfidenceMergeData")
    public Object[][] makeReferenceConfidenceMergeData() {
        final List<Object[]> tests = new ArrayList<>();
        final int start = 10;
        final SimpleInterval loc = new SimpleInterval("20", start, start);
        final VariantContext VCbase = new VariantContextBuilder("test", "20", start, start, Arrays.asList(Aref)).make();
        final VariantContext VCbase2 = new VariantContextBuilder("test2", "20", start, start, Arrays.asList(Aref)).make();
        final VariantContext VCprevBase = new VariantContextBuilder("test", "20", start-1, start-1, Arrays.asList(Aref)).make();

        final int[] standardPLs = new int[]{30, 20, 10, 71, 72, 73};
        final int[] reorderedSecondAllelePLs = new int[]{30, 71, 73, 20, 72, 10};

        final List<Allele> noCalls = new ArrayList<>(2);
        noCalls.add(Allele.NO_CALL);
        noCalls.add(Allele.NO_CALL);

        final List<Allele> A_ALT = Arrays.asList(Aref, Allele.NON_REF_ALLELE);
        final Genotype gA_ALT = new GenotypeBuilder("A").PL(new int[]{0, 100, 1000}).alleles(noCalls).make();
        final VariantContext vcA_ALT = new VariantContextBuilder(VCbase).alleles(A_ALT).genotypes(gA_ALT).make();

        final Allele AAref = Allele.create("AA", true);
        final List<Allele> AA_ALT = Arrays.asList(AAref, Allele.NON_REF_ALLELE);
        final Genotype gAA_ALT = new GenotypeBuilder("AA").PL(new int[]{0, 80, 800}).alleles(noCalls).make();
        final VariantContext vcAA_ALT = new VariantContextBuilder(VCprevBase).alleles(AA_ALT).genotypes(gAA_ALT).make();

        final List<Allele> A_C = Arrays.asList(Aref, C);
        final Genotype gA_C = new GenotypeBuilder("A_C").PL(new int[]{30, 20, 10}).alleles(noCalls).make();
        final List<Allele> A_C_ALT = Arrays.asList(Aref, C, Allele.NON_REF_ALLELE);
        final Genotype gA_C_ALT = new GenotypeBuilder("A_C").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcA_C = new VariantContextBuilder(VCbase2).alleles(A_C_ALT).genotypes(gA_C).make();
        final VariantContext vcA_C_ALT = new VariantContextBuilder(VCbase).alleles(A_C_ALT).genotypes(gA_C_ALT).make();

        final List<Allele> A_G_ALT = Arrays.asList(Aref, G, Allele.NON_REF_ALLELE);
        final Genotype gA_G_ALT = new GenotypeBuilder("A_G").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcA_G_ALT = new VariantContextBuilder(VCbase).alleles(A_G_ALT).genotypes(gA_G_ALT).make();

        final List<Allele> A_C_G = Arrays.asList(Aref, C, G);
        final Genotype gA_C_G = new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30}).alleles(noCalls).make();
        final List<Allele> A_C_G_ALT = Arrays.asList(Aref, C, G, Allele.NON_REF_ALLELE);
        final Genotype gA_C_G_ALT = new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30, 71, 72, 73, 74}).alleles(noCalls).make();
        final VariantContext vcA_C_G = new VariantContextBuilder(VCbase2).alleles(A_C_G_ALT).genotypes(gA_C_G).make();
        final VariantContext vcA_C_G_ALT = new VariantContextBuilder(VCbase).alleles(A_C_G_ALT).genotypes(gA_C_G_ALT).make();

        final List<Allele> A_ATC_ALT = Arrays.asList(Aref, ATC, Allele.NON_REF_ALLELE);
        final Genotype gA_ATC_ALT = new GenotypeBuilder("A_ATC").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcA_ATC_ALT = new VariantContextBuilder(VCbase).alleles(A_ATC_ALT).genotypes(gA_ATC_ALT).make();

        final Allele A = Allele.create("A", false);
        final List<Allele> AA_A_ALT = Arrays.asList(AAref, A, Allele.NON_REF_ALLELE);
        final Genotype gAA_A_ALT = new GenotypeBuilder("AA_A").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcAA_A_ALT = new VariantContextBuilder(VCprevBase).alleles(AA_A_ALT).genotypes(gAA_A_ALT).make();
        final List<Allele> A_C_del = Arrays.asList(Aref, C, del);


        // first test the case of a single record
        tests.add(new Object[]{"test00",Arrays.asList(vcA_C_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(gA_C).make()});

        // now, test pairs:
        // a SNP with another SNP
        tests.add(new Object[]{"test01",Arrays.asList(vcA_C_ALT, vcA_G_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(gA_C_ALT, new GenotypeBuilder("A_G").PL(reorderedSecondAllelePLs).alleles(noCalls).make()).make()});
        // a SNP with an indel
        tests.add(new Object[]{"test02",Arrays.asList(vcA_C_ALT, vcA_ATC_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(Arrays.asList(Aref, C, ATC)).genotypes(gA_C_ALT, new GenotypeBuilder("A_ATC").PL(reorderedSecondAllelePLs).alleles(noCalls).make()).make()});
        // a SNP with 2 SNPs
        tests.add(new Object[]{"test03",Arrays.asList(vcA_C_ALT, vcA_C_G_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(gA_C_ALT, gA_C_G).make()});
        // a SNP with a ref record
        tests.add(new Object[]{"test04",Arrays.asList(vcA_C_ALT, vcA_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(gA_C, gA_ALT).make()});

        // spanning records:
        // a SNP with a spanning ref record
        tests.add(new Object[]{"test05",Arrays.asList(vcA_C_ALT, vcAA_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(gA_C, gAA_ALT).make()});
        // a SNP with a spanning deletion
        tests.add(new Object[]{"test06",Arrays.asList(vcA_C_ALT, vcAA_A_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_del).genotypes(new GenotypeBuilder("A_C").PL(new int[]{30, 20, 10, 71, 72, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("AA_A").PL(new int[]{30, 71, 73, 20, 72, 10}).alleles(noCalls).make()).make()});

        // combination of all
        tests.add(new Object[]{"test07",Arrays.asList(vcA_C_ALT, vcA_G_ALT, vcA_ATC_ALT, vcA_C_G_ALT, vcA_ALT, vcAA_ALT, vcAA_A_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(Arrays.asList(Aref, C, G, ATC, del)).genotypes(
                        new GenotypeBuilder("A_C").PL(new int[]{30, 20, 10, 71, 72, 73, 71, 72, 73, 73, 71, 72, 73, 73, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_G").PL(new int[]{30, 71, 73, 20, 72, 10, 71, 73, 72, 73, 71, 73, 72, 73, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_ATC").PL(new int[]{30, 71, 73, 71, 73, 73, 20, 72, 72, 10, 71, 73, 73, 72, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30, 71, 72, 73, 74, 71, 72, 73, 74, 74}).alleles(noCalls).make(),
                        new GenotypeBuilder("A").PL(new int[]{0, 100, 1000, 100, 1000, 1000, 100, 1000, 1000, 1000, 100, 1000, 1000, 1000, 1000}).alleles(noCalls).make(),
                        new GenotypeBuilder("AA").PL(new int[]{0, 80, 800, 80, 800, 800, 80, 800, 800, 800, 80, 800, 800, 800, 800}).alleles(noCalls).make(),
                        new GenotypeBuilder("AA_A").PL(new int[]{30, 71, 73, 71, 73, 73, 71, 73, 73, 73, 20, 72, 72, 72, 10}).alleles(noCalls).make()).make()});

        // just spanning ref contexts, trying both instances where we want/do not want ref-only contexts
        tests.add(new Object[]{"test08",Arrays.asList(vcAA_ALT),
                loc, false, false,
                null});
        tests.add(new Object[]{"test09", Arrays.asList(vcAA_ALT),
                loc, true, false,
                new VariantContextBuilder(VCbase).alleles(Arrays.asList(Allele.create("A", true))).genotypes(new GenotypeBuilder("AA").PL(new int[]{0}).alleles(noCalls).make()).make()});

        // test uniquification of sample names
        tests.add(new Object[]{"test10",Arrays.asList(vcA_C, vcA_C_ALT), loc, false, true,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(
                        new GenotypeBuilder("A_C.test2").PL(new int[]{30, 20, 10}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C.test").PL(new int[]{30, 20, 10}).alleles(noCalls).make()).make()});

        tests.add(new Object[]{"test11",Arrays.asList(vcA_C_G, vcA_C_G_ALT), loc, false, true,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(
                        new GenotypeBuilder("A_C_G.test2").PL(new int[]{40, 20, 30, 20, 10, 30}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C_G.test").PL(new int[]{40, 20, 30, 20, 10, 30}).alleles(noCalls).make()).make()});

        // test creation of AD with proper allele indexing with and without PLs
        // Note: this test hard codes the filtering out of Allele.NON_REF_ALLELE alleles, so no AD value is expected
        final Genotype gA_ATC_ALT_AD_and_PLs = new GenotypeBuilder("A_ATC").PL(new int[]{30, 20, 10, 71, 72, 73}).AD(new int[]{20,10}).alleles(noCalls).make();
        final VariantContext vcA_ATC_ALT_AD_and_PLs = new VariantContextBuilder(VCbase).alleles(A_ATC_ALT).genotypes(gA_ATC_ALT_AD_and_PLs).make();
        final Genotype gA_C_G_ALT_AD_and_PLs = new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30, 71, 72, 73, 74}).AD(new int[]{30,0,8}).alleles(noCalls).make();
        final VariantContext vcA_C_G_ALT_AD_and_PLs = new VariantContextBuilder(VCbase).alleles(A_C_G_ALT).genotypes(gA_C_G_ALT_AD_and_PLs).make();
        final List<Allele> A_C_G_ATC = Arrays.asList(Aref, ATC, C, G);

        // 1 and 2 alt alleles (excluding Allele.NON_REF_ALLELEs) w no overlaps should give 4 AD values (1 ref, and 3 distinct alt alleles)
        tests.add(new Object[]{"test12",Arrays.asList(vcA_ATC_ALT_AD_and_PLs, vcA_C_G_ALT_AD_and_PLs), loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G_ATC).genotypes(
                        new GenotypeBuilder("A_ATC").AD(new int[]{20,10,0,0}).PL(new int[]{30,20,10,71,72,73,71,72,73,73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C_G").AD(new int[]{30,0,0,8}).PL(new int[]{40,71,74,20,72,30,20,73,10,30}).alleles(noCalls).make()).make()});


        final Genotype gA_C_G_ALT_noPLs = new GenotypeBuilder("A_C_G").AD(new int[]{60,9,20}).alleles(noCalls).make();
        final VariantContext vcA_C_G_ALT_noPLs = new VariantContextBuilder(VCbase).alleles(A_C_G_ALT).genotypes(gA_C_G_ALT_noPLs).make();

        // 2 and 2 alt alleles (excluding Allele.NON_REF_ALLELEs) w all overlaps should give 3 AD values (1 ref, and 2 distinct alt alleles)
        tests.add(new Object[]{"test13",Arrays.asList(vcA_C_G_ALT_noPLs, vcA_C_G_ALT_noPLs), loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(
                        new GenotypeBuilder("A_C_G").AD(new int[]{60,9,20}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C_G").AD(new int[]{60,9,20}).alleles(noCalls).make()).make()});

        final Genotype gA_ATC_AA_ALT_noPLs = new GenotypeBuilder("A_ATC_AA").AD(new int[]{30,8,40}).alleles(noCalls).make();
        final Allele AA = Allele.create("AA");
        final List<Allele> A_ATC_AA_ALT = Arrays.asList(Aref, ATC, AA, Allele.NON_REF_ALLELE);
        final VariantContext vcA_ATC_AA_ALT_noPLs = new VariantContextBuilder(VCbase).alleles(A_ATC_AA_ALT).genotypes(gA_ATC_AA_ALT_noPLs).make();
        final List<Allele> A_C_G_ATC_AA = Arrays.asList(Aref, C, G, ATC, AA);
        // 2 and 2 alt alleles (excluding Allele.NON_REF_ALLELEs) w no overlaps should give 5 AD values (1 ref, and 4 distinct alt alleles)
        tests.add(new Object[]{"test14",Arrays.asList(vcA_C_G_ALT_noPLs, vcA_ATC_AA_ALT_noPLs), loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G_ATC_AA).genotypes(
                        new GenotypeBuilder("A_C_G").AD(new int[]{60,9,20,0,0}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_ATC_AA").AD(new int[]{30,0,0,8,40}).alleles(noCalls).make()).make()});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider
    public Object[][] allelesToRemap(){
        VariantContextBuilder builder = new VariantContextBuilder().loc("1",1,1);
        final Allele CTC =  Allele.create("CTC");
        return new Object[][]{
                {builder.alleles(Arrays.asList(Aref, C, Allele.NON_REF_ALLELE)).make(), Aref,
                        Arrays.asList(Aref, C, Allele.NON_REF_ALLELE)},
                {builder.alleles(Arrays.asList(Aref, C)).make(), ATCref,
                        Arrays.asList(ATCref, CTC)},
                {builder.alleles(Arrays.asList(Aref, C, del)).make(), Aref,
                        Arrays.asList( Aref, C, del) },
                {builder.alleles(Arrays.asList(Aref, del)).make(), Aref,
                        Arrays.asList( Aref, del)},
        };
    }

    @Test(dataProvider = "allelesToRemap")
    public void testRemapAlleles(final VariantContext vc, final Allele refAllele, final List<Allele> expected){
        final List<Allele> result = ReferenceConfidenceVariantContextMerger.remapAlleles(vc, refAllele);
        Assert.assertEquals(result, expected, "result");
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testBadAlleleRemap(){
        VariantContext vc = new VariantContextBuilder().loc("1",1,1).alleles(Arrays.asList(ATCref, C)).make();
        ReferenceConfidenceVariantContextMerger.remapAlleles(vc, Aref);
    }

    @DataProvider
    public Object[][] getSpanningDeletionCases(){
        final VariantContext mixedVC = new VariantContextBuilder().loc("1", 2, 4)
                .alleles(Arrays.asList(ATCref, C, G, Allele.NON_REF_ALLELE))
                .make();

        final VariantContext spanningDelAllele = new VariantContextBuilder().loc("1", 2, 2)
                .alleles(Arrays.asList(Aref, Allele.SPAN_DEL))
                .make();

        final VariantContext deprecatedSpanningDel = new VariantContextBuilder().loc("1", 2, 2)
                .alleles(Arrays.asList(Aref, GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED))
                .make();

        final VariantContext notSpanning = new VariantContextBuilder().loc("1", 2,2).alleles(Arrays.asList(Aref, C)).make();

        final List<Allele> irrelevant = Collections.emptyList();

        return new Object[][] {
                { new ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(mixedVC, irrelevant, true), true },
                { new ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(mixedVC, irrelevant, false), false},
                { new ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(notSpanning, irrelevant, true), false},
                { new ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(notSpanning, irrelevant, false), false},
                { new  ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(spanningDelAllele, irrelevant, false), true},
                { new ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(deprecatedSpanningDel, irrelevant, false), true},
                { new  ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(spanningDelAllele, irrelevant, true), true},
                { new ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(deprecatedSpanningDel, irrelevant, true), true},
                { new  ReferenceConfidenceVariantContextMerger.VCWithNewAlleles(spanningDelAllele, irrelevant, false), true},
        };
    }

    @Test(dataProvider = "getSpanningDeletionCases")
    public void testVCWithNewAllelesIsSpanningDeletion(ReferenceConfidenceVariantContextMerger.VCWithNewAlleles vcWithNewAlleles, boolean expected){
        Assert.assertEquals(vcWithNewAlleles.isSpanningDeletion(), expected);
    }

    @Test
    public void testMedianCalculationOnMixedSerializedTypes() {
        // Merging attributes by median calculation requires sorting the values, which in turn requires a list
        // of values with homogeneous boxed representations. Make sure that FLOAT attributes with a serialized
        // representation that looks like an integer (with no decimal point, i.e. "0") get boxed into the same
        // type as other floating point values for that attribute to ensure successful sorting.
        final double medianRankSum = 1.46;
        final VCFHeader vcfHeader = new VCFHeader();
        vcfHeader.addMetaDataLine(GATKVCFHeaderLines.getInfoLine(MAP_QUAL_RANK_SUM_KEY));

        final VariantContextBuilder vcBuilder = new VariantContextBuilder("vc1", "20", 10, 10, Arrays.asList(Aref));

        // create 3 VCs with one each of a small value, the median value, and a large value for MQ_RankSum
        final List<VariantContext> toMergeVCs = new ArrayList<>(3);

        // use a literal string for this one to ensure that we have at least one test value that
        // has no embedded decimal point to emulate conditions found in the wild
        Map<String, Object> attributes1 = new HashMap<>();
        attributes1.put(MAP_QUAL_RANK_SUM_KEY, "0");
        toMergeVCs.add(vcBuilder.attributes(attributes1).make());

        Map<String, Object> attributes2 = new HashMap<>();
        attributes2.put(MAP_QUAL_RANK_SUM_KEY, Double.toString(medianRankSum));
        toMergeVCs.add(vcBuilder.attributes(attributes2).make());

        Map<String, Object> attributes3 = new HashMap<>();
        attributes3.put(MAP_QUAL_RANK_SUM_KEY, "2.46");
        toMergeVCs.add(vcBuilder.attributes(attributes3).make());

        // merge and make sure we get the median value
        final ReferenceConfidenceVariantContextMerger merger = new ReferenceConfidenceVariantContextMerger(getAnnotationEngine(), vcfHeader);
        final VariantContext mergedVC = merger.merge(
                toMergeVCs,
                new SimpleInterval("20", 10, 10),
                (byte) 'A',
                true,
                false);

        Assert.assertEquals(mergedVC.getAttributeAsDouble(MAP_QUAL_RANK_SUM_KEY,-1.0), medianRankSum);
    }
}
