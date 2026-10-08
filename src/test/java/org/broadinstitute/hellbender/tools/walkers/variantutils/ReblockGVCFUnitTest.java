package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.Coverage;
import org.broadinstitute.hellbender.tools.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_ReadPosRankSumTest;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriterUnitTest;
import org.broadinstitute.hellbender.utils.variant.writers.MockVcfWriter;
import org.broadinstitute.hellbender.utils.variant.writers.ReblockingGVCFWriter;
import org.broadinstitute.hellbender.utils.variant.writers.ReblockingOptions;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class ReblockGVCFUnitTest extends CommandLineProgramTest {
    private final static Allele LONG_REF = Allele.create("ACTG", true);
    private final static Allele DELETION = Allele.create("A", false);
    private final static Allele SHORT_REF = Allele.create("A", true);
    private final static Allele LONG_SNP = Allele.create("TCTA", false);
    private final static Allele SHORT_INS = Allele.create("AT", false);
    private final static Allele LONG_INS = Allele.create("ATT", false);
    private final static int EXAMPLE_DP = 18;
    public static final int DEFAULT_START = 10;

    @Test
    public void testCleanUpHighQualityVariant() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        //We need an annotation engine for cleanUpHighQualityVariant(), but this is just a dummy; annotations won't initialize properly without runCommandLine()
        reblocker.createAnnotationEngine();
        //...and a vcfwriter
        reblocker.vcfWriter = new ReblockingGVCFWriter(new MockVcfWriter(), Arrays.asList(20, 100), true, null, new ReblockingOptions());
        reblocker.dropLowQuals = true;
        reblocker.doQualApprox = true;

        final Genotype g0 = VariantContextTestUtils.makeG("sample1", LONG_REF, DELETION, 41, 0, 37, 200, 100, 200, 400, 600, 800, 1200);
        final Genotype g = addAD(g0,13,17,0,0);
        final VariantContext extraAlt0 = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, LONG_SNP, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final Map<String, Object> attr = new HashMap<>();
        attr.put(VCFConstants.DEPTH_KEY, 32);
        final VariantContext extraAlt = addAttributes(extraAlt0, attr);
        //we'll call this with the same VC again under the assumption that STAND_CALL_CONF is zero so no alleles/GTs change
        final VariantContext cleaned1 = reblocker.cleanUpHighQualityVariant(extraAlt);
        Assert.assertEquals(cleaned1.getAlleles().size(), 3);
        Assert.assertTrue(cleaned1.getAlleles().contains(LONG_REF));
        Assert.assertTrue(cleaned1.getAlleles().contains(DELETION));
        Assert.assertTrue(cleaned1.getAlleles().contains(Allele.NON_REF_ALLELE));
        Assert.assertTrue(cleaned1.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        Assert.assertEquals(cleaned1.getAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY), 41);
        Assert.assertTrue(cleaned1.hasAttribute(GATKVCFConstants.VARIANT_DEPTH_KEY));
        Assert.assertEquals(cleaned1.getAttribute(GATKVCFConstants.VARIANT_DEPTH_KEY), 30);
        Assert.assertTrue(cleaned1.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        Assert.assertEquals(cleaned1.getAttributeAsString(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,"").split(",")[1], "32");

        final Genotype hetNonRef = VariantContextTestUtils.makeG("sample2", DELETION, LONG_SNP, 891,879,1128,84,0,30,891,879,84,891);
        final VariantContext keepAlts = makeDeletionVC("keepAllAlts", Arrays.asList(LONG_REF, DELETION, LONG_SNP, Allele.NON_REF_ALLELE), LONG_REF.length(), hetNonRef);
        final VariantContext cleaned2 = reblocker.cleanUpHighQualityVariant(keepAlts);
        Assert.assertEquals(cleaned2.getAlleles().size(), 4);
        Assert.assertTrue(cleaned2.getAlleles().contains(LONG_REF));
        Assert.assertTrue(cleaned2.getAlleles().contains(DELETION));
        Assert.assertTrue(cleaned2.getAlleles().contains(LONG_SNP));
        Assert.assertTrue(cleaned2.getAlleles().contains(Allele.NON_REF_ALLELE));

        //if a "high quality" variant has a called * that gets dropped, it might turn into a low quality variant
        reblocker.dropLowQuals = false;
        final Genotype withStar = VariantContextTestUtils.makeG("sample3", SHORT_REF, Allele.SPAN_DEL, 502,65,98,347,0,370,404,129,341,535,398,123,335,462,456);
        final VariantContext trickyLowQual = makeDeletionVC("withStar", Arrays.asList(SHORT_REF, Allele.SPAN_DEL, SHORT_INS, LONG_INS, Allele.NON_REF_ALLELE), SHORT_REF.length(), withStar);
        final VariantContext cleaned3 = reblocker.cleanUpHighQualityVariant(trickyLowQual);
        Assert.assertNull(cleaned3);
    }

    @Test
    public void testLowQualVariantToGQ0HomRef() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        reblocker.vcfWriter = new ReblockingGVCFWriter(new MockVcfWriter(), Arrays.asList(20, 100), true, null, new ReblockingOptions());

        reblocker.dropLowQuals = true;
        final Genotype g = VariantContextTestUtils.makeG("sample1", 11, LONG_REF, Allele.NON_REF_ALLELE, 200, 100, 200, 11, 0, 37);
        final VariantContext toBeNoCalled = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final VariantContextBuilder dropped = reblocker.lowQualVariantToGQ0HomRef(toBeNoCalled);
        Assert.assertNull(dropped);

        reblocker.dropLowQuals = false;
        final VariantContext modified = reblocker.lowQualVariantToGQ0HomRef(toBeNoCalled).make();
        Assert.assertTrue(modified.getAttributes().containsKey(VCFConstants.END_KEY));
        Assert.assertEquals(modified.getAttributes().get(VCFConstants.END_KEY), 13);
        Assert.assertEquals(modified.getReference(), SHORT_REF);
        Assert.assertEquals(modified.getAlternateAllele(0), Allele.NON_REF_ALLELE);
        Assert.assertTrue(!modified.filtersWereApplied());
        Assert.assertEquals(modified.getLog10PError(), VariantContext.NO_LOG10_PERROR);

        final Genotype longPls = VariantContextTestUtils.makeG("sample1", Allele.NO_CALL, Allele.NO_CALL, 0,0,0,0,0,0);
        final VariantContext lotsOfZeroPls = makeNoDepthVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), longPls);
        final VariantContext properlySubset = reblocker.lowQualVariantToGQ0HomRef(lotsOfZeroPls).make();
        Assert.assertEquals(properlySubset.getGenotype(0).getPL().length, 3);
        Assert.assertEquals(properlySubset.getGenotype(0).getPL(), new int[]{0,0,0});

        //No-calls were throwing NPEs.  Now they're not.
        final Genotype g2 = VariantContextTestUtils.makeG("sample1", Allele.NO_CALL,Allele.NO_CALL);
        final VariantContext noData = makeDeletionVC("noData", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g2);
        final VariantContext notCrashing = reblocker.lowQualVariantToGQ0HomRef(noData).make();
        final Genotype outGenotype = notCrashing.getGenotype(0);
        Assert.assertTrue(outGenotype.isHomRef());
        Assert.assertEquals(outGenotype.getGQ(), 0);
        Assert.assertTrue(Arrays.stream(outGenotype.getPL()).allMatch(x -> x == 0));
        Assert.assertTrue(notCrashing.getGenotype(0).isHomRef());  //get rid of no-calls -- GQ0 hom ref instead

        //haploid hom ref call
        final int[] pls = {0, 35, 72};
        final int gq = 35;
        final GenotypeBuilder gb = new GenotypeBuilder("male_sample", Collections.singletonList(LONG_REF)).PL(pls).GQ(gq);
        final VariantContextBuilder vb = new VariantContextBuilder();
        vb.chr("20").start(10001).stop(10004).alleles(Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE)).log10PError(-3.0).genotypes(gb.make());
        final VariantContext vc = vb.make();

        final VariantContext haploidRefBlock = reblocker.lowQualVariantToGQ0HomRef(vc).make();
        final Genotype newG = haploidRefBlock.getGenotype("male_sample");

        Assert.assertEquals(newG.getPloidy(), 1);
        Assert.assertTrue(newG.hasGQ());
        Assert.assertEquals(newG.getGQ(), 35);
    }

    @Test
    public void testCalledHomRefGetsAltGQ() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final Genotype g3 = VariantContextTestUtils.makeG("sample1", 11, LONG_REF, LONG_REF, 0, 11, 37, 100, 200, 400);
        final VariantContext twoAltsHomRef = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g3);
        final GenotypeBuilder takeGoodAltGQ = reblocker.changeCallToHomRefVersusNonRef(twoAltsHomRef, new HashMap<>());
        final Genotype nowRefBlock = takeGoodAltGQ.make();
        Assert.assertTrue(nowRefBlock.hasGQ());
        Assert.assertEquals(nowRefBlock.getGQ(), 11);
        Assert.assertEquals(nowRefBlock.getDP(), 18);
        Assert.assertEquals((int)nowRefBlock.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY), 18);
    }

    @Test
    public void testChangeCallToGQ0HomRef() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final Genotype g = VariantContextTestUtils.makeG("sample1", LONG_REF, Allele.NON_REF_ALLELE, 200, 100, 200, 11, 0, 37);
        final VariantContext toBeNoCalled = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final Map<String, Object> noAttributesMap = new HashMap<>();
        final GenotypeBuilder noCalled = reblocker.changeCallToHomRefVersusNonRef(toBeNoCalled, noAttributesMap);
        final Genotype newG = noCalled.make();
        Assert.assertTrue(noAttributesMap.containsKey(VCFConstants.END_KEY));
        Assert.assertEquals(noAttributesMap.get(VCFConstants.END_KEY), 13);
        Assert.assertEquals(newG.getAllele(0), SHORT_REF);
        Assert.assertEquals(newG.getAllele(1), SHORT_REF);
        Assert.assertTrue(!newG.hasAD());
    }

    @Test  //no-calls can be dropped or reblocked just like hom-refs, i.e. we don't have to preserve them like variants
    public void testBadCalls() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final Genotype g2 = VariantContextTestUtils.makeG("sample1", Allele.NO_CALL,Allele.NO_CALL);
        final VariantContext noData = makeDeletionVC("noData", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g2);
        Assert.assertTrue(reblocker.shouldBeReblocked(noData));

        final Genotype g3 = VariantContextTestUtils.makeG("sample1", LONG_REF, Allele.NON_REF_ALLELE);
        final VariantContext nonRefCall = makeDeletionVC("nonRefCall", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g3);
        Assert.assertTrue(reblocker.shouldBeReblocked(nonRefCall));
    }

    @Test
    public void testPosteriors() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        reblocker.vcfWriter = new ReblockingGVCFWriter(new MockVcfWriter(), Arrays.asList(20, 100), true, null, new ReblockingOptions());
        reblocker.posteriorsKey = "GP";

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(LONG_REF, LONG_REF));
        final double[] posteriors = {0,2,5.01,2,4,5.01,2,4,4,5.01};
        final int[] pls = {2,0,50,2,289,52,38,325,407,88};
        gb.attribute("GP", posteriors).PL(pls).GQ(-2); //an older version of GATK-DRAGEN output a -2 GQ
        final VariantContext vc = makeDeletionVC("DRAGEN", Arrays.asList(LONG_REF, DELETION, LONG_SNP, Allele.NON_REF_ALLELE), LONG_REF.length(), gb.make());
        Assert.assertTrue(reblocker.shouldBeReblocked(vc));
        final VariantContext out = reblocker.lowQualVariantToGQ0HomRef(vc).make();
        Assert.assertTrue(out.getGenotype(0).isHomRef());
        Assert.assertEquals(out.getGenotype(0).getGQ(), 2);

        final GenotypeBuilder gb2 = new GenotypeBuilder("sample1", Arrays.asList(LONG_REF, LONG_REF));
        final double gqForNoPLs = 34.77;
        final int inputGQ = 32;
        Assert.assertNotEquals(gqForNoPLs, inputGQ);
        final double[] posteriors2 = {0,gqForNoPLs,37.78,39.03,73.8,42.04};
        gb2.attribute("GP", posteriors2).GQ(inputGQ);
        final VariantContext vc2 = makeDeletionVC("DRAGEN", Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), gb2.make());
        final VariantContext out2 = reblocker.lowQualVariantToGQ0HomRef(vc2).make();
        final Genotype gOut2 = out2.getGenotype(0);
        Assert.assertTrue(gOut2.isHomRef());
        Assert.assertEquals(gOut2.getGQ(), (int)Math.round(gqForNoPLs));
        Assert.assertTrue(gOut2.hasPL());
        Assert.assertEquals(gOut2.getPL().length, 3);
    }

    @DataProvider(name = "overlappingDeletionCases")
    public Object[][] createOverlappingDeletionCases() {
        return new Object[][] {
                //{100000, 10, 100005, 10, 99, 99, 2},
                //{100000, 10, 100005, 10, 99, 5, 2},
                //{100000, 10, 100005, 10, 5, 99, 2},
                //{100000, 10, 100005, 10, 5, 5, 1},
                //{100000, 15, 100010, 5, 99, 99, 2},
                {100000, 15, 100010, 5, 99, 5, 1},
                {100000, 15, 100005, 5, 5, 99, 3},
                {100000, 15, 100005, 5, 5, 5, 1}
        };
    }

    @Test(dataProvider = "overlappingDeletionCases")
    public void testOverlappingDeletions(final int del1start, final int del1length,
                                         final int del2start, final int del2length,
                                         final int del1qual, final int del2qual, final int numExpected) throws IOException {
        final String inputPrefix = "overlappingDeletions";
        final String inputSuffix = ".g.vcf";
        final File inputFile = File.createTempFile(inputPrefix, inputSuffix);
        final GVCFWriter gvcfWriter= setUpWriter(inputFile, new File(GATKBaseTest.FULL_HG19_DICT));

        final ReferenceSequenceFile ref = ReferenceUtils.createReferenceReader(new GATKPath(GATKBaseTest.b37Reference));
        final Allele del1Ref = Allele.create(ReferenceUtils.getRefBasesAtPosition(ref, "20", del1start, del1length), true);
        final Allele del1Alt = Allele.create(ReferenceUtils.getRefBaseAtPosition(ref, "20", del1start), false);
        final Allele del2Ref = Allele.create(ReferenceUtils.getRefBasesAtPosition(ref, "20", del2start, del2length), true);
        final Allele del2Alt = Allele.create(ReferenceUtils.getRefBaseAtPosition(ref, "20", del2start), false);
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.chr("20").start(del1start).stop(del1start+del1length-1).attribute(VCFConstants.DEPTH_KEY, 10).alleles(Arrays.asList(del1Ref, del1Alt, Allele.NON_REF_ALLELE));
        final VariantContext del1 = VariantContextTestUtils.makeGVCFVariantContext(variantContextBuilder, Arrays.asList(del1Ref, del1Alt), del1qual);

        variantContextBuilder.chr("20").start(del2start).stop(del2start+del2length-1).attribute(VCFConstants.DEPTH_KEY, 10).alleles(Arrays.asList(del2Ref, del2Alt, Allele.NON_REF_ALLELE));
        final VariantContext del2 = VariantContextTestUtils.makeGVCFVariantContext(variantContextBuilder, Arrays.asList(del2Ref, del2Alt), del2qual);

        gvcfWriter.add(del1);
        gvcfWriter.add(del2);
        gvcfWriter.close();

        final File outputFile = File.createTempFile(inputPrefix,".reblocked" + inputSuffix);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", inputFile)
            .add(ReblockGVCF.RGQ_THRESHOLD_SHORT_NAME, 10.0)
            .addReference(b37_reference_20_21)
            .addOutput(outputFile);
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        Assert.assertEquals(outVCs.getRight().size(), numExpected);
    }

    @Test
    public void testComplicatedOverlaps() throws IOException {
        final String inputPrefix = "overlappingDeletions";
        final String inputSuffix = ".g.vcf";
        final File inputFile = File.createTempFile(inputPrefix, inputSuffix);
        final GVCFWriter gvcfWriter= setUpWriter(inputFile, new File(GATKBaseTest.FULL_HG19_DICT));

        final Genotype withStar = VariantContextTestUtils.makeG("sample3", SHORT_REF, Allele.SPAN_DEL, 502,65,98,347,0,370,404,129,341,535,398,123,335,462,456);
        final VariantContext trickyLowQual = makeDeletionVC("withStar", Arrays.asList(SHORT_REF, Allele.SPAN_DEL, SHORT_INS, LONG_INS, Allele.NON_REF_ALLELE), SHORT_REF.length(), withStar);

        gvcfWriter.add(trickyLowQual);
        gvcfWriter.close();

        final File outputFile = File.createTempFile(inputPrefix,".reblocked" + inputSuffix);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", inputFile)
                .add(ReblockGVCF.RGQ_THRESHOLD_SHORT_NAME, 10.0)
                .addReference(b37_reference_20_21)
                .addOutput(outputFile);
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        Assert.assertEquals(outVCs.getRight().size(), 1);
        final VariantContext outputVC = outVCs.getRight().get(0);
        Assert.assertTrue(outputVC.getGenotype(0).isHomRef());
        Assert.assertTrue(outputVC.hasAttribute(VCFConstants.END_KEY));
        Assert.assertEquals(outputVC.getAttributeAsInt(VCFConstants.END_KEY, 0), DEFAULT_START);

        //
    }

    @Test
    public void testIndelTrimming() throws IOException {
        final String inputPrefix = "altToTrim";
        final String inputSuffix = ".g.vcf";
        final File inputFile = new File(inputPrefix+inputSuffix);
        final GVCFWriter gvcfWriter= setUpWriter(inputFile, new File(GATKBaseTest.FULL_HG19_DICT));

        final int longestDelLength = 20;
        final int del1start = 200001;  //chr20:200001 is a 'g'
        final int goodDelLength = 10;
        final ReferenceSequenceFile ref = ReferenceUtils.createReferenceReader(new GATKPath(GATKBaseTest.b37Reference));
        final Allele del1Ref = Allele.create(ReferenceUtils.getRefBasesAtPosition(ref, "20", del1start, longestDelLength), true);  //+1 for anchor base
        final Allele del1Alt2 = Allele.create(ReferenceUtils.getRefBaseAtPosition(ref, "20", del1start), false);
        final Allele del1Alt1 = Allele.extend(Allele.create(del1Alt2, true), ReferenceUtils.getRefBasesAtPosition(ref, "20", del1start+goodDelLength, goodDelLength));
        //add a SNP that lexicagraphically precedes the (shorter) del alleles
        final byte[] snpBases = ReferenceUtils.getRefBasesAtPosition(ref, "20", del1start, longestDelLength);
        snpBases[0] = (byte)'a';
        final Allele del1Alt3 = Allele.create(snpBases, false);
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.chr("20").start(del1start).stop(del1start+longestDelLength-1).attribute(VCFConstants.DEPTH_KEY, 10)
                .alleles(Arrays.asList(del1Ref, del1Alt1, del1Alt2, del1Alt3, Allele.NON_REF_ALLELE));
        final GenotypeBuilder gb = new GenotypeBuilder(VariantContextTestUtils.SAMPLE_NAME, Arrays.asList(del1Ref, del1Alt1));
        gb.PL(new int[]{50, 0, 100, 150, 200, 300, 400, 500, 600, 1000, 5000, 5000, 5000, 5000, 5000});
        variantContextBuilder.genotypes(gb.make());
        final VariantContext del1 = variantContextBuilder.make();

        final int goodStart = del1start + longestDelLength;
        final Allele goodRef = Allele.create(ReferenceUtils.getRefBasesAtPosition(ref, "20", goodStart, 1), true);
        final Allele goodSNP = VariantContextTestUtils.makeAnySNPAlt(goodRef);  //generate valid data, but make it agnostic to position and reference genome
        variantContextBuilder.start(goodStart).stop(goodStart).alleles(Arrays.asList(goodRef, goodSNP, Allele.NON_REF_ALLELE));
        final GenotypeBuilder gb2 = new GenotypeBuilder(VariantContextTestUtils.SAMPLE_NAME, Arrays.asList(goodRef, goodSNP));
        gb2.PL(new int[]{50, 0, 100, 150, 200, 300});
        gb2.GQ(50);
        variantContextBuilder.genotypes(gb2.make());
        final VariantContext keepVar = variantContextBuilder.make();

        gvcfWriter.add(del1);
        gvcfWriter.add(keepVar);
        gvcfWriter.close();

        final File outputFile = File.createTempFile(inputPrefix,".reblocked" + inputSuffix);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", inputFile)
                .addReference(b37_reference_20_21)
                .addOutput(outputFile);
        runCommandLine(args);

        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath()).getRight();
        Assert.assertEquals(outVCs.size(), 3);
        //make sure vc0 had uncalled del allele dropped and remaining allele trimmed
        final VariantContext vc0 = outVCs.get(0);
        Assert.assertTrue(vc0.isVariant());
        Assert.assertEquals(vc0.getAlleles().size(), 3);
        final Genotype g0 = vc0.getGenotype(0);
        Assert.assertTrue(g0.getAllele(0).isReference());
        Assert.assertEquals(g0.getAllele(0).length(), 10);
        Assert.assertFalse(g0.getAllele(1).isReference());
        Assert.assertEquals(g0.getAllele(1).length(), 1);
        final VariantContext vc1 = outVCs.get(1);
        Assert.assertTrue(vc1.isReferenceBlock());
        Assert.assertEquals(vc1.getStart(), del1start + goodDelLength);
        Assert.assertEquals(vc1.getGenotype(0).getLikelihoods().getAsPLs()[1], 100);  //should take ref block likelihoods from del not SNP
        Assert.assertTrue(outVCs.get(2).isVariant());
    }

    @Test
    public void testAnnotationSubsetting() {
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Arrays.asList(new Coverage(), new AS_RMSMappingQuality(),
                new RMSMappingQuality(), new AS_ReadPosRankSumTest(), new AS_FisherStrand()), null, Collections.emptyList(), false, false);

        final Genotype g0 = VariantContextTestUtils.makeG("sample1", LONG_REF, LONG_SNP, 41, 0, 37, 200, 100, 200, 400, 600, 800, 1200, 5000, 5000, 5000, 5000, 5000);
        final Genotype g = addAD(g0,13,17,0,0,1);

        final VariantContext originalVCbase = makeDeletionVC("", Arrays.asList(LONG_REF, LONG_SNP, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), g);
        final VariantContextBuilder originalBuilder = new VariantContextBuilder(originalVCbase);
        final Map<String, Object> origAttributes = new LinkedHashMap<>();
        origAttributes.put(VCFConstants.DEPTH_KEY, 93);
        origAttributes.put(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, "329810.0");
        origAttributes.put(GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY, "|-2.1,1|NaN|NaN");
        origAttributes.put(GATKVCFConstants.AS_SB_TABLE_KEY, "17,18|8,5|0,0|0,1");
        origAttributes.put(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "123769.00|46800.00|0.00|0.00");
        originalBuilder.attributes(origAttributes);
        final VariantContext originalVC = originalBuilder.make();
        final Genotype newG = VariantContextTestUtils.makeG("sample1", LONG_REF, LONG_SNP, 41, 0, 37, 200, 100, 200);
        final VariantContext regenotypedVC = makeDeletionVC("", Arrays.asList(LONG_REF, LONG_SNP, Allele.NON_REF_ALLELE), LONG_REF.length(), newG);

        final Map<String, Object> subsetAnnotations = ReblockGVCF.subsetAnnotationsIfNecessary(annotationEngine, true, null, originalVC, regenotypedVC, Collections.emptyList());
        Assert.assertTrue(subsetAnnotations.containsKey(VCFConstants.DEPTH_KEY));
        Assert.assertEquals(subsetAnnotations.get(VCFConstants.DEPTH_KEY), 93);

        Assert.assertEquals(subsetAnnotations.get(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY), "329810,93");
        Assert.assertEquals(subsetAnnotations.get(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED), 329810.0);
        Assert.assertEquals(subsetAnnotations.get(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED), 93);

        Assert.assertEquals(subsetAnnotations.get(GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY), "|-2.1,1|");
        Assert.assertEquals(subsetAnnotations.get(GATKVCFConstants.AS_SB_TABLE_KEY), "17,18|8,5|0,0");
        Assert.assertEquals(subsetAnnotations.get(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY), "123769.00|46800.00|0.00");
    }

    private GVCFWriter setUpWriter(final File outputFile, final File dictionary) throws IOException {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputPath(outputFile.toPath());
        final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(dictionary.toPath());
        builder.setReferenceDictionary(dict);
        final VariantContextWriter vcfWriter = builder.build();
        final GVCFWriter gvcfWriter= new GVCFWriter(vcfWriter, Arrays.asList(20,100), true);
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(VariantContextTestUtils.SAMPLE_NAME));
        result.setSequenceDictionary(dict);
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.String,  "genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS, VCFHeaderLineCount.R,
                VCFHeaderLineType.Integer, "Allele depth"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1,
                VCFHeaderLineType.Integer, " depth"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1,
                VCFHeaderLineType.Integer, " depth"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1,
                VCFHeaderLineType.Integer, "Genotype quality"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_PL_KEY, VCFHeaderLineCount.G,
                VCFHeaderLineType.Integer, "Phred-scaled likelihoods"));
        gvcfWriter.writeHeader(result);
        return gvcfWriter;
    }

    @Test
    public void testLowQualityAfterSubsetting() {

    }

    private static final Allele CHR_M_REF = Allele.create("A", true);
    private static final Allele CHR_M_ALT = Allele.create("G", false);

    /**
     * Wire up a reblocker with a MockVcfWriter we can read back, using the same GQ bands
     * (20, 100) and floorBlocks setting as the other unit tests in this class.
     */
    private static MockVcfWriter attachMockWriter(final ReblockGVCF reblocker) {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        reblocker.createAnnotationEngine();
        reblocker.vcfWriter = new ReblockingGVCFWriter(mockWriter, Arrays.asList(20, 100), true, null, new ReblockingOptions());
        return mockWriter;
    }

    /**
     * A DRAGEN-style somatic reference block: SQ, but no GQ and no PL.
     */
    private static VariantContext makeSomaticRefBlock(final int start, final int end, final String sq, final int minDP) {
        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        gb.attribute(GATKVCFConstants.SOMATIC_QUALITY_KEY, sq)
          .attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, minDP)
          .AD(new int[]{76, 0})
          .DP(76)
          .noGQ()
          .noPL();
        return new VariantContextBuilder("test", "chrM", start, end,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, end)
                .genotypes(gb.make()).unfiltered().make();
    }

    /**
     * A DRAGEN-style somatic ALT call: SQ, but no GQ and no PL.
     */
    private static VariantContext makeSomaticVariant(final int position) {
        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_ALT, CHR_M_ALT));
        gb.attribute(GATKVCFConstants.SOMATIC_QUALITY_KEY, "97.96,0.00")
          .AD(new int[]{0, 88, 0})
          .DP(88)
          .noGQ()
          .noPL();
        return new VariantContextBuilder("test", "chrM", position, position,
                Arrays.asList(CHR_M_REF, CHR_M_ALT, Allele.NON_REF_ALLELE))
                .genotypes(gb.make()).unfiltered().make();
    }

    /**
     * Somatic-style mitochondrial records from DRAGEN carry SQ instead of GQ/PL. Both the reference
     * blocks and the ALT calls should reach the output with their genotypes untouched -- in particular
     * still carrying SQ and still lacking GQ, rather than being banded or having quality synthesized.
     *
     * Note the submission order: the reference block covers 2-72 and so must be submitted before the
     * call at 73. Submitting the call first advances the writer's output end past the block, and the
     * block is then correctly discarded as already-covered, which is a different code path.
     */
    @Test
    public void testSomaticSQRecordsPassThroughUntouched() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        reblocker.regenotypeVC(makeSomaticRefBlock(2, 72, "10", 45));
        reblocker.regenotypeVC(makeSomaticVariant(73));
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 2, "both the somatic ref block and the somatic call should be emitted");

        final VariantContext outRefBlock = emitted.get(0);
        Assert.assertEquals(outRefBlock.getStart(), 2, "ref block start should be unchanged");
        Assert.assertEquals(outRefBlock.getEnd(), 72, "ref block end should be unchanged");
        final Genotype refBlockGenotype = outRefBlock.getGenotype(0);
        Assert.assertEquals(refBlockGenotype.getExtendedAttribute(GATKVCFConstants.SOMATIC_QUALITY_KEY), "10", "SQ should survive on the ref block");
        Assert.assertFalse(refBlockGenotype.hasGQ(), "no GQ should be synthesized for the somatic ref block");
        Assert.assertFalse(refBlockGenotype.hasPL(), "no PL should be synthesized for the somatic ref block");

        final VariantContext outVariant = emitted.get(1);
        Assert.assertEquals(outVariant.getStart(), 73, "call should be emitted at its original position");
        final Genotype variantGenotype = outVariant.getGenotype(0);
        Assert.assertEquals(variantGenotype.getExtendedAttribute(GATKVCFConstants.SOMATIC_QUALITY_KEY), "97.96,0.00", "SQ should survive on the call");
        Assert.assertFalse(variantGenotype.hasGQ(), "no GQ should be synthesized for the somatic call");
        Assert.assertFalse(variantGenotype.hasPL(), "no PL should be synthesized for the somatic call");
    }

    /**
     * Adjacent somatic reference blocks must not be merged into a single banded block the way ordinary
     * hom-ref blocks are -- there is no GQ to band on, and their SQ values differ. This is what exercises
     * the somatic branch in ReblockingGVCFBlockCombiner.addHomRefSite: without it these blocks reach the
     * "hom-ref genotypes must contain GQ or PL" check and the tool fails on valid DRAGEN input.
     */
    @Test
    public void testAdjacentSomaticRefBlocksAreNotBanded() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        reblocker.regenotypeVC(makeSomaticRefBlock(2, 72, "10", 45));
        reblocker.regenotypeVC(makeSomaticRefBlock(73, 126, "6", 82));
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 2, "adjacent somatic ref blocks should stay separate, not be banded together");
        Assert.assertEquals(emitted.get(0).getStart(), 2);
        Assert.assertEquals(emitted.get(0).getEnd(), 72);
        Assert.assertEquals(emitted.get(1).getStart(), 73);
        Assert.assertEquals(emitted.get(1).getEnd(), 126);
        Assert.assertEquals(emitted.get(0).getGenotype(0).getExtendedAttribute(GATKVCFConstants.SOMATIC_QUALITY_KEY), "10");
        Assert.assertEquals(emitted.get(1).getGenotype(0).getExtendedAttribute(GATKVCFConstants.SOMATIC_QUALITY_KEY), "6");
    }

    /**
     * A record with both SQ and GQ/PL should NOT take the somatic passthrough path; it should be
     * reblocked normally, as seen in some mixed DRAGEN outputs.
     *
     * The record here is a low-quality variant, chosen deliberately: for a hom-ref block the two paths
     * both hand the same VariantContext to the writer and so cannot be told apart from the output. A
     * low-quality variant diverges -- the normal path collapses it to a GQ0 hom-ref block carrying only
     * <NON_REF>, whereas passthrough would emit it unchanged with its real ALT allele still attached.
     * That surviving ALT is what this test watches for.
     */
    @Test
    public void testRecordWithSQAndGQIsNotPassedThrough() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);
        reblocker.dropLowQuals = false;

        final Genotype lowQualG = VariantContextTestUtils.makeG("sample1", 11, LONG_REF, Allele.NON_REF_ALLELE, 200, 100, 200, 11, 0, 37);
        final Genotype withSQ = new GenotypeBuilder(lowQualG).attribute(GATKVCFConstants.SOMATIC_QUALITY_KEY, "10").make();
        final VariantContext lowQualVariantWithSQ = makeDeletionVC("lowQualVarWithSQ",
                Arrays.asList(LONG_REF, DELETION, Allele.NON_REF_ALLELE), LONG_REF.length(), withSQ);

        // sanity check on the fixture: the input really does carry SQ alongside GQ and PL
        Assert.assertTrue(lowQualVariantWithSQ.getGenotype(0).hasExtendedAttribute(GATKVCFConstants.SOMATIC_QUALITY_KEY));
        Assert.assertTrue(lowQualVariantWithSQ.getGenotype(0).hasGQ());
        Assert.assertTrue(lowQualVariantWithSQ.getGenotype(0).hasPL());

        reblocker.regenotypeVC(lowQualVariantWithSQ);
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 1, "the record should be emitted once");
        Assert.assertEquals(emitted.get(0).getAlternateAlleles(), Collections.singletonList(Allele.NON_REF_ALLELE),
                "the low quality variant should have been reblocked to a hom-ref block; a surviving ALT allele means "
                        + "it was passed through as somatic instead");
        Assert.assertFalse(emitted.get(0).hasAllele(DELETION), "the real ALT allele should not survive reblocking");
    }

    /**
     * A somatic record bypasses the block machinery in ReblockingGVCFBlockCombiner.addHomRefSite and is emitted
     * directly. Any band left open by preceding ordinary hom-ref records covers earlier positions, so it has to be
     * closed out first -- otherwise it is flushed at end-of-input and lands *after* the somatic record, leaving the
     * output out of position order.
     *
     * The blocks here are deliberately contiguous (2-50 then 51-100) so that nothing else triggers a flush in
     * between: the only thing that can emit the first band before the second record is the somatic branch itself.
     */
    @Test
    public void testOpenBandIsFlushedBeforeSomaticRecordIsEmitted() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        // an ordinary hom-ref block, which opens a band
        final GenotypeBuilder normalGB = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        normalGB.GQ(40).PL(new int[]{0, 40, 400}).AD(new int[]{50, 0}).DP(50);
        final VariantContext normalBlock = new VariantContextBuilder("test", "chrM", 2, 50,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, 50)
                .genotypes(normalGB.make()).unfiltered().make();

        reblocker.regenotypeVC(normalBlock);
        reblocker.regenotypeVC(makeSomaticRefBlock(51, 100, "6", 82));
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 2, "both blocks should reach the output");
        Assert.assertEquals(emitted.get(0).getStart(), 2,
                "the banded hom-ref block must be emitted first; emitting the somatic record ahead of it puts the "
                        + "output out of position order");
        Assert.assertEquals(emitted.get(1).getStart(), 51, "the somatic block should follow the band it came after");

        // and the records should still be the ones we expect, not swapped in content
        Assert.assertTrue(emitted.get(0).getGenotype(0).hasGQ(), "the first record should be the banded one");
        Assert.assertEquals(emitted.get(1).getGenotype(0).getExtendedAttribute(GATKVCFConstants.SOMATIC_QUALITY_KEY), "6",
                "the second record should be the somatic one");
    }

    /**
     * The somatic guard in ReblockingGVCFBlockCombiner.addHomRefSite carries the same "no GQ and no PL"
     * requirement as the one in regenotypeVC, and needs its own coverage: a hom-ref block carrying SQ
     * alongside a real GQ must still be banded by the combiner rather than passed straight through.
     *
     * The discriminator is the GQ. Banding with floorBlocks set floors it to the band minimum of 20;
     * a passthrough would leave the original 42 intact.
     */
    @Test
    public void testCombinerBandsHomRefBlockThatHasBothSQAndGQ() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        gb.attribute(GATKVCFConstants.SOMATIC_QUALITY_KEY, "10")
          .GQ(42)
          .PL(new int[]{0, 42, 420})
          .AD(new int[]{76, 0})
          .DP(76);
        final VariantContext sqAndGQBlock = new VariantContextBuilder("test", "chrM", 2, 72,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, 72)
                .genotypes(gb.make()).unfiltered().make();

        reblocker.regenotypeVC(sqAndGQBlock);
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 1, "the block should be emitted once");
        final Genotype emittedGenotype = emitted.get(0).getGenotype(0);
        Assert.assertTrue(emittedGenotype.hasGQ(), "a block with GQ should keep a GQ");
        Assert.assertEquals(emittedGenotype.getGQ(), 20,
                "GQ should be floored to the band minimum, proving the combiner banded the block rather than "
                        + "treating it as somatic and passing it through");
    }

    /**
     * DRAGEN's targeted caller emits records with no <NON_REF> allele, which would otherwise trip the
     * "not a GVCF" UserException. Flagged targeted records are dropped instead of failing the run.
     *
     * This pins the current drop behavior deliberately -- see the CAVEAT comment on the bypass in
     * ReblockGVCF#apply about dropping not being equivalent to asserting hom-ref. If we later decide to
     * emit a reference block over these positions instead, this test should be updated, not deleted.
     */
    @Test
    public void testTargetedCallWithoutNonRefIsDropped() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_ALT));
        gb.GQ(50).PL(new int[]{500, 0, 600}).AD(new int[]{10, 20}).DP(30);
        final VariantContext targetedNoNonRef = new VariantContextBuilder("test", "chrM", 150, 150,
                Arrays.asList(CHR_M_REF, CHR_M_ALT))
                .attribute(GATKVCFConstants.TARGETED_KEY, true)
                .genotypes(gb.make()).unfiltered().make();

        reblocker.apply(targetedNoNonRef, null, null, null);
        reblocker.vcfWriter.close();

        Assert.assertEquals(mockWriter.getEmitted().size(), 0,
                "a flagged targeted record with no <NON_REF> allele should be dropped, not emitted");
        Assert.assertEquals(reblocker.droppedTargetedRecordCount, 1,
                "the dropped record should be counted so closeTool() can report it");
    }

    /**
     * Dropped targeted records are counted rather than logged individually, so the count has to accumulate across
     * the traversal. Records that are not dropped must not be counted.
     */
    @Test
    public void testDroppedTargetedRecordsAreCounted() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        attachMockWriter(reblocker);

        Assert.assertEquals(reblocker.droppedTargetedRecordCount, 0, "nothing dropped before any records are seen");

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_ALT));
        gb.GQ(50).PL(new int[]{500, 0, 600}).AD(new int[]{10, 20}).DP(30);
        for (final int position : new int[]{150, 151, 152}) {
            reblocker.apply(new VariantContextBuilder("test", "chrM", position, position,
                    Arrays.asList(CHR_M_REF, CHR_M_ALT))
                    .attribute(GATKVCFConstants.TARGETED_KEY, true)
                    .genotypes(gb.make()).unfiltered().make(), null, null, null);
        }
        Assert.assertEquals(reblocker.droppedTargetedRecordCount, 3, "each dropped record should be counted");

        // a targeted record that still has <NON_REF> is processed, not dropped, so it must not be counted
        reblocker.apply(new VariantContextBuilder(makeSomaticVariant(73))
                .attribute(GATKVCFConstants.TARGETED_KEY, true).make(), null, null, null);
        Assert.assertEquals(reblocker.droppedTargetedRecordCount, 3, "records that are not dropped should not be counted");
    }

    /**
     * The targeted bypass must be narrow: a record missing <NON_REF> that is NOT flagged as targeted is
     * still a malformed GVCF record and must still fail loudly rather than being silently discarded.
     */
    @Test
    public void testUnflaggedRecordWithoutNonRefStillThrows() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        attachMockWriter(reblocker);

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_ALT));
        gb.GQ(50).PL(new int[]{500, 0, 600}).AD(new int[]{10, 20}).DP(30);
        final VariantContext noNonRef = new VariantContextBuilder("test", "chrM", 150, 150,
                Arrays.asList(CHR_M_REF, CHR_M_ALT))
                .genotypes(gb.make()).unfiltered().make();

        Assert.assertThrows(UserException.class, () -> reblocker.apply(noNonRef, null, null, null));
    }

    /**
     * An explicit TARGETED=false is equivalent to the flag being absent, and must not open the bypass.
     */
    @Test
    public void testExplicitlyFalseTargetedFlagStillThrows() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        attachMockWriter(reblocker);

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_ALT));
        gb.GQ(50).PL(new int[]{500, 0, 600}).AD(new int[]{10, 20}).DP(30);
        final VariantContext notTargeted = new VariantContextBuilder("test", "chrM", 150, 150,
                Arrays.asList(CHR_M_REF, CHR_M_ALT))
                .attribute(GATKVCFConstants.TARGETED_KEY, false)
                .genotypes(gb.make()).unfiltered().make();

        Assert.assertThrows(UserException.class, () -> reblocker.apply(notTargeted, null, null, null));
    }

    /**
     * The bypass keys on the missing <NON_REF> allele, not on the TARGETED flag alone. A targeted record
     * that does carry <NON_REF> is a well-formed GVCF record and must still be processed and emitted.
     */
    @Test
    public void testTargetedCallWithNonRefIsNotDropped() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        final VariantContext somatic = makeSomaticVariant(73);
        final VariantContext targetedWithNonRef = new VariantContextBuilder(somatic)
                .attribute(GATKVCFConstants.TARGETED_KEY, true).make();

        reblocker.apply(targetedWithNonRef, null, null, null);
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 1,
                "a targeted record that still has <NON_REF> should be processed, not dropped");
        Assert.assertEquals(emitted.get(0).getStart(), 73);
    }

    private VariantContext makeDeletionVC(final String source, final List<Allele> alleles, final int refLength, final Genotype... genotypes) {
        final int start = DEFAULT_START;
        final int stop = start+refLength-1;
        return new VariantContextBuilder(source, "1", start, stop, alleles)
                .genotypes(Arrays.asList(genotypes)).unfiltered().log10PError(-3.0).attribute(VCFConstants.DEPTH_KEY, EXAMPLE_DP).make();
    }

    private VariantContext makeNoDepthVC(final String source, final List<Allele> alleles, final int refLength, final Genotype... genotypes) {
        final int start = DEFAULT_START;
        final int stop = start+refLength-1;
        return new VariantContextBuilder(source, "1", start, stop, alleles)
                .genotypes(Arrays.asList(genotypes)).unfiltered().log10PError(-3.0).make();
    }

    private Genotype addAD(final Genotype g, final int... ads) {
        return new GenotypeBuilder(g).AD(ads).make();
    }

    /**
     * A reference block carrying a negative GQ and PL -- malformed, but DRAGEN 4.4.6 emits them. Reblocking should
     * raise both to zero rather than failing the run. Before this, such a record reached the GQ banding code and
     * threw "GQ ... didn't fit into any partition", killing a whole-genome run over a single bad block.
     *
     * Values here mirror the real record: a zero-depth block at the chrY PAR1 boundary with GQ -1240 / PL 0,0,-1240.
     */
    @Test
    public void testNegativeGqAndPlAreFlooredToZero() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        gb.GQ(-1240).PL(new int[]{0, 0, -1240}).DP(0).AD(new int[]{0, 0});
        final VariantContext malformedBlock = new VariantContextBuilder("test", "chrM", 100, 148,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, 148)
                .genotypes(gb.make()).unfiltered().make();

        final VariantContext floored = reblocker.floorNegativeQualities(malformedBlock);
        final Genotype flooredGenotype = floored.getGenotype(0);

        Assert.assertEquals(flooredGenotype.getGQ(), 0, "a negative GQ should be raised to zero");
        Assert.assertEquals(flooredGenotype.getPL(), new int[]{0, 0, 0}, "negative PLs should be raised to zero");
        Assert.assertEquals(floored.getStart(), 100, "flooring must not disturb the block's span");
        Assert.assertEquals(floored.getAttributeAsInt(VCFConstants.END_KEY, -1), 148);
    }

    /**
     * The floor must not touch records whose qualities are already valid -- including GQ 0, which is legitimate and
     * must not be confused with a malformed value. A well-formed record should come back as the same object so the
     * check costs nothing on the overwhelming majority of records.
     */
    @Test
    public void testValidQualitiesAreLeftAlone() {
        final ReblockGVCF reblocker = new ReblockGVCF();

        final GenotypeBuilder healthy = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        healthy.GQ(45).PL(new int[]{0, 45, 450}).DP(12);
        final VariantContext healthyBlock = new VariantContextBuilder("test", "chrM", 100, 148,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, 148)
                .genotypes(healthy.make()).unfiltered().make();
        Assert.assertSame(reblocker.floorNegativeQualities(healthyBlock), healthyBlock,
                "a valid record should be returned untouched, not rebuilt");

        final GenotypeBuilder zeroGq = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        zeroGq.GQ(0).PL(new int[]{0, 0, 0}).DP(0);
        final VariantContext zeroGqBlock = new VariantContextBuilder("test", "chrM", 200, 248,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, 248)
                .genotypes(zeroGq.make()).unfiltered().make();
        Assert.assertSame(reblocker.floorNegativeQualities(zeroGqBlock), zeroGqBlock,
                "GQ 0 is valid and must not be treated as malformed");
    }

    /**
     * A malformed block must actually survive reblocking end to end, not merely pass the floor in isolation. This is
     * the behaviour that was broken: GQ banding rejected the negative value and the run died.
     */
    @Test
    public void testBlockWithNegativeQualitiesStillReblocks() {
        final ReblockGVCF reblocker = new ReblockGVCF();
        final MockVcfWriter mockWriter = attachMockWriter(reblocker);

        final GenotypeBuilder gb = new GenotypeBuilder("sample1", Arrays.asList(CHR_M_REF, CHR_M_REF));
        gb.GQ(-1240).PL(new int[]{0, 0, -1240}).DP(0).AD(new int[]{0, 0});
        final VariantContext malformedBlock = new VariantContextBuilder("test", "chrM", 100, 148,
                Arrays.asList(CHR_M_REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, 148)
                .genotypes(gb.make()).unfiltered().make();

        reblocker.apply(malformedBlock, null, null, null);
        reblocker.vcfWriter.close();

        final List<VariantContext> emitted = mockWriter.getEmitted();
        Assert.assertEquals(emitted.size(), 1, "the block should be reblocked and emitted, not rejected");
        Assert.assertEquals(emitted.get(0).getGenotype(0).getGQ(), 0,
                "the emitted block should carry the floored GQ");
    }

    private VariantContext addAttributes(final VariantContext vc, final Map<String, Object> attributes) {
        return new VariantContextBuilder(vc).attributes(attributes).make();
    }

}
