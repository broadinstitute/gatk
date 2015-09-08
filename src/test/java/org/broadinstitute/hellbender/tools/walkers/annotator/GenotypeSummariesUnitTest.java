package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Map;

public final class GenotypeSummariesUnitTest {

    @Test
    public void testBasicGenotypeSummaries() {
        final String sample1 = "NA1";
        final String sample2 = "NA2";
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> ./. with missing GQ
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(noCallAllele, noCallAllele)).make());
        final VariantContext testVC = (new VariantContextBuilder())
                                        .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();

        final GenotypeSummaries GS = new GenotypeSummaries();
        final Map<String,Object> resultMap = GS.annotate(null, testVC, null);

        Assert.assertEquals(resultMap.get(GATKVCFConstants.NOCALL_CHROM_KEY), 1); // 1 no-called called sample
        Assert.assertEquals(Double.parseDouble((String) resultMap.get(GATKVCFConstants.GQ_MEAN_KEY)), 30.0, 1E-4); // mean GQ is 30
        Assert.assertFalse(resultMap.containsKey(GATKVCFConstants.GQ_STDEV_KEY)); // no stddev with only one data point

        Assert.assertEquals(GS.getKeyNames(), Arrays.asList(GATKVCFConstants.NOCALL_CHROM_KEY, GATKVCFConstants.GQ_MEAN_KEY, GATKVCFConstants.GQ_STDEV_KEY));
    }

    @Test
    public void testGenotypeSummaries() {
        final String sample1 = "NA1";
        final String sample2 = "NA2";
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());
        final VariantContext testVC = (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();

        final GenotypeSummaries GS = new GenotypeSummaries();
        final Map<String,Object> resultMap = GS.annotate(null, testVC, null);

        Assert.assertEquals(resultMap.get(GATKVCFConstants.NOCALL_CHROM_KEY), 0); // 1 no-called called sample
        Assert.assertEquals(Double.parseDouble((String) resultMap.get(GATKVCFConstants.GQ_MEAN_KEY)), 35.0, 1E-4); // mean GQ is 35
        Assert.assertEquals(Double.parseDouble((String) resultMap.get(GATKVCFConstants.GQ_STDEV_KEY)), Math.sqrt(Math.pow(30-35, 2)+Math.pow(40-35, 2)), 1E-2);
        Assert.assertEquals(GS.getKeyNames(), Arrays.asList(GATKVCFConstants.NOCALL_CHROM_KEY, GATKVCFConstants.GQ_MEAN_KEY, GATKVCFConstants.GQ_STDEV_KEY));
    }

    @Test
    public void testNoGenotypes() {
        final String sample1 = "NA1";
        final String sample2 = "NA2";
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final GenotypesContext testGC = GenotypesContext.create(2);
        final Allele noCallAllele = Allele.NO_CALL;

        // sample1 -> ./. with missing GQ
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(noCallAllele, noCallAllele)).make());
        // sample2 -> ./. with missing GQ
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(noCallAllele, noCallAllele)).make());
        final VariantContext testVC = (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).make();

        final GenotypeSummaries GS = new GenotypeSummaries();
        final Map<String,Object> resultMap = GS.annotate(null, testVC, null);
        Assert.assertNull(resultMap);
    }

}
