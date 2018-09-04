package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 7/28/16.
 */
public class AlleleFrequencyCalculatorUnitTest extends GATKBaseTest {
    private static final double EPS = 1.0e-3;
    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();

    private static final Allele A = Allele.create("A", true);
    private static final Allele B = Allele.create("C");
    private static final Allele C = Allele.create("G");

    private static final int DIPLOID = 2;
    private static final int TRIPLOID = 3;

    private static final int BIALLELIC = 2;
    private static final int TRIALLELIC = 3;

    private static final int EXTREMELY_CONFIDENT_PL = 1000;
    private static final int FAIRLY_CONFIDENT_PL = 20;

    private static final int DEFAULT_PLOIDY = 2;

    private static int sampleNameCounter = 0;

    @Test
    public void testSymmetries() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B,C);
        final Genotype AA = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,2}, FAIRLY_CONFIDENT_PL);
        final Genotype BB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {1,2}, FAIRLY_CONFIDENT_PL);
        final Genotype CC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {2,2}, FAIRLY_CONFIDENT_PL);
        final Genotype AB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,1,1}, FAIRLY_CONFIDENT_PL);
        final Genotype AC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,2,1}, FAIRLY_CONFIDENT_PL);

        final Genotype BBB = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {1,3}, FAIRLY_CONFIDENT_PL);
        final Genotype CCC = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {2,3}, FAIRLY_CONFIDENT_PL);

        // make pairs of VCs tht differ only by B <--> C
        final List<Pair<VariantContext, VariantContext>> switchBWithCPairs = Arrays.asList(
                new ImmutablePair<>(makeVC(alleles, AA, BB), makeVC(alleles, AA, CC)),
                new ImmutablePair<>(makeVC(alleles, AA, AB), makeVC(alleles, AA, AC)),
                new ImmutablePair<>(makeVC(alleles, AB, AB), makeVC(alleles, AC, AC)),
                new ImmutablePair<>(makeVC(alleles, AA, AA, BB), makeVC(alleles, AA, AA, CC)),
                new ImmutablePair<>(makeVC(alleles, AA, AB, AB), makeVC(alleles, AA, AC, AC)),
                new ImmutablePair<>(makeVC(alleles, AA, BBB), makeVC(alleles, AA, CCC))
        );
        for (final Pair<VariantContext, VariantContext> pair : switchBWithCPairs) {
            final VariantContext vc1 = pair.getLeft();
            final VariantContext vc2 = pair.getRight();
            final AFCalculationResult result1 = afCalc.getLog10PNonRef(vc1);
            final AFCalculationResult result2 = afCalc.getLog10PNonRef(vc2);
            Assert.assertEquals(result1.getLog10PosteriorOfAFEq0(), result2.getLog10PosteriorOfAFEq0(), EPS);
            Assert.assertEquals(result1.getLog10PosteriorOfAFEq0ForAllele(B), result2.getLog10PosteriorOfAFEq0ForAllele(C), EPS);
            Assert.assertEquals(result1.getLog10PosteriorOfAFEq0ForAllele(C), result2.getLog10PosteriorOfAFEq0ForAllele(B), EPS);
        }
    }

    @Test
    public void testMLECounts() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 1, 1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B,C);
        final Genotype AA = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,2}, FAIRLY_CONFIDENT_PL);
        final Genotype BB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {1,2}, FAIRLY_CONFIDENT_PL);
        final Genotype AB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,1,1}, FAIRLY_CONFIDENT_PL);
        final Genotype AC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,2,1}, FAIRLY_CONFIDENT_PL);

        final Genotype BBB = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {1,3}, FAIRLY_CONFIDENT_PL);
        final Genotype CCC = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {2,3}, FAIRLY_CONFIDENT_PL);

        final List<Pair<VariantContext, int[]>> vcWithExpectedCounts = Arrays.asList(
                new ImmutablePair<>(makeVC(alleles, AA, BB), new int[] {2,0}),
                new ImmutablePair<>(makeVC(alleles, AA, AB), new int[] {1,0}),
                new ImmutablePair<>(makeVC(alleles, AB, AB), new int[] {2,0}),
                new ImmutablePair<>(makeVC(alleles, AA, AA, BB), new int[] {2,0}),
                new ImmutablePair<>(makeVC(alleles, AA, AB, AB), new int[] {2,0}),
                new ImmutablePair<>(makeVC(alleles, AA, BBB), new int[] {3,0}),
                new ImmutablePair<>(makeVC(alleles, AA, BBB, CCC), new int[] {3,3}),
                new ImmutablePair<>(makeVC(alleles, AA, AB, AC), new int[] {1,1}),
                new ImmutablePair<>(makeVC(alleles, AA, AB, AC, BBB, CCC), new int[] {4,4})

        );
        for (final Pair<VariantContext, int[]> pair : vcWithExpectedCounts) {
            final VariantContext vc = pair.getLeft();
            final int[] expected = pair.getRight();
            final int[] actual = afCalc.getLog10PNonRef(vc).getAlleleCountsOfMLE();
            Assert.assertEquals(actual, expected);
        }
    }

    // many samples with low confidence should yield a non-zero MLE, in contrast to the old exact model
    @Test
    public void testManySamplesWithLowConfidence() {
        // prior corresponding to 1000 observations of ref, 1 of a SNP
        // for this test, we want many pseudocounts in the prior because the new AF calculator learns the allele frequency
        // and we don't want the complication of the posterior being differetn from the prior
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1000, 1, 1, DEFAULT_PLOIDY);    //prior corresponding to 1000 observations of ref, 1 of a SNP
        final List<Allele> alleles = Arrays.asList(A,B);

        // for FAIRLY_CONFIDENT_PL = 20, this genotype has about 100 times greater likelihood to be het than hom ref
        // with our prior giving 1000 times as much weight to ref, this implies a 1 in 5 chance of each sample having a copy of the alt allele
        // (that is, 100/1000 times the combinatorial factor of 2).  Thus the MLE for up to 2 samples should be zero
        // for five samples we should have one
        // for ten samples we will have more than twice as many as for five since the counts fromt he samples start to influence
        // the estimated allele frequency
        final Genotype AB = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,1,1,1}, FAIRLY_CONFIDENT_PL);

        final List<VariantContext> vcsWithDifferentNumbersOfSamples = IntStream.range(1, 11)
                .mapToObj(n -> makeVC(alleles, Collections.nCopies(n, AB))).collect(Collectors.toList());
        final int[] counts = vcsWithDifferentNumbersOfSamples.stream().mapToInt(vc -> afCalc.getLog10PNonRef(vc).getAlleleCountAtMLE(B)).toArray();
        Assert.assertEquals(counts[0],0); // one sample
        Assert.assertEquals(counts[1],0); // two samples
        Assert.assertEquals(counts[4],2); // five samples
        Assert.assertTrue(counts[8] >= 3); // ten samples
    }

    @Test
    public void testApproximateMultiplicativeConfidence() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 1, 1, DEFAULT_PLOIDY);    //flat prior -- we will choose genotypes such that the posterior remains flat
        final List<Allele> alleles = Arrays.asList(A,B);

        final Genotype AA = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,2}, FAIRLY_CONFIDENT_PL);
        final Genotype BB = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {1,2}, FAIRLY_CONFIDENT_PL);

        final List<VariantContext> vcsWithDifferentNumbersOfSamples = new ArrayList<>();
        final List<Genotype> genotypeList = new ArrayList<>();

        for (int n = 0; n < 10; n++) {
            genotypeList.add(AA);
            genotypeList.add(BB);   //adding both keeps the flat prior.  Thus the posterior will equal the likelihood
            vcsWithDifferentNumbersOfSamples.add(makeVC(alleles, genotypeList));
        }

        // since we maintain a flat allele frequency distribution, the probability of being ref as each successive sample is added
        // is multiplied by the probability of any one.  Thus we get an arithmetic series in log space
        final double[] log10PRefs = vcsWithDifferentNumbersOfSamples.stream()
                .mapToDouble(vc -> afCalc.getLog10PNonRef(vc).getLog10LikelihoodOfAFEq0()).toArray();

        for (int n = 0; n < 9; n++) {
            Assert.assertEquals(log10PRefs[n+1] - log10PRefs[n], log10PRefs[0], 0.01);
        }
    }

    @Test
    public void testManyRefSamplesDontKillGoodVariant() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B);
        final Genotype AA = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,2}, EXTREMELY_CONFIDENT_PL);
        final Genotype AB = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,1,1,1}, EXTREMELY_CONFIDENT_PL);
        for (final int numRef : new int[]{1, 10, 100, 1000, 10000, 100000}) {
            final List<Genotype> genotypeList = new ArrayList<>(Collections.nCopies(numRef, AA));
            genotypeList.add(AB);
            final VariantContext vc = makeVC(alleles, genotypeList);
            final double log10PRef = afCalc.getLog10PNonRef(vc).getLog10LikelihoodOfAFEq0();
            Assert.assertTrue(log10PRef < (-EXTREMELY_CONFIDENT_PL/10) + Math.log10(numRef) + 1);
        }
    }

    // spanning deletion alleles should be treated as non-variant, so we make a few random PLs and test that the
    // qual is invariant to swapping the ref/ref PL with the ref/SPAN_DEL PL
    @Test
    public void testSpanningDeletionIsNotConsideredVariant() {
        final int ploidy = 2;
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, ploidy);
        final List<Allele> alleles = Arrays.asList(A, B, Allele.SPAN_DEL);

        // some pls that have high likelihood for span del allele but not for the SNP (B)
        final int[] spanDelPls = new int[] {50, 100, 100, 0, 100, 100};

        // some pls that weakly support the SNP
        final int[] lowQualSnpPls = new int[] {10,0,40,100,70,300};


        final Genotype spanDel = makeGenotype(ploidy, spanDelPls);
        final Genotype lowQualSNP = makeGenotype(ploidy, lowQualSnpPls);

        // first test the span del genotype alone.  Its best PL containing the SNP is 100, so we expect a variant probability
        // of about 10^(-100/10) -- a bit less due to the prior bias in favor of the reference
        final VariantContext vcSpanDel = makeVC(alleles, Arrays.asList(spanDel));
        final double log10PVariant = afCalc.getLog10PNonRef(vcSpanDel).getLog10LikelihoodOfAFGT0();
        Assert.assertTrue(log10PVariant < - 10);

        // now test a realistic situation of two samples, one with a low-quality SNP and one with the spanning deletion
        // we want to find that the spanning deletion has little effect on the qual
        // In fact, we also want to check that it *decreases* the qual, because it's essentially one more hom ref sample
        // Furthermore, to be precise it should be really behave almost identically to a hom ref *haploid* sample,
        // so we check that, too
        final VariantContext vcLowQualSnp = makeVC(alleles, Arrays.asList(lowQualSNP));
        final double lowQualSNPQualScore = afCalc.getLog10PNonRef(vcLowQualSnp).getLog10LikelihoodOfAFGT0();
        final VariantContext vcBoth = makeVC(alleles, Arrays.asList(lowQualSNP, spanDel));
        final double bothQualScore = afCalc.getLog10PNonRef(vcBoth).getLog10LikelihoodOfAFGT0();
        Assert.assertEquals(lowQualSNPQualScore, bothQualScore, 0.1);
        Assert.assertTrue(bothQualScore < lowQualSNPQualScore);

        final int[] haploidRefPls = new int[] {0, 100, 100};
        final Genotype haploidRef = makeGenotype(1, haploidRefPls);

        final VariantContext vcLowQualSnpAndHaploidRef = makeVC(alleles, Arrays.asList(lowQualSNP, haploidRef));
        final double lowQualSNPAndHaplpidRefQualScore = afCalc.getLog10PNonRef(vcLowQualSnpAndHaploidRef).getLog10LikelihoodOfAFGT0();
        Assert.assertEquals(bothQualScore, lowQualSNPAndHaplpidRefQualScore, 1e-5);

        // as a final test, we check that getting rid of the spanning deletion allele, in the sense that
        // REF / SPAN_DEL --> haploid REF; REF / SNP --> REF / SNP
        // does not affect the qual score

        final int[] haploidRefPlsWithoutSpanDel = new int[] {0, 100};
        final int[] snpPlsWithoutSpanDel = new int[] {10, 0, 40};
        final VariantContext vcNoSpanDel = makeVC(Arrays.asList(A,B), Arrays.asList(makeGenotype(ploidy, snpPlsWithoutSpanDel),
                makeGenotype(1, haploidRefPlsWithoutSpanDel)));
        final double noSpanDelQualScore = afCalc.getLog10PNonRef(vcNoSpanDel).getLog10LikelihoodOfAFGT0();
        Assert.assertEquals(bothQualScore, noSpanDelQualScore, 1e-6);
    }

    @Test
    public void testPresenceOfUnlikelySpanningDeletionDoesntAffectResults() {
        final int ploidy = 2;
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, ploidy);
        final List<Allele> allelesWithoutSpanDel = Arrays.asList(A, B);
        final List<Allele> allelesWithSpanDel = Arrays.asList(A, B, Allele.SPAN_DEL);
        
        // make PLs that support an A/B genotype
        final int[] plsWithoutSpanDel = new int[] {50, 0, 50};
        final int[] plsWithSpanDel = new int[] {50, 0, 50, 100, 100, 100};
        final Genotype genotypeWithoutSpanDel = makeGenotype(ploidy, plsWithoutSpanDel);
        final Genotype genotypeWithSpanDel = makeGenotype(ploidy, plsWithSpanDel);
        final VariantContext vcWithoutSpanDel = makeVC(allelesWithoutSpanDel, Arrays.asList(genotypeWithoutSpanDel));
        final VariantContext vcWithSpanDel = makeVC(allelesWithSpanDel, Arrays.asList(genotypeWithSpanDel));
        final double log10PVariantWithoutSpanDel = afCalc.getLog10PNonRef(vcWithoutSpanDel).getLog10LikelihoodOfAFGT0();
        final double log10PVariantWithSpanDel = afCalc.getLog10PNonRef(vcWithSpanDel).getLog10LikelihoodOfAFGT0();
        Assert.assertEquals(log10PVariantWithoutSpanDel, log10PVariantWithSpanDel, 0.0001);
    }

    // test that a finite precision bug for span del sites with a very unlikely alt allele doesn't occur
    @Test
    public void testSpanningDeletionWithVeryUnlikelyAltAllele() {
        final int ploidy = 4;
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, ploidy);
        final List<Allele> alleles = Arrays.asList(A, Allele.SPAN_DEL, B);

        // make PLs that don't support the alt allele
        final List<int[]> pls = Arrays.asList(new int[] {0,10000,10000,10000,10000, 10000,10000,10000,10000,10000,10000,10000,10000,10000,10000});
        final VariantContext vc = makeVC(alleles, pls.stream().map(pl -> makeGenotype(ploidy, pl)).collect(Collectors.toList()));
        final double log10PVariant = afCalc.getLog10PNonRef(vc).getLog10LikelihoodOfAFGT0();
    }

    // make PLs that correspond to an obvious call i.e. one PL is relatively big and the rest are zero
    // alleleCounts is the GenotypeAlleleCounts format for the obvious genotype, with repeats but in no particular order
    private static int[] PLsForObviousCall(final int ploidy, final int numAlleles, final int[] alleleCounts, final int PL)   {
        final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, numAlleles);
        final int[] result = Collections.nCopies(glCalc.genotypeCount(), PL).stream().mapToInt(n->n).toArray();
        result[glCalc.alleleCountsToIndex(alleleCounts)] = 0;
        return result;
    }

    private static Genotype genotypeWithObviousCall(final int ploidy, final int numAlleles, final int[] alleles, final int PL) {
        return makeGenotype(ploidy, PLsForObviousCall(ploidy, numAlleles, alleles, PL));
    }
    //note the call is irrelevant to the AFCalculator, which only looks at PLs
    private static Genotype makeGenotype(final int ploidy, int ... pls) {
        return new GenotypeBuilder("sample" + sampleNameCounter++).alleles(Collections.nCopies(ploidy, Allele.NO_CALL)).PL(pls).make();
    }

    private static VariantContext makeVC(final List<Allele> alleles, final Genotype... genotypes) {
        return new VariantContextBuilder().chr("chr1").alleles(alleles).genotypes(genotypes).make();
    }

    private static VariantContext makeVC(final List<Allele> alleles, final Collection<Genotype> genotypes) {
        return new VariantContextBuilder().chr("chr1").alleles(alleles).genotypes(genotypes).make();
    }
}