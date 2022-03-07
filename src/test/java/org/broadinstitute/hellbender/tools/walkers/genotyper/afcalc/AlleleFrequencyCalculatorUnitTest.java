package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeIndexCalculator;
import org.broadinstitute.hellbender.utils.MathUtils;
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
            final AFCalculationResult result1 = afCalc.calculate(vc1);
            final AFCalculationResult result2 = afCalc.calculate(vc2);
            Assert.assertEquals(result1.log10ProbOnlyRefAlleleExists(), result2.log10ProbOnlyRefAlleleExists(), EPS);
            Assert.assertEquals(result1.getLog10PosteriorOfAlleleAbsent(B), result2.getLog10PosteriorOfAlleleAbsent(C), EPS);
            Assert.assertEquals(result1.getLog10PosteriorOfAlleleAbsent(C), result2.getLog10PosteriorOfAlleleAbsent(B), EPS);
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
            final int[] actual = afCalc.calculate(vc).getAlleleCountsOfMLE();
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
        final int[] counts = vcsWithDifferentNumbersOfSamples.stream().mapToInt(vc -> afCalc.calculate(vc).getAlleleCountAtMLE(B)).toArray();
        Assert.assertEquals(counts[0],0); // one sample
        Assert.assertEquals(counts[1],0); // two samples
        Assert.assertEquals(counts[4],2); // five samples
        Assert.assertTrue(counts[8] >= 3); // ten samples
    }

    // a previous implementation of {@link AlleleFrequencyCalculator} had finite precision errors that manifested
    // when there were very many confident samples.  This is a regression test for that bug.
    // See https://github.com/broadinstitute/gatk/issues/4833
    @Test
    public void testManyVeryConfidentSamples() {
        // flat prior to simplify back-of-the-envelope calculations
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 1, 1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B,C);


        final Genotype AC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,2,1}, EXTREMELY_CONFIDENT_PL);
        for (final int numSamples : new int[] {100, 1000}) {

            final VariantContext vc = makeVC(alleles, Collections.nCopies(numSamples, AC));
            final AFCalculationResult result = afCalc.calculate(vc);
            Assert.assertEquals(result.getAlleleCountAtMLE(B), 0);
            Assert.assertEquals(result.getAlleleCountAtMLE(C), numSamples);

            Assert.assertEquals(result.log10ProbOnlyRefAlleleExists(), result.getLog10PosteriorOfAlleleAbsent(C), numSamples * 0.01);

            // with a large number of samples all with the AC genotype, the calculator will learn that the frequencies of the A and C alleles
            // are 1/2, while the frequency of the B allele is 0.  Thus the only genotypes with appreciable priors are AA, AC, and CC
            // with priors of 1/4, 1/2, and 1/4 and relative likelihoods of 1, 10^(PL/10), and 1
            // the posterior probability of each sample having the C allele is thus
            // (1 + 2*10^(PL/10))/(1 + 2*10^(PL/10) + 1) = (1 + x/2)/(1 + x), where x = 10^(-PL/10)

            // to first-order in x, which is an extremely good approximation, this is 1 - x/2
            // thus the probability that N identical samples don't have the C allele is (x/2)^N, and the log-10 probability of this is
            // N * [log_10(1/2) - PL/10]
            final double expectedLog10ProbabilityOfNoCAllele = numSamples * (MathUtils.LOG10_ONE_HALF - EXTREMELY_CONFIDENT_PL / 10);
            Assert.assertEquals(result.getLog10PosteriorOfAlleleAbsent(C), expectedLog10ProbabilityOfNoCAllele, numSamples * 0.01);
        }
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
                .mapToDouble(vc -> afCalc.calculate(vc).log10ProbOnlyRefAlleleExists()).toArray();

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
            final double log10PRef = afCalc.calculate(vc).log10ProbOnlyRefAlleleExists();
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
        final double log10PVariant = afCalc.calculate(vcSpanDel).log10ProbVariantPresent();
        Assert.assertTrue(log10PVariant < - 10);

        // now test a realistic situation of two samples, one with a low-quality SNP and one with the spanning deletion
        // we want to find that the spanning deletion has little effect on the qual
        // In fact, we also want to check that it *decreases* the qual, because it's essentially one more hom ref sample
        // Furthermore, to be precise it should be really behave almost identically to a hom ref *haploid* sample,
        // so we check that, too
        final VariantContext vcLowQualSnp = makeVC(alleles, Arrays.asList(lowQualSNP));
        final double lowQualSNPQualScore = afCalc.calculate(vcLowQualSnp).log10ProbVariantPresent();
        final VariantContext vcBoth = makeVC(alleles, Arrays.asList(lowQualSNP, spanDel));
        final double bothQualScore = afCalc.calculate(vcBoth).log10ProbVariantPresent();
        Assert.assertEquals(lowQualSNPQualScore, bothQualScore, 0.1);
        Assert.assertTrue(bothQualScore < lowQualSNPQualScore);

        final int[] haploidRefPls = new int[] {0, 100, 100};
        final Genotype haploidRef = makeGenotype(1, haploidRefPls);

        final VariantContext vcLowQualSnpAndHaploidRef = makeVC(alleles, Arrays.asList(lowQualSNP, haploidRef));
        final double lowQualSNPAndHaplpidRefQualScore = afCalc.calculate(vcLowQualSnpAndHaploidRef).log10ProbVariantPresent();
        Assert.assertEquals(bothQualScore, lowQualSNPAndHaplpidRefQualScore, 1e-5);

        // as a final test, we check that getting rid of the spanning deletion allele, in the sense that
        // REF / SPAN_DEL --> haploid REF; REF / SNP --> REF / SNP
        // does not affect the qual score

        final int[] haploidRefPlsWithoutSpanDel = new int[] {0, 100};
        final int[] snpPlsWithoutSpanDel = new int[] {10, 0, 40};
        final VariantContext vcNoSpanDel = makeVC(Arrays.asList(A,B), Arrays.asList(makeGenotype(ploidy, snpPlsWithoutSpanDel),
                makeGenotype(1, haploidRefPlsWithoutSpanDel)));
        final double noSpanDelQualScore = afCalc.calculate(vcNoSpanDel).log10ProbVariantPresent();
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
        final double log10PVariantWithoutSpanDel = afCalc.calculate(vcWithoutSpanDel).log10ProbVariantPresent();
        final double log10PVariantWithSpanDel = afCalc.calculate(vcWithSpanDel).log10ProbVariantPresent();
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
        final double log10PVariant = afCalc.calculate(vc).log10ProbVariantPresent();
    }

    //check the exact integration calculation
    @Test
    public void testSingleSampleBiallelicShortcut() {
        // in the haploid case, if the AF calc has equal pseudocounts of ref and alt, the posterior is proportional to the likelihoods:
        for (final double pseudocount : new double[] { 1, 5, 10} ) {
            final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(pseudocount, pseudocount, pseudocount, DEFAULT_PLOIDY);
            for (int pl : new int[]{10, 100, 1000}) {
                final double[] log10Likelihoods = new double[]{0, pl / 10.0};
                final double result = afCalc.calculateSingleSampleBiallelicNonRefPosterior(log10Likelihoods, false);
                final double expected = MathUtils.normalizeFromLog10ToLinearSpace(log10Likelihoods)[1];
                Assert.assertEquals(result, expected, 1.0e-10);
            }
        }

        // in the diploid case, we roughly multiply the prior by the likelihoods -- it's not exact because the allele frequency is a random variable and not a single value
        for (final double heterozygosity : new double[] {0.1, 0.01, 0.001}) {
            final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(100, 100*heterozygosity, 100*heterozygosity, DEFAULT_PLOIDY);
            for (int pl : new int[]{10, 100, 1000}) {
                final double log10Likelihood = pl / 10.0;
                final double[] log10Likelihoods = new double[]{0, log10Likelihood, -100}; //assume hom alt is very unlikely
                final double[] log10Priors = new double[] {Math.log10(MathUtils.square(1 - heterozygosity)), Math.log10(2*heterozygosity*(1-heterozygosity)), MathUtils.square(heterozygosity)};
                final double result = afCalc.calculateSingleSampleBiallelicNonRefPosterior(log10Likelihoods, false);
                final double expected = 1 - MathUtils.normalizeFromLog10ToLinearSpace(MathArrays.ebeAdd(log10Likelihoods, log10Priors))[0];
                Assert.assertEquals(result, expected, 0.03);
            }
        }


    }

    // make PLs that correspond to an obvious call i.e. one PL is relatively big and the rest are zero
    // alleleCounts is the GenotypeAlleleCounts format for the obvious genotype, with repeats but in no particular order
    private static int[] PLsForObviousCall(final int ploidy, final int numAlleles, final int[] alleleCounts, final int PL)   {
        final int[] result = Collections.nCopies(GenotypeIndexCalculator.genotypeCount(ploidy, numAlleles), PL).stream().mapToInt(n->n).toArray();
        result[GenotypeIndexCalculator.alleleCountsToIndex(alleleCounts)] = 0;
        return result;
    }

    private static Genotype genotypeWithObviousCall(final int ploidy, final int numAlleles, final int[] alleles, final int PL) {
        return makeGenotype(ploidy, PLsForObviousCall(ploidy, numAlleles, alleles, PL));
    }
    //note the call is irrelevant to the AlleleFrequencyCalculator, which only looks at PLs
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