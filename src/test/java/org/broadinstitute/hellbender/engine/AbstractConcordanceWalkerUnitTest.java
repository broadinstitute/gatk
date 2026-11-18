package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AbstractConcordanceWalkerUnitTest extends GATKBaseTest {

    /**
     * Matches records with identical alleles
     */
    private static class TestConcordanceWalker extends AbstractConcordanceWalker {
        private final boolean allowTruthToMatchMultipleEvals;

        TestConcordanceWalker(final boolean allowTruthToMatchMultipleEvals) {
            this.allowTruthToMatchMultipleEvals = allowTruthToMatchMultipleEvals;
        }

        @Override
        protected boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval) {
            return truth.getAlleles().equals(eval.getAlleles());
        }

        @Override
        protected boolean allowTruthToMatchMultipleEvals() {
            return allowTruthToMatchMultipleEvals;
        }

        @Override
        protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) { }
    }

    private static VariantContext makeVariant(final String id, final Allele... alleles) {
        return new VariantContextBuilder("test", "chr1", 100, 100, Arrays.asList(alleles)).id(id).make();
    }

    private static List<AbstractConcordanceWalker.TruthVersusEval> match(final boolean allowTruthToMatchMultipleEvals,
                                                                         final List<VariantContext> truth, final List<VariantContext> eval) {
        return new ArrayList<>(new TestConcordanceWalker(allowTruthToMatchMultipleEvals).matchVariants(truth, eval));
    }

    private static void assertTruthVersusEval(final AbstractConcordanceWalker.TruthVersusEval actual, final VariantContext expectedTruth, final VariantContext expectedEval) {
        Assert.assertSame(actual.getTruth(), expectedTruth);
        Assert.assertSame(actual.getEval(), expectedEval);
    }

    @Test
    public void testMatchesCompatibleRecords() {
        final VariantContext truthC = makeVariant("truthC", Allele.REF_A, Allele.ALT_C);
        final VariantContext truthG = makeVariant("truthG", Allele.REF_A, Allele.ALT_G);
        final VariantContext evalG = makeVariant("evalG", Allele.REF_A, Allele.ALT_G);
        final VariantContext evalC = makeVariant("evalC", Allele.REF_A, Allele.ALT_C);

        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(false, Arrays.asList(truthC, truthG), Arrays.asList(evalG, evalC));
        Assert.assertEquals(result.size(), 2);
        assertTruthVersusEval(result.get(0), truthG, evalG);
        assertTruthVersusEval(result.get(1), truthC, evalC);
    }

    @Test
    public void testUnmatchedRecords() {
        final VariantContext truthC = makeVariant("truthC", Allele.REF_A, Allele.ALT_C);
        final VariantContext truthG = makeVariant("truthG", Allele.REF_A, Allele.ALT_G);
        final VariantContext truthT = makeVariant("truthT", Allele.REF_A, Allele.ALT_T);
        final VariantContext evalG = makeVariant("evalG", Allele.REF_A, Allele.ALT_G);
        final VariantContext evalA = makeVariant("evalA", Allele.REF_C, Allele.ALT_A);

        //eval records come first in input order, followed by unmatched truth records in input order
        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(false, Arrays.asList(truthT, truthG, truthC), Arrays.asList(evalA, evalG));
        Assert.assertEquals(result.size(), 4);
        assertTruthVersusEval(result.get(0), null, evalA);
        assertTruthVersusEval(result.get(1), truthG, evalG);
        assertTruthVersusEval(result.get(2), truthT, null);
        assertTruthVersusEval(result.get(3), truthC, null);
    }

    @Test
    public void testUnmatchedTruthOrderIsDeterministic() {
        final List<VariantContext> truth = new ArrayList<>();
        for (int i = 0; i < 50; i++) {
            truth.add(makeVariant("truth" + i, Allele.REF_A, Allele.ALT_C));
        }

        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(false, truth, Collections.emptyList());
        Assert.assertEquals(result.size(), truth.size());
        for (int i = 0; i < truth.size(); i++) {
            assertTruthVersusEval(result.get(i), truth.get(i), null);
        }
    }

    @Test
    public void testFirstCompatibleTruthIsMatchedRegardlessOfId() {
        final VariantContext truth1 = makeVariant("rs1", Allele.REF_A, Allele.ALT_C);
        final VariantContext truth2 = makeVariant("rs2", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval = makeVariant("rs2", Allele.REF_A, Allele.ALT_C);

        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(false, Arrays.asList(truth1, truth2), Collections.singletonList(eval));
        Assert.assertEquals(result.size(), 2);
        assertTruthVersusEval(result.get(0), truth1, eval);
        assertTruthVersusEval(result.get(1), truth2, null);
    }

    @Test
    public void testTruthMatchedOnlyOnceByDefault() {
        final VariantContext truth = makeVariant("truth", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval1 = makeVariant("eval1", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval2 = makeVariant("eval2", Allele.REF_A, Allele.ALT_C);

        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(false, Collections.singletonList(truth), Arrays.asList(eval1, eval2));
        Assert.assertEquals(result.size(), 2);
        assertTruthVersusEval(result.get(0), truth, eval1);
        assertTruthVersusEval(result.get(1), null, eval2);
    }

    @Test
    public void testOneToOneMatchingUsesNextCompatibleTruth() {
        final VariantContext truth1 = makeVariant("truth1", Allele.REF_A, Allele.ALT_C);
        final VariantContext truth2 = makeVariant("truth2", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval1 = makeVariant("eval1", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval2 = makeVariant("eval2", Allele.REF_A, Allele.ALT_C);

        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(false, Arrays.asList(truth1, truth2), Arrays.asList(eval1, eval2));
        Assert.assertEquals(result.size(), 2);
        assertTruthVersusEval(result.get(0), truth1, eval1);
        assertTruthVersusEval(result.get(1), truth2, eval2);
    }

    @Test
    public void testTruthMatchedMultipleTimesWhenAllowed() {
        final VariantContext truth = makeVariant("truth", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval1 = makeVariant("eval1", Allele.REF_A, Allele.ALT_C);
        final VariantContext eval2 = makeVariant("eval2", Allele.REF_A, Allele.ALT_C);

        final List<AbstractConcordanceWalker.TruthVersusEval> result = match(true, Collections.singletonList(truth), Arrays.asList(eval1, eval2));
        Assert.assertEquals(result.size(), 2);
        assertTruthVersusEval(result.get(0), truth, eval1);
        assertTruthVersusEval(result.get(1), truth, eval2);
    }
}
