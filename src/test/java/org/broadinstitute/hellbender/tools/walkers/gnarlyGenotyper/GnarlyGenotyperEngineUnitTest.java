package org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Created by gauthier on 6/1/18.
 */
public class GnarlyGenotyperEngineUnitTest {
    private final static Allele Aref = Allele.create("A", true);
    private final static Allele oneInserted = Allele.create("AA");
    private final static Allele twoInserted = Allele.create("AAA");
    private final static Allele threeInserted = Allele.create("AAAA");
    private final static Allele fourRepeats = Allele.create("AAAAA");
    private final static Allele fiveRepeats = Allele.create("AAAAAA");

    //[Sample1 C/CTCTTCTTCTTCTTCTTCTTCTTCTTCT GQ 71 DP 13 AD 0,9,2,0,0,0,0 PL 491,129,305,364,0,351,435,71,357,422,435,71,357,422,422,435,71,357,422,422,422,435,71,357,422,422,422,422 {SB=[0, 0, 2, 11]}]
    //[Sample2 CTCTTCTTCTTCTTCTTCT*/C GQ 99 DP 15 AD 7,8,0,0,0,0,0 PL 289,0,205,310,231,541,310,231,541,541,310,231,541,541,541,310,231,541,541,541,541,310,231,541,541,541,541,541 {SB=[5, 2, 1, 7]}]

    private final static int[] sample1pls = {491,129,305,364,0,351,435,71,357,422,435,71,357,422,422,435,71,357,422,422,422,435,71,357,422,422,422,422};
    private final static int[] sample2pls = {289,0,205,310,231,541,310,231,541,541,310,231,541,541,541,310,231,541,541,541,541,310,231,541,541,541,541,541};

    //use more alts than the maxAltAllelesToOutput for the engine, forcing on-the-fly generation of PL counts not in the cache
    @Test
    public void testLotsOfAlts() {
        final GnarlyGenotyperEngine engine = new GnarlyGenotyperEngine(false, 4, false, true);

        final Genotype g1 = VariantContextTestUtils.makeG("g1", oneInserted, twoInserted, sample1pls);
        final Genotype g2 = VariantContextTestUtils.makeG("g1", Aref, oneInserted, sample2pls);
        final VariantContext vc = VariantContextTestUtils.makeVC("test", Arrays.asList(Aref, oneInserted, twoInserted, threeInserted, fourRepeats, fiveRepeats), g1, g2);

        final Map<Allele, Integer> alleleCounts = new HashMap<>();
        alleleCounts.put(Aref, 1);
        alleleCounts.put(oneInserted, 1);
        alleleCounts.put(twoInserted,2);
        final int[] sbSum = {10,20,30,40};
        final int[] rawGenotypeCounts = {0, 1, 1};

        final GenotypesContext genotypes = engine.iterateOnGenotypes(vc, Arrays.asList(Aref, oneInserted, twoInserted, threeInserted, fourRepeats),
                alleleCounts, sbSum, false, true, rawGenotypeCounts);

        Assert.assertTrue(genotypes.get(0).hasPL() && genotypes.get(0).getPL().length == 15);
        Assert.assertTrue(genotypes.get(1).hasPL() && genotypes.get(1).getPL().length == 15);

        // repeat, but this time request that PLs are not returned in the GenotypesContext and verify
        final GenotypesContext genotypesNoPl = engine.iterateOnGenotypes(vc, Arrays.asList(Aref, oneInserted, twoInserted, threeInserted, fourRepeats),
                alleleCounts, sbSum, false, false, rawGenotypeCounts);

        Assert.assertFalse(genotypesNoPl.get(0).hasPL());
        Assert.assertFalse(genotypesNoPl.get(1).hasPL());
    }

    //use more alts than the maxAltAllelesToOutput for the engine, forcing on-the-fly generation of GLCalculator not in the cache
    @Test
    public void testGenotypeCallForLotsOfAlts() {
        final GnarlyGenotyperEngine engine = new GnarlyGenotyperEngine(false, 4, true);

        final Genotype g1 = VariantContextTestUtils.makeG("g1", oneInserted, twoInserted, sample1pls);
        final Genotype g2 = VariantContextTestUtils.makeG("g1", Aref, oneInserted, sample2pls);
        final VariantContext vc = VariantContextTestUtils.makeVC("test", Arrays.asList(Aref, oneInserted, twoInserted, threeInserted, fourRepeats, fiveRepeats), g1, g2);

        final GenotypeBuilder builder1 = new GenotypeBuilder(g1);
        engine.makeGenotypeCall(g1, builder1, GenotypeLikelihoods.fromPLs(sample1pls).getAsVector(), Arrays.asList(Aref, oneInserted, twoInserted, threeInserted, fourRepeats, fiveRepeats));
        final List<Allele> calledAlleles1 = builder1.make().getAlleles();
        Assert.assertTrue(calledAlleles1.size() == 2 && calledAlleles1.contains(oneInserted) && calledAlleles1.contains(twoInserted));

        final GenotypeBuilder builder2 = new GenotypeBuilder(g2);
        engine.makeGenotypeCall(g2, builder2, GenotypeLikelihoods.fromPLs(sample2pls).getAsVector(), Arrays.asList(Aref, oneInserted, twoInserted, threeInserted, fourRepeats, fiveRepeats));
        final List<Allele> calledAlleles2 = builder2.make().getAlleles();
        Assert.assertTrue(calledAlleles2.size() == 2 && calledAlleles2.contains(Aref) && calledAlleles2.contains(oneInserted));
    }
}