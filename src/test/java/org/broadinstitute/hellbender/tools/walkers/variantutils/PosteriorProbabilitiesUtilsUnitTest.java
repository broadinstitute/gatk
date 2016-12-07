package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@SuppressWarnings("unchecked")
public final class PosteriorProbabilitiesUtilsUnitTest extends BaseTest {

    final Allele Aref = Allele.create("A", true);
    final Allele T = Allele.create("T");
    final Allele C = Allele.create("C");
    final Allele ATC = Allele.create("ATC");
    final Allele ATCATC = Allele.create("ATCATC");

    private String arraysEq(final int[] a, final int[] b) {
        if ( a.length != b.length ) {
            return String.format("NEQ: %s | %s", Arrays.toString(a), Arrays.toString(b));
        }
        for ( int idx = 0; idx < a.length; idx++) {
            if ( a[idx] - b[idx] > 1 || b[idx] - a[idx] > 1) {
                return String.format("NEQ: %s | %s", Arrays.toString(a), Arrays.toString(b));
            }
        }

        return "";
    }

    private int[] _mleparse(final List<Integer> s) {
        final int[] mle = new int[s.size()];
        for ( int idx = 0; idx < mle.length; idx ++) {
            mle[idx] = s.get(idx);
        }

        return mle;
    }

    private Genotype makeGwithPLs(final String sample, final Allele a1, final Allele a2, final double[] pls) {
        final Genotype gt = new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
        if ( pls != null && pls.length > 0 ) {
            Assert.assertNotNull(gt.getPL());
            Assert.assertTrue(gt.getPL().length > 0);
            for ( final int i : gt.getPL() ) {
                Assert.assertTrue(i >= 0);
            }
            Assert.assertNotEquals(Arrays.toString(gt.getPL()),"[0]");
        }
        return gt;
    }

    private Genotype makeG(final String sample, final Allele a1, final Allele a2) {
        return GenotypeBuilder.create(sample, Arrays.asList(a1, a2));
    }

    private Genotype makeG(final String sample, final Allele a1, final Allele a2, final int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private VariantContext makeVC(final String source, final List<Allele> alleles, final Genotype... genotypes) {
        final int start = 10;
        final int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(Arrays.asList(genotypes)).filters((String)null).make();
    }

    @Test
    public void testCalculatePosteriorNoExternalData() {
        VariantContext test1 = makeVC("1", Arrays.asList(Aref,T), makeG("s1",Aref,T,20,0,10),
                makeG("s2",T,T,60,40,0),
                makeG("s3",Aref,Aref,0,30,90),
                makeG("s4",Aref,T,20,0,10),
                makeG("s5",T,T,60,40,0),
                makeG("s6",Aref,Aref,0,30,90),
                makeG("s7",Aref,T,20,0,10),
                makeG("s8",T,T,60,40,0),
                makeG("s9",Aref,Aref,0,30,90),
                makeG("s10",Aref,T,20,0,10),
                makeG("s11",T,T,60,40,0),
                makeG("s12",Aref,Aref,0,30,90));
        test1 = new VariantContextBuilder(test1).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,12).make();
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(test1, new ArrayList<>(), 0, 0.001, true, false, false);
        final Genotype test1exp1 = makeGwithPLs("s1",Aref,T,new double[]{-2.26110257, -0.02700903, -1.26110257});
        Assert.assertTrue(test1exp1.hasPL());
        final Genotype test1exp2 = makeGwithPLs("s2",T,T,new double[]{-6.000075e+00, -3.765981e+00, -7.488009e-05});
        final Genotype test1exp3 = makeGwithPLs("s3",Aref,Aref,new double[]{-0.0007438855, -2.7666503408, -9.0007438855});
        Assert.assertTrue(List.class.isAssignableFrom(test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY).getClass()));
        Assert.assertEquals(arraysEq(test1exp1.getPL(), _mleparse((List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test1exp2.getPL(),_mleparse((List<Integer>)test1result.getGenotype(1).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test1exp3.getPL(),_mleparse((List<Integer>)test1result.getGenotype(2).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");

        // AA AB BB AC BC CC
        // AA AC CC AT CT TT
        VariantContext test2 = makeVC("2", Arrays.asList(Aref,C,T),
                makeG("s1",Aref,T,30,10,60,0,15,90),
                makeG("s2",Aref,C,40,0,10,30,40,80),
                makeG("s3",Aref,Aref,0,5,8,15,20,40),
                makeG("s4",C,T,80,40,12,20,0,10),
                makeG("s5",Aref,T,30,10,60,0,15,90),
                makeG("s6",Aref,C,40,0,10,30,40,80),
                makeG("s7",Aref,Aref,0,5,8,15,20,40),
                makeG("s8",C,T,80,40,12,20,0,10),
                makeG("s9",Aref,T,30,10,60,0,15,90),
                makeG("s10",Aref,C,40,0,10,30,40,80),
                makeG("s11",Aref,Aref,0,5,8,15,20,40),
                makeG("s12",C,T,80,40,12,20,0,10));
        test2 = new VariantContextBuilder(test2).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, new ArrayList<>(Arrays.asList(6, 6))).make();
        final VariantContext test2result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(test2, new ArrayList<>(),5,0.001,true,false, false);
        final Genotype test2exp1 = makeGwithPLs("s1",Aref,T,new double[]{-2.823957, -1.000000, -6.686344,  0.000000, -1.952251, -9.686344});
        final Genotype test2exp2 = makeGwithPLs("s2",Aref,C,new double[]{-3.823957,  0.000000, -1.686344, -3.000000, -4.452251, -8.686344});
        final Genotype test2exp3 = makeGwithPLs("s3",Aref,Aref,new double[] {0.000000, -0.676043, -1.662387, -1.676043, -2.628294, -4.862387});
        final Genotype test2exp4 = makeGwithPLs("s4",C,T,new double[]{-7.371706, -3.547749, -1.434094, -1.547749,  0.000000, -1.234094});
        Assert.assertEquals(arraysEq(test2exp1.getPL(), _mleparse((List<Integer>)test2result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test2exp2.getPL(), _mleparse((List<Integer>)test2result.getGenotype(1).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test2exp3.getPL(), _mleparse((List<Integer>)test2result.getGenotype(2).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test2exp4.getPL(), _mleparse((List<Integer>)test2result.getGenotype(3).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
    }

    @Test
    public void testCalculatePosteriorSamplePlusExternal() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(3);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,2).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        supplTest1.add(new VariantContextBuilder(makeVC("3", Arrays.asList(Aref,T))).attribute(VCFConstants.ALLELE_COUNT_KEY,4).attribute(VCFConstants.ALLELE_NUMBER_KEY,22).make());
        supplTest1.add(makeVC("4", Arrays.asList(Aref,T),
                makeG("s_1",T,T),
                makeG("s_2",Aref,T)));
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false, false);
        // the counts here are ref=30, alt=14
        final Genotype test1exp1 = makeGwithPLs("t1",T,T,new double[]{-3.370985, -1.415172, -0.01721766});
        final Genotype test1exp2 = makeGwithPLs("t2",Aref,T,new double[]{-1.763792, -0.007978791, -3.010024});
        final Genotype test1exp3 = makeGwithPLs("t3",Aref,T,new double[]{-2.165587, -0.009773643, -1.811819});
        Assert.assertEquals(arraysEq(test1exp1.getPL(),_mleparse((List<Integer>) test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test1exp2.getPL(),_mleparse((List<Integer>) test1result.getGenotype(1).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
        Assert.assertEquals(arraysEq(test1exp3.getPL(),_mleparse((List<Integer>) test1result.getGenotype(2).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");

        final VariantContext testNonOverlapping = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,3,1,0));
        final List<VariantContext> other = Arrays.asList(makeVC("2", Arrays.asList(Aref,C),makeG("s2",C,C,10,2,0)));
        final VariantContext test2result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testNonOverlapping,other,0,0.001,true,false,false);
        final Genotype test2exp1 = makeGwithPLs("SGV",T,T,new double[]{-4.078345, -3.276502, -0.0002661066});
        Assert.assertEquals(arraysEq(test2exp1.getPL(),_mleparse((List<Integer>) test2result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY))), "");
    }

    @Test
     public void testCalculatePosteriorHOM_VARtoHET() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,1,0));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,500).attribute(VCFConstants.ALLELE_NUMBER_KEY,1000).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);

        final int[] GP = _mleparse( (List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        Assert.assertTrue(GP[2] > GP[1]);
    }

    @Test
    public void testCalculatePosteriorHETtoHOM_VAR() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,0,1));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,900).attribute(VCFConstants.ALLELE_NUMBER_KEY,1000).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);

        final int[] GP = _mleparse( (List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        Assert.assertTrue(GP[2] < GP[1]);
    }

    @Test
    public void testCalculatePosteriorHOM_REFtoHET() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,0,1,40));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,500).attribute(VCFConstants.ALLELE_NUMBER_KEY,1000).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);

        final int[] GP = _mleparse( (List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        Assert.assertTrue(GP[0] > GP[1]);
    }

    @Test
    public void testCalculatePosteriorHETtoHOM_REF() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,1,0,40));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,100).attribute(VCFConstants.ALLELE_NUMBER_KEY,1000).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);

        final int[] GP = _mleparse( (List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        Assert.assertTrue(GP[0] < GP[1]);
    }

    @Test
    public void testMLEACgreaterThanAN() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,11).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);
    }

    @Test(expectedExceptions = {UserException.class})
    public void testWrongNumberACvalues() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T,C))).attribute(VCFConstants.ALLELE_COUNT_KEY,5).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());

        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);
    }

    @Test(expectedExceptions = {UserException.class})
    public void testWrongNumberMLEACvalues() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T,C))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,5).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);
    }

    @Test
    public void testMultipleACvalues() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T,C))).attribute(VCFConstants.ALLELE_COUNT_KEY, Arrays.asList(5,4)).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);
    }

    @Test
    public void testMultipleMLEACvalues() {
        final VariantContext testOverlappingBase = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T,C))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, Arrays.asList(5,4)).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(testOverlappingBase,supplTest1,0,0.001,true,false,false);
    }

    @Test
    public void testInputIndel() {
        final VariantContext inputIndel = makeVC("1", Arrays.asList(Aref, ATC), makeG("s1",ATC,ATC,40,20,0),
                makeG("s2",Aref,ATC,18,0,24),
                makeG("s3",Aref,ATC,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,T,C))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, Arrays.asList(5,4)).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(inputIndel,supplTest1,0,0.001,true,false,false);

        final int[] GPs = _mleparse( (List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        final int[] PLs = test1result.getGenotype(0).getPL();
        Assert.assertEquals(PLs,GPs);
    }

    @Test
    public void testPriorIndel() {
        final VariantContext inputIndel = makeVC("1", Arrays.asList(Aref,T), makeG("s1",T,T,40,20,0),
                makeG("s2",Aref,T,18,0,24),
                makeG("s3",Aref,T,22,0,12));
        final List<VariantContext> supplTest1 = new ArrayList<>(1);
        supplTest1.add(new VariantContextBuilder(makeVC("2", Arrays.asList(Aref,ATC,ATCATC))).attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, Arrays.asList(5,4)).attribute(VCFConstants.ALLELE_NUMBER_KEY,10).make());
        final VariantContext test1result = PosteriorProbabilitiesUtils.calculatePosteriorProbs(inputIndel,supplTest1,0,0.001,true,false,false);

        final int[] GPs = _mleparse( (List<Integer>)test1result.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        final int[] PLs = test1result.getGenotype(0).getPL();
        Assert.assertEquals(PLs,GPs);
    }

    private double[] pl2gl(final int[] pl) {
        final double[] gl = new double[pl.length];
        for ( int idx = 0; idx < gl.length; idx++ ) {
            gl[idx] = pl[idx]/(-10.0);
        }

        return MathUtils.normalizeLog10(gl);
    }

    @Test
    public void testCalculatePosterior() {
        final int[][] likelihood_PLs  = {
                new int[]{3,0,3},
                new int[]{99,0,99},
                new int[]{50,20,0},
                new int[]{10,0,50},
                new int[]{80,60,0},
                new int[]{0,42,44}};

        final int[] altCounts = {10,40,90};
        final int[] altAlleleNum = {100,500,1000};

        final double[] expected_post_10_100 = {
                9.250326e-03, 3.020208e-01, 6.887289e-01,
                7.693433e-12, 1.000000e+00, 5.728111e-10,
                1.340156e-07, 2.192982e-03, 9.978069e-01,
                6.073718e-03, 9.938811e-01, 4.522159e-05,
                1.343101e-10, 2.197802e-07, 9.999998e-01,
                9.960193e-01, 1.028366e-03, 2.952290e-03
        };

        final double[] expected_post_10_500 = {
                4.226647e-04, 7.513277e-02, 9.244446e-01,
                1.413080e-12, 1.000000e+00, 3.090662e-09,
                4.570232e-09, 4.071661e-04, 9.995928e-01,
                1.120916e-03, 9.986339e-01, 2.451646e-04,
                4.572093e-12, 4.073320e-08, 1.000000e+00,
                9.151689e-01, 5.144399e-03, 7.968675e-02
        };

        final double[] expected_post_10_1000 = {
                1.077685e-04, 3.870477e-02, 9.611875e-01,
                6.994030e-13, 1.000000e+00, 6.237975e-09,
                1.120976e-09, 2.017756e-04, 9.997982e-01,
                5.549722e-04, 9.989500e-01, 4.949797e-04,
                1.121202e-12, 2.018163e-08, 1.000000e+00,
                7.318346e-01, 8.311615e-03, 2.598538e-01
        };

        final double[] expected_post_40_100 = {
                1.102354e-01, 6.437516e-01, 2.460131e-01,
                4.301328e-11, 1.000000e+00, 9.599306e-11,
                4.422850e-06, 1.294493e-02, 9.870507e-01,
                3.303763e-02, 9.669550e-01, 7.373032e-06,
                4.480868e-09, 1.311474e-06, 9.999987e-01,
                9.997266e-01, 1.846199e-04, 8.882157e-05
        };

        final double[] expected_post_40_500 = {
                5.711785e-03, 2.557266e-01, 7.385617e-01,
                5.610428e-12, 1.000000e+00, 7.254558e-10,
                7.720262e-08, 1.732352e-03, 9.982676e-01,
                4.436495e-03, 9.955061e-01, 5.736604e-05,
                7.733659e-11, 1.735358e-07, 9.999998e-01,
                9.934793e-01, 1.406575e-03, 5.114153e-03
        };

        final double[] expected_post_40_1000 = {
                1.522132e-03, 1.422229e-01, 8.562549e-01,
                2.688330e-12, 1.000000e+00, 1.512284e-09,
                1.776184e-08, 8.317737e-04, 9.991682e-01,
                2.130611e-03, 9.977495e-01, 1.198547e-04,
                1.777662e-11, 8.324661e-08, 9.999999e-01,
                9.752770e-01, 2.881677e-03, 2.184131e-02
        };

        final double[] expected_post_90_100 = {
                6.887289e-01, 3.020208e-01, 9.250326e-03,
                5.728111e-10, 1.000000e+00, 7.693433e-12,
                6.394346e-04, 1.405351e-01, 8.588255e-01,
                3.127146e-01, 6.872849e-01, 4.200075e-07,
                7.445327e-07, 1.636336e-05, 9.999829e-01,
                9.999856e-01, 1.386699e-05, 5.346906e-07
        };

        final double[] expected_post_90_500 = {
                2.528165e-02, 4.545461e-01, 5.201723e-01,
                1.397100e-11, 1.000000e+00, 2.874546e-10,
                4.839050e-07, 4.360463e-03, 9.956391e-01,
                1.097551e-02, 9.890019e-01, 2.258221e-05,
                4.860244e-10, 4.379560e-07, 9.999996e-01,
                9.986143e-01, 5.677671e-04, 8.179741e-04
        };

        final double[] expected_post_90_1000 = {
                7.035938e-03, 2.807708e-01, 7.121932e-01,
                6.294627e-12, 1.000000e+00, 6.371561e-10,
                9.859771e-08, 1.971954e-03, 9.980279e-01,
                4.974874e-03, 9.949748e-01, 5.035678e-05,
                9.879252e-11, 1.975850e-07, 9.999998e-01,
                9.947362e-01, 1.255272e-03, 4.008518e-03
        };

        final double[][] expectations = {
                expected_post_10_100,
                expected_post_10_500,
                expected_post_10_1000,
                expected_post_40_100,
                expected_post_40_500,
                expected_post_40_1000,
                expected_post_90_100,
                expected_post_90_500,
                expected_post_90_1000
        };

        int testIndex = 0;
        for ( final int altCount : altCounts ) {
            for ( final int numAlt : altAlleleNum ) {
                final double[] knownCounts = new double[2];
                knownCounts[0] = altCount;
                knownCounts[1] = numAlt-altCount;
                int expected_index = 0;
                for ( int gl_index = 0; gl_index < likelihood_PLs.length; gl_index++ ) {
                    final double[] post = PosteriorProbabilitiesUtils.calculatePosteriorProbs(pl2gl(likelihood_PLs[gl_index]), knownCounts, 2,false);
                    for ( int i = 0; i < post.length; i++ ) {
                        final double expected = expectations[testIndex][expected_index++];
                        final double observed = Math.pow(10.0,post[i]);
                        final double err = Math.abs( (expected-observed)/expected );
                        Assert.assertTrue(err < 1e-4, String.format("Counts: %s | Expected: %e | Observed: %e | pre %s | prior %s | post %s",
                                Arrays.toString(knownCounts), expected,observed, Arrays.toString(pl2gl(likelihood_PLs[gl_index])),
                                Arrays.toString(PosteriorProbabilitiesUtils.getDirichletPrior(knownCounts,2,false)), Arrays.toString(post)));
                    }
                }
                testIndex++;
            }
        }
    }

    private boolean arraysApproxEqual(final double[] a, final double[] b, final double tol) {
        if ( a.length != b.length ) {
            return false;
        }

        for ( int idx = 0; idx < a.length; idx++ ) {
            if ( Math.abs(a[idx]-b[idx]) > tol ) {
                return false;
            }
        }

        return true;
    }

    private String errMsgArray(final double[] a, final double[] b) {
        return String.format("Expected %s, Observed %s", Arrays.toString(a), Arrays.toString(b));
    }

    @Test
    public void testPosteriorMultiAllelic() {
        // AA AB BB AC BC CC AD BD CD DD
        final int[] PL_one = {40,20,30,0,15,25};
        final int[] PL_two = {0,20,10,99,99,99};
        final int[] PL_three = {50,40,0,30,30,10,20,40,80,50};
        final int[] PL_four  = {99,90,85,10,5,30,40,20,40,30,0,12,20,14,5};
        final int[] PL_five = {60,20,30,0,40,10,8,12,18,22,40,12,80,60,20};
        final double[] counts_one = {100.001,40.001,2.001};
        final double[] counts_two = {2504.001,16.001,218.001};
        final double[] counts_three = {10000.001,500.001,25.001,0.001};
        final double[] counts_four = {4140.001,812.001,32.001,104.001,12.001};
        final double[] counts_five = {80.001,40.001,8970.001,200.001,1922.001};

        final double[] expected_one = { -2.684035, -0.7852596, -2.4735, -0.08608339, -1.984017, -4.409852 };
        final double[] expected_two = { -5.736189e-05, -3.893688, -5.362878, -10.65938, -12.85386, -12.0186};
        final double[] expected_three = {-2.403234, -2.403276, -0.004467802, -2.70429, -4.005319, -3.59033, -6.102247, -9.403276, -14.70429, -13.40284};
        final double[] expected_four = {-7.828677, -7.335196, -7.843136, -0.7395892, -0.947033, -5.139092, -3.227715,
                -1.935159, -5.339552, -4.124552, -0.1655353, -2.072979, -4.277372, -3.165498, -3.469589 };
        final double[] expected_five = { -9.170334, -5.175724, -6.767055, -0.8250021, -5.126027, -0.07628661, -3.276762,
                -3.977787, -2.227065, -4.57769, -5.494041, -2.995066, -7.444344, -7.096104, -2.414187};

        final double[] post1 = PosteriorProbabilitiesUtils.calculatePosteriorProbs(pl2gl(PL_one),counts_one,2,false);
        final double[] post2 = PosteriorProbabilitiesUtils.calculatePosteriorProbs(pl2gl(PL_two),counts_two,2,false);
        final double[] post3 = PosteriorProbabilitiesUtils.calculatePosteriorProbs(pl2gl(PL_three),counts_three,2,false);
        final double[] post4 = PosteriorProbabilitiesUtils.calculatePosteriorProbs(pl2gl(PL_four),counts_four,2,false);
        final double[] post5 = PosteriorProbabilitiesUtils.calculatePosteriorProbs(pl2gl(PL_five),counts_five,2,false);

        final double[] expecPrior5 = {-4.2878195, -4.2932090, -4.8845400, -1.9424874, -2.2435120, -0.1937719, -3.5942477,
                -3.8952723, -1.5445506, -3.4951749, -2.6115263, -2.9125508, -0.5618292, -2.2135895,
                -1.5316722};

        Assert.assertTrue(arraysApproxEqual(expecPrior5, PosteriorProbabilitiesUtils.getDirichletPrior(counts_five,2,false),1e-5),errMsgArray(expecPrior5, PosteriorProbabilitiesUtils.getDirichletPrior(counts_five,2,false)));

        Assert.assertTrue(arraysApproxEqual(expected_one,post1,1e-6),errMsgArray(expected_one,post1));
        Assert.assertTrue(arraysApproxEqual(expected_two,post2,1e-5),errMsgArray(expected_two,post2));
        Assert.assertTrue(arraysApproxEqual(expected_three,post3,1e-5),errMsgArray(expected_three,post3));
        Assert.assertTrue(arraysApproxEqual(expected_four,post4,1e-5),errMsgArray(expected_four,post4));
        Assert.assertTrue(arraysApproxEqual(expected_five,post5,1e-5),errMsgArray(expected_five,post5));
    }
}
