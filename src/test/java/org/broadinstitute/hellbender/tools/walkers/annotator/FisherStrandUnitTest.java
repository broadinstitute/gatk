package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;


public final class FisherStrandUnitTest {
    private static final double DELTA_PRECISION = 10e-8;
    private static final double DELTA_PRECISION_FOR_PHRED = 0.001;

    @DataProvider(name = "UsingTable")
    public Object[][] makeUsingTable() {

        /* NOTE: the expected P values were computed in R after `normalizing` the table, as follows
             fisher = function(v) {
                return(fisher.test(matrix(v, nrow=2, ncol=2))$p.value)
             }
             normalize = function(v) {
                s=sum(as.numeric(v))
                if (s < 2 * 200){ return(v) }
                return(v/(s/200))
            }
        */
        //> fisher(floor(normalize(c(2068, 6796, 1133, 0))))
        //[1] 5.919091e-13

        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{0, 0, 0, 0,                     1.0});
        tests.add(new Object[]{100000, 100000, 100000, 100000, 1.0}  );
        tests.add(new Object[]{0, 0, 100000, 100000,           1.0});
        tests.add(new Object[]{100000,100000,100000,0,         1.312515e-15});

        tests.add(new Object[]{100,100,100,0,                  3.730187e-23});
        tests.add(new Object[]{13736,9047,41,1433,             6.162592e-05});
        tests.add(new Object[]{66, 14, 64, 4,                  4.243330e-02});
        tests.add(new Object[]{351169, 306836, 153739, 2379,   2.193607e-09});
        tests.add(new Object[]{116449, 131216, 289, 16957,     1.340052e-03});
        tests.add(new Object[]{137, 159, 9, 23,                6.088506e-02});
        tests.add(new Object[]{129, 90, 21, 20,                3.919603e-01});
        tests.add(new Object[]{14054, 9160, 16, 7827,          7.466277e-17});
        tests.add(new Object[]{32803, 9184, 32117, 3283,       1.795855e-02});
        tests.add(new Object[]{2068, 6796, 1133, 0,            5.919091e-13});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "UsingTable")
    public void testUsingTableData(final int refpos, final int refneg, final int altpos, final int altneg, double expectedPvalue) {
        final int[][] contingencyTable = new int[2][2];
        contingencyTable[0][0] = refpos;
        contingencyTable[0][1] = refneg;
        contingencyTable[1][0] = altpos;
        contingencyTable[1][1] = altneg;
        final double pvalue = FisherStrand.pValueForContingencyTable(contingencyTable);
        Assert.assertEquals(pvalue, expectedPvalue, DELTA_PRECISION, "Pvalues");

        final double actualPhredPval = Double.parseDouble((String) new FisherStrand().annotationForOneTable(pvalue).get(GATKVCFConstants.FISHER_STRAND_KEY));
        final double expectedPhredPvalue =  QualityUtils.phredScaleErrorRate(Math.max(expectedPvalue, FisherStrand.MIN_PVALUE));
        Assert.assertEquals(actualPhredPval, expectedPhredPvalue, DELTA_PRECISION_FOR_PHRED, "phred pvalues");

        Assert.assertEquals(new FisherStrand().getKeyNames(), Collections.singletonList(GATKVCFConstants.FISHER_STRAND_KEY));
        Assert.assertEquals(new FisherStrand().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.FISHER_STRAND_KEY)));
    }

    @DataProvider(name = "UsingOverflowTable")
    public Object[][] makeUsingOverflowTable() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE});
        tests.add(new Object[]{0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE});
        tests.add(new Object[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 0});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "UsingOverflowTable", expectedExceptions = ArithmeticException.class)
    public void testUsingOverflowTable(final int refpos, final int refneg, final int altpos, final int altneg) {
        final int[][] contingencyTable = new int[2][2];
        contingencyTable[0][0] = refpos;
        contingencyTable[0][1] = refneg;
        contingencyTable[1][0] = altpos;
        contingencyTable[1][1] = altneg;
        final Double p = FisherStrand.pValueForContingencyTable(contingencyTable);//this call will overflow
    }

    @Test
    public void testUsingGT() {
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
        final int depth = 1+2+3+4;

        final String sbbs = "1,2,3,4";
        final int[][] table = {{1, 2}, {3, 4}};

        final Genotype gAC = new GenotypeBuilder("1", AC).DP(depth).attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, sbbs).make();

        final double log10PError = -5;

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError)
                .genotypes(Arrays.asList(gAC))
                .make();

        final Map<String, Object> annotatedMap = new FisherStrand().annotate(null, vc, null);
        Assert.assertNotNull(annotatedMap, vc.toString());
        final String pPhredStr = (String)annotatedMap.get(GATKVCFConstants.FISHER_STRAND_KEY);

        final double actualPhredPval = Double.parseDouble(pPhredStr);

        final double expectedPvalue = FisherStrand.pValueForContingencyTable(table);
        final double expectedPhredPvalue =  QualityUtils.phredScaleErrorRate(Math.max(expectedPvalue, FisherStrand.MIN_PVALUE));
        Assert.assertEquals(actualPhredPval, expectedPhredPvalue, DELTA_PRECISION_FOR_PHRED, "phred pvalues");

        Assert.assertEquals(Math.pow(10, actualPhredPval/-10.0), expectedPvalue, DELTA_PRECISION);
    }

    private GATKRead makeRead(final boolean forward) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setIsReverseStrand(!forward);
        return read;
    }

    private final String sample1 = "NA1";
    private VariantContext makeVC( final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    @Test
    public void testUsingMap(){

        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final int[][] table= {{1, 1},  // alt: one read in each direction,
                              {2, 0}}; //ref: 2 reads fwd, 0 reads back

        final GATKRead read1 = makeRead(true);
        final GATKRead read2 = makeRead(true);
        final GATKRead read3 = makeRead(false);
        final GATKRead read4 = makeRead(true);
        map.add(read1, alleleAlt, -1.0);
        map.add(read1, alleleRef, -100.0);

        map.add(read2, alleleAlt, -1.0);
        map.add(read2, alleleRef, -100.0);

        map.add(read3, alleleAlt, -100.0);
        map.add(read3, alleleRef, -1.0);

        map.add(read4, alleleAlt, -100.0);
        map.add(read4, alleleRef, -1.0);

        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap = Collections.singletonMap(sample1, map);


        final VariantContext vc= makeVC(alleleRef, alleleAlt);

        final Map<String, Object> annotatedMap = new FisherStrand().annotate(null, vc, stratifiedPerReadAlleleLikelihoodMap);
        Assert.assertNotNull(annotatedMap, vc.toString());
        final String pPhredStr = (String)annotatedMap.get(GATKVCFConstants.FISHER_STRAND_KEY);

        final double actualPhredPval = Double.parseDouble(pPhredStr);

        final double expectedPvalue = FisherStrand.pValueForContingencyTable(table);
        final double expectedPhredPvalue =  QualityUtils.phredScaleErrorRate(Math.max(expectedPvalue, FisherStrand.MIN_PVALUE));
        Assert.assertEquals(actualPhredPval, expectedPhredPvalue, DELTA_PRECISION_FOR_PHRED, "phred pvalues");

        Assert.assertEquals(Math.pow(10, actualPhredPval / -10.0), expectedPvalue, DELTA_PRECISION);
    }

    @Test
    public void testNullIfNonVariant() {
        FisherStrand fs= new FisherStrand();
        Assert.assertNull(fs.annotate(null, when(mock(VariantContext.class).isVariant()).thenReturn(false).getMock(), null));

    }
}
