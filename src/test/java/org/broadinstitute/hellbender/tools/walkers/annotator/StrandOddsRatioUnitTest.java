package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandOddsRatio;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class StrandOddsRatioUnitTest {
    private static double DELTA_PRECISION = 0.001;

    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    @DataProvider(name = "UsingSOR")
    public Object[][] makeUsingSORData() {
        final double LOG_OF_TWO =  0.6931472;
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{0, 0, 0, 0, LOG_OF_TWO});
        tests.add(new Object[]{100000, 100000, 100000, 100000, LOG_OF_TWO}  );
        tests.add(new Object[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, LOG_OF_TWO}  );

        tests.add(new Object[]{0, 0, 100000, 100000, LOG_OF_TWO});
        tests.add(new Object[]{0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE, LOG_OF_TWO});

        tests.add(new Object[]{100000,100000,100000,0, 23.02587});
        tests.add(new Object[]{100,100,100,0, 9.230339});
        tests.add(new Object[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE,0, 42.97513});

        tests.add(new Object[]{13736,9047,41,1433, 7.061479});
        tests.add(new Object[]{66, 14, 64, 4, 2.248203});
        tests.add(new Object[]{351169, 306836, 153739, 2379, 8.066731});
        tests.add(new Object[]{116449, 131216, 289, 16957, 7.898818});
        tests.add(new Object[]{137, 159, 9, 23, 1.664854});
        tests.add(new Object[]{129, 90, 21, 20, 0.4303384});
        tests.add(new Object[]{14054, 9160, 16, 7827, 12.2645});
        tests.add(new Object[]{32803, 9184, 32117, 3283, 2.139932});
        tests.add(new Object[]{2068, 6796, 1133, 0, 14.06701});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "UsingSOR")
    public void testUsingSOR(final int refpos, final int refneg, final int altpos, final int altneg, double expectedOddsRatio ) {
        int[][] contingencyTable = new int[2][2];
        contingencyTable[0][0] = refpos;
        contingencyTable[0][1] = refneg;
        contingencyTable[1][0] = altpos;
        contingencyTable[1][1] = altneg;
        final double ratio = StrandOddsRatio.calculateSOR(contingencyTable);
        Assert.assertEquals(ratio, expectedOddsRatio, DELTA_PRECISION, "Pass");
        Assert.assertEquals(Double.parseDouble((String) new StrandOddsRatio().annotationForOneTable(ratio).get(GATKVCFConstants.STRAND_ODDS_RATIO_KEY)), expectedOddsRatio, DELTA_PRECISION, "Pass");
    }

    @Test
    public void testLabels() {
        Assert.assertEquals(new StrandOddsRatio().getKeyNames(), Collections.singletonList(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));
        Assert.assertEquals(new StrandOddsRatio().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STRAND_ODDS_RATIO_KEY)));
    }

    @Test
    public void testUsingGT() {
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
        final int depth = 20;
        final String sbbs = "1,2,3,4";
        final Genotype gAC = new GenotypeBuilder("1", AC).DP(depth).attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, sbbs).make();

        final double log10PError = -5;

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError)
                .genotypes(Arrays.asList(gAC))
                .make();

        final Map<String, Object> annotatedMap = new StrandOddsRatio().annotate(null, vc, null);
        Assert.assertNotNull(annotatedMap, vc.toString());
        final String sor = (String)annotatedMap.get(GATKVCFConstants.STRAND_ODDS_RATIO_KEY);

        final double expectedSOR = StrandOddsRatio.calculateSOR(new int[][]{{1,2},{3,4}});
        Assert.assertEquals(Double.valueOf(sor), expectedSOR, 0.001);
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
    public void testUsingLikelihoods(){
        final InfoFieldAnnotation ann = new StrandOddsRatio();
        final String key = GATKVCFConstants.STRAND_ODDS_RATIO_KEY;

        final int[][] table= {{2, 0},  //ref: 2 reads fwd, 0 reads back
                              {1, 1}}; // alt: one read in each direction

        final List<GATKRead> refReads = Arrays.asList(makeRead(true), makeRead(true));
        final List<GATKRead> altReads = Arrays.asList(makeRead(false), makeRead(true));
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc= makeVC(REF, ALT);

        final Map<String, Object> annotatedMap = ann.annotate(null, vc, likelihoods);
        final String sor = (String)annotatedMap.get(key);
        final double expectedSOR = StrandOddsRatio.calculateSOR(table);

        Assert.assertEquals(Double.valueOf(sor), expectedSOR, 0.001);
    }

    @Test
    public void testUsingLikelihoods_Raw(){
        final AS_StrandBiasTest ann= new AS_StrandOddsRatio();
        final String key = GATKVCFConstants.AS_SB_TABLE_KEY;

        final int[][] table= {{2, 0},  //ref: 2 reads fwd, 0 reads back
                {1, 1}}; // alt: one read in each direction

        final List<GATKRead> refReads = Arrays.asList(makeRead(true), makeRead(true));
        final List<GATKRead> altReads = Arrays.asList(makeRead(false), makeRead(true));
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);
        final VariantContext vc= makeVC(REF, ALT);

        final Map<String, Object> annotatedMapRaw = ann.annotateRawData(null, vc, likelihoods);
        final Map<String, Object> annotatedMapNonRaw = ann.annotateRawData(null, vc, likelihoods);
        final String sorStringRaw = (String)annotatedMapRaw.get(key);
        final String sorStringNonRaw = (String)annotatedMapNonRaw.get(key);

        final String expectedString = AS_StrandBiasTest.rawValueAsString(table);

        Assert.assertEquals(sorStringRaw, expectedString);
        Assert.assertEquals(sorStringNonRaw, expectedString);
    }
}
