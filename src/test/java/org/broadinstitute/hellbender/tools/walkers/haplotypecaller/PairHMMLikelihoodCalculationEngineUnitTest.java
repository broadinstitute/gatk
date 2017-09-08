package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.base.Strings;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Unit tests for PairHMMLikelihoodCalculationEngine
 */
public final class PairHMMLikelihoodCalculationEngineUnitTest extends BaseTest {

    Allele Aref, T, C, G, Cref, ATC, ATCATC;

    @BeforeClass
    public void setup() {
        // alleles
        Aref = Allele.create("A", true);
        Cref = Allele.create("C", true);
        T = Allele.create("T");
        C = Allele.create("C");
        G = Allele.create("G");
        ATC = Allele.create("ATC");
        ATCATC = Allele.create("ATCATC");
    }

    @Test
    public void testNormalizeDiploidLikelihoodMatrixFromLog10() {
        double[][] log10_likelihoodMatrix = {
            {-90.2,     0,      0},
            {-190.1, -2.1,      0},
            {-7.0,  -17.5,  -35.9}
        };
        double[][] normalizedMatrix = {
            {-88.1,     0,      0},
            {-188.0,  0.0,      0},
            {-4.9,  -15.4,  -33.8}
        };


        Assert.assertTrue(compareDoubleArrays(normalizeDiploidLikelihoodMatrixFromLog10(log10_likelihoodMatrix), normalizedMatrix));

        double[][] log10_likelihoodMatrix2 = {
                {-90.2,     0,      0,        0},
                {-190.1, -2.1,      0,        0},
                {-7.0,  -17.5,  -35.9,        0},
                {-7.0,  -17.5,  -35.9,  -1000.0},
        };
        double[][] normalizedMatrix2 = {
                {-88.1,     0,      0,        0},
                {-188.0,  0.0,      0,        0},
                {-4.9,  -15.4,  -33.8,        0},
                {-4.9,  -15.4,  -33.8,   -997.9},
        };
        Assert.assertTrue(compareDoubleArrays(normalizeDiploidLikelihoodMatrixFromLog10(log10_likelihoodMatrix2), normalizedMatrix2));
    }

    static double[][] normalizeDiploidLikelihoodMatrixFromLog10( final double[][] log10_likelihoodMatrix ) {
        final double[] genotypeLog10Likelihoods = copyLowerTriangleToArray(log10_likelihoodMatrix);
        final double[] normalizedGenotypeLikelihoods = MathUtils.scaleLogSpaceArrayForNumericalStability(genotypeLog10Likelihoods);
        copyArrayToLowerTriangle(normalizedGenotypeLikelihoods, log10_likelihoodMatrix);
        return log10_likelihoodMatrix;
    }

    private static double[] copyLowerTriangleToArray(final double[][] from) {
        final int n = from.length;
        final double[] to = new double[n*(n+1)/2]; //one triangle
        int index = 0;
        for( int i = 0; i < n; i++ ) {
            for( int j = 0; j <= i; j++ ){
                to[index++] = from[i][j];
            }
        }
        return to;
    }

    private static void copyArrayToLowerTriangle(final double[] from, final double[][] to){
        int index = 0;
        for( int i = 0; i < to.length; i++ ) {
            for( int j = 0; j <= i; j++ ){
                to[i][j] = from[index++];
            }
        }
    }


    @DataProvider(name = "PcrErrorModelTestProvider")
    public Object[][] createPcrErrorModelTestData() {
        List<Object[]> tests = new ArrayList<>();

        for ( final String repeat : Arrays.asList("A", "AC", "ACG", "ACGT") ) {
            for ( final int repeatLength : Arrays.asList(1, 2, 3, 5, 10, 15) ) {
                tests.add(new Object[]{repeat, repeatLength});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PcrErrorModelTestProvider", enabled = true)
    public void createPcrErrorModelTest(final String repeat, final int repeatLength) {

        final PairHMMLikelihoodCalculationEngine engine = new PairHMMLikelihoodCalculationEngine((byte)0, new PairHMMNativeArguments(),
                PairHMM.Implementation.ORIGINAL, 0.0,
                PairHMMLikelihoodCalculationEngine.PCRErrorModel.CONSERVATIVE);

        final String readString = Strings.repeat(repeat, repeatLength);
        final byte[] insQuals = new byte[readString.length()];
        final byte[] delQuals = new byte[readString.length()];
        Arrays.fill(insQuals, (byte) PairHMMLikelihoodCalculationEngine.INITIAL_QSCORE);
        Arrays.fill(delQuals, (byte) PairHMMLikelihoodCalculationEngine.INITIAL_QSCORE);

        engine.applyPCRErrorModel(readString.getBytes(), insQuals, delQuals);

        for ( int i = 1; i < insQuals.length; i++ ) {

            final int repeatLengthFromCovariate = PairHMMLikelihoodCalculationEngine.findTandemRepeatUnits(readString.getBytes(), i-1).getRight();
            final byte adjustedScore = PairHMMLikelihoodCalculationEngine.getErrorModelAdjustedQual(repeatLengthFromCovariate, 3.0);

            Assert.assertEquals(insQuals[i - 1], adjustedScore);
            Assert.assertEquals(delQuals[i - 1], adjustedScore);
        }
    }

    //Private function to compare 2d arrays
    private static boolean compareDoubleArrays(final double[][] b1, final double[][] b2) {
        if( b1.length != b2.length ) {
            return false; // sanity check
        }

        for( int i=0; i < b1.length; i++ ){
            if( b1[i].length != b2[i].length) {
                return false; // sanity check
            }
            for( int j=0; j < b1.length; j++ ){
                if ( MathUtils.compareDoubles(b1[i][j], b2[i][j]) != 0 && !Double.isInfinite(b1[i][j]) && !Double.isInfinite(b2[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }

    @Test
    public void testComputeLikelihoods(){
        final LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();

        PairHMMLikelihoodCalculationEngine.writeLikelihoodsToFile = true;

        final ReadLikelihoodCalculationEngine lce = new PairHMMLikelihoodCalculationEngine((byte) SAMUtils.MAX_PHRED_SCORE, new PairHMMNativeArguments(),
                PairHMM.Implementation.LOGLESS_CACHING, MathUtils.logToLog10(QualityUtils.qualToErrorProbLog10(LEAC.phredScaledGlobalReadMismappingRate)),
                PairHMMLikelihoodCalculationEngine.PCRErrorModel.CONSERVATIVE);

        final Map<String, List<GATKRead>> perSampleReadList= new HashMap<>();
        final int n = 10 ;
        final GATKRead read1= ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(n + "M"));
        read1.setMappingQuality(60);
        final String sample1 = "sample1";
        perSampleReadList.put(sample1, Arrays.asList(read1));

        final SampleList samples = new IndexedSampleList(sample1);

        final AssemblyResultSet assemblyResultSet = new AssemblyResultSet();
        final byte[] bases = Strings.repeat("A", n+1).getBytes();
        final Haplotype hap1 = new Haplotype(bases, true);
        hap1.setGenomeLocation(read1);
        assemblyResultSet.add(hap1);

        final byte[] basesModified= bases;
        basesModified[5] = 'C';//different bases
        final Haplotype hap2 = new Haplotype(basesModified, false);
        hap2.setGenomeLocation(read1);//use same loc
        assemblyResultSet.add(hap2);


        final ReadLikelihoods<Haplotype> likes = lce.computeReadLikelihoods(assemblyResultSet, samples, perSampleReadList);
        final LikelihoodMatrix<Haplotype> mtx = likes.sampleMatrix(0);

        Assert.assertEquals(mtx.numberOfAlleles(), 2);
        Assert.assertEquals(mtx.numberOfReads(), 1);
        final double v1 = mtx.get(0, 0);
        final double v2 = mtx.get(1, 0);

        Assert.assertTrue(v1 > v2, "matching haplotype should have a higher likelihood");
        lce.close();
        new File(PairHMMLikelihoodCalculationEngine.LIKELIHOODS_FILENAME).delete();
    }
}
