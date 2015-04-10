package org.broadinstitute.hellbender.tools.walkers.vqsr;


import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.UnvalidatingGenomeLoc;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static java.lang.Double.*;
import static java.util.Collections.singletonList;

public class VariantDataManagerUnitTest extends BaseTest {

    @Test
    public final void testCalculateSortOrder() {
        final double passingQual = 400.0;
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();

        VariantDataManager vdm = new VariantDataManager(new String[0], VRAC);

        final List<VariantDatum> theData = new ArrayList<>();
        final VariantDatum datum1 = new VariantDatum();
        datum1.atTrainingSite = true;
        datum1.failingSTDThreshold = false;
        datum1.originalQual = passingQual;
        datum1.annotations = new double[]{0.0,-10.0,10.0};
        datum1.isNull = new boolean[]{false, false, false};
        theData.add(datum1);

        final VariantDatum datum2 = new VariantDatum();
        datum2.atTrainingSite = true;
        datum2.failingSTDThreshold = false;
        datum2.originalQual = passingQual;
        datum2.annotations = new double[]{0.0,-9.0,15.0};
        datum2.isNull = new boolean[]{false, false, false};
        theData.add(datum2);

        final VariantDatum datum3 = new VariantDatum();
        datum3.atTrainingSite = false;
        datum3.failingSTDThreshold = false;
        datum3.originalQual = passingQual;
        datum3.annotations = new double[]{0.0,1.0,999.0};
        datum3.isNull = new boolean[]{false, false, false};
        theData.add(datum3);

        final VariantDatum datum4 = new VariantDatum();
        datum4.atTrainingSite = false;
        datum4.failingSTDThreshold = false;
        datum4.originalQual = passingQual;
        datum4.annotations = new double[]{0.015,2.0,1001.11};
        datum4.isNull = new boolean[]{false, false, false};
        theData.add(datum4);

        vdm.addData(theData);

        final double[] meanVector = new double[3];
        for( int iii = 0; iii < meanVector.length; iii++ ) {
            meanVector[iii] = vdm.mean(iii, true);
        }
        final int[] order = vdm.calculateSortOrder(meanVector);
        Assert.assertEquals(new int[]{2, 1, 0}, order);
    }

    @Test
    public final void testDownSamplingTrainingData() {
        final int MAX_NUM_TRAINING_DATA = 5000;
        final double passingQual = 400.0;
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.MAX_NUM_TRAINING_DATA = MAX_NUM_TRAINING_DATA;

        VariantDataManager vdm = new VariantDataManager(new String[0], VRAC);
        final List<VariantDatum> theData = new ArrayList<>();
        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 10; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = true;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 2; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = false;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        vdm.addData(theData);
        final List<VariantDatum> trainingData = vdm.getTrainingData();

        Assert.assertTrue(trainingData.size() == MAX_NUM_TRAINING_DATA);
    }

    @Test
    public final void testDropAggregateData() {
        final int MAX_NUM_TRAINING_DATA = 5000;
        final double passingQual = 400.0;
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.MAX_NUM_TRAINING_DATA = MAX_NUM_TRAINING_DATA;

        VariantDataManager vdm = new VariantDataManager(new String[0], VRAC);
        final List<VariantDatum> theData = new ArrayList<>();
        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 10; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = true;
            datum.isAggregate = false;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 2; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = false;
            datum.isAggregate = true;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        vdm.addData(theData);
        vdm.dropAggregateData();

        Assert.assertTrue(vdm.getData().stream().noneMatch(d -> d.isAggregate));
    }

    @Test
    public void testDataStatistics(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD", "HRun"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double mean1 = 2.0;
        double sd1= 1.5;
        NormalDistribution normal1 = new NormalDistribution(rng, mean1, sd1);
        double mean2= 1.0;
        double sd2= 1.0;
        NormalDistribution normal2 = new NormalDistribution(rng, mean2, sd2);

        int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.annotations = new double[]{normal1.sample(),normal2.sample()};
            vd.atTrainingSite = false;
            vd.isNull = new boolean[2];
            theData.add(vd);
        }
        vdm.addData(theData);
        final double sampleMean1 = vdm.mean(0, false);
        final double sampleMean2 = vdm.mean(1, false);

        Assert.assertEquals(sampleMean1, mean1, sd1 / Math.sqrt(N), "wrong mean 1");
        Assert.assertEquals(sampleMean2, mean2, sd2 / Math.sqrt(N), "wrong mean 2");

        //TODO not sure what better value to put instead of 1/100
        Assert.assertEquals(vdm.standardDeviation(sampleMean1, 0, false), sd1, 1.0 / 10, "wrong stddev 1");
        Assert.assertEquals(vdm.standardDeviation(sampleMean2, 1, false), sd2, 1.0 / 10, "wrong stddev 2");
    }

    @Test
    public void testDataStatisticsAtTrainingOrNonTrainingSites() {
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMeanTrainingSites = 10.0;
        double realSDTrainingSites= 5.0;
        NormalDistribution normalTraining = new NormalDistribution(rng, realMeanTrainingSites, realSDTrainingSites);

        double realMeanNonTrainingSites= 1.0;
        double readSDNonTrainingSites= 1.0;
        NormalDistribution normalNonTraining = new NormalDistribution(rng, realMeanNonTrainingSites, readSDNonTrainingSites);

        final int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;

            final boolean trainingSite = rng.nextBoolean();
            vd.atTrainingSite = trainingSite;
            if (trainingSite) {
                vd.annotations = new double[]{normalTraining.sample()};
            } else {
                vd.annotations = new double[]{normalNonTraining.sample()};
            }
            vd.isNull = new boolean[1];
            theData.add(vd);
        }
        vdm.addData(theData);
        final double sampleMeanTrainingSites = vdm.mean(0, true);
        final double sampleMeanNonTrainingSites = vdm.mean(0, false);

        Assert.assertEquals(sampleMeanTrainingSites, realMeanTrainingSites, realSDTrainingSites / Math.sqrt(N), "wrong mean training sites");
        Assert.assertEquals(sampleMeanNonTrainingSites, realMeanNonTrainingSites, readSDNonTrainingSites / Math.sqrt(N), "wrong mean non training sites");

        //TODO not sure what better value to put instead of 1/100
        Assert.assertEquals(vdm.standardDeviation(sampleMeanTrainingSites, 0, true), realSDTrainingSites, 1.0 / 10, "wrong stddev training sites");
        Assert.assertEquals(vdm.standardDeviation(sampleMeanNonTrainingSites, 0, false), readSDNonTrainingSites, 1.0 / 10, "wrong stddev non training sites");
    }

    @Test
    public void testDataStatisticsSkipNullAnnotations() {
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMeanTrainingSites = 10.0;
        double realSDTrainingSites= 5.0;
        NormalDistribution normalTraining = new NormalDistribution(rng, realMeanTrainingSites, realSDTrainingSites);

        double realMeanNullSites= 1.0;
        double readSDNonNullSites= 1.0;
        NormalDistribution normalNullSites = new NormalDistribution(rng, realMeanNullSites, readSDNonNullSites);

        int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.atTrainingSite = true;
            vd.isNull = new boolean[1];
            final boolean isNull = rng.nextBoolean();

            vd.isNull[0] = isNull;
            if (isNull) {
                vd.annotations = new double[]{normalNullSites.sample()};
            } else {
                vd.annotations = new double[]{normalTraining.sample()};
            }
            theData.add(vd);
        }
        vdm.addData(theData);
        final double sampleMeanTrainingSites = vdm.mean(0, true);
        Assert.assertEquals(sampleMeanTrainingSites, realMeanTrainingSites, realSDTrainingSites / Math.sqrt(N), "wrong mean training sites");

        //TODO not sure what better value to put instead of 1/100
        Assert.assertEquals(vdm.standardDeviation(sampleMeanTrainingSites, 0, true), realSDTrainingSites, 1.0 / 10, "wrong stddev training sites");
    }


    @Test
    public void testNormalizeDataRemovesDataBySTDThreshold(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMeanTrainingSites = 2.0;
        double realSDTrainingSites= 2.0;
        NormalDistribution normalTraining = new NormalDistribution(rng, realMeanTrainingSites, realSDTrainingSites);

        final Map<VariantDatum, Double> hugeValues = new HashMap<>();
        final int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;

            final boolean randomHugeValue = rng.nextFloat() < 0.001; //0.1% of bad values
            if (randomHugeValue){
                vd.annotations = new double[]{VRAC.STD_THRESHOLD * 2 * realSDTrainingSites + normalTraining.sample()};
                hugeValues.put(vd, vd.annotations[0]);
            } else {
                vd.annotations = new double[]{normalTraining.sample()};
            }
            vd.isNull = new boolean[1];
            vd.atTrainingSite = true;
            theData.add(vd);
        }
        vdm.addData(theData);

        vdm.normalizeData();

        //All data with huge values should fail
        Assert.assertTrue(hugeValues.keySet().stream().allMatch(vd -> vd.failingSTDThreshold), "some huge values did not fail STD threshold");

        //no data with not huge values should fail
        Assert.assertTrue(vdm.getData().stream().filter(vd -> !hugeValues.containsKey(vd)).noneMatch(vd -> vd.failingSTDThreshold), "some normal values failed STD threshold");
    }

    @Test
    public void testNormalizeData(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMeanTrainingSites = 10.0;
        double realSDTrainingSites= 5.0;
        NormalDistribution normalTraining = new NormalDistribution(rng, realMeanTrainingSites, realSDTrainingSites);

        final int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.annotations = new double[]{normalTraining.sample()};
            vd.isNull = new boolean[1];
            vd.atTrainingSite = true;
            theData.add(vd);
        }
        vdm.addData(theData);

        final double sampleMeanBeforeNorm = vdm.mean(0, true);
        Assert.assertEquals(sampleMeanBeforeNorm, realMeanTrainingSites, realSDTrainingSites / Math.sqrt(N), "wrong mean training sites");
        Assert.assertEquals(vdm.standardDeviation(sampleMeanBeforeNorm, 0, true), realSDTrainingSites, 1.0 / 10, "wrong stddev training sites");

        vdm.normalizeData();

        //After normalization, mean should be 0.0 and std should be 1.0
        final double sampleMeanAfterNorm = vdm.mean(0, true);
        Assert.assertEquals(sampleMeanAfterNorm, 0.0, realSDTrainingSites / Math.sqrt(N), "wrong mean after normalization");
        Assert.assertEquals(vdm.standardDeviation(sampleMeanAfterNorm, 0, true), 1.0, 1.0 / 10, "wrong stddev after normalization");
    }

    @Test
    public void testNormalizeDataAndDenormalizeDataItems(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMeanTrainingSites = 10.0;
        double realSDTrainingSites= 5.0;
        NormalDistribution normalTraining = new NormalDistribution(rng, realMeanTrainingSites, realSDTrainingSites);

        final int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.annotations = new double[]{normalTraining.sample()};
            vd.isNull = new boolean[1];
            vd.atTrainingSite = true;
            theData.add(vd);
        }
        vdm.addData(theData);
        vdm.normalizeData();

        Assert.assertEquals(vdm.denormalizeDatum(0.0, 0), realMeanTrainingSites, realSDTrainingSites / 10); //the mean
        Assert.assertEquals(vdm.denormalizeDatum(1.0, 0), realMeanTrainingSites + 1.0 * realSDTrainingSites, realSDTrainingSites / 10); //the mean  + 1.0 * stdDev
        Assert.assertEquals(vdm.denormalizeDatum(-1.0, 0), realMeanTrainingSites - 1.0 * realSDTrainingSites, realSDTrainingSites / 10); //the mean - 1.0 * stdDev
    }

    @Test
    public void testNormalizeDataFindWorstItemsItems(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMeanTrainingSites = 2.0;
        double realSDTrainingSites= 2.0;
        NormalDistribution normalDataDistro = new NormalDistribution(rng, realMeanTrainingSites, realSDTrainingSites);

        final int NgoodTraining = 10;
        for (int i = 0; i < NgoodTraining; i++){
            VariantDatum vd = new VariantDatum();
            vd.annotations = new double[]{normalDataDistro.sample()};
            vd.isNull = new boolean[1];
            vd.atTrainingSite = true;
            vd.lod = 0.0;
            theData.add(vd);
        }

        final int Ngood = 10;
        final List<VariantDatum> knownGoods = new ArrayList<>(Ngood);
        for (int i = 0; i < Ngood; i++){
            VariantDatum vd = new VariantDatum();
            vd.annotations = new double[]{normalDataDistro.sample()};
            vd.isNull = new boolean[1];
            vd.atTrainingSite = false;
            vd.lod = 0.0;
            theData.add(vd);
            knownGoods.add(vd);
        }
        final int Nbad = 1;
        final List<VariantDatum> knownBads = new ArrayList<>(Nbad);
        for (int i = 0; i < Nbad; i++){
            VariantDatum vd = new VariantDatum();
            vd.annotations = new double[]{normalDataDistro.sample()};
            vd.isNull = new boolean[1];
            vd.atTrainingSite = false;
            vd.lod = VRAC.BAD_LOD_CUTOFF - 1;
            theData.add(vd);
            knownBads.add(vd);
        }
        final int NbadButExcluded = 1;
        for (int i = 0; i < NbadButExcluded; i++){
            VariantDatum vd = new VariantDatum();
            //make outliers
            vd.annotations = new double[]{ realMeanTrainingSites  + normalDataDistro.sample() *realSDTrainingSites*VRAC.STD_THRESHOLD*2};
            vd.isNull = new boolean[1];
            vd.atTrainingSite = false;
            vd.lod = VRAC.BAD_LOD_CUTOFF - 1;
            theData.add(vd);
        }

        vdm.addData(theData);
        vdm.normalizeData();

        final List<VariantDatum> worst = vdm.selectWorstVariants();
        Assert.assertEquals(new HashSet<>(worst), new HashSet<>(knownBads));
        Assert.assertTrue(worst.stream().allMatch(d -> d.atAntiTrainingSite));
        Assert.assertTrue(vdm.getData().stream().filter(d -> !worst.contains(d)).noneMatch(d -> d.atAntiTrainingSite));

        final List<VariantDatum> eval = vdm.getEvaluationData();
        Assert.assertEquals(new HashSet<>(eval), new HashSet<>(knownGoods));
    }

    @Test
    public void testNormalizeDataOrdersAnnotationsByVariance(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"FOO", "BAR", "BAZ"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        int SEED = 1337;
        RandomGenerator rng= new Well19937c(SEED);

        double realMean1 = 1.0;
        double realStdDev1= 1.0;
        NormalDistribution normal1 = new NormalDistribution(rng, realMean1, realStdDev1);

        double realMean2 = 3.0;
        double realStdDev2= 3.0;
        NormalDistribution normal2 = new NormalDistribution(rng, realMean2, realStdDev2);

        double realMean3 = 2.0;
        double realStdDev3= 2.0;
        NormalDistribution normal3 = new NormalDistribution(rng, realMean3, realStdDev3);

        final int N = 10000;
        for (int i = 0; i < N; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.annotations = new double[]{normal1.sample(), normal2.sample(), normal3.sample()};
            vd.isNull = new boolean[3];
            vd.atTrainingSite = i % 10 == 0;
            theData.add(vd);
        }
        vdm.addData(theData);

        //before normalization, nothing changes
        Assert.assertEquals(vdm.getAnnotationKeys(), new String[]{"FOO", "BAR", "BAZ"});
        vdm.normalizeData();

        //after normalization, order is by decreasing variance of the original data
        Assert.assertEquals(vdm.getAnnotationKeys(), new String[]{"BAR", "BAZ", "FOO"});
        final double[] trainingMeans = vdm.getOriginalTrainingMeans();
        //the means are the same as variances in this test so we test the means here
        Assert.assertTrue(trainingMeans[0] > trainingMeans[1] && trainingMeans[1] > trainingMeans[2], Arrays.toString(trainingMeans));

        final double[] trainingStdDevs = vdm.getOriginalStdDevs();
        Assert.assertTrue(trainingStdDevs[0] > trainingStdDevs[1] && trainingStdDevs[1] > trainingStdDevs[2], Arrays.toString(trainingStdDevs));

    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNormalizeDataFailWhenNoDataInTrainingSet(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        for (int i = 0; i < 100; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.annotations = new double[]{2.0};
            vd.isNull = new boolean[1];
            theData.add(vd);
        }
        vdm.addData(theData);
        vdm.normalizeData(); //will fail - no data in training set
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNormalizeDataFailWhenZeroVariance(){
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
        final List<VariantDatum> theData = new ArrayList<>();

        for (int i = 0; i < 100; i++){
            VariantDatum vd = new VariantDatum();
            vd.failingSTDThreshold = false;
            vd.annotations = new double[]{2.0};
            vd.atTrainingSite= true;
            vd.isNull = new boolean[1];
            theData.add(vd);
        }
        vdm.addData(theData);
        vdm.normalizeData(); //will fail - all data is the same (zero variance)
    }

    @Test
    public void testAnnotationDecoding() {
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD", "FS", "HaplotypeScore", "InbreedingCoeff", "SOR"}, VRAC);

        final Map<String, Object> attrs = new HashMap<>();
        attrs.put("QD", "1.0");              //QD requires no jitter,
        attrs.put("FS", "1.1");              //bottoms out at 0.0
        attrs.put("HaplotypeScore", "1.2");  //bottoms out at 0.0
        attrs.put("InbreedingCoeff", "1.3"); //bottoms out at 0.0
        attrs.put("SOR", "4.0");              //bottoms out at ln 2
        VariantContext vc = makeVC(attrs);
        VariantDatum vd = new VariantDatum();
        vdm.decodeAnnotations(vd, vc);

        final double[] anns = vd.annotations;
        Assert.assertEquals(anns, new double[]{1.0, 1.1, 1.2, 1.3, 4.0}, Arrays.toString(anns));
        Assert.assertEquals(vd.isNull, new boolean[5], Arrays.toString(vd.isNull));
    }

    @Test
    public void testAnnotationDecodingWithJitter() {
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD", "FS", "HaplotypeScore", "InbreedingCoeff", "SOR"}, VRAC);

        final Map<String, Object> attrs = new HashMap<>();
        attrs.put("QD", "0.0");              //QD requires no jitter,
        attrs.put("FS", "0.0");              //bottoms out at 0.0
        attrs.put("HaplotypeScore", "0.0");  //bottoms out at 0.0
        attrs.put("InbreedingCoeff", "0.0"); //bottoms out at 0.0
        attrs.put("SOR", String.valueOf(Math.log(2.0)));              //bottoms out at ln 2
        VariantContext vc = makeVC(attrs);
        VariantDatum vd = new VariantDatum();
        vdm.decodeAnnotations(vd, vc);

        final double[] anns = vd.annotations;
        Assert.assertEquals(anns[0], 0.0, "QD");

        //FS not equal but approximately equal
        Assert.assertNotEquals(anns[1], 0.0, "FS");
        Assert.assertEquals(anns[1], 0.0, 6 * 0.01, "FS"); //6 std deviations away from the mean

        //HaplotypeScore
        Assert.assertNotEquals(anns[2], 0.0, "HaplotypeScore");
        Assert.assertEquals(anns[2], 0.0, 6 * 0.01, "HaplotypeScore"); //6 std deviations away from the mean

        //InbreedingCoeff
        Assert.assertNotEquals(anns[3], 0.0, "InbreedingCoeff");
        Assert.assertEquals(anns[3], 0.0, 6*0.01, "InbreedingCoeff"); //6 std deviations away from the mean

        //SOR
        Assert.assertNotEquals(anns[4], Math.log(2.0), "SOR");
        Assert.assertEquals(anns[4], Math.log(2.0), 6 * 0.01, "SOR"); //6 std deviations away from the mean

        Assert.assertEquals(vd.isNull, new boolean[5], Arrays.toString(vd.isNull));
    }

    @Test
    public void testAnnotationDecodingInfinite() {
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VariantDataManager vdm = new VariantDataManager(new String[]{"QD1", "QD2"}, VRAC);

        final Map<String, Object> attrs = new HashMap<>();
        attrs.put("QD1", String.valueOf(NEGATIVE_INFINITY));
        attrs.put("QD2", String.valueOf(POSITIVE_INFINITY));
        VariantContext vc = makeVC(attrs);
        VariantDatum vd = new VariantDatum();
        vdm.decodeAnnotations(vd, vc);

        final double[] anns = vd.annotations;
        Assert.assertEquals(anns, new double[]{NaN, NaN}, Arrays.toString(anns));
        Assert.assertEquals(vd.isNull, new boolean[]{true, true}, Arrays.toString(vd.isNull));
    }

    private static VariantContext makeVC(Map<String, Object> attrs) {
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final List<Allele> AC = Arrays.asList(A, C);
        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(50).make();
        final double[] hetPL = MathUtils.normalizeFromRealSpace(new double[]{0.09, 0.9, 0.01});
        final Genotype acGT = new GenotypeBuilder(base).alleles(AC).AD(new int[]{10,2}).PL(hetPL).GQ(30).make();
        final Set<String> filters = Collections.singleton("FAIL");

        return makeVC("source", AC, singletonList(acGT), filters, attrs);
    }

    private static VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes, Set<String> filters, Map<String, Object> attrs) {
        int start = 10;
        int stop = start + alleles.get(0).length() - 1; // alleles.contains(ATC) ? start + 3 : start;

        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(genotypes).filters(filters).attributes(attrs).make();
    }

    @Test
    public void testRecalTable() throws IOException {

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOptions(VariantContextWriterBuilder.NO_OPTIONS);

        final File recalFile = createTempFile("recalTable", ".vcf");
        final File expectedRecalFile = new File(getTestDataDir() + "/" + "expected.recalTable.vcf");

        //auto close
        try (VariantContextWriter recalWriter = builder
                .setOutputFile(recalFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .build()) {

            final Set<VCFHeaderLine> hInfo = new HashSet<>();
            VQSRUtils.addVQSRStandardHeaderLines(hInfo);
            recalWriter.writeHeader(new VCFHeader(hInfo));

            final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
            VariantDataManager vdm = new VariantDataManager(new String[]{"QD"}, VRAC);
            final List<VariantDatum> theData = new ArrayList<>();

            final double[] qds = {1.0, 2.0, 3.0, 4.0, 5.0};
            final double[] lods = {10.0, 20.0, 30.0, 40.0, 50.0};
            final boolean[] training = {true, true, false, true, true};
            for (int i = 0; i < qds.length; i++) {
                VariantDatum vd = new VariantDatum();
                vd.failingSTDThreshold = false;
                vd.annotations = new double[]{qds[i]};
                vd.isNull = new boolean[1];
                vd.atTrainingSite = training[i];
                vd.lod = lods[i];

                int start = 10 * i;
                vd.loc = new UnvalidatingGenomeLoc("1", 0, start, start);
                theData.add(vd);
            }
            vdm.addData(theData);
            vdm.normalizeData();
            vdm.writeOutRecalibrationTable(recalWriter);
        }
        final List<String> expectedLines = new XReadLines(expectedRecalFile).readLines();
        final List<String> actualLines = new XReadLines(recalFile).readLines();
        Assert.assertEquals(expectedLines, actualLines);
    }
}

