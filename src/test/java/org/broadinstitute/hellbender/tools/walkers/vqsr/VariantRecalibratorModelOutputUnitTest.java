package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.collect.Lists;
import org.apache.log4j.Logger;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class VariantRecalibratorModelOutputUnitTest extends GATKBaseTest {
    protected final static Logger logger = Logger.getLogger(VariantRecalibratorModelOutputUnitTest.class);
    private final boolean printTables = false;
    private final int numAnnotations = 6;
    private final double shrinkage = 1.0;
    private final double dirichlet = 0.001;
    private final double priorCounts = 20.0;
    private final double epsilon = 1e-6;
    private final File modelReportFile = createTempFile("vqsr_model", "report");
    private final int dataSize = 1000; //use this a placeholder so we can test model reading without data

    @Test
    public void testVQSRModelOutput() {
        GaussianMixtureModel goodModel = getGoodGMM();
        GaussianMixtureModel badModel = getBadGMM();

        if (printTables) {
            System.out.println("Good model mean matrix:");
            System.out.println(vectorToString(goodModel.getModelGaussians().get(0).mu));
            System.out.println(vectorToString(goodModel.getModelGaussians().get(1).mu));
            System.out.println("\n\n");

            System.out.println("Good model covariance matrices:");
            goodModel.getModelGaussians().get(0).sigma.print(10, 3);
            goodModel.getModelGaussians().get(1).sigma.print(10, 3);
            System.out.println("\n\n");

            System.out.println("Bad model mean matrix:\n");
            System.out.println(vectorToString(badModel.getModelGaussians().get(0).mu));
            System.out.println("\n\n");

            System.out.println("Bad model covariance matrix:");
            badModel.getModelGaussians().get(0).sigma.print(10, 3);
        }

        VariantRecalibrator vqsr = new VariantRecalibrator();
        List<String> annotationList = getAnnotationList();

        GATKReport report = vqsr.writeModelReport(goodModel, badModel, annotationList);
        //this generates input data for testVQSRModelInput
        try(PrintStream modelReporter = new PrintStream(modelReportFile)) {
            report.print(modelReporter);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        //Check values for Gaussian means
        GATKReportTable goodMus = report.getTable("PositiveModelMeans");
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(goodModel.getModelGaussians().get(0).mu[i], (Double)goodMus.get(0,annotationList.get(i)), epsilon);
        }
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(goodModel.getModelGaussians().get(1).mu[i], (Double)goodMus.get(1,annotationList.get(i)), epsilon);
        }

        GATKReportTable badMus = report.getTable("NegativeModelMeans");
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(badModel.getModelGaussians().get(0).mu[i], (Double)badMus.get(0,annotationList.get(i)), epsilon);
        }

        //Check values for Gaussian covariances
        GATKReportTable goodSigma = report.getTable("PositiveModelCovariances");
        for(int i = 0; i < annotationList.size(); i++) {
            for(int j = 0; j < annotationList.size(); j++) {
                Assert.assertEquals(goodModel.getModelGaussians().get(0).sigma.get(i,j), (Double)goodSigma.get(i,annotationList.get(j)), epsilon);
            }
        }

        //add annotationList.size() to row indexes for second Gaussian because the matrices are concatenated by row in the report
        for(int i = 0; i < annotationList.size(); i++) {
            for(int j = 0; j < annotationList.size(); j++) {
                Assert.assertEquals(goodModel.getModelGaussians().get(1).sigma.get(i,j), (Double)goodSigma.get(annotationList.size()+i,annotationList.get(j)), epsilon);
            }
        }

        GATKReportTable badSigma = report.getTable("NegativeModelCovariances");
        for(int i = 0; i < annotationList.size(); i++) {
            for(int j = 0; j < annotationList.size(); j++) {
                Assert.assertEquals(badModel.getModelGaussians().get(0).sigma.get(i,j), (Double)badSigma.get(i,annotationList.get(j)), epsilon);
            }
        }
    }


    @Test (dependsOnMethods = {"testVQSRModelOutput"})
    public void testVQSRModelInput(){
        final GATKReport report = new GATKReport(modelReportFile);

        // Now test model report reading
        // Read all the tables
        final GATKReportTable badMus = report.getTable("NegativeModelMeans");
        final GATKReportTable badSigma = report.getTable("NegativeModelCovariances");
        final GATKReportTable nPMixTable = report.getTable("BadGaussianPMix");

        final GATKReportTable goodMus = report.getTable("PositiveModelMeans");
        final GATKReportTable goodSigma = report.getTable("PositiveModelCovariances");
        final GATKReportTable pPMixTable = report.getTable("GoodGaussianPMix");

        List<String> annotationList = getAnnotationList();
        VariantRecalibrator vqsr = new VariantRecalibrator();

        GaussianMixtureModel goodModelFromFile = vqsr.GMMFromTables(goodMus, goodSigma, pPMixTable, annotationList.size(), dataSize);
        GaussianMixtureModel badModelFromFile = vqsr.GMMFromTables(badMus, badSigma, nPMixTable, annotationList.size(), dataSize);

        testGMMsForEquality(getGoodGMM(), goodModelFromFile, epsilon);
        testGMMsForEquality(getBadGMM(), badModelFromFile, epsilon);
    }

    @Test
    //This is tested separately to avoid setting up a VariantDataManager and populating it with fake data
    public void testAnnotationNormalizationOutput() {
        final VariantRecalibrator vqsr = new VariantRecalibrator();
        final List<String> annotationList = new ArrayList<>();
        annotationList.add("QD");
        annotationList.add("FS");
        annotationList.add("ReadPosRankSum");
        annotationList.add("MQ");
        annotationList.add("MQRankSum");
        annotationList.add("SOR");

        final double epsilon = 1e-6;

        double[] meanVector = {16.13, 2.45, 0.37, 59.08, 0.14, 0.91};
        final String columnName = "Mean";
        final String formatString = "%.3f";
        GATKReportTable vectorTable = vqsr.makeVectorTable("AnnotationMeans", "Mean for each annotation, used to normalize data", annotationList, meanVector, columnName, formatString);
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(meanVector[i], (Double)vectorTable.get(i, columnName), epsilon);
        }

        if (printTables) {
            final GATKReport report = new GATKReport();
            report.addTable(vectorTable);
            report.print(System.out);
        }
    }

    private String vectorToString(double[] meanVec) {
        String returnString = "";
        for (int j = 0; j < meanVec.length; j++) {
            returnString += String.format("%.3f", meanVec[j]);
            if (j < meanVec.length-1)
                returnString += ",";
        }
        return returnString;
    }

    private void testGMMsForEquality(GaussianMixtureModel gmm1, GaussianMixtureModel gmm2, double epsilon){
        Assert.assertEquals(gmm1.getModelGaussians().size(), gmm2.getModelGaussians().size(), 0);

        for(int k = 0; k < gmm1.getModelGaussians().size(); k++) {
            final MultivariateGaussian g = gmm1.getModelGaussians().get(k);
            final MultivariateGaussian gFile = gmm2.getModelGaussians().get(k);

            Assert.assertEquals(g.pMixtureLog10, gFile.pMixtureLog10);

            for(int i = 0; i < g.mu.length; i++){
                Assert.assertEquals(g.mu[i], gFile.mu[i], epsilon);
            }

            for(int i = 0; i < g.sigma.getRowDimension(); i++) {
                for (int j = 0; j < g.sigma.getColumnDimension(); j++) {
                    Assert.assertEquals(g.sigma.get(i, j), gFile.sigma.get(i, j), epsilon);
                }
            }
        }
    }

    private List<String> getAnnotationList(){
        List<String> annotationList = new ArrayList<>();
        annotationList.add("QD");
        annotationList.add("MQ");
        annotationList.add("FS");
        annotationList.add("SOR");
        annotationList.add("ReadPosRankSum");
        annotationList.add("MQRankSum");
        return annotationList;
    }

    private GaussianMixtureModel getGoodGMM(){
        Random rand = new Random(12878);
        MultivariateGaussian goodGaussian1 = new MultivariateGaussian(dataSize, numAnnotations);
        goodGaussian1.initializeRandomMu(rand);
        goodGaussian1.initializeRandomSigma(rand);

        MultivariateGaussian goodGaussian2 = new MultivariateGaussian(dataSize, numAnnotations);
        goodGaussian2.initializeRandomMu(rand);
        goodGaussian2.initializeRandomSigma(rand);

        List<MultivariateGaussian> goodGaussianList = new ArrayList<>();
        goodGaussianList.add(goodGaussian1);
        goodGaussianList.add(goodGaussian2);

        return new GaussianMixtureModel(goodGaussianList, shrinkage, dirichlet, priorCounts);
    }

    private GaussianMixtureModel getBadGMM(){
        Random rand = new Random(12878);
        MultivariateGaussian badGaussian1 = new MultivariateGaussian(dataSize, numAnnotations);

        badGaussian1.initializeRandomMu(rand);
        badGaussian1.initializeRandomSigma(rand);

        List<MultivariateGaussian> badGaussianList = new ArrayList<>();
        badGaussianList.add(badGaussian1);

        return new GaussianMixtureModel(badGaussianList, shrinkage, dirichlet, priorCounts);
    }

    @Test
    public void testAnnotationOrderAndValidate() {
        final VariantRecalibrator vqsr = new VariantRecalibrator();
        final List<String> annotationList = new ArrayList<>();
        annotationList.add("QD");
        annotationList.add("FS");
        annotationList.add("ReadPosRankSum");
        annotationList.add("MQ");
        annotationList.add("MQRankSum");
        annotationList.add("SOR");

        double[] meanVector = {16.13, 2.45, 0.37, 59.08, 0.14, 0.91};
        final String columnName = "Mean";
        final String formatString = "%.3f";
        GATKReportTable annotationTable = vqsr.makeVectorTable("AnnotationMeans", "Mean for each annotation, used to normalize data", annotationList, meanVector, columnName, formatString);
        vqsr.orderAndValidateAnnotations(annotationTable, annotationList);

        double epsilon = 1e-7;
        for (int i = 0; i < vqsr.annotationOrder.size(); i++){
            Assert.assertEquals(i, vqsr.annotationOrder.get(i), epsilon);
        }

        final List<String> reversed = Lists.reverse(annotationList);
        vqsr.orderAndValidateAnnotations(annotationTable, reversed);
        for (int i = 0; i < vqsr.annotationOrder.size(); i++){
            Assert.assertEquals(reversed.size()-i-1, vqsr.annotationOrder.get(i), epsilon);
        }

        // Now break things...
        // we should throw an error if there are too many or too few annotations on the command line.
        annotationList.add("BaseQRankSum");
        Assert.assertThrows(CommandLineException.class, () -> vqsr.orderAndValidateAnnotations(annotationTable, annotationList));

        annotationList.remove(0);
        annotationList.remove(annotationList.size()-1);
        Assert.assertThrows(CommandLineException.class, () -> vqsr.orderAndValidateAnnotations(annotationTable, annotationList));
    }

}