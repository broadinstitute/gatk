package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.apache.log4j.Logger;

import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class VariantRecalibratorModelOutputUnitTest {
    protected final static Logger logger = Logger.getLogger(VariantRecalibratorModelOutputUnitTest.class);
    private final boolean printTables = true;

    @Test
    public void testVQSRModelOutput() {
        final int numAnnotations = 6;
        final double shrinkage = 1.0;
        final double dirichlet = 0.001;
        final double priorCounts = 20.0;
        final double epsilon = 1e-6;

        Random rand = new Random(12878);
        MultivariateGaussian goodGaussian1 = new MultivariateGaussian(1, numAnnotations);
        goodGaussian1.initializeRandomMu(rand);
        goodGaussian1.initializeRandomSigma(rand);

        MultivariateGaussian goodGaussian2 = new MultivariateGaussian(1, numAnnotations);
        goodGaussian2.initializeRandomMu(rand);
        goodGaussian2.initializeRandomSigma(rand);

        MultivariateGaussian badGaussian1 = new MultivariateGaussian(1, numAnnotations);
        badGaussian1.initializeRandomMu(rand);
        badGaussian1.initializeRandomSigma(rand);

        List<MultivariateGaussian> goodGaussianList = new ArrayList<>();
        goodGaussianList.add(goodGaussian1);
        goodGaussianList.add(goodGaussian2);

        List<MultivariateGaussian> badGaussianList = new ArrayList<>();
        badGaussianList.add(badGaussian1);

        GaussianMixtureModel goodModel = new GaussianMixtureModel(goodGaussianList, shrinkage, dirichlet, priorCounts);
        GaussianMixtureModel badModel = new GaussianMixtureModel(badGaussianList, shrinkage, dirichlet, priorCounts);

        if (printTables) {
            logger.info("Good model mean matrix:");
            logger.info(vectorToString(goodGaussian1.mu));
            logger.info(vectorToString(goodGaussian2.mu));
            logger.info("\n\n");

            logger.info("Good model covariance matrices:");
            goodGaussian1.sigma.print(10, 3);
            goodGaussian2.sigma.print(10, 3);
            logger.info("\n\n");

            logger.info("Bad model mean matrix:\n");
            logger.info(vectorToString(badGaussian1.mu));
            logger.info("\n\n");

            logger.info("Bad model covariance matrix:");
            badGaussian1.sigma.print(10, 3);
        }

        VariantRecalibrator vqsr = new VariantRecalibrator();
        List<String> annotationList = new ArrayList<>();
        annotationList.add("QD");
        annotationList.add("MQ");
        annotationList.add("FS");
        annotationList.add("SOR");
        annotationList.add("ReadPosRankSum");
        annotationList.add("MQRankSum");


        GATKReport report = vqsr.writeModelReport(goodModel, badModel, annotationList);
        if(printTables)
            report.print(System.out);

        //Check values for Gaussian means
        GATKReportTable goodMus = report.getTable("PositiveModelMeans");
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(goodGaussian1.mu[i], (Double)goodMus.get(0,annotationList.get(i)), epsilon);
        }
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(goodGaussian2.mu[i], (Double)goodMus.get(1,annotationList.get(i)), epsilon);
        }

        GATKReportTable badMus = report.getTable("NegativeModelMeans");
        for(int i = 0; i < annotationList.size(); i++) {
            Assert.assertEquals(badGaussian1.mu[i], (Double)badMus.get(0,annotationList.get(i)), epsilon);
        }

        //Check values for Gaussian covariances
        GATKReportTable goodSigma = report.getTable("PositiveModelCovariances");
        for(int i = 0; i < annotationList.size(); i++) {
            for(int j = 0; j < annotationList.size(); j++) {
                Assert.assertEquals(goodGaussian1.sigma.get(i,j), (Double)goodSigma.get(i,annotationList.get(j)), epsilon);
            }
        }

        //add annotationList.size() to row indexes for second Gaussian because the matrices are concatenated by row in the report
        for(int i = 0; i < annotationList.size(); i++) {
            for(int j = 0; j < annotationList.size(); j++) {
                Assert.assertEquals(goodGaussian2.sigma.get(i,j), (Double)goodSigma.get(annotationList.size()+i,annotationList.get(j)), epsilon);
            }
        }

        GATKReportTable badSigma = report.getTable("NegativeModelCovariances");
        for(int i = 0; i < annotationList.size(); i++) {
            for(int j = 0; j < annotationList.size(); j++) {
                Assert.assertEquals(badGaussian1.sigma.get(i,j), (Double)badSigma.get(i,annotationList.get(j)), epsilon);
            }
        }
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

}