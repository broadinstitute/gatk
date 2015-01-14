/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.hellbender.tools.recalibration;


import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;


public class RecalDatumUnitTest extends BaseTest {

    // --------------------------------------------------------------------------------
    //
    // merge case Provider
    //
    // --------------------------------------------------------------------------------

    private class RecalDatumTestProvider extends TestDataProvider {
        int exError, exTotal, reportedQual;

        private RecalDatumTestProvider(int E, int N, int reportedQual) {
            super(RecalDatumTestProvider.class);

            this.exError = E;
            this.exTotal = N;
            this.reportedQual = reportedQual;
        }

        public double getErrorRate() {
            return (exError + 1) / (1.0 * (exTotal + 2));
        }

        public double getErrorRatePhredScaled() {
            return QualityUtils.phredScaleErrorRate(getErrorRate());
        }

        public int getReportedQual() {
            return reportedQual;
        }

        public RecalDatum makeRecalDatum() {
            return new RecalDatum((long)exTotal, (double)exError, (byte)getReportedQual());
        }

        @Override
        public String toString() {
            return String.format("exError=%d, exTotal=%d, reportedQual=%d", exError, exTotal, reportedQual);
        }
    }

    private static boolean createdDatumTestProviders = false;

    @DataProvider(name = "RecalDatumTestProvider")
    public Object[][] makeRecalDatumTestProvider() {
        if ( !createdDatumTestProviders ) {
            for ( int E : Arrays.asList(1, 10, 100, 1000, 10000) )
                for ( int N : Arrays.asList(10, 100, 1000, 10000, 100000, 1000000) )
                    for ( int reportedQual : Arrays.asList(10, 20) )
                        if ( E <= N )
                            new RecalDatumTestProvider(E, N, reportedQual);
            createdDatumTestProviders = true;
        }

        return RecalDatumTestProvider.getTests(RecalDatumTestProvider.class);
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumBasics(RecalDatumTestProvider cfg) {
        final RecalDatum datum = cfg.makeRecalDatum();
        assertBasicFeaturesOfRecalDatum(datum, cfg);
    }

    private static void assertBasicFeaturesOfRecalDatum(final RecalDatum datum, final RecalDatumTestProvider cfg) {
        Assert.assertEquals(datum.getNumMismatches(), cfg.exError, 1E-6);
        Assert.assertEquals(datum.getNumObservations(), cfg.exTotal, 1E-6);
        if ( cfg.getReportedQual() != -1 )
            Assert.assertEquals(datum.getEstimatedQReportedAsByte(), cfg.getReportedQual());
        assertEqualsDoubleSmart(datum.getEmpiricalErrorRate(), cfg.getErrorRate());

        final double e = datum.getEmpiricalQuality();
        Assert.assertTrue(datum.getEmpiricalQualityAsByte() >= Math.floor(e));
        Assert.assertTrue(datum.getEmpiricalQualityAsByte() <= Math.ceil(e));
        Assert.assertNotNull(datum.toString());
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumCopyAndCombine(RecalDatumTestProvider cfg) {
        final RecalDatum datum = cfg.makeRecalDatum();
        final RecalDatum copy = new RecalDatum(datum);
        assertBasicFeaturesOfRecalDatum(copy, cfg);

        RecalDatumTestProvider combinedCfg = new RecalDatumTestProvider(cfg.exError * 2, cfg.exTotal * 2, cfg.reportedQual);
        copy.combine(datum);
        assertBasicFeaturesOfRecalDatum(copy, combinedCfg);
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumModification(RecalDatumTestProvider cfg) {
        RecalDatum datum = cfg.makeRecalDatum();
        datum.setEmpiricalQuality(10.1);
        Assert.assertEquals(datum.getEmpiricalQuality(), 10.1);

        datum.setEstimatedQReported(10.1);
        Assert.assertEquals(datum.getEstimatedQReported(), 10.1);
        Assert.assertEquals(datum.getEstimatedQReportedAsByte(), 10);

        datum = cfg.makeRecalDatum();
        cfg.exTotal = 100000;
        datum.setNumObservations(cfg.exTotal);
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        cfg.exError = 1000;
        datum.setNumMismatches(cfg.exError);
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.increment(true);
        cfg.exError++;
        cfg.exTotal++;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.increment(false);
        cfg.exTotal++;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.incrementNumObservations(2);
        cfg.exTotal += 2;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.incrementNumMismatches(2);
        cfg.exError += 2;
        assertBasicFeaturesOfRecalDatum(datum, cfg);


        datum = cfg.makeRecalDatum();
        datum.increment(10, 5);
        cfg.exError += 5;
        cfg.exTotal += 10;
        assertBasicFeaturesOfRecalDatum(datum, cfg);
    }

    @Test
    public void testNoObs() {
        final RecalDatum rd = new RecalDatum(0L, 0.0, (byte)10);
        Assert.assertEquals(rd.getEmpiricalErrorRate(), 0.0);
    }

    @Test
    public void testlog10QempPrior() {
        for ( int Qemp = 0; Qemp <= QualityUtils.MAX_SAM_QUAL_SCORE; Qemp++ ) {
            for ( int Qrep = 0; Qrep <= QualityUtils.MAX_SAM_QUAL_SCORE; Qrep++ ) {
                final double log10prior = RecalDatum.log10QempPrior(Qemp, Qrep);
                Assert.assertTrue(log10prior < 0.0);
                Assert.assertFalse(Double.isInfinite(log10prior));
                Assert.assertFalse(Double.isNaN(log10prior));
            }
        }

        final int Qrep = 20;
        int maxQemp = -1;
        double maxQempValue = -Double.MAX_VALUE;
        for ( int Qemp = 0; Qemp <= QualityUtils.MAX_SAM_QUAL_SCORE; Qemp++ ) {
            final double log10prior = RecalDatum.log10QempPrior(Qemp, Qrep);
            if ( log10prior > maxQempValue ) {
                maxQemp = Qemp;
                maxQempValue = log10prior;
            }
        }
        Assert.assertEquals(maxQemp, Qrep);
    }

    @Test
    public void testBayesianEstimateOfEmpiricalQuality() {

        final int Qrep = 20;

        // test no shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(0, 0, Qrep), (double) Qrep);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10, 0, Qrep), (double) Qrep);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(1000, 10, Qrep), (double) Qrep);

        // test small shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10, 10, Qrep), Qrep - 1.0);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(1000, 0, Qrep), Qrep + 1.0);

        // test medium shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10000, 0, Qrep), Qrep + 3.0);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10000, 10, Qrep), Qrep + 3.0);

        // test large shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(100000, 10, Qrep), Qrep + 8.0);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(1000000, 10, Qrep), Qrep + 16.0);
    }

    @Test
    public void testlog10QempLikelihood() {

        final double[] Qemps = new double[] { 0.0, 10.0, 20.0, 30.0 };
        final int[] observations = new int[] {0, 10, 1000, 1000000};
        final int[] errors = new int[] {0, 10, 1000, 1000000};

        for ( double Qemp : Qemps ) {
            for ( int observation : observations ) {
                for ( int error : errors ) {
                    if ( error > observation )
                        continue;

                    final double log10likelihood = RecalDatum.log10QempLikelihood(Qemp, observation, error);
                    Assert.assertTrue(observation == 0 ? MathUtils.compareDoubles(log10likelihood, 0.0) == 0 : log10likelihood < 0.0);
                    Assert.assertFalse(Double.isInfinite(log10likelihood));
                    Assert.assertFalse(Double.isNaN(log10likelihood));
                }
            }
        }

        long bigNum = new Long((long) Integer.MAX_VALUE);
        bigNum *= 2L;
        final double log10likelihood = RecalDatum.log10QempLikelihood(30, bigNum, 100000);
        Assert.assertTrue(log10likelihood < 0.0);
        Assert.assertFalse(Double.isInfinite(log10likelihood));
        Assert.assertFalse(Double.isNaN(log10likelihood));
    }

    @Test
    public void basicHierarchicalBayesianQualityEstimateTest() {

        for( double epsilon = 15.0; epsilon <= 60.0; epsilon += 2.0 ) {
            double RG_Q = 45.0;
            RecalDatum RG = new RecalDatum( (long)100000000, (long) (100000000 * 1.0 / (Math.pow(10.0, RG_Q / 10.0))), (byte)RG_Q);
            double Q = 30.0;
            RecalDatum QS = new RecalDatum( (long)100000000, (long) (100000000 * 1.0 / (Math.pow(10.0, Q / 10.0))), (byte)Q);
            RecalDatum COV = new RecalDatum( (long)15, (long) 1, (byte)45.0); // no data here so Bayesian prior has a huge effect on the empirical quality

            // initial epsilon condition shouldn't matter when there are a lot of observations
            Assert.assertEquals(BaseRecalibration.hierarchicalBayesianQualityEstimate(epsilon, RG, QS, Collections.singletonList(COV)), Q, 1E-4);
        }

        for( double epsilon = 15.0; epsilon <= 60.0; epsilon += 2.0 ) {
            double RG_Q = 45.0;
            RecalDatum RG = new RecalDatum( (long)10, (long) (10 * 1.0 / (Math.pow(10.0, RG_Q / 10.0))), (byte)RG_Q);
            double Q = 30.0;
            RecalDatum QS = new RecalDatum( (long)10, (long) (10 * 1.0 / (Math.pow(10.0, Q / 10.0))), (byte)Q);
            RecalDatum COV = new RecalDatum( (long)15, (long) 1, (byte)45.0); // no data here so Bayesian prior has a huge effect on the empirical quality

            // initial epsilon condition dominates when there is no data
            Assert.assertEquals(BaseRecalibration.hierarchicalBayesianQualityEstimate(epsilon, RG, QS, Collections.singletonList(COV)), epsilon, 1E-4);
        }

    }
}