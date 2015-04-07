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

package org.broadinstitute.hellbender.tools.walkers.vqsr;

import Jama.Matrix;
import org.apache.commons.math.special.Gamma;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.collections.ExpandingArrayList;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class MultivariateGaussian {
    public double pMixtureLog10;
    public double sumProb;
    final public double[] mu;
    final public Matrix sigma;
    public double hyperParameter_a;
    public double hyperParameter_b;
    public double hyperParameter_lambda;
    private double cachedDenomLog10;
    private Matrix cachedSigmaInverse;
    final private ExpandingArrayList<Double> pVarInGaussian;

    public MultivariateGaussian( final int numAnnotations ) {
        mu = new double[numAnnotations];
        sigma = new Matrix(numAnnotations, numAnnotations);
        pVarInGaussian = new ExpandingArrayList<>();
    }

    public void zeroOutMu() {
        Arrays.fill(mu, 0.0);
    }

    public void zeroOutSigma() {
        final double[][] zeroSigma = new double[mu.length][mu.length];
        for( final double[] row : zeroSigma ) {
            Arrays.fill(row, 0);
        }
        final Matrix tmp = new Matrix(zeroSigma);
        sigma.setMatrix(0, mu.length - 1, 0, mu.length - 1, tmp);
    }

    public void initializeRandomMu( final Random rand ) {
        for( int jjj = 0; jjj < mu.length; jjj++ ) {
            mu[jjj] = -4.0 + 8.0 * rand.nextDouble();
        }
    }

    public void initializeRandomSigma( final Random rand ) {
        final double[][] randSigma = new double[mu.length][mu.length];
        for( int iii = 0; iii < mu.length; iii++ ) {
            for( int jjj = iii; jjj < mu.length; jjj++ ) {
                randSigma[jjj][iii] = 0.55 + 1.25 * rand.nextDouble();
                if( rand.nextBoolean() ) {
                    randSigma[jjj][iii] *= -1.0;
                }
                if( iii != jjj ) { randSigma[iii][jjj] = 0.0; } // Sigma is a symmetric, positive-definite matrix created by taking a lower diagonal matrix and multiplying it by its transpose
            }
        }
        Matrix tmp = new Matrix( randSigma );
        tmp = tmp.times(tmp.transpose());
        sigma.setMatrix(0, mu.length - 1, 0, mu.length - 1, tmp);
    }

    public double calculateDistanceFromMeanSquared( final VariantDatum datum ) {
        return MathUtils.distanceSquared( datum.annotations, mu );
    }

    public void incrementMu( final VariantDatum datum ) {
        incrementMu( datum, 1.0 );
    }
    
    public void incrementMu( final VariantDatum datum, final double prob ) {
        for( int jjj = 0; jjj < mu.length; jjj++ ) {
            mu[jjj] += prob * datum.annotations[jjj];
        }
    }

    public void divideEqualsMu( final double x ) {
        for( int jjj = 0; jjj < mu.length; jjj++ ) {
            mu[jjj] /= x;
        }
    }

    private void precomputeInverse() {
        try {
            cachedSigmaInverse = sigma.inverse();
        } catch( Exception e ) {
            throw new UserException("Error during clustering. Most likely there are too few variants used during Gaussian mixture modeling. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example).");
        }
    }


    public void precomputeDenominatorForEvaluation() {
        precomputeInverse();
        cachedDenomLog10 = Math.log10(Math.pow(2.0 * Math.PI, -1.0 * ((double) mu.length) / 2.0)) + Math.log10(Math.pow(sigma.det(), -0.5)) ;
    }

    public void precomputeDenominatorForVariationalBayes( final double sumHyperParameterLambda ) {

        // Variational Bayes calculations from Bishop
        precomputeInverse();
        cachedSigmaInverse.timesEquals( hyperParameter_a );
        double sum = 0.0;
        for(int jjj = 1; jjj <= mu.length; jjj++) {
            sum += Gamma.digamma( (hyperParameter_a + 1.0 - jjj) / 2.0 );
        }
        sum -= Math.log(sigma.det());
        sum += Math.log(2.0) * mu.length;
        final double lambda = 0.5 * sum;
        final double pi = Gamma.digamma( hyperParameter_lambda ) - Gamma.digamma( sumHyperParameterLambda );
        final double beta = (-1.0 * mu.length) / (2.0 * hyperParameter_b);
        cachedDenomLog10 = (pi / Math.log(10.0)) + (lambda / Math.log(10.0)) + (beta / Math.log(10.0));
    }

    public double evaluateDatumLog10( final VariantDatum datum ) {
        double sumKernel = 0.0;
        final double[] crossProdTmp = new double[mu.length];
        Arrays.fill(crossProdTmp, 0.0);
        for( int iii = 0; iii < mu.length; iii++ ) {
            for( int jjj = 0; jjj < mu.length; jjj++ ) {
                crossProdTmp[iii] += (datum.annotations[jjj] - mu[jjj]) * cachedSigmaInverse.get(jjj, iii);
            }
        }
        for( int iii = 0; iii < mu.length; iii++ ) {
            sumKernel += crossProdTmp[iii] * (datum.annotations[iii] - mu[iii]);
        }
        
        return (( -0.5 * sumKernel ) / Math.log(10.0)) + cachedDenomLog10; // This is the definition of a Gaussian PDF Log10
    }

    public void assignPVarInGaussian( final double pVar ) {
        pVarInGaussian.add( pVar );
    }

    public void resetPVarInGaussian() {
        pVarInGaussian.clear();
    }

    public void maximizeGaussian( final List<VariantDatum> data, final double[] empiricalMu, final Matrix empiricalSigma,
                                  final double SHRINKAGE, final double DIRICHLET_PARAMETER, final double DEGREES_OF_FREEDOM ) {
        sumProb = 1E-10;
        final Matrix wishart = new Matrix(mu.length, mu.length);
        zeroOutMu();
        zeroOutSigma();
        
        int datumIndex = 0;
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            sumProb += prob;
            incrementMu( datum, prob );
        }
        divideEqualsMu( sumProb );

        final double shrinkageFactor = (SHRINKAGE * sumProb) / (SHRINKAGE + sumProb);
        for( int iii = 0; iii < mu.length; iii++ ) {
            for( int jjj = 0; jjj < mu.length; jjj++ ) {
                wishart.set(iii, jjj, shrinkageFactor * (mu[iii] - empiricalMu[iii]) * (mu[jjj] - empiricalMu[jjj]));
            }
        }

        datumIndex = 0;
        final Matrix pVarSigma = new Matrix(mu.length, mu.length);
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            for( int iii = 0; iii < mu.length; iii++ ) {
                for( int jjj = 0; jjj < mu.length; jjj++ ) {
                    pVarSigma.set(iii, jjj, prob * (datum.annotations[iii]-mu[iii]) * (datum.annotations[jjj]-mu[jjj]));
                }
            }
            sigma.plusEquals( pVarSigma );
        }

        sigma.plusEquals( empiricalSigma );
        sigma.plusEquals( wishart );

        for( int iii = 0; iii < mu.length; iii++ ) {
            mu[iii] = (sumProb * mu[iii] + SHRINKAGE * empiricalMu[iii]) / (sumProb + SHRINKAGE);
        }

        hyperParameter_a = sumProb + DEGREES_OF_FREEDOM;
        hyperParameter_b = sumProb + SHRINKAGE;
        hyperParameter_lambda = sumProb + DIRICHLET_PARAMETER;

        resetPVarInGaussian(); // clean up some memory
    }

    public void evaluateFinalModelParameters( final List<VariantDatum> data ) {
        sumProb = 0.0;
        zeroOutMu();
        zeroOutSigma();

        int datumIndex = 0;
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            sumProb += prob;
            incrementMu( datum, prob );
        }
        divideEqualsMu( sumProb );

        datumIndex = 0;
        final Matrix pVarSigma = new Matrix(mu.length, mu.length);
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            for( int iii = 0; iii < mu.length; iii++ ) {
                for( int jjj = 0; jjj < mu.length; jjj++ ) {
                    pVarSigma.set(iii, jjj, prob * (datum.annotations[iii]-mu[iii]) * (datum.annotations[jjj]-mu[jjj]));
                }
            }
            sigma.plusEquals( pVarSigma );
        }
        sigma.timesEquals( 1.0 / sumProb );

        resetPVarInGaussian(); // clean up some memory
    }
}