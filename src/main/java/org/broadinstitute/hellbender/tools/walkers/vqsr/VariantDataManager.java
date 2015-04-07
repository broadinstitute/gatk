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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.ExpandingArrayList;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDataManager {
    private List<VariantDatum> data = Collections.emptyList();
    private double[] meanVector;
    private double[] varianceVector; // this is really the standard deviation
    public List<String> annotationKeys;
    private final VariantRecalibratorArgumentCollection VRAC;
    protected final static Logger logger = Logger.getLogger(VariantDataManager.class);
    protected final List<TrainingSet> trainingSets;

    public VariantDataManager( final List<String> annotationKeys, final VariantRecalibratorArgumentCollection VRAC ) {
        this.data = Collections.emptyList();
        this.annotationKeys = new ArrayList<>( annotationKeys );
        this.VRAC = VRAC;
        meanVector = new double[this.annotationKeys.size()];
        varianceVector = new double[this.annotationKeys.size()];
        trainingSets = new ArrayList<>();
    }

    public void setData( final List<VariantDatum> data ) {
        this.data = data;
    }

    public List<VariantDatum> getData() {
        return data;
    }

    public void normalizeData() {
        boolean foundZeroVarianceAnnotation = false;
        for( int iii = 0; iii < meanVector.length; iii++ ) {
            final double theMean = mean(iii, true);
            final double theSTD = standardDeviation(theMean, iii, true);
            logger.info( annotationKeys.get(iii) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );
            if( Double.isNaN(theMean) ) {
                throw new UserException.BadInput("Values for " + annotationKeys.get(iii) + " annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations. See " + HelpConstants.forumPost("discussion/49/using-variant-annotator"));
            }

            foundZeroVarianceAnnotation = foundZeroVarianceAnnotation || (theSTD < 1E-5);
            meanVector[iii] = theMean;
            varianceVector[iii] = theSTD;
            for( final VariantDatum datum : data ) {
                // Transform each data point via: (x - mean) / standard deviation
                datum.annotations[iii] = ( datum.isNull[iii] ? 0.1 * Utils.getRandomGenerator().nextGaussian() : ( datum.annotations[iii] - theMean ) / theSTD );
            }
        }
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput( "Found annotations with zero variance. They must be excluded before proceeding." );
        }

        // trim data by standard deviation threshold and mark failing data for exclusion later
        for( final VariantDatum datum : data ) {
            boolean remove = false;
            for( final double val : datum.annotations ) {
                remove = remove || (Math.abs(val) > VRAC.STD_THRESHOLD);
            }
            datum.failingSTDThreshold = remove;
        }

        // re-order the data by increasing standard deviation so that the results don't depend on the order things were specified on the command line
        // standard deviation over the training points is used as a simple proxy for information content, perhaps there is a better thing to use here
        final List<Integer> theOrder = calculateSortOrder(meanVector);
        annotationKeys = reorderList(annotationKeys, theOrder);
        varianceVector = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(varianceVector), theOrder));
        meanVector = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(meanVector), theOrder));
        for( final VariantDatum datum : data ) {
            datum.annotations = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(datum.annotations), theOrder));
            datum.isNull = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(datum.isNull), theOrder));
        }
        logger.info("Annotations are now ordered by their information content: " + annotationKeys.toString());
    }

    /**
     * Get a list of indices which give the ascending sort order of the data array
     * @param inputVector the data to consider
     * @return a non-null list of integers with length matching the length of the input array
     */
    protected List<Integer> calculateSortOrder(final double[] inputVector) {
        final List<Integer> theOrder = new ArrayList<>(inputVector.length);
        final List<MyDoubleForSorting> toBeSorted = new ArrayList<>(inputVector.length);
        int count = 0;
        for( int iii = 0; iii < inputVector.length; iii++ ) {
            toBeSorted.add(new MyDoubleForSorting(-1.0 * Math.abs(inputVector[iii] - mean(iii, false)), count++));
        }
        Collections.sort(toBeSorted);
        for( final MyDoubleForSorting d : toBeSorted ) {
            theOrder.add(d.originalIndex); // read off the sort order by looking at the index field
        }
        return theOrder;
    }

    // small private class to assist in reading off the new ordering of the annotation array
    private class MyDoubleForSorting implements Comparable<MyDoubleForSorting> {
        final Double myData;
        final int originalIndex;

        public MyDoubleForSorting(final double myData, final int originalIndex) {
            this.myData = myData;
            this.originalIndex = originalIndex;
        }

        @Override
        public int compareTo(final MyDoubleForSorting other) {
            return myData.compareTo(other.myData);
        }
    }

    /**
     * Convenience connector method to work with arrays instead of lists. See ##reorderList##
     */
    private <T> T[] reorderArray(final T[] data, final List<Integer> order) {
        return reorderList(Arrays.asList(data), order).toArray(data);
    }

    /**
     * Reorder the given data list to be in the specified order
     * @param data the data to reorder
     * @param order the new order to use
     * @return a reordered list of data
     */
    private <T> List<T> reorderList(final List<T> data, final List<Integer> order) {
        final List<T> returnList = new ArrayList<>(data.size());
        for( final int index : order ) {
            returnList.add( data.get(index) );
        }
        return returnList;
    }

    /**
     * Convert a normalized point to it's original annotation value
     *
     * norm = (orig - mu) / sigma
     * orig = norm * sigma + mu
     *
     * @param normalizedValue the normalized value of the ith annotation
     * @param annI the index of the annotation value
     * @return the denormalized value for the annotation
     */
    public double denormalizeDatum(final double normalizedValue, final int annI) {
        final double mu = meanVector[annI];
        final double sigma = varianceVector[annI];
        return normalizedValue * sigma + mu;
    }

    public void addTrainingSet( final TrainingSet trainingSet ) {
        trainingSets.add( trainingSet );
    }

    public List<String> getAnnotationKeys() {
        return annotationKeys;
    }

    public boolean checkHasTrainingSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isTraining ) { return true; }
        }
        return false;
    }

    public boolean checkHasTruthSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isTruth ) { return true; }
        }
        return false;
    }

    public List<VariantDatum> getTrainingData() {
        final List<VariantDatum> trainingData = new ExpandingArrayList<>();
        for( final VariantDatum datum : data ) {
            if( datum.atTrainingSite && !datum.failingSTDThreshold ) {
                trainingData.add( datum );
            }
        }
        logger.info( "Training with " + trainingData.size() + " variants after standard deviation thresholding." );
        if( trainingData.size() < VRAC.MIN_NUM_BAD_VARIANTS ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
        } else if( trainingData.size() > VRAC.MAX_NUM_TRAINING_DATA ) {
            logger.warn( "WARNING: Very large training set detected. Downsampling to " + VRAC.MAX_NUM_TRAINING_DATA + " training variants." );
            Collections.shuffle(trainingData, Utils.getRandomGenerator());
            return trainingData.subList(0, VRAC.MAX_NUM_TRAINING_DATA);
        }
        return trainingData;
    }

    public List<VariantDatum> selectWorstVariants() {
        final List<VariantDatum> trainingData = new ExpandingArrayList<>();

        for( final VariantDatum datum : data ) {
            if( datum != null && !datum.failingSTDThreshold && !Double.isInfinite(datum.lod) && datum.lod < VRAC.BAD_LOD_CUTOFF ) {
                datum.atAntiTrainingSite = true;
                trainingData.add( datum );
            }
        }

        logger.info( "Training with worst " + trainingData.size() + " scoring variants --> variants with LOD <= " + String.format("%.4f", VRAC.BAD_LOD_CUTOFF) + "." );

        return trainingData;
    }

    public List<VariantDatum> getEvaluationData() {
        final List<VariantDatum> evaluationData = new ExpandingArrayList<>();

        for( final VariantDatum datum : data ) {
            if( datum != null && !datum.failingSTDThreshold && !datum.atTrainingSite && !datum.atAntiTrainingSite ) {
                evaluationData.add( datum );
            }
        }

        return evaluationData;
    }

    /**
     * Remove all VariantDatum's from the data list which are marked as aggregate data
     */
    public void dropAggregateData() {
        final Iterator<VariantDatum> iter = data.iterator();
        while (iter.hasNext()) {
            final VariantDatum datum = iter.next();
            if( datum.isAggregate ) {
                iter.remove();
            }
        }
    }

    public List<VariantDatum> getRandomDataForPlotting( final int numToAdd, final List<VariantDatum> trainingData, final List<VariantDatum> antiTrainingData, final List<VariantDatum> evaluationData ) {
        final List<VariantDatum> returnData = new ExpandingArrayList<>();
        Collections.shuffle(trainingData, Utils.getRandomGenerator());
        Collections.shuffle(antiTrainingData, Utils.getRandomGenerator());
        Collections.shuffle(evaluationData, Utils.getRandomGenerator());
        returnData.addAll(trainingData.subList(0, Math.min(numToAdd, trainingData.size())));
        returnData.addAll(antiTrainingData.subList(0, Math.min(numToAdd, antiTrainingData.size())));
        returnData.addAll(evaluationData.subList(0, Math.min(numToAdd, evaluationData.size())));
        Collections.shuffle(returnData, Utils.getRandomGenerator());
        return returnData;
    }

    protected double mean( final int index, final boolean trainingData ) {
        double sum = 0.0;
        int numNonNull = 0;
        for( final VariantDatum datum : data ) {
            if( (trainingData == datum.atTrainingSite) && !datum.isNull[index] ) { sum += datum.annotations[index]; numNonNull++; }
        }
        return sum / ((double) numNonNull);
    }

    protected double standardDeviation( final double mean, final int index, final boolean trainingData ) {
        double sum = 0.0;
        int numNonNull = 0;
        for( final VariantDatum datum : data ) {
            if( (trainingData == datum.atTrainingSite) && !datum.isNull[index] ) { sum += ((datum.annotations[index] - mean)*(datum.annotations[index] - mean)); numNonNull++; }
        }
        return Math.sqrt(sum / ((double) numNonNull));
    }

    public void decodeAnnotations( final VariantDatum datum, final VariantContext vc, final boolean jitter ) {
        final double[] annotations = new double[annotationKeys.size()];
        final boolean[] isNull = new boolean[annotationKeys.size()];
        int iii = 0;
        for( final String key : annotationKeys ) {
            isNull[iii] = false;
            annotations[iii] = decodeAnnotation( key, vc, jitter );
            if( Double.isNaN(annotations[iii]) ) { isNull[iii] = true; }
            iii++;
        }
        datum.annotations = annotations;
        datum.isNull = isNull;
    }

    private static double decodeAnnotation( final String annotationKey, final VariantContext vc, final boolean jitter ) {
        double value;

        final double LOG_OF_TWO = 0.6931472;

        try {
            value = vc.getAttributeAsDouble( annotationKey, Double.NaN );
            if( Double.isInfinite(value) ) { value = Double.NaN; }
            if( jitter && annotationKey.equalsIgnoreCase("HaplotypeScore") && MathUtils.compareDoubles(value, 0.0, 0.01) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }
            if( jitter && annotationKey.equalsIgnoreCase("FS") && MathUtils.compareDoubles(value, 0.0, 0.01) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }
            if( jitter && annotationKey.equalsIgnoreCase("InbreedingCoeff") && MathUtils.compareDoubles(value, 0.0, 0.01) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }
            if( jitter && annotationKey.equalsIgnoreCase("SOR") && MathUtils.compareDoubles(value, LOG_OF_TWO, 0.01) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }   //min SOR is 2.0, then we take ln
        } catch( Exception e ) {
            value = Double.NaN; // The VQSR works with missing data by marginalizing over the missing dimension when evaluating the Gaussian mixture model
        }

        return value;
    }

    public void parseTrainingSets( final RefMetaDataTracker tracker, final GenomeLoc genomeLoc, final VariantContext evalVC, final VariantDatum datum, final boolean TRUST_ALL_POLYMORPHIC ) {
        datum.isKnown = false;
        datum.atTruthSite = false;
        datum.atTrainingSite = false;
        datum.atAntiTrainingSite = false;
        datum.prior = 2.0;

        for( final TrainingSet trainingSet : trainingSets ) {
            for( final VariantContext trainVC : tracker.getValues(trainingSet.rodBinding, genomeLoc) ) {
                if( isValidVariant( evalVC, trainVC, TRUST_ALL_POLYMORPHIC ) ) {
                    datum.isKnown = datum.isKnown || trainingSet.isKnown;
                    datum.atTruthSite = datum.atTruthSite || trainingSet.isTruth;
                    datum.atTrainingSite = datum.atTrainingSite || trainingSet.isTraining;
                    datum.prior = Math.max(datum.prior, trainingSet.prior);
                    datum.consensusCount += ( trainingSet.isConsensus ? 1 : 0 );
                }
                if( trainVC != null ) {
                    datum.atAntiTrainingSite = datum.atAntiTrainingSite || trainingSet.isAntiTraining;
                }
            }
        }
    }

    private boolean isValidVariant( final VariantContext evalVC, final VariantContext trainVC, final boolean TRUST_ALL_POLYMORPHIC) {
        return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() && checkVariationClass( evalVC, trainVC ) &&
                (TRUST_ALL_POLYMORPHIC || !trainVC.hasGenotypes() || trainVC.isPolymorphicInSamples());
    }

    protected static boolean checkVariationClass( final VariantContext evalVC, final VariantContext trainVC ) {
        switch( trainVC.getType() ) {
            case SNP:
            case MNP:
                return checkVariationClass( evalVC, VariantRecalibratorArgumentCollection.Mode.SNP );
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return checkVariationClass( evalVC, VariantRecalibratorArgumentCollection.Mode.INDEL );
            default:
                return false;
        }
    }

    protected static boolean checkVariationClass( final VariantContext evalVC, final VariantRecalibratorArgumentCollection.Mode mode ) {
        switch( mode ) {
            case SNP:
                return evalVC.isSNP() || evalVC.isMNP();
            case INDEL:
                return evalVC.isStructuralIndel() || evalVC.isIndel() || evalVC.isMixed() || evalVC.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new ReviewedGATKException( "Encountered unknown recal mode: " + mode );
        }
    }

    public void writeOutRecalibrationTable( final VariantContextWriter recalWriter ) {
        // we need to sort in coordinate order in order to produce a valid VCF
        Collections.sort(data, new Comparator<VariantDatum>() {
            public int compare(VariantDatum vd1, VariantDatum vd2) {
                return vd1.loc.compareTo(vd2.loc);
            }
        });

        // create dummy alleles to be used
        final List<Allele> alleles = Arrays.asList(Allele.create("N", true), Allele.create("<VQSR>", false));

        for( final VariantDatum datum : data ) {
            VariantContextBuilder builder = new VariantContextBuilder("VQSR", datum.loc.getContig(), datum.loc.getStart(), datum.loc.getStop(), alleles);
            builder.attribute(VCFConstants.END_KEY, datum.loc.getStop());
            builder.attribute(GATKVCFConstants.VQS_LOD_KEY, String.format("%.4f", datum.lod));
            builder.attribute(GATKVCFConstants.CULPRIT_KEY, (datum.worstAnnotation != -1 ? annotationKeys.get(datum.worstAnnotation) : "NULL"));

            if ( datum.atTrainingSite ) builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
            if ( datum.atAntiTrainingSite ) builder.attribute(GATKVCFConstants.NEGATIVE_LABEL_KEY, true);

            recalWriter.add(builder.make());
        }
    }
}
