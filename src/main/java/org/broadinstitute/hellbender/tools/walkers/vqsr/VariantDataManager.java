package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.DoubleStream;

import static java.lang.Math.*;
import static java.util.Collections.shuffle;
import static java.util.stream.Collectors.toList;
import static org.broadinstitute.hellbender.utils.MathUtils.equalDoubles;
import static org.broadinstitute.hellbender.utils.Utils.getRandomGenerator;

/*
 * Helper class for dealing with data for VQSR.
 * Package-private because it's not usable outside of VQSR.
 *
 * The way to interact with this data structure is that you
 * create it,
 * addData,
 * normalizeData,
 * decodeAnnotations on  data items (can be different data items than those held by the manager)
 *
 */
final class VariantDataManager {
    protected final static Logger logger = LogManager.getLogger(VariantDataManager.class);

    private final List<VariantDatum> data;     //the data
    private final double[] meanVector;         // means of the annotations - used to denormalize values
    private final double[] stdDevVector;       // standard deviations of the annotations - used to denormalize values
    private final String[] annotationKeys;     //names of the annotations
    private final VariantRecalibratorArgumentCollection VRAC;

    public VariantDataManager(final String[] annotationKeys, final VariantRecalibratorArgumentCollection VRAC) {
        this.data = new ArrayList<>(1000);
        this.annotationKeys = Arrays.copyOf(annotationKeys, annotationKeys.length);
        this.VRAC = VRAC;
        this.meanVector = new double[this.annotationKeys.length];
        this.stdDevVector = new double[this.annotationKeys.length];
    }

    /**
     * Copies the data items to the internal list.
     */
    public void addData(final Collection<VariantDatum> data) {
        this.data.addAll(data);
    }

    /**
     * Returns an unmodifiable view of the data list.
     */
    public List<VariantDatum> getData() {
        return Collections.unmodifiableList(data);
    }

    /**
     * Normalizes data so that each annotation has a mean 0 and standard deviation 1.
     * Additionally, it marks data items as outliers (setting failingSTDThreshold to true) if the value is too far off (STD_THRESHOLD)).
     *
     * @throws  UserException.BadInput is any of the annotations has values with zero or near-zero variance
     * (ie. all values are the same or very nearly the same).
     */
    public void normalizeData() {
        //transform each 'dimension' by shifting each data item to zero  mean and dividing by the std dev to set variance to 1.
        for( int i = 0; i < meanVector.length; i++ ) {
            final double theMean = mean(i, true);
            final double theSTD = standardDeviation(theMean, i, true);
            logger.info( annotationKeys[i] + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );

            if( Double.isNaN(theMean) ) {
                throw new UserException.BadInput("Values for " + annotationKeys[i] + " annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations. See " + HelpConstants.forumPost("discussion/49/using-variant-annotator"));
            }
            if( theSTD < 1E-5 ) {
                throw new UserException.BadInput( "Annotation " + annotationKeys[i] + " has zero or nearly zero variance. It needs to be removed." );
            }

            meanVector[i] = theMean;
            stdDevVector[i] = theSTD;
            for( final VariantDatum datum : data ) {
                // Transform each data point via: (x - mean) / standard deviation
                if (datum.isNull[i]) {
                    datum.annotations[i] = 0.1 * getRandomGenerator().nextGaussian(); //XXX why the 0.1 ?
                } else {
                    datum.annotations[i] = (datum.annotations[i] - theMean) / theSTD;
                }
            }
        }

        // trim data by standard deviation threshold and mark failing data for exclusion later
        data.forEach(d -> d.failingSTDThreshold = DoubleStream.of(d.annotations).anyMatch(val -> (abs(val) > VRAC.STD_THRESHOLD)));

        // re-order the data by increasing standard deviation so that the results don't depend on the order things were specified on the command line
        // standard deviation over the training points is used as a simple proxy for information content, perhaps there is a better thing to use here
        //TODO - the code uses non-training data but comment says training.
        final int[] theOrder = calculateSortOrder(meanVector);
        reorderStringArray(annotationKeys, theOrder);
        reorderDoublesArray(stdDevVector, theOrder);
        reorderDoublesArray(meanVector, theOrder);
        for( final VariantDatum datum : data ) {
            reorderDoublesArray(datum.annotations, theOrder);
            reorderBooleanArray(datum.isNull, theOrder);
        }
        logger.info("Annotations are now ordered by their information content: " + Arrays.toString(annotationKeys));
    }


    /*
    * Reorders the array of Strings in place, according to the permutation given as the order.
    */
    private static void reorderStringArray(String[] args, int[] theOrder) {
        final String[] copy = Arrays.copyOf(args, args.length);
        for(int i = 0; i < args.length; i++) {
            args[i] = copy[theOrder[i]];
        }
    }

    /*
     * Reorders the array of booleans in place, according to the permutation given as the order.
     */
    private static void reorderBooleanArray(boolean[] args, int[] theOrder) {
        final boolean[] copy = Arrays.copyOf(args, args.length);
        for(int i = 0; i < args.length; i++) {
            args[i] = copy[theOrder[i]];
        }
    }

    /*
     * Reorders the array of doubles in place, according to the permutation given as the order.
     */
    private static void reorderDoublesArray(double[] args, int[] theOrder) {
        final double[] copy = Arrays.copyOf(args, args.length);
        for(int i = 0; i < args.length; i++) {
            args[i] = copy[theOrder[i]];
        }
    }

    @VisibleForTesting
    protected double[] getOriginalTrainingMeans() {
        return Arrays.copyOf(meanVector, meanVector.length);
    }

    @VisibleForTesting
    protected double[] getOriginalStdDevs() {
        return Arrays.copyOf(stdDevVector, stdDevVector.length);
    }

    /**
     * Returns the copy of the annotation name array.
     */
    public String[] getAnnotationKeys() {
        return Arrays.copyOf(annotationKeys, annotationKeys.length);
    }

    /**
     * Get a list of indices which give the ascending sort order of the data array
     * @param originalMeans the data to consider
     * @return a non-null list of integers with length matching the length of the input array
     */
    @VisibleForTesting
    protected int[] calculateSortOrder(final double[] originalMeans) {
        final List<MyDoubleForSorting> toBeSorted = new ArrayList<>(originalMeans.length);
        for( int i = 0; i < originalMeans.length; i++ ) {
            toBeSorted.add(new MyDoubleForSorting(abs(originalMeans[i] - mean(i, false)), i));
        }
        Collections.sort(toBeSorted);
        Collections.reverse(toBeSorted);

        return toBeSorted.stream().mapToInt(d -> d.originalIndex).toArray();
    }

    // small private class to assist in reading off the new ordering of the annotation array
    private static final class MyDoubleForSorting implements Comparable<MyDoubleForSorting> {
        final double myData;
        final int originalIndex;

        public MyDoubleForSorting(final double myData, final int originalIndex) {
            this.myData = myData;
            this.originalIndex = originalIndex;
        }

        @Override
        public int compareTo(final MyDoubleForSorting other) {
            return Double.compare(myData, other.myData);
        }
    }

    /**
     * Convert a normalized point to its original annotation value.
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
        final double sigma = stdDevVector[annI];
        return normalizedValue * sigma + mu;
    }

    /**
     * Returns the list of data items that are marked as "training" and are not marked as outliers.
     * If the list is too long ( MAX_NUM_TRAINING_DATA ) then it gets trimmed.
     */
    public List<VariantDatum> getTrainingData() {
        final List<VariantDatum> trainingData = data.stream()
                .filter(d -> d.atTrainingSite)
                .filter(d -> !d.failingSTDThreshold).collect(toList());

        logger.info( "Training with " + trainingData.size() + " variants after standard deviation thresholding." );
        if( trainingData.size() < VRAC.MIN_NUM_BAD_VARIANTS ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
        } else if( trainingData.size() > VRAC.MAX_NUM_TRAINING_DATA ) {
            logger.warn( "WARNING: Very large training set detected. Downsampling to " + VRAC.MAX_NUM_TRAINING_DATA + " training variants." );
            shuffle(trainingData, getRandomGenerator());
            return trainingData.subList(0, VRAC.MAX_NUM_TRAINING_DATA);
        }
        return trainingData;
    }

    /**
     * Returns the list of variants with LOD < BAD_LOD_CUTOFF, for building the negative model.
     * Returned items are those that have normalized values within (-STD_THRESHOLD , STD_THRESHOLD).
     * Returned data items have the atAntiTrainingSite set to true.
     */
    public List<VariantDatum> selectWorstVariants() {
        final List<VariantDatum> result = data.stream()
                .filter(vd -> !vd.failingSTDThreshold)
                .filter(vd -> !Double.isInfinite(vd.lod))
                .filter(vd -> vd.lod < VRAC.BAD_LOD_CUTOFF).collect(toList());

        result.forEach(vd -> vd.atAntiTrainingSite = true);

        logger.info( "Training with worst " + result.size() + " scoring variants --> variants with LOD <= " + String.format("%.4f", VRAC.BAD_LOD_CUTOFF) + "." );

        return result;
    }

    /**
     * Returns the list of all data items that:
     * are not outliers,
     * are not in the training set,
     * are not part of the 'worst variants'
     */
    public List<VariantDatum> getEvaluationData() {
        return data.stream()
                .filter(vd -> !vd.failingSTDThreshold)
                .filter(vd -> !vd.atTrainingSite)
                .filter(vd -> !vd.atAntiTrainingSite)
                .collect(toList());
    }

    /**
     * Remove all VariantDatum objects from the list which are marked as aggregate data
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

    /**
     * Computes the mean of the given annotation, either at training sites or at non-training sites.
     * Skips over annotations marked as 'null'.
     * Returns the mean or NaN when there are no data to compute the mean over.
     */
    @VisibleForTesting
    protected double mean(final int index, final boolean trainingData) {
        return data.stream()
                .filter(d -> trainingData == d.atTrainingSite)
                .filter(d -> !d.isNull[index])
                .mapToDouble(d -> d.annotations[index])
                .average().orElse(Double.NaN);
    }

    /**
     * Computes the std dev of the given annotation, either at training sites or at non-training sites.
     * Skips over annotations marked as 'null'.
     */
    @VisibleForTesting
    protected double standardDeviation( final double mean, final int index, final boolean trainingData ) {
        final double variance = data.stream()
                .filter(d -> trainingData == d.atTrainingSite)
                .filter(d -> !d.isNull[index])
                .mapToDouble(d -> (d.annotations[index] - mean))
                .map(d -> d*d)
                .average().orElse(Double.NaN);

        return sqrt(variance);
    }

    /**
     * Sets the annotations field for the VariantDatum based on the values of annotations from the VariantContext.
     */
    public void decodeAnnotations( final VariantDatum datum, final VariantContext vc) {
        final double[] annotations = new double[annotationKeys.length];
        final boolean[] isNull = new boolean[annotationKeys.length];
        int i = 0;
        for( final String key : annotationKeys ) {
            annotations[i] = decodeAnnotation( key, vc );
            if( Double.isNaN(annotations[i]) ) { isNull[i] = true; }
            i++;
        }
        datum.annotations = annotations;
        datum.isNull = isNull;
    }

    private static double decodeAnnotation( final String annotationKey, final VariantContext vc) {
        double value = vc.getAttributeAsDouble( annotationKey, Double.NaN );

        if( Double.isInfinite(value) ) {
            value = Double.NaN;
        }

        //XXX: some annotations have values that either have too little variance or have a hard threshold.
        //XXX: In those cases, we add a little bit of gaussian noise to make the data more suitable for the gaussian mixture modeling.
        //XXX: This is a HACK, to be clear. And the list of those annotations is hard-wired, which is also a HACK.
        if (needsJitter(annotationKey, value)){
            value += 0.01 * getRandomGenerator().nextGaussian();
        }

        return value;
    }

    private static boolean needsJitter(String annotationKey, double value) {
        switch(annotationKey.toLowerCase()){
            case "haplotypescore":  return equalDoubles(value, 0.0, 0.01);
            case "fs":              return equalDoubles(value, 0.0, 0.01);
            case "inbreedingcoeff": return equalDoubles(value, 0.0, 0.01);
            case "sor":             return equalDoubles(value, log(2.0), 0.01);  //min SOR is 2.0, then we take ln
            default: return false;
        }
    }

    /**
     * Print our the recalibration table from this manager to the provided writer.
     */
    public void writeOutRecalibrationTable( final VariantContextWriter recalWriter ) {
        // we need to sort in coordinate order in order to produce a valid VCF
        data.sort( (vd1, vd2) -> vd1.loc.compareTo(vd2.loc));

        // create dummy alleles to be used
        final List<Allele> alleles = Arrays.asList(Allele.create("N", true), Allele.create("<VQSR>", false));

        for( final VariantDatum datum : data ) {
            VariantContextBuilder builder = new VariantContextBuilder("VQSR", datum.loc.getContig(), datum.loc.getStart(), datum.loc.getStop(), alleles);
            builder.attribute(VCFConstants.END_KEY, datum.loc.getStop());
            builder.attribute(GATKVCFConstants.VQS_LOD_KEY, String.format("%.4f", datum.lod));
            builder.attribute(GATKVCFConstants.CULPRIT_KEY, (datum.worstAnnotation != -1 ? annotationKeys[datum.worstAnnotation] : "NULL"));

            if ( datum.atTrainingSite ) builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
            if ( datum.atAntiTrainingSite ) builder.attribute(GATKVCFConstants.NEGATIVE_LABEL_KEY, true);

            recalWriter.add(builder.make());
        }
    }
}
