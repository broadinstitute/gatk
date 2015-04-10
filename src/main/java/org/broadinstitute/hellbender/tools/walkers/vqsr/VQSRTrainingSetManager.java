package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;

/**
 * This class is used to manage access to training sets in VQSR.
 */
public final class VQSRTrainingSetManager {

    private final List<DataSet> dataSets;

    public VQSRTrainingSetManager() {
        this.dataSets = new ArrayList<>();
    }

    /*
     * Helper class to represent a data set used in VQSR.
     */
    private final static class DataSet {
        private final FeatureInput<VariantContext> fi;
        private final boolean isKnown;
        private final boolean isTraining;
        private final boolean isAntiTraining;
        private final boolean isTruth;
        private final boolean isConsensus;
        private final double prior;

        private final static Logger logger = LogManager.getLogger(DataSet.class);

        public DataSet(final FeatureInput<VariantContext> fi) {
            this.fi = fi;

            final String name = fi.getName();

            // Parse the tags to decide which tracks have which properties
            isKnown = Objects.equals(fi.getAttribute("known"), "true");
            isTraining = Objects.equals(fi.getAttribute("training"), "true");
            isAntiTraining = Objects.equals(fi.getAttribute("bad"), "true");
            isTruth = Objects.equals(fi.getAttribute("truth"), "true");
            isConsensus = Objects.equals(fi.getAttribute("consensus"), "true");
            prior = ( fi.getAttribute("prior") != null ? Double.parseDouble(fi.getAttribute("prior")) : 0.0 );

            // Report back to the user which tracks were found and the properties that were detected
            if( !isConsensus && !isAntiTraining ) {
                logger.info( String.format("Found %s track: \tKnown = %s \tTraining = %s \tTruth = %s \tPrior = Q%.1f", name, isKnown, isTraining, isTruth, prior) );
            } else if( isConsensus ) {
                logger.info( String.format("Found consensus track: %s", name) );
            } else {
                logger.info( String.format("Found bad sites training track: %s", name) );
            }
        }

        /**
         * Annotates the VariantDatum with information from this data set.
         */
        public void annotateDatum(FeatureContext tracker, VariantContext evalVC, VariantDatum datum, boolean trust_all_polymorphic) {
            for( final VariantContext trainVC : tracker.getValues(this.getFeatureInput()) ) {
                if( isValidVariant( evalVC, trainVC, trust_all_polymorphic ) ) {
                    datum.isKnown = datum.isKnown || this.isKnown;
                    datum.atTruthSite = datum.atTruthSite || this.isTruth;
                    datum.atTrainingSite = datum.atTrainingSite || this.isTraining;
                    datum.prior = Math.max( datum.prior, this.prior );
                    datum.consensusCount += ( this.isConsensus ? 1 : 0 );
                }
                if( trainVC != null ) {
                    datum.atAntiTrainingSite = datum.atAntiTrainingSite || this.isAntiTraining;
                }
            }
        }

        private static boolean isValidVariant( final VariantContext evalVC, final VariantContext trainVC, final boolean trustAllPolymorphic) {
            return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() && checkVariationClass( evalVC, trainVC ) &&
                    (trustAllPolymorphic || !trainVC.hasGenotypes() || trainVC.isPolymorphicInSamples());
        }

        private static boolean checkVariationClass( final VariantContext evalVC, final VariantContext trainVC ) {
            switch( trainVC.getType() ) {
                case SNP:
                case MNP:
                    return evalVC.isSNP() || evalVC.isMNP();
                case INDEL:
                case MIXED:
                case SYMBOLIC:
                    return evalVC.isStructuralIndel() || evalVC.isIndel() || evalVC.isMixed() || evalVC.isSymbolic();
                default:
                    return false;
            }
        }

        public FeatureInput<VariantContext> getFeatureInput(){
            return fi;
        }
    }

    /**
     * Adds data sets to this manager.
     */
    public void addDataSets(final Collection<FeatureInput<VariantContext>> dataSetsArgs) {
        dataSetsArgs.forEach(this::addDataSet);
    }

    /**
     * Adds a data set to this manager.
     */
    public void addDataSet(final FeatureInput<VariantContext> dataset) {
        dataSets.add(new DataSet(dataset));
    }

    /**
     * Returns whether the manager has a data set that is marked as 'training'.
     */
    public boolean hasTrainingSet() {
        return dataSets.stream().anyMatch(ts -> ts.isTraining);
    }

    /**
     * Returns whether the manager has a data set that is marked as 'truth'.
     */
    public boolean hasTruthSet() {
        return dataSets.stream().anyMatch(ts -> ts.isTruth);
    }

    public void annotateDatum(final FeatureContext tracker, final VariantContext evalVC, final VariantDatum datum, final boolean trust_all_polymorphic) {
        resetDatum(datum);
        dataSets.forEach(ts -> ts.annotateDatum(tracker, evalVC, datum, trust_all_polymorphic));
    }

    private void resetDatum(VariantDatum datum) {
        datum.isKnown = false;
        datum.atTruthSite = false;
        datum.atTrainingSite = false;
        datum.atAntiTrainingSite = false;
        datum.prior = 2.0;
    }
}
