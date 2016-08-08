package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * A set of training variants for use with VQSR.
 */
public class TrainingSet {

    public FeatureInput<VariantContext> variantSource;
    public boolean isKnown = false;
    public boolean isTraining = false;
    public boolean isAntiTraining = false;
    public boolean isTruth = false;
    public boolean isConsensus = false;
    public double prior = 0.0;

    protected final Logger logger = LogManager.getLogger(TrainingSet.class);

    /**
     * A set of training variants for use with VQSR.
     * @param variantContextFeatureSource source of training variants. May not be  null.
     */
    public TrainingSet(final FeatureInput<VariantContext> variantContextFeatureSource) {
        Utils.nonNull(variantContextFeatureSource);
        this.variantSource = variantContextFeatureSource;

        final String name = variantContextFeatureSource.getName();

        // Parse the tags to decide which tracks have which properties
        isKnown = getBooleanAttribute(variantSource, "known");
        isTraining = getBooleanAttribute(variantSource, "training");
        isAntiTraining = getBooleanAttribute(variantSource, "bad");
        isTruth = getBooleanAttribute(variantSource, "truth");
        isConsensus = getBooleanAttribute(variantSource, "consensus");
        prior = getDoubleAttributeOrElse(variantSource, "prior", prior);

        // Report back to the user which tracks were found and the properties that were detected
        if( !isConsensus && !isAntiTraining ) {
            logger.info( String.format( "Found %s track: \tKnown = %s \tTraining = %s \tTruth = %s \tPrior = Q%.1f", name, isKnown, isTraining, isTruth, prior) );
        } else if( isConsensus ) {
            logger.info( String.format( "Found consensus track: %s", name) );
        } else {
            logger.info( String.format( "Found bad sites training track: %s", name) );
        }
    }

    private boolean getBooleanAttribute(final FeatureInput<VariantContext> variantSource, final String key) {
        String attributeValue = variantSource.getAttribute(key);
        return attributeValue != null && attributeValue.equals("true");
    }

    private double getDoubleAttributeOrElse(
            final FeatureInput<VariantContext> variantSource, final String key, final double defaultValue) {
        String attributeValue = variantSource.getAttribute(key);
        try {
            return attributeValue == null ?
                    defaultValue :
                    Double.valueOf(attributeValue);
        }
        catch (NumberFormatException e) {
            throw new UserException.MalformedFile("Malformed floating point value" + key);
        }
    }

}
