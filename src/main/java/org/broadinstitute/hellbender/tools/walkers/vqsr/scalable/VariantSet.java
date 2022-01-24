package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

final class VariantSet {

    final FeatureInput<VariantContext> variantSource;
    final boolean isTraining;
    final boolean isTruth;

    private static final Logger logger = LogManager.getLogger(VariantSet.class);

    /**
     * @param variantContextFeatureSource source of variants. May not be null.
     */
    VariantSet(final FeatureInput<VariantContext> variantContextFeatureSource) {
        Utils.nonNull(variantContextFeatureSource);
        this.variantSource = variantContextFeatureSource;

        final String name = variantContextFeatureSource.getName();

        // Parse the tags to decide which tracks have which properties
        isTraining = getBooleanAttribute(variantSource, "training");
        isTruth = getBooleanAttribute(variantSource, "truth");

        // Report back to the user which tracks were found and the properties that were detected
        logger.info( String.format( "Found %s track: \tTraining = %s \tTruth = %s", name, isTraining, isTruth) );
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
