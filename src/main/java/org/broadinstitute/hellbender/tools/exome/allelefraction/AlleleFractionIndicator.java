package org.broadinstitute.hellbender.tools.exome.allelefraction;

/**
 * Possible hidden states of hets in allele fraction model.
 *
 * NOTE: these are not part of AlleleFractionState because AlleleFractionState pertains to a collapsed model in which latent variables
 * (this indicator and allelic biases) have been marginalized out.  These hidden states are, however, necessary
 * for initializing the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public enum AlleleFractionIndicator {
    ALT_MINOR, REF_MINOR, OUTLIER
}
