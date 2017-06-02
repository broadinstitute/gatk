package org.broadinstitute.hellbender.utils.hmm.interfaces;

/**
 * An interface for hidden states that can provide a scalar value.
 *
 * For example, {@link org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState}
 * naturally provides the copy ratio.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface ScalarProducer {

    /**
     * Convert the hidden state to a double type scalar
     * @return a double value
     */
    double getScalar();

}
