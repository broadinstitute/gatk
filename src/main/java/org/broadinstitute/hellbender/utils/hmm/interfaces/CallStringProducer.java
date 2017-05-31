package org.broadinstitute.hellbender.utils.hmm.interfaces;

/**
 * An interface for hidden states that provide a call string (i.e. a unique string identifier)
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface CallStringProducer {

    String getCallString();
}
