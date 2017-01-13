package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.List;

/**
 * This class stores the results of forward-backward and Viterbi algorithms for HMM-based
 * copy ratio (or copy number) models.
 *
 * @param <D> copy ratio emission data type
 * @param <S> type of hidden state
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CopyRatioHiddenMarkovModelResults<D, S> implements Serializable {

    private static final long serialVersionUID = 1891158919985229044L;

    private final ForwardBackwardAlgorithm.Result<D, Target, S> fbResult;
    private final List<S> viterbiResult;

    /**
     * Public constructor.
     *
     * @param fbResult the result of forward-backward algorithm
     * @param viterbiResult the result of Viterbi algorithm
     */
    public CopyRatioHiddenMarkovModelResults(@Nonnull final ForwardBackwardAlgorithm.Result<D, Target, S> fbResult,
                                             @Nonnull final List<S> viterbiResult) {
        this.fbResult = Utils.nonNull(fbResult, "The forward-backward result must be non-null");
        this.viterbiResult = Utils.nonNull(viterbiResult, "The viterbi result must be non-null");
        Utils.validateArg(viterbiResult.size() == fbResult.positions().size(), "The target list used in the" +
                " forward-backward algorithm has a different length that the Viterbi chain of states. This is" +
                " an inconsistency.");
    }

    public ForwardBackwardAlgorithm.Result<D, Target, S> getForwardBackwardResult() {
        return fbResult;
    }

    public List<S> getViterbiResult() {
        return viterbiResult;
    }
}
