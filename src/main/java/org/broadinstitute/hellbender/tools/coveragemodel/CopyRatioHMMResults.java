package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collections;
import java.util.List;

/**
 * This class stores the results of forward-backward and Viterbi algorithms for HMM-based
 * copy ratio (or copy number) models.
 *
 * @param <DATA> copy ratio emission data type
 * @param <STATE> type of hidden state
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CopyRatioHMMResults<DATA, STATE> implements Serializable {

    private static final long serialVersionUID = 1891158919985229044L;
    private final CopyRatioCallingMetadata metaData;
    private final ForwardBackwardAlgorithm.Result<DATA, Target, STATE> fbResult;
    private final List<STATE> viterbiResult;

    /**
     * Public constructor.
     *
     * @param fbResult the result of forward-backward algorithm
     * @param viterbiResult the result of Viterbi algorithm
     */
    public CopyRatioHMMResults(@Nonnull final CopyRatioCallingMetadata metaData,
                               @Nonnull final ForwardBackwardAlgorithm.Result<DATA, Target, STATE> fbResult,
                               @Nonnull final List<STATE> viterbiResult) {
        this.metaData = Utils.nonNull(metaData, "Copy ratio calling metadata must be non-null");
        this.fbResult = Utils.nonNull(fbResult, "The forward-backward result must be non-null");
        this.viterbiResult = Collections.unmodifiableList(Utils.nonNull(viterbiResult, "The viterbi result must be non-null"));
        Utils.validateArg(viterbiResult.size() == fbResult.positions().size(), "The target list used in the" +
                " forward-backward algorithm has a different length that the Viterbi chain of states. This is" +
                " an inconsistency.");
    }

    public ForwardBackwardAlgorithm.Result<DATA, Target, STATE> getForwardBackwardResult() {
        return fbResult;
    }

    public List<STATE> getViterbiResult() {
        return viterbiResult;
    }

    public CopyRatioCallingMetadata getMetaData() {
        return metaData;
    }
}
