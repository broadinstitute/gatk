package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;
import java.util.List;

/**
 * Set of arguments for the {@link FlowFeatureMapper}
 */
public class AddFlowSNVQualityArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    public enum SnvqModeEnum {
        Legacy,
        Optimistic,
        Pessimistic,
        Geometric
    };

    /**
     *  maximum value for
     *  delta in score
     **/
    @Argument(fullName = "limit-phred-score", doc = "Limit value for phred scores", optional = true)
    public double limitPhredScore = Double.NaN;

    /**
     *  include duplicate read?
     **/
    @Argument(fullName = "include-dup-reads", doc = "include duplicate reads", optional = true)
    public boolean includeDupReads = false;

    /**
     *  keep supplementary alignments?
     **/
    @Argument(fullName = "keep-supplementary-alignments", doc = "keep supplementary alignments ?", optional = true)
    public boolean keepSupplementaryAlignments = false;

    /**
     * snvq computation mode
     */
    @Argument(fullName = "snvq-mode", doc = "ksnvq computation mode", optional = true)
    public SnvqModeEnum snvMode = SnvqModeEnum.Geometric;

    /**
     *  debug read names?
     **/
    @Hidden
    @Argument(fullName = "debug-read-name", doc = "debug specific reads?", optional = true)
    public List<String> debugReadName = null;
}
