package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;
import java.util.List;

/**
 * Set of arguments for the {@link FlowFeatureMapper}
 */
public class ApplySNVQRArgumentCollection implements Serializable{
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

    @Hidden
    @Argument(fullName = "debug-read-folder", doc = "folder to write debugged reads", optional = true)
    public String debugReadFolder = "/tmp";

    @Hidden
    @Argument(fullName = "keep-negatives", doc = "keep nevative scores?", optional = true)
    public boolean keepNegatives = false;

    @Hidden
    @Argument(fullName = "debug-negatives", doc = "debug negative scores?", optional = true)
    public boolean debugNegatives = false;

    @Argument(fullName = "negative-score-override", doc = "value to override negative scores with", optional = true)
    public double negativeScoreOverride = 0;
    /**
     *  snv-qual model file  - if not specified, this feature is off
     **/
    @Argument(fullName = "model", doc = "XGBoost model file. in this mode ... <need to insert explanation here>", optional = true)
    public String model;

    /**
     *  svn-qual conf file
     **/
    @Argument(fullName = "conf", doc = "Snvconf (json) configuration file. <need to insert explanation here>", optional = true)
    public String conf;

    @Argument(fullName = "replace-quality-mode", doc = "<doc>", optional = true)
    public boolean replaceQualityMode = false;

    @Argument(fullName = "stats-path-prefix", doc = "<doc>", optional = true)
    public String statsPathPrefix = null;

}
