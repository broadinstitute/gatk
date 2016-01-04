package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

/**
 * Calls an R script to create plots of target data
 */
public final class CopyRatioSegmentedPlotter {
    private static final String R_SCRIPT = "TangentResultsPlotting.R";

    private CopyRatioSegmentedPlotter() {
    }
    /**
     * @param sample_name Name of the sample being run through the segmenter
     * @param tnFile Tangent-normalized targets file
     * @param preTnFile Targets file before tangent normalization
     * @param segFile Segmented tangent file
     * @param outputDir Full path to the outputted segment file
     * @param log Input tangent file has had a log2 transform applied
     */
    public static void writeSegmentedCopyRatioPlot(final String sample_name, final String tnFile, final String preTnFile, final String segFile, final String outputDir, final Boolean log, final Boolean sexChrs) {
        String logArg = log ? "TRUE" : "FALSE";
        String schr = sexChrs ? "TRUE" : "FALSE";
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, CopyRatioSegmentedPlotter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name="+sample_name, "--targets_file="+tnFile, "--pre_tn_file="+preTnFile, "--seg_file="+segFile,
                "--output_dir="+outputDir, "--log2_input="+logArg, "--sex_chrs="+schr);
        executor.exec();
    }
}