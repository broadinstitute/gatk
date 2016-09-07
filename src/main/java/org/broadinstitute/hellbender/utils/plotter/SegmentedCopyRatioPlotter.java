package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;

/**
 * Calls an R script to create plots of target data
 */
public final class SegmentedCopyRatioPlotter {
    private static final String CNV_PLOTTING_R_LIBRARY = "CNV_plotting_library.R";
    private static final String COPY_RATIO_PLOTTING_R_SCRIPT = "TangentResultsPlotting.R";

    private SegmentedCopyRatioPlotter() {
    }

    /**
     * @param sampleName Name of the sample being run through the segmenter
     * @param tnFile Tangent-normalized targets file
     * @param preTnFile Targets file before tangent normalization
     * @param segmentsFile Segmented tangent file
     * @param outputDir Full path to the output plot directory
     * @param outputPrefix Prefix to use for output plot files
     * @param log Input tangent file has had a log2 transform applied
     * @param useSexChromosomes plot results for X and Y chromosomes
     */
    public static void writeSegmentedCopyRatioPlot(final String sampleName, final File tnFile, final File preTnFile,
                                                   final File segmentsFile, final String outputDir, final String outputPrefix,
                                                   final Boolean log, final Boolean useSexChromosomes) {
        final String logArg = log ? "TRUE" : "FALSE";
        final String schr = useSexChromosomes ? "TRUE" : "FALSE";
        final RScriptExecutor executor = new RScriptExecutor();
        //This leads to the R statement source("CNV_plotting_library.R") before the main script runs
        executor.addScript(new Resource(CNV_PLOTTING_R_LIBRARY, SegmentedCopyRatioPlotter.class));
        executor.addScript(new Resource(COPY_RATIO_PLOTTING_R_SCRIPT, SegmentedCopyRatioPlotter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name=" + sampleName, "--targets_file=" + tnFile, "--pre_tn_file=" + preTnFile,
                "--segments_file=" + segmentsFile, "--output_dir=" + outputDir, "--output_prefix=" + outputPrefix,
                "--log2_input=" + logArg, "--sex_chrs=" + schr);
        executor.exec();
    }
}