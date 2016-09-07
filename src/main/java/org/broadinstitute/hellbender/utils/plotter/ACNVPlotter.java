package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;

/**
 * Calls an R script to create plots of target data
 */
public final class ACNVPlotter {
    private static final String CNV_PLOTTING_R_LIBRARY = "CNV_plotting_library.R";
    private static final String ACNV_PLOTTING_R_SCRIPT = "ACNVResultsPlotting.R";
    private static final String ACNV_PLOTTING_PER_SEG_R_SCRIPT = "ACNVResultsPlottingPerSeg.R";

    private ACNVPlotter() {}

    /**
     * @param sampleName Name of the sample being run through the segmenter
     * @param snpCountsFile File containing het SNP positions, ref counts, and alt counts
     * @param coverageFile File containing tangent normalized coverage of targets
     * @param segmentsFile ACS modelled segment file
     * @param outputDir Full path to the output directory
     * @param outputPrefix Prefix to use for output plot files
     * @param useSexChromosomes plot results for X and Y chromosomes
     */
    public static void writeSegmentedAlleleFractionPlot(final String sampleName, final File snpCountsFile, final File coverageFile,
                                                        final File segmentsFile, final String outputDir, final String outputPrefix, final Boolean useSexChromosomes) {
        final RScriptExecutor executor = new RScriptExecutor();
        final String schr = useSexChromosomes ? "TRUE" : "FALSE";
        //This leads to the R statement source("CNV_plotting_library.R") before the main script runs
        executor.addScript(new Resource(CNV_PLOTTING_R_LIBRARY, ACNVPlotter.class));
        executor.addScript(new Resource(ACNV_PLOTTING_R_SCRIPT, ACNVPlotter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name=" + sampleName, "--snp_counts_file=" + snpCountsFile,  "--coverage_file=" + coverageFile,
                "--segments_file=" + segmentsFile, "--output_dir=" + outputDir, "--output_prefix=" + outputPrefix, "--sex_chrs=" + schr);
        executor.exec();
    }

    /**
     * @param sampleName Name of the sample being run through the segmenter
     * @param snpCountsFile File containing het SNP positions, ref counts, and alt counts
     * @param coverageFile File containing tangent normalized coverage of targets
     * @param segmentsFile ACS modelled segment file
     * @param outputDir Full path to the output directory
     * @param outputPrefix Prefix to use for output plot files
     * @param useSexChromosomes plot results for X and Y chromosomes
     */
    public static void writeSegmentedAlleleFractionPlotPerSeg(final String sampleName, final File snpCountsFile, final File coverageFile,
                                                              final File segmentsFile, final String outputDir, final String outputPrefix, final Boolean useSexChromosomes) {
        final RScriptExecutor executor = new RScriptExecutor();
        final String schr = useSexChromosomes ? "TRUE" : "FALSE";
        //This leads to the R statement source("CNV_plotting_library.R") before the main script runs
        executor.addScript(new Resource(CNV_PLOTTING_R_LIBRARY, ACNVPlotter.class));
        executor.addScript(new Resource(ACNV_PLOTTING_PER_SEG_R_SCRIPT, ACNVPlotter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name=" + sampleName, "--snp_counts_file=" + snpCountsFile,  "--coverage_file=" + coverageFile,
                "--segments_file=" + segmentsFile, "--output_dir=" + outputDir, "--output_prefix=" + outputPrefix, "--sex_chrs=" + schr);
        executor.exec();
    }
}