package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

/**
 * Calls an R script to create plots of target data
 */
public final class AlleleFractionSegmentedPlotter {
    private static final String R_SCRIPT = "AlleleFractionResultsPlotting.R";

    private AlleleFractionSegmentedPlotter() {}

    /**
     * @param sample_name Name of the sample being run through the segmenter
     * @param snpCountsFile File containing het SNP positions, ref counts, and alt counts
     * @param segmentsFile ACS modelled segment file
     * @param outputDir Full path to the outputted plot file
     * @param useSexChromosomes plot results for X and Y chromosomes
     */
    public static void writeSegmentedAlleleFractionPlot(final String sample_name, final String snpCountsFile,
            final String segmentsFile, final String outputDir, final Boolean useSexChromosomes) {
        final RScriptExecutor executor = new RScriptExecutor();
        String schr = useSexChromosomes ? "TRUE" : "FALSE";
        executor.addScript(new Resource(R_SCRIPT, AlleleFractionSegmentedPlotter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name="+sample_name, "--snp_counts_file="+snpCountsFile,  "--segments_file="+segmentsFile,
                "--output_dir="+outputDir, "--sex_chrs="+schr);
        executor.exec();
    }
}