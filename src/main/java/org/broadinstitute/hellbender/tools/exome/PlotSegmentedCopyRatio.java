package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.plotter.CopyRatioSegmentedPlotter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Create plots of copy number variant data.  Please note that this tool is only supported for hg19 and b37 references.  All other references may fail.",
        oneLineSummary = "Create plots of copy number variant data.",
        programGroup = CopyNumberProgramGroup.class
)
public final class PlotSegmentedCopyRatio extends CommandLineProgram {

    @Argument(
            doc = "Genomic targets file after tangent normalization has been applied, produced by NormalizeSomaticReadCounts: tn",
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            optional = false
    )
    protected String tangentFile;

    @Argument(
            doc = "Genomic targets before tangent normalization file, produced by NormalizeSomaticReadCounts: preTN",
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected String preTangentFile;

    @Argument(
            doc = "File of segmented regions of the genome",
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            optional = false
    )
    protected String segmentFile;

    @Argument(
            doc = "Directory to write plots",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected String plotDir;

    @Argument(
            doc = "If input data has had a log2 transform applied",
            shortName = ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.LOG2_LONG_NAME,
            optional = true
    )
    protected Boolean log = false;

    @Argument(
            doc = "Plot sex chromosomes",
            shortName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_LONG_NAME,
            optional = true
    )
    protected Boolean sexChrs = false;

    @Override
    protected Object doWork() {
        final String sampleName = TargetCoverageUtils.getSampleNameForCLIsFromTargetCoverageFile(new File(tangentFile));
        createPlot(sampleName, tangentFile, preTangentFile, segmentFile, plotDir, log, sexChrs);
        return "Success";
    }

    private void createPlot(String sampleName, String tangentFile, String preTangentFile, String segmentFile, String outFile, boolean log, boolean sexChrs){
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(sampleName, tangentFile, preTangentFile, segmentFile, outFile, log, sexChrs);
    }
}
