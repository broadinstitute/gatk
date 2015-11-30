package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.utils.plotter.CopyRatioSegmentedPlotter;

@CommandLineProgramProperties(
        summary = "Create plots of copy number variant data.  Please note that this tool is only supported for hg19 and b37 references.  All other references may fail.",
        oneLineSummary = "Create plots of copy number variant data.",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class PlotSegmentedCopyRatio extends CommandLineProgram {

    @Argument(
            doc = "Name of the sample we are plotting",
            fullName = ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME,
            optional = false
    )
    protected String sampleName;

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

    @Override
    protected Object doWork() {
        createPlot(sampleName, tangentFile, preTangentFile, segmentFile, plotDir, log);
        return "Success";
    }

    private void createPlot(String sampleName, String tangentFile, String preTangentFile, String segmentFile, String outFile, boolean log){
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(sampleName, tangentFile, preTangentFile, segmentFile, outFile, log);
    }
}
