package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.utils.plotter.CopyRatioSegmentedPlotter;

@CommandLineProgramProperties(
        summary = "Create plots of copy number variant data.  Please note that this tool is only supported for hg19 and b37 references.  All other references may fail.",
        oneLineSummary = "Create plots of copy number variant data.",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class PlotSegmentedCopyRatio extends CommandLineProgram {

    public static final String SAMPLE_NAME_LONG_NAME = "sampleName";
    public static final String SAMPLE_NAME_SHORT_NAME = "S";

    public static final String TARGETS_FILE_LONG_NAME = "targets";
    public static final String TARGETS_FILE_SHORT_NAME = "T";

    public static final String PRE_TANGENT_TARGETS_FILE_LONG_NAME = "preTangent";
    public static final String PRE_TANGENT_TARGETS_FILE_SHORT_NAME = "P";

    public static final String SEGMENT_FILE_LONG_NAME = "segments";
    public static final String SEGMENT_FILE_SHORT_NAME = "seg";

    public static final String LOG2_LONG_NAME= "log2Input";
    public static final String LOG2_SHORT_NAME = "log";

    @Argument(
            doc = "Name of the sample we are plotting",
            shortName = SAMPLE_NAME_SHORT_NAME,
            fullName = SAMPLE_NAME_LONG_NAME,
            optional = false
    )
    protected String sampleName;

    @Argument(
            doc = "Genomic targets file after tangent normalization has been applied, produced by NormalizeSomaticReadCounts: tn",
            shortName = TARGETS_FILE_SHORT_NAME,
            fullName =  TARGETS_FILE_LONG_NAME,
            optional = false
    )
    protected String tangentFile;

    @Argument(
            doc = "Genomic targets before tangent normalization file, produced by NormalizeSomaticReadCounts: preTN",
            shortName = PRE_TANGENT_TARGETS_FILE_SHORT_NAME,
            fullName =  PRE_TANGENT_TARGETS_FILE_LONG_NAME,
            optional = false
    )
    protected String preTangentFile;

    @Argument(
            doc = "File of segmented regions of the genome based on copy-ratio produced by CBS",
            shortName = SEGMENT_FILE_SHORT_NAME,
            fullName =  SEGMENT_FILE_LONG_NAME,
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
            shortName = LOG2_SHORT_NAME,
            fullName = LOG2_LONG_NAME,
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
