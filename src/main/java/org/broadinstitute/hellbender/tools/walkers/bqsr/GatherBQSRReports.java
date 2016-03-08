package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Gathers scattered BQSR recalibration reports into a single file",
        oneLineSummary = "Gathers scattered BQSR recalibration reports into a single file",
        programGroup = ReadProgramGroup.class
)
public final class GatherBQSRReports extends CommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc="List of scattered BQSR report files")
    public final List<File> inputReports = new ArrayList<>();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="File to output the gathered file to")
    public File outputReport;


    @Override
    protected Object doWork() {
        inputReports.forEach(IOUtil::assertFileIsReadable);
        IOUtil.assertFileIsWritable(outputReport);

        RecalibrationReport.gatherReportsIntoOneFile(inputReports, outputReport);

        return 0;
    }
}
