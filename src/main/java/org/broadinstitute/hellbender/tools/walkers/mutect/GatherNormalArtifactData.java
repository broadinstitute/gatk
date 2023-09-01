package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary="Combine output files from GetNormalArtifactData in the order defined by a sequence dictionary",
        oneLineSummary = "Combine output files from GetNormalArtifactData in the order defined by a sequence dictionary",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class GatherNormalArtifactData extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "an output of GetNormalArtifactData")
    final List<File> input = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "output")
    final File output = null;

    @Override
    protected Object doWork() {

        try ( NormalArtifactRecord.NormalArtifactWriter writer = new NormalArtifactRecord.NormalArtifactWriter(IOUtils.fileToPath(output)) ) {
            for (final File inputFile : input) {
                writer.writeAllRecords(NormalArtifactRecord.readFromFile(inputFile));
            }
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", output));
        }

        return "SUCCESS";
    }
}
