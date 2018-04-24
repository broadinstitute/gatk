package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Gathers scattered VQSLOD tranches into a single file. For use when running VariantRecalibrator on " +
                "scattered input using the -scatterTranches mode.",
        oneLineSummary = "Gathers scattered VQSLOD tranches into a single file",
        programGroup = OtherProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class GatherTranches extends CommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc="List of scattered tranches files")
    public final List<File> inputReports = new ArrayList<>();

    /**
     * Add truth sensitivity slices through the call set at the given values. The default values are 100.0, 99.9, 99.0, and 90.0
     * which will result in 4 estimated tranches in the final call set: the full set of calls (100% sensitivity at the accessible
     * sites in the truth set), a 99.9% truth sensitivity tranche, along with progressively smaller tranches at 99% and 90%.
     */
    @Argument(fullName="truth-sensitivity-tranche",
            shortName="tranche",
            doc="The levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)",
            optional=true)
    private List<Double> TS_TRANCHES = new ArrayList<>(Arrays.asList(100.0, 99.9, 99.0, 90.0));

    /**
     * Use either SNP for recalibrating only SNPs (emitting indels untouched in the output VCF) or INDEL for indels (emitting SNPs untouched in the output VCF). There is also a BOTH option for recalibrating both SNPs and indels simultaneously, but this is meant for testing purposes only and should not be used in actual analyses.
     */
    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ", optional = false)
    public VariantRecalibratorArgumentCollection.Mode MODE = VariantRecalibratorArgumentCollection.Mode.SNP;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="File to output the gathered tranches file to")
    public File outputReport;

    private PrintStream tranchesStream;

    @Override
    protected Object doWork() {
        inputReports.forEach(IOUtil::assertFileIsReadable);

        try {
            tranchesStream = new PrintStream(outputReport);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputReport, e);
        }

        //use a data structure to hold the tranches from each scatter shard in a format that's easy to merge
        final TreeMap<Double, List<VQSLODTranche>> scatteredTranches = new TreeMap<>();
        for (final File trancheFile : inputReports) {
            try {
                for (final VQSLODTranche currentTranche : VQSLODTranche.readTranches(trancheFile)) {
                    if (scatteredTranches.containsKey(currentTranche.minVQSLod)) {
                        scatteredTranches.get(currentTranche.minVQSLod).add(currentTranche);
                    }
                    else {
                        scatteredTranches.put(currentTranche.minVQSLod, new ArrayList<>(Arrays.asList(currentTranche)));
                    }
                }
            } catch (IOException e) {
                throw new UserException.CouldNotReadInputFile(trancheFile, e);
            }
        }

        tranchesStream.print(TruthSensitivityTranche.printHeader());
        tranchesStream.print(Tranche.tranchesString(VQSLODTranche.mergeAndConvertTranches(scatteredTranches, TS_TRANCHES, MODE)));

        return 0;
    }
}
