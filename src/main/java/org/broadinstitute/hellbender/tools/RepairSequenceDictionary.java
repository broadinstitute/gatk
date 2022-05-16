package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

@DocumentedFeature
@CommandLineProgramProperties(
        summary="Repair a sequence dictionary",
        oneLineSummary = "Repair a sequence dictionary",
        programGroup = DiagnosticsAndQCProgramGroup.class)
public class RepairSequenceDictionary extends CommandLineProgram {
    @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME)
    public final GATKPath readsInput = null;

    public final ReferenceInputArgumentCollection referenceArgumentCollection =
            new OptionalReferenceInputArgumentCollection();

    @Override
    protected Object doWork() {
//        try (final SAMReader samReader = new S)
//            final SAMSequenceDictionary sequenceDictionary;
//            checkCRAM();
        return null;
    }

    protected String[] customCommandLineValidation() {
        if (!readsInput.isCram()) {
            throw new CommandLineException(String.format(
                    "The provided in put must be a CRAM file:", readsInput));
        }
        return null;
    }
}
