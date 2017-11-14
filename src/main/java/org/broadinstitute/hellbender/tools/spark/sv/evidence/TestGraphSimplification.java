package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FermiLiteAssemblyHandler;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) total junk", summary = "a complete mess",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class TestGraphSimplification extends CommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "input fastq",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME)
    private String fastqFile;

    @Override
    protected Object doWork() {
        // read the reads
        final List<SVFastqUtils.FastqRead> reads = SVFastqUtils.readFastqFile(fastqFile);

        // assemble them
        final FermiLiteAssembler assembler = new FermiLiteAssembler();
        assembler.setCleaningFlag(0x60);
        final FermiLiteAssembly assembly1 = assembler.createAssembly(reads);
        final FermiLiteAssembly assembly2 = FermiLiteAssemblyHandler.removeShadowedContigs(assembly1);
        final FermiLiteAssembly assembly3 = FermiLiteAssemblyHandler.removeUnbranchedConnections(assembly2);
        final FermiLiteAssembly assembly4 = FermiLiteAssemblyHandler.expandAssemblyGraph(assembly3);
        try {
            assembly4.writeGFA(new FileOutputStream("/dev/stdout"));
        } catch ( final IOException ioe ) {
            throw new GATKException("can't write output", ioe);
        }
        return null;
    }
}
