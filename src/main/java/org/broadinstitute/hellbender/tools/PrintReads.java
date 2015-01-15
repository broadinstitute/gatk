package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.io.File;

@CommandLineProgramProperties(
	usage = "Prints reads from the input to the output.",
    usageShort = "Print reads",
    programGroup = ReadProgramGroup.class
)
public class PrintReads extends ReadWalker {

    @Option(fullName = "output", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private SAMFileWriter outputWriter;

    @Override
    public void onTraversalStart() {
        final SAMFileHeader outputHeader = getHeaderForReads().clone();
        outputWriter = new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, REFERENCE_FILE);
    }

    @Override
    public void apply( SAMRecord read, ReferenceContext referenceContext ) {
        outputWriter.addAlignment(read);
    }

    @Override
    public Object onTraversalDone() {
        CloserUtil.close(outputWriter);
        return null;
    }
}
