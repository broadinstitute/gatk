package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;

@CommandLineProgramProperties(
	usage = "Prints reads from the input to the output.",
        usageShort = "Print reads",
        programGroup = ReadProgramGroup.class
)
public class PrintReads extends CommandLineProgram {

    @Option(fullName = "input", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The file to print reads from.")
    public File INPUT;

    @Option(fullName = "output", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    /**
     * Does the work of the tool, ie prints the reads from the inout to the output.
     * Returns null.
     */
    @Override
    public Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileHeader outputHeader = in.getFileHeader().clone();
        System.out.print(outputHeader);

        final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, REFERENCE_SEQUENCE);
        in.forEach(outputWriter::addAlignment);
        CloserUtil.close(in);
        CloserUtil.close(outputWriter);
        return null;
    }

}
