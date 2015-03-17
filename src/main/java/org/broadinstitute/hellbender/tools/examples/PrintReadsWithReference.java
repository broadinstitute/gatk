package org.broadinstitute.hellbender.tools.examples;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Example/toy program that prints reads from the provided file or files with corresponding reference bases
 * (if a reference is provided). Intended to show how to implement the ReadWalker interface.
 */
@CommandLineProgramProperties(
        usage = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        usageShort = "Print reads with reference context",
        programGroup = ReadProgramGroup.class
)
public class PrintReadsWithReference extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    public File OUTPUT_FILE = null;

    private PrintStream outputStream = null;

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public void onTraversalStart() {
        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }
    }

    @Override
    public void apply( SAMRecord read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputStream.printf("Read at %s:%d-%d:\n%s\n", read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getReadString());
        if ( referenceContext.hasBackingDataSource() )
            outputStream.println("Reference Context:\n" + new String(referenceContext.getBases()));
        outputStream.println();
    }

    @Override
    public Object onTraversalDone() {
        if ( outputStream != null )
            outputStream.close();

        return null;
    }
}
