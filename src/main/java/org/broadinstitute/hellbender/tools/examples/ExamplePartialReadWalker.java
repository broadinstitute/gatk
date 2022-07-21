package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;

/**
 * Example/toy program that prints reads from the provided file or files with corresponding reference bases
 * (if a reference is provided). Intended to show how to implement the ReadWalker interface.
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) until a read by a given name is encountered (which is not printed)",
        oneLineSummary = "Print reads until a specific read is encountered",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExamplePartialReadWalker extends PartialReadWalker {

    public static final String FULL_NAME_STOP_ON_READ_NAME = "stop-on-read-name";


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private GATKPath OUTPUT_FILE = null;

    @Argument(fullName = FULL_NAME_STOP_ON_READ_NAME, doc = "stop when a read by this name is encountered")
    public String stopOnReadName;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE.getOutputStream()) : System.out;
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputStream.printf("Read %s at %s:%d-%d:\n%s\n", read.getName(), read.getContig(), read.getStart(), read.getEnd(), read.getBasesString());
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }

    @Override
    protected boolean shouldExitEarly(GATKRead read) {
        return read.getName().equals(stopOnReadName);
    }
}
