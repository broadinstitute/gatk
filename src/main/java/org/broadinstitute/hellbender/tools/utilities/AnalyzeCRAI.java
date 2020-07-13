package org.broadinstitute.hellbender.tools.utilities;

import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.cram.CRAIIndex;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;

/**
 * Analyzer for CRAM (.crai) index files.
 */
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Print information about a .crai index",
        oneLineSummary = "Print information about a .crai index",
        programGroup = OtherProgramGroup.class
)
public class AnalyzeCRAI extends CommandLineProgram {

    @Argument(
            shortName=StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName=StandardArgumentDefinitions.INPUT_LONG_NAME,
            doc="File to be analyzed",
            optional=false)
    private GATKPath inputPath;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "File to which to write file analysis information (if not specified, prints to standard output)",
            optional = true)
    private GATKPath outputPath;

    /**
     * Run the analyzer for the file.
     * @return
     */
    protected Object doWork() {
        try (final InputStream fis = new BufferedInputStream(inputPath.getInputStream());
             final PrintStream outputStream = new PrintStream(new BufferedOutputStream(
                    outputPath == null ?
                            System.out :
                            outputPath.getOutputStream()))) {

                final CRAIIndex craiIndex = CRAMCRAIIndexer.readIndex(fis);

            System.out.println("\nSeqId AlignmentStart AlignmentSpan ContainerOffset SliceOffset SliceSize\n");
            craiIndex.getCRAIEntries().stream().forEach(outputStream::println);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return null;
    }

}

