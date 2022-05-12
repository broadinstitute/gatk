package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.CramHeader;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.regex.Pattern;

/***
 * SplitCRAM - split a cram file into smaller cram files (shards) containing a minimal number of records
 * while still respecting container boundaries.
 *
 * The tool operates on a CRAM container level and therefore is efficient but not exact in the number of
 * records on each output file (container boundaries are maintained)
 *
 * Note that CRAM files have relative record counters embedded in each container. These are not reset by
 * this tool. Therefore, the resulting files may not contain correct record counter values.
 *
 * <h3>Usage</h3>
 *
 * ./gatk SplitCRAM \
 *  -I
 *  input.cram
 *  -O
 *  output_%04d.cram
 *  --shard-records
 *  5000000
 *
 *
 * Notes:
 * 1. shard-records is optional. defaults to 10M
 * 2. output filename should contain a %d formatter pattern
 */

@CommandLineProgramProperties(
        summary = "Splits CRAM files efficiently by taking advantage of their container based structure",
        oneLineSummary = "Split CRAM files to smaller files efficiently",
        programGroup = FlowBasedProgramGroup.class
)
@WorkflowProperties
@ExperimentalFeature
public class SplitCRAM extends CommandLineProgram {

    public static final int DEFAULT_SHARD_RECORDS = 10000000;
    public static final String SHARD_RECORDS_FULL_NAME = "shard-records";
    public static final Pattern numeratorFormat = Pattern.compile("%[0-9]*d");

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "input cram file to split")
    private GATKPath cramInput = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "output cram file template. should contain %d, which will be replaced by shard index", optional = true)
    private String cramOutputTemplate = "output_%04d.cram";

    @Argument(fullName = SHARD_RECORDS_FULL_NAME, doc = "minimum threshold for number of records per shard.", optional = true)
    private long shardRecords = DEFAULT_SHARD_RECORDS;

    // locals
    CramContainerIterator cramContainerIterator;
    int shard;

    @Override
    protected void onStartup() {
        super.onStartup();

        // check that output template contains a %d formatter
        if ( !numeratorFormat.matcher(cramOutputTemplate).find() ) {
            throw new IllegalArgumentException("output template missing a %d enumerator formatter: " + cramOutputTemplate);
        }
    }

    @Override
    protected Object doWork() {

        try (final CramContainerIterator cramContainerIterator = new CramContainerIterator(new BufferedInputStream(cramInput.getInputStream())) ){

            // get header
            final CramHeader cramHeader = cramContainerIterator.getCramHeader();

            // iterate
            while (cramContainerIterator.hasNext()) {

                try (final OutputStream os = nextOutputStream()) {

                    // write headers
                    CramIO.writeCramHeader(cramContainerIterator.getCramHeader(), os);
                    Container.writeSAMFileHeaderContainer(cramContainerIterator.getCramHeader().getCRAMVersion(), cramContainerIterator.getSamFileHeader(), os);

                    // iterate
                    long records = 0;
                    while (cramContainerIterator.hasNext() && (records < shardRecords)) {

                        // get next container
                        final Container container = cramContainerIterator.next();

                        // write container to output stream
                        container.write(cramHeader.getCRAMVersion(), os);

                        // update record count
                        records += container.getContainerHeader().getNumberOfRecords();
                    }

                    CramIO.writeCramEOF(cramContainerIterator.getCramHeader().getCRAMVersion(), os);
                }
            }
        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return null;
    }

    private OutputStream nextOutputStream() {

        final String filename = String.format(cramOutputTemplate, shard++);
        final GATKPath path = new GATKPath(filename);

        return new BufferedOutputStream(path.getOutputStream());
    }
}
