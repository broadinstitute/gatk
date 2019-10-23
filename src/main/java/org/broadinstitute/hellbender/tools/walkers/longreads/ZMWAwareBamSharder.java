package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

@CommandLineProgramProperties(
        summary = "TODO",
        oneLineSummary = "TODO",
        programGroup = ReadDataManipulationProgramGroup.class
)
public class ZMWAwareBamSharder extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String outputBAM;

    @Argument(fullName = "output-index",
            doc="Write output index to this file")
    public String outputIndexFile;

    @Argument(fullName = "target-shard-size",
              doc = "each shard should have approximately this many reads")
    public int targetShardSize;

    private SAMFileGATKReadWriter outputWriter;

    private PrintWriter indexWriter;

    private String currentZMW;
    private String previousZMW;
    private int numReadsInCurrentShard = 0;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(outputBAM), true);
        try {
            indexWriter = new PrintWriter(outputIndexFile);
        } catch ( IOException e ) {
            throw new UserException("Oops", e);
        }
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        previousZMW = currentZMW;
        currentZMW = getZMWForRead(read.getName());
        outputWriter.addRead(read);
        ++numReadsInCurrentShard;

        if ( previousZMW == null || ( ! currentZMW.equals(previousZMW) && numReadsInCurrentShard >= targetShardSize) ) {
            writeIndexEntryForRead(read);
        }
    }

    @Override
    public void closeTool() {
        if ( indexWriter != null ) {
            indexWriter.close();
        }

        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }

    private void writeIndexEntryForRead(final GATKRead read) {
        final SAMRecord unwrapped = read.convertToSAMRecord(getHeaderForReads());
        final BAMFileSpan fileSource = (BAMFileSpan)unwrapped.getFileSource().getFilePointer();

        final List<Chunk> chunks = fileSource.getChunks();
        if ( chunks.size() > 1 ) {
            throw new UserException("Found more than 1 chunk for read " + read);
        }

        indexWriter.printf("%s\t%d\t%d\n", read.getName(), chunks.get(0).getChunkStart(), chunks.get(0).getChunkEnd());
    }

    private static String getZMWForRead(final String readName) {
        final String[] readNameTokens = readName.split("/", -1);
        if ( readNameTokens.length != 3 ) {
            throw new UserException("Read name " + readName + " does not appear to be a standard PacBio read name");
        }

        return readNameTokens[1];
    }
}
