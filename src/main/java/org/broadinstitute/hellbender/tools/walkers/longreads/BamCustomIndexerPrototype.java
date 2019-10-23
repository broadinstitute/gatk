package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.SAMFileSource;
import htsjdk.samtools.SAMRecord;
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
public class BamCustomIndexerPrototype extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String outputBAM;

    @Argument(fullName = "output-index",
            doc="Write output index to this file")
    public String outputIndexFile;

    public static final String JONN_TEST_BAM = "src/test/resources/chrM_and_chr20_subset.corrected.bam";

    private SAMFileGATKReadWriter outputWriter;

    private PrintWriter indexWriter;

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
        outputWriter.addRead(read);

        final SAMRecord unwrapped = read.convertToSAMRecord(getHeaderForReads());
        final BAMFileSpan fileSource = (BAMFileSpan)unwrapped.getFileSource().getFilePointer();

        final List<Chunk> chunks = fileSource.getChunks();
        if ( chunks.size() > 1 ) {
            throw new UserException("Found more than 1 chunk for read " + read);
        }

        indexWriter.printf("%s\t%d\t%d\n", read.getName(), chunks.get(0).getChunkStart(), chunks.get(0).getChunkEnd());

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
}
