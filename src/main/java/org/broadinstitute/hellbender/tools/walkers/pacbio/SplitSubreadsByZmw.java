package org.broadinstitute.hellbender.tools.walkers.pacbio;

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

/**
 * Shards PacBio subread BAMs, making sure to keep subreads from the same ZMW in the same file",
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A collection of BAM files where care has been keep reads from the same ZMW in the same file</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Split reads in BAM file by sample name, read group and library name</h4>
 * <pre>
 *   gatk SplitSubreadsByZmw \
 *     -I input.bam \
 *     -O outputDirectory \
 *     --num-reads-per-split 10000
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Shards PacBio subread BAMs, making sure to keep subreads from the same ZMW in the same file",
        oneLineSummary = "Shards PacBio subread BAMs, making sure to keep subreads from the same ZMW in the same file",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class SplitSubreadsByZmw extends ReadWalker {
    public static final String READS_PER_SPLIT_LONG_NAME = "num-reads-per-split";
    public static final String READS_PER_SPLIT_SHORT_NAME = "nr";

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The directory to output SAM/BAM/CRAM files."
    )
    public File OUTPUT_DIRECTORY = new File("");

    @Argument(
            fullName = READS_PER_SPLIT_LONG_NAME,
            shortName = READS_PER_SPLIT_SHORT_NAME,
            doc = "Approximate number of reads per split"
    )
    public int NUM_READS_PER_SPLIT = 10000;

    private SAMFileGATKReadWriter writer;
    private int shardIndex = 0;

    private GATKRead lastRead = null;
    private int numWrittenToShard = 0;

    @Override
    public void onTraversalStart() {
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIRECTORY);
        if ( readArguments.getReadFiles().size() != 1 ) {
            throw new UserException("This tool only accepts a single SAM/BAM/CRAM as input");
        }
        writer = createWriter(shardIndex);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        if (numWrittenToShard < NUM_READS_PER_SPLIT) {
            writer.addRead(read);
        } else {
            int lastZmw = lastRead.getAttributeAsInteger("zm");
            int thisZmw = read.getAttributeAsInteger("zm");

            if (lastZmw != thisZmw) {
                writer.close();

                shardIndex++;
                numWrittenToShard = 0;
                writer = createWriter(shardIndex);
            }

            writer.addRead(read);
        }

        lastRead = read;
        numWrittenToShard++;
    }

    @Override
    public void closeTool() {
        if ( writer != null ) {
            writer.close();
        }
    }

    /**
     * Creates SAMFileWriter instances for the read shards
     * @param shardIndex index of shard
     * @return A SAMFileWriter.
     */
    private SAMFileGATKReadWriter createWriter(final int shardIndex) {
        final String base = FilenameUtils.getBaseName(readArguments.getReadFiles().get(0).getName());
        final String extension = FilenameUtils.getExtension(readArguments.getReadFiles().get(0).getName());
        final String key = String.format("%06d", shardIndex);
        final File outFile = new File(OUTPUT_DIRECTORY, base + "." + key + "." + extension);

        return createSAMWriter(outFile, true);
    }
}
