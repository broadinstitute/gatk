package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.TestSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;

import java.io.*;


/**
 * Spark troubleshooting utility that converts a headerless Hadoop bam shard (eg., a part0000, part0001, etc. file
 * produced by {@link ReadsSparkSink}) into a readable bam file by adding a header and a BGZF terminator.
 *
 * This tool is not intended for use with Hadoop bam shards that already have a header -- these shards are
 * already readable using samtools. Currently {@link ReadsSparkSink} saves the "shards" with a header for the
 * {@link ReadsWriteFormat#SHARDED} case, and without a header for the {@link ReadsWriteFormat#SINGLE} case.
 */
@CommandLineProgramProperties(
        summary = "Convert a headerless hadoop bam shard into a readable bam",
        oneLineSummary = "Convert a headerless hadoop bam shard into a readable bam",
        programGroup = TestSparkProgramGroup.class
)
public final class ConvertHeaderlessHadoopBamShardToBam extends CommandLineProgram {

    @Argument(shortName = "bamShard", fullName = "bamShard", doc = "Headerless Hadoop BAM shard to be converted into a readable BAM", optional = false)
    private File bamShard = null;

    @Argument(shortName = "bamWithHeader", fullName = "bamWithHeader", doc = "Well-formed BAM whose header to use for the converted fragment", optional = false)
    private File bamWithHeader = null;

    @Argument(shortName = "O", fullName = "outputBam", doc = "Location to write the converted BAM shard", optional = false)
    private File outputBam = null;

    @Override
    protected Object doWork(){
        SAMFileHeader header = null;
        try ( final SamReader headerReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamWithHeader) ) {
            header = headerReader.getFileHeader();
        }
        catch ( IOException e ) {
            throw new UserException("Error reading header from " + bamWithHeader.getAbsolutePath(), e);
        }

        SparkUtils.convertHeaderlessHadoopBamShardToBam(bamShard, header, outputBam);
        return null;
    }
}
