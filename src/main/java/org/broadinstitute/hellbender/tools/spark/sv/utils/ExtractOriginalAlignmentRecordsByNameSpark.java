package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Subsets reads by name (basically a parallel version of "grep -f", or "grep -vf")
 *
 * <p>Reads a file of read (i.e., template) names, and searches a SAM/BAM/CRAM to find names that match.
 * The matching reads are copied to an output file.</p>
 * <p>Unlike FilterSamReads (Picard), this tool can take input reads in any order
 * (e.g., unsorted or coordinate-sorted).</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of read names.</li>
 *     <li>An input file of reads.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A file containing matching reads.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk ExtractOriginalAlignmentRecordsByNameSpark \
 *     -I input_reads.bam \
 *     --read-name-file read_names.txt \
 *     -O output_reads.bam
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Subsets reads by name",
        summary =
            "Reads a file of read (i.e., template) names, and searches a SAM/BAM/CRAM to find names that match." +
            " The matching reads are copied to an output file." +
            " Unlike FilterSamReads (Picard), this tool can take input reads in any order" +
            " (e.g., unsorted or coordinate-sorted).",
        programGroup = ReadDataManipulationProgramGroup.class)
public final class ExtractOriginalAlignmentRecordsByNameSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "file containing list of read names", shortName = "f", fullName = "read-name-file")
    private String readNameFile;

    @Argument(doc = "invert the list, i.e. filter out reads whose name appear in the given file",
            shortName = "v", fullName = "invert-match", optional = true)
    private Boolean invertFilter = false;

    @Argument(doc = "file to write reads to", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputSAM;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final Broadcast<Set<String>> namesToLookForBroadcast = ctx.broadcast(parseReadNames());

        final Function<GATKRead, Boolean> predicate = getGatkReadBooleanFunction(namesToLookForBroadcast, invertFilter);

        final JavaRDD<GATKRead> reads = getUnfilteredReads().filter(predicate).cache();
        writeReads(ctx, outputSAM, reads, getHeaderForReads());

        logger.info("Found " + reads.count() + " alignment records for " +
                    namesToLookForBroadcast.getValue().size() + " unique read names.");
    }

    private static Function<GATKRead, Boolean> getGatkReadBooleanFunction(final Broadcast<Set<String>> namesToLookForBroadcast,
                                                                          final boolean invertFilter) {
        return invertFilter
                ? read -> !namesToLookForBroadcast.getValue().contains(read.getName())
                : read -> namesToLookForBroadcast.getValue().contains(read.getName());
    }

    private Set<String> parseReadNames() {

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(readNameFile))) ) {
            return rdr.lines().map(s -> s.replaceAll("^@", "")
                                         .replaceAll("/1$", "")
                                         .replaceAll("/2$", ""))
                    .collect(Collectors.toCollection(HashSet::new));
        } catch ( final IOException ioe ) {
            throw new GATKException("Unable to read names file from " + readNameFile, ioe);
        }
    }
}
