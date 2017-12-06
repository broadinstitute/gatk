package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.pipelines.PrintReadsSpark;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Find reads by name.
 *
 * <p>Reads a file of read (i.e., template) names, and searches an input file of reads to find names that match.
 * The matching reads are copied to an output file.
 * Unlike {@link PrintReadsSpark}, this tool does not require reads to be coordinate-sorted.</p>
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
 */
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Find reads by name.",
        summary =
            "Reads a file of read (i.e., template) names, and searches an input file of reads to find names that match." +
            " The matching reads are copied to an output file. Unlike PrintReadsSpark, this tool does not require" +
            " reads to be coordinate-sorted.",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class ExtractOriginalAlignmentRecordsByNameSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "file containing list of read names", shortName = "names",
            fullName = "read-name-file")
    private String readNameFile;

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

        final JavaRDD<GATKRead> reads =
                getUnfilteredReads().filter(read -> namesToLookForBroadcast.getValue().contains(read.getName())).cache();
        writeReads(ctx, outputSAM, reads);

        logger.info("Found " + reads.count() + " alignment records for " +
                    namesToLookForBroadcast.getValue().size() + " unique read names.");
    }

    private Set<String> parseReadNames() {

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(readNameFile))) ) {
            return rdr.lines().map(s -> s.replace("^@", "")
                                         .replace("/1$", "")
                                         .replace("/2$", ""))
                    .collect(Collectors.toCollection(HashSet::new));
        } catch ( final IOException ioe ) {
            throw new GATKException("Unable to read names file from " + readNameFile, ioe);
        }
    }
}
