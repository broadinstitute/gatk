package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import org.apache.logging.log4j.LogManager;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

/**
 * (Internal) Collects read metrics relevant to structural variant discovery
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It executes the first step
 * in the workflow that the StructuralVariationDiscoveryPipelineSpark tool undertakes, but is packaged as a separately
 * runnable tool for the convenience of developers.</p>
 * <p>This tool takes a SAM/BAM/CRAM as input and calculates metrics about the reads:
 * fragment length statistics by read group, mean length, coverage, partition statistics, etc.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of reads aligned to reference.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file describing the collected metrics.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk CalcMetadataSpark \
 *     -I input_reads.bam \
 *     -O statistics.txt
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Collects read metrics relevant to structural variant discovery",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It executes the first step" +
        " in the workflow that the StructuralVariationDiscoveryPipelineSpark tool undertakes, but is packaged as a separately" +
        " runnable tool for the convenience of developers." +
        " This tool takes a SAM/BAM/CRAM as input and calculates metrics about the reads:" +
        " fragment length statistics by read group, mean length, coverage, partition statistics, etc.",
        programGroup = DiagnosticsAndQCProgramGroup.class)
public class CalcMetadataSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "output file for metadata", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "write metadata as serialized binary data, rather than as human-readable text",
            fullName = "write-as-binary")
    private boolean writeAsBinary = false;

    @Override public boolean requiresReads() { return true; }

    @Override protected void runTool( final JavaSparkContext ctx ) {
        final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params =
                new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();
        final ReadMetadata readMetadata =
                new ReadMetadata(Collections.emptySet(), getHeaderForReads(),
                        10000, getUnfilteredReads(),
                        new SVReadFilter(params), LogManager.getLogger(CalcMetadataSpark.class));
        if ( !writeAsBinary ) {
            ReadMetadata.writeMetadata(readMetadata, outputFile);
        } else {
            try ( final OutputStream os = new FileOutputStream(outputFile) ) {
                new ReadMetadata.Serializer().write(new Kryo(), new Output(os), readMetadata);
            } catch ( final IOException ioe ) {
                throw new UserException("Can't create output file " + outputFile, ioe);
            }
        }
    }
}
