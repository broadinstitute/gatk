package org.broadinstitute.hellbender.tools.spark.longread;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.List;

/**
 * Collect some metrics on alignments produced from long read alignments.
 *
 * <ul>
 * <li>Total clipped bases</li>
 * <li>Total count of unmapped queries</li>
 * <li>Total count of supplementary (SAM flag 2048) alignment records</li>
 * <li>Total count of secondary (SAM flag 256) alignment records</li>
 * <li>average number CIGAR operations per 100 base pair (mapped queries)</li>
 * </ul>
 * <p>
 * In addition to summary metrics, we also collect a metric that aims to measure contiguity of the alignments.
 * That is, for one alignment, there typically are several operations in the CIGAR.
 * We construct a vector of lengths of the alignment blocks, then normalize it against the total non-clipped bases.
 * This is output as a vector of 2-digit integers (that is, the fraction with precision up to 2 decimal places then * 100)
 * for each SAM record.
 * The output format is
 * "read_name\tclipped_start\tclipped_end\tSAM_flag\tvector_separated_by_,"
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Collect some metrics on alignments produced from long read alignments",
        summary = "Collect some metrics on alignments produced from long read alignments",
        programGroup = LongReadAnalysisProgramGroup.class)
public class ExtractAlignmentsWithLongCigarsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String output;

    @Argument(doc = "cigar operation threshold",
            shortName = "c", fullName = "cigar-opcnt")
    private int cigarOptCnt;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This is overriden because the alignment might be not "WellFormed".
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(new ReadFilterLibrary.AllowAllReadsReadFilter());
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final SAMFileHeader headerForReads = getHeaderForReads();
        final int cigarOperationsThreshold = cigarOptCnt;
        final JavaRDD<GATKRead> longCigarAlignments = getReads().filter(read -> read.getCigarElements().size() > cigarOperationsThreshold).cache();

        writeReads(ctx,output,longCigarAlignments,headerForReads,true);

        final List<String> longCigarReadNames = longCigarAlignments
                .map(read -> {
                    final String cnt = String.valueOf(read.getCigarElements().size());
                    return read.getName() + "\t" + cnt + "\t" + read.getCigar().toString();
                }).collect();
        final String readNamesOut = output.replace(".bam", ".readnames.txt");
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(readNamesOut))) ) {
            for (final String s : longCigarReadNames) { writer.write(s + "\n"); }
        } catch ( final IOException ioe ) {
            throw new UserException.CouldNotCreateOutputFile("Can't write intervals file " + readNamesOut, ioe);
        }
    }
}
