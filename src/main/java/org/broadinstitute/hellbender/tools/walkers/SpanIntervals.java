package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.IntStream;

/**
 * This tool takes in intervals via the standard arguments of
 * {@link IntervalArgumentCollection} and splits them into interval files for scattering. The resulting files contain
 * equal number of bases.
 *
 * <p>Standard GATK engine arguments include -L and -XL, interval padding, and interval set rule etc.
 * For example, for the -L argument, the tool accepts GATK-style intervals (.list or .intervals), BED files
 * and VCF files.  See --subdivision-mode parameter for more options.</p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk SplitIntervals \
 *   -R ref_fasta.fa \
 *   -L intervals.list \
 *   --scatter-count 10 \
 *   -O interval-files-folder
 * </pre>
 *
 * <p>
 *    The -O argument specifies a directory name for the scatter intervals files. Each file will be named, e.g 0000-scattered.interval_list,
 *    0001-scattered.interval_list, 0002-scattered.interval_list and so on.
 *    The default --scatter-count is 1 and so this value should be changed to utilize the tool's functionality.
 *    Specify --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION to avoid splitting input intervals -- that is, the set
 *    of input intervals is split, but individual intervals are left intact.  This may affect results when using assembly-based callers downstream.
 * </p>
 *
 * */
@CommandLineProgramProperties(
        summary = "Output a Picard interval list containing one interval per contig spanning all the input intervals.",
        oneLineSummary = "Output a Picard interval list containing one interval per contig spanning all the input intervals.",
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class SpanIntervals extends GATKTool {

    public static final String INTERVAL_FILE_EXTENSION_FULL_NAME = "extension";

    public static final String PICARD_INTERVAL_FILE_EXTENSION = "interval_list";
    public static final String DEFAULT_EXTENSION = "-scattered." + PICARD_INTERVAL_FILE_EXTENSION;

    @Argument(doc = "The file into which to write the spanning interval(s).",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File outputFilename;

    @Override
    public void onTraversalStart() {

        //TODO: check the extension and suggest using a Picard interval list file extension

        // in general dictionary will be from the reference, but using -I reads.bam or -F variants.vcf
        // to use the sequence dict from a bam or vcf is also supported
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        if (!hasUserSuppliedIntervals()) {
            throw new IllegalArgumentException("Input intervals must be supplied with the -L argument");
        }
        final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);

        final List<SimpleInterval> spanningIntervals = IntervalUtils.getSpanningIntervals(intervals, sequenceDictionary);

        final IntervalList intervalList = new IntervalList(sequenceDictionary);
        spanningIntervals.stream().map(si -> new Interval(si.getContig(), si.getStart(), si.getEnd())).forEach(intervalList::add);

        intervalList.write(outputFilename);
    }

    @Override
    public void traverse() { }  // no traversal for this tool!
}
