package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import picard.util.IntervalListScatterer;

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
 *    The -O argument specifies a directory name for the scatter intervals files. Each file will be named, e.g 0000-scattered.intervals,
 *    0001-scattered.intervals, 0002-scattered.intervals and so on.
 *    The default --scatter_count is 1 and so this value should be changed to utilize the tool's functionality.
 *    Specify --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION to avoid splitting input intervals -- that is, the set
 *    of input intervals is split, but individual intervals are left intact.  This may affect results when using assembly-based callers downstream.
 * </p>
 *
 * */
@CommandLineProgramProperties(
        summary = "Split intervals into sub-interval files.",
        oneLineSummary = "Split intervals into sub-interval files.",
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class SplitIntervals extends GATKTool {

    public static final String SCATTER_COUNT_SHORT_NAME = "scatter";
    public static final String SCATTER_COUNT_LONG_NAME = "scatter-count";

    public static final String SUBDIVISION_MODE_SHORT_NAME = "mode";
    public static final String SUBDIVISION_MODE_lONG_NAME = "subdivision-mode";


    @Argument(fullName = SCATTER_COUNT_LONG_NAME, shortName = SCATTER_COUNT_SHORT_NAME,
            doc = "scatter count: number of output interval files to split into", optional = true)
    private int scatterCount = 1;

    @Argument(fullName = SUBDIVISION_MODE_lONG_NAME, shortName = SUBDIVISION_MODE_SHORT_NAME, doc = "How to divide intervals.")
    private IntervalListScatterer.Mode subdivisionMode = IntervalListScatterer.Mode.INTERVAL_SUBDIVISION;

    @Argument(doc = "The directory into which to write the scattered interval sub-directories.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File outputDir;

    @Override
    public void onTraversalStart() {
        ParamUtils.isPositive(scatterCount, "scatter-count must be > 0.");

        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        // in general dictionary will be from the reference, but using -I reads.bam or -F variants.vcf
        // to use the sequence dict from a bam or vcf is also supported
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? intervalArgumentCollection.getIntervals(sequenceDictionary)
                : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);

        final IntervalList intervalList = new IntervalList(sequenceDictionary);
        intervals.stream().map(si -> new Interval(si.getContig(), si.getStart(), si.getEnd())).forEach(intervalList::add);
        final IntervalListScatterer scatterer = new IntervalListScatterer(subdivisionMode);
        final List<IntervalList> scattered = scatterer.scatter(intervalList, scatterCount, false);

        final DecimalFormat formatter = new DecimalFormat("0000");
        IntStream.range(0, scattered.size()).forEach(n -> scattered.get(n).write(new File(outputDir, formatter.format(n) + "-scattered.intervals")));
    }

    @Override
    public void traverse() { }  // no traversal for this tool!
}
