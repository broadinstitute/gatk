package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.tools.picard.interval.IntervalListScatterer;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.IntStream;

/**
 * This tool takes in intervals via the standard arguments (-L and -XL, but also interval padding, interval set rule etc) of
 * {@link IntervalArgumentCollection} and splits them into equal sub-intervals for scattering.  Except for the actual splitting
 * the GATK engine handles all of the logic.
 *
 * Created by David Benjamin on 4/25/17.
 */
@CommandLineProgramProperties(
        summary = "Split intervals into sub-interval files.",
        oneLineSummary = "Split intervals into sub-interval files.",
        programGroup = VariantProgramGroup.class
)
public class SplitIntervals extends GATKTool {

    public static final String SCATTER_COUNT_SHORT_NAME = "scatter";
    public static final String SCATTER_COUNT_LONG_NAME = "scatter_count";

    public static final String SUBDIVISION_MODE_SHORT_NAME = "mode";
    public static final String SUBDIVISION_MODE_lONG_NAME = "subdivision_mode";


    @Argument(fullName = SCATTER_COUNT_LONG_NAME, shortName = SCATTER_COUNT_SHORT_NAME,
            doc = "scatter count: number of sub-intervals to split into", optional = true)
    private int scatterCount = 1;

    @Argument(fullName = SUBDIVISION_MODE_lONG_NAME, shortName = SUBDIVISION_MODE_SHORT_NAME, doc = "How to divide intervals.")
    private IntervalListScatterer.Mode subdivisionMode = IntervalListScatterer.Mode.INTERVAL_SUBDIVISION;

    @Argument(doc = "The directory into which to write the scattered interval sub-directories.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File outputDir;

    @Override
    public void onTraversalStart() {
        ParamUtils.isPositive(scatterCount, "scatter count must be > 0.");

        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        // in general dictionary will be from the reference, but using -I reads.bam or -F variants.vcf
        // to use the sequence dict from a bam or vcf is also supported
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        final List<SimpleInterval> intervals = hasIntervals() ? intervalArgumentCollection.getIntervals(sequenceDictionary)
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
