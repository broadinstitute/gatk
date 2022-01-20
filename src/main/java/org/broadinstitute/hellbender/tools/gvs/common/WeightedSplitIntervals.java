package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.collections4.iterators.PushbackIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


@CommandLineProgramProperties(
        summary = "Split intervals into equally weight sub-interval files",
        oneLineSummary = "Split intervals into equally weight sub-interval files",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)

public class WeightedSplitIntervals extends GATKTool {

    public static final String SCATTER_COUNT_SHORT_NAME = "scatter";
    public static final String SCATTER_COUNT_LONG_NAME = "scatter-count";

    public static final String INTERVAL_FILE_PREFIX_FULL_NAME = "interval-file-prefix";
    public static final String INTERVAL_FILE_EXTENSION_FULL_NAME = "extension";

    public static final String DEFAULT_PREFIX = "";
    public static final String PICARD_INTERVAL_FILE_EXTENSION = "interval_list";
    public static final String DEFAULT_EXTENSION = "-scattered." + PICARD_INTERVAL_FILE_EXTENSION;

    public static final String INTERVAL_NUMBER_OF_DIGITS_FULL_NAME = "interval-file-num-digits";
    public static final int DEFAULT_NUMBER_OF_DIGITS = 4;  //to preserve backward compatibility

    public static final String WEIGHTS_BED_FILE_FULL_NAME = "weight-bed-file";
    public static final String DEFAULT_WEIGHT_FULL_NAME = "default-weight-per-base";

    @Argument(fullName = SCATTER_COUNT_LONG_NAME, shortName = SCATTER_COUNT_SHORT_NAME,
            doc = "scatter count: number of output interval files to split into", optional = true)
    private int scatterCount = 1;

    @Argument(doc = "The directory into which to write the scattered interval sub-directories.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File outputDir;

    @Argument(doc = "Prefix to use when writing interval files", fullName = INTERVAL_FILE_PREFIX_FULL_NAME, optional = true)
    public String prefix = DEFAULT_PREFIX;

    @Argument(doc = "Extension to use when writing interval files", fullName = INTERVAL_FILE_EXTENSION_FULL_NAME, optional = true)
    public String extension = DEFAULT_EXTENSION;

    @Argument(doc = "Number of digits to use when writing interval files", fullName = INTERVAL_NUMBER_OF_DIGITS_FULL_NAME, minValue = 1, optional = true)
    public int numDigits = DEFAULT_NUMBER_OF_DIGITS;

    @Argument(doc = "Scattered interval files do not contain intervals from multiple contigs.  This is applied after the initial scatter, so that the requested scatter count is a lower bound on the number of actual scattered files.", fullName = "dont-mix-contigs", optional = true)
    public boolean dontMixContigs = false;

    @Argument(doc = "BED file of genomic weights, represented by the score field", fullName = WEIGHTS_BED_FILE_FULL_NAME, optional = false)
    public File weightsBedFile;

    @Argument(doc = "Default weight (per base) to use if weight not found in BED file", fullName = DEFAULT_WEIGHT_FULL_NAME, optional = true)
    public long defaultWeightPerBase = 0L;

    @Override
    public void onTraversalStart() {
        ParamUtils.isPositive(scatterCount, "scatter-count must be > 0.");

        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        // in general dictionary will be from the reference, but using -I reads.bam or -F variants.vcf
        // to use the sequence dict from a bam or vcf is also supported
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        final List<SimpleInterval> inputIntervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        final IntervalList intervalsWithWeights = preprocessIntervalsWithWeights(sequenceDictionary, inputIntervals);

        // calculate total weights and target weight per shard
        float totalWeight = 0;
        for (Interval i : intervalsWithWeights) {
            totalWeight += ((WeightedInterval) i).getWeight();
        }
        float targetWeightPerScatter = totalWeight / (float) scatterCount;

        System.out.println("Total Weight:"+totalWeight + " scatterCount: " + scatterCount + " target:" + targetWeightPerScatter);

        final int maxNumberOfPlaces = Math.max((int)Math.floor(Math.log10(scatterCount-1))+1, numDigits);
        final String formatString = "%0" + maxNumberOfPlaces + "d";

        int scatterPiece = 0;
        float cumulativeWeight = 0;
        IntervalList currentList = new IntervalList(sequenceDictionary);
        String lastContig = null;

        PushbackIterator<Interval> iter = new PushbackIterator<>(intervalsWithWeights.iterator());
        while (iter.hasNext()) {
            WeightedInterval wi = (WeightedInterval) iter.next();

            // if we're not mixing contigs, but we switched contigs, emit the list
            if (dontMixContigs && lastContig != null && lastContig != wi.getContig() ) {
                // write out the current list (uniqued and sorted) and start a new one
                currentList.uniqued().sorted().write(new File(outputDir, prefix + String.format(formatString, scatterPiece++) + extension));
                currentList = new IntervalList(sequenceDictionary);
                lastContig = wi.getContig();
                cumulativeWeight = 0;
            }

            // if the interval fits completely, just add it
            if (cumulativeWeight + wi.getWeight() <= targetWeightPerScatter ) {
                cumulativeWeight += wi.getWeight();
                currentList.add(wi);

            // if it would push us over the edge
            } else {
                // add a piece of it
                float remainingSpace = targetWeightPerScatter - cumulativeWeight;

                // how many bases can we take?
                int basesToTake = (int) Math.floor(remainingSpace /  wi.getWeightPerBase());

                // split and add the first part into this list
                WeightedInterval[] pair = wi.split(basesToTake);
                currentList.add(pair[0]);

                // push the remainder back onto the iterator
                iter.pushback(pair[1]);

                // add uniqued, sorted output list and reset
                currentList.uniqued().sorted().write(new File(outputDir, prefix + String.format(formatString, scatterPiece++) + extension));
                currentList = new IntervalList(sequenceDictionary);
                lastContig = wi.getContig();
                cumulativeWeight = 0;
            }
        }
        // write the final list
        currentList.uniqued().sorted().write(new File(outputDir, prefix + String.format(formatString, scatterPiece++) + extension));
    }

    @Override
    public void traverse() { }  // no traversal for this tool!

    protected IntervalList preprocessIntervalsWithWeights(SAMSequenceDictionary sequenceDictionary, Collection<SimpleInterval> intervals) {
        OverlapDetector<WeightedInterval> od;
        try {
            // read the BED of weights
            FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(weightsBedFile.toPath().toUri().toString(), new BEDCodec(), false);
            List<WeightedInterval> weights = new ArrayList<>();
            bedReader.iterator().stream().map(f -> new WeightedInterval(f.getContig(), f.getStart(), f.getEnd(), f.getScore())).forEach(weights::add);

            // weights should be entirely disjoint sets of intervals
            IntervalUtils.validateNoOverlappingIntervals(weights);

            // create the overlap detector
            od = OverlapDetector.create(weights);
        } catch (IOException e) {
            throw new GATKException("Error reading BED file", e);
        }

        final IntervalList intervalList = new IntervalList(sequenceDictionary);
        for (SimpleInterval si : intervals) {
            List<WeightedInterval> l = applyWeightsToInterval(sequenceDictionary, new Interval(si.getContig(), si.getStart(), si.getEnd()), od, defaultWeightPerBase);
            intervalList.addall(Collections.unmodifiableList(l));
        }

        return intervalList.sorted();
    }

    protected static List<WeightedInterval> applyWeightsToInterval(SAMSequenceDictionary sequenceDictionary, Interval interval, OverlapDetector<WeightedInterval> weightsOverlapDetector, long defaultWeightPerBase) {
        List<WeightedInterval> outputIntervals = new ArrayList<>();

        for( WeightedInterval w : weightsOverlapDetector.getOverlaps(interval) ) {
            float weightPerBase = w.getWeight() / (float) w.length();
            // System.out.println("Weight per base of " + weightPerBase + " from length " + w.length() + " from interval " + w );

            // get the weighted piece and add it to the output
            Interval piece = interval.intersect(w);
            outputIntervals.add(new WeightedInterval(piece, (float) piece.length() * weightPerBase ));
        }

        // now we need to emit any pieces of the interval that were not covered by a weight
        IntervalList original = new IntervalList(sequenceDictionary);
        original.add(interval);

        IntervalList weighted = new IntervalList(sequenceDictionary);
        weighted.addall(Collections.unmodifiableList(outputIntervals));

        IntervalList uncovered = IntervalList.subtract(original, weighted);
        for( Interval piece : uncovered ) {
            outputIntervals.add(new WeightedInterval(piece, (float) piece.length() * defaultWeightPerBase ));
        }

        IntervalUtils.validateNoOverlappingIntervals(outputIntervals);

        return outputIntervals;
    }
}