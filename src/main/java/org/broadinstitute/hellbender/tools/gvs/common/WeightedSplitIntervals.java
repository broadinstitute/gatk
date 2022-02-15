package org.broadinstitute.hellbender.tools.gvs.common;

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
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.SplitIntervals.*;


@CommandLineProgramProperties(
        summary = "Split intervals into equally weighted sub-interval files",
        oneLineSummary = "Split intervals into equally weighted sub-interval files",
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class WeightedSplitIntervals extends GATKTool {
    public static final String INTERVAL_NUMBER_OF_DIGITS_FULL_NAME = "interval-file-num-digits";
    public static final String WEIGHTS_BED_FILE_FULL_NAME = "weight-bed-file";
    public static final String DEFAULT_WEIGHT_FULL_NAME = "default-weight-per-base";

    @Argument(fullName = SCATTER_COUNT_LONG_NAME, shortName = SCATTER_COUNT_SHORT_NAME,
            doc = "scatter count: number of output interval files to split into", optional = true, minValue = 1)
    private int scatterCount = 1;

    @Argument(doc = "The directory into which to write the scattered interval sub-directories.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File outputDir;

    @Argument(doc = "Extension to use when writing interval files", fullName = INTERVAL_FILE_EXTENSION_FULL_NAME, optional = true)
    public String extension = DEFAULT_EXTENSION;

    @Argument(doc = "Number of digits to use when writing interval files", fullName = INTERVAL_NUMBER_OF_DIGITS_FULL_NAME, minValue = 1, optional = true)
    public int numDigits = DEFAULT_NUMBER_OF_DIGITS;

    @Argument(doc = "Scattered interval files do not contain intervals from multiple contigs.  This is applied after the initial scatter, so that the requested scatter count is a lower bound on the number of actual scattered files.", fullName = DONT_MIX_CONTIGS_LONG_NAME, optional = true)
    public boolean dontMixContigs = false;

    @Argument(doc = "BED file of genomic weights, represented by the score field", fullName = WEIGHTS_BED_FILE_FULL_NAME, optional = false)
    public GATKPath weightsBedFile;

    @Argument(doc = "Default weight (per base) to use if weight not found in BED file", fullName = DEFAULT_WEIGHT_FULL_NAME, optional = true)
    public float defaultWeightPerBase = 0;

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        ParamUtils.isPositive(scatterCount, "scatter-count must be > 0.");

        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        final List<SimpleInterval> inputIntervals = getTraversalIntervals();
        final IntervalList intervalsWithWeights = preprocessIntervalsWithWeights(sequenceDictionary, inputIntervals);

        // calculate total weights and target weight per shard
        float totalWeight = 0;
        for (Interval i : intervalsWithWeights) {
            totalWeight += ((WeightedInterval) i).getWeight();
        }
        final float targetWeightPerScatter = totalWeight / (float) scatterCount;

        logger.info("Total Weight:"+totalWeight + " scatterCount: " + scatterCount + " target:" + targetWeightPerScatter);

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
            if (dontMixContigs && lastContig != null && !lastContig.equals(wi.getContig()) ) {
                // write out the current list (uniqued and sorted) and start a new one
                writeIntervalList(formatString, scatterPiece++, currentList);
                currentList = new IntervalList(sequenceDictionary);
                cumulativeWeight = 0;
            }

            // in either case, set the last seen contig
            lastContig = wi.getContig();

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
                writeIntervalList(formatString, scatterPiece++, currentList);
                currentList = new IntervalList(sequenceDictionary);
                cumulativeWeight = 0;
            }
        }
        // write the final list
        writeIntervalList(formatString, scatterPiece++, currentList);
    }

    private void writeIntervalList(String formatString, int scatterPiece, IntervalList currentList) {
        currentList.uniqued().sorted().write(new File(outputDir, String.format(formatString, scatterPiece) + extension));
    }

    @Override
    public void traverse() { }  // no traversal for this tool!

    private OverlapDetector<WeightedInterval> constructWeightsOverlapDetector() {

        // read the BED of weights
        try (FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(weightsBedFile.toPath().toUri().toString(), new BEDCodec(), false)) {
            List<WeightedInterval> weights = new ArrayList<>();
            bedReader.iterator().stream().map(f -> new WeightedInterval(f, f.getScore())).forEach(weights::add);

            // weights should be entirely disjoint sets of intervals
            IntervalUtils.validateNoOverlappingIntervals(weights);

            // create the overlap detector
            return OverlapDetector.create(weights);
        } catch (IOException e) {
            throw new UserException("Error reading BED file", e);
        }
    }

    private IntervalList preprocessIntervalsWithWeights(SAMSequenceDictionary sequenceDictionary, Collection<SimpleInterval> intervals) {
        OverlapDetector<WeightedInterval> od = constructWeightsOverlapDetector();

        final IntervalList intervalList = new IntervalList(sequenceDictionary);
        for (SimpleInterval si : intervals) {
            List<WeightedInterval> l = applyWeightsToInterval(sequenceDictionary, new Interval(si), od, defaultWeightPerBase);
            intervalList.addall(Collections.unmodifiableList(l));
        }

        return intervalList.sorted();
    }

    static List<WeightedInterval> applyWeightsToInterval(SAMSequenceDictionary sequenceDictionary, Interval interval, OverlapDetector<WeightedInterval> weightsOverlapDetector, float defaultWeightPerBase) {
        List<WeightedInterval> outputIntervals = new ArrayList<>();

        for( WeightedInterval w : weightsOverlapDetector.getOverlaps(interval) ) {
            float weightPerBase = w.getWeight() / (float) w.length();

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
