package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTree;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.*;
import scala.Tuple2;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Preprocess SV event intervals for CNV calling
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@ExperimentalFeature
@CommandLineProgramProperties(
        oneLineSummary = "Collates event intervals for cohort CNV calling",
        summary = "This tool uses the output from StructuralVariantIntervalsForCNV from multiple samples to create an intervals file for CNV calling",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class PreprocessSVCohortIntervals extends GATKTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Intervals over which to generate bins", fullName = "target-intervals", optional = true)
    private String eventIntervalsPath;

    @Argument(doc = "Sample counts files", fullName = "counts", optional = true)
    private List<String> countsPathList;

    @Argument(doc = "Output directory", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPath;

    @ArgumentCollection
    private final StructuralVariationDiscoveryArgumentCollection.PreprocessSVCohortIntervals arguments = new StructuralVariationDiscoveryArgumentCollection.PreprocessSVCohortIntervals();

    private SAMSequenceDictionary dictionary;

    @Override
    public void traverse() {
        dictionary = getBestAvailableSequenceDictionary();
        //Note we cannot use -L because it automatically merges the intervals
        final List<SimpleInterval> intervals;
        if (eventIntervalsPath != null) {
            intervals = getIntervals(eventIntervalsPath, dictionary);
        } else {
            intervals = dictionary.getSequences().stream()
                    .map(sequence -> new SimpleInterval(sequence.getSequenceName(), 1, sequence.getSequenceLength()))
                    .collect(Collectors.toList());
        }
        final List<SimpleInterval> binIntervals = stratifyAndBinIntervals(intervals);
        logger.info("Generated " + binIntervals.size() + " intervals at bin size " + arguments.binSize);

        final IntervalList allRebinnedIntervals = new IntervalList(dictionary);
        allRebinnedIntervals.addall(binIntervals.stream()
                .map(interval -> new Interval(interval.getContig(), interval.getStart(), interval.getEnd()))
                .collect(Collectors.toList()));
        final Path allIntervalsOutput = Paths.get(outputPath, "SV-intervals.interval_list");
        allRebinnedIntervals.write(allIntervalsOutput.toFile());

        final Map<String, List<SimpleInterval>> contigToIntervalMap = binIntervals.stream().collect(Collectors.groupingBy(SimpleInterval::getContig));
        for (final String key : contigToIntervalMap.keySet()) {
            final IntervalList rebinnedIntervals = new IntervalList(dictionary);
            rebinnedIntervals.addall(contigToIntervalMap.get(key).stream()
                    .map(interval -> new Interval(interval.getContig(), interval.getStart(), interval.getEnd()))
                    .collect(Collectors.toList()));
            final Path sampleIntervalOutput = Paths.get(outputPath, "SV-intervals-" + key + ".interval_list");
            rebinnedIntervals.write(sampleIntervalOutput.toFile());
        }

        for (int i = 0; i < countsPathList.size(); i++) {
            final String countsPath = countsPathList.get(i);
            logger.info("Loading " + countsPath);
            final SimpleCountCollection counts = SimpleCountCollection.read(new File(countsPath));
            final String sampleName = counts.getMetadata().getSampleName();
            validateCountCollection(counts);
            logger.info("Rebinning counts...");
            final List<SimpleCount> rebinnedCounts = rebinCounts(counts, binIntervals);
            final SimpleCountCollection rebinnedCountCollection = new SimpleCountCollection(counts.getMetadata(), rebinnedCounts);
            final Path sampleCountsOutput = Paths.get(outputPath, sampleName + "-counts.tsv");
            rebinnedCountCollection.write(sampleCountsOutput.toFile());
        }
    }

    private static List<SimpleInterval> getIntervals(final String path, final SAMSequenceDictionary dictionary) {
        if (path != null) {
            final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
            return IntervalUtils.parseIntervalArguments(genomeLocParser, path).stream()
                    .map(genomeLoc -> new SimpleInterval(genomeLoc.getContig(), genomeLoc.getStart(), genomeLoc.getEnd()))
                    .collect(Collectors.toList());
        }
        return Collections.emptyList();
    }

    private List<SimpleCount> rebinCounts(final SimpleCountCollection counts, final List<SimpleInterval> bins) {
        final Map<String, IntervalTree<Tuple2<Integer,Integer>>> binTrees = new HashMap<>(dictionary.size());
        for (final SimpleInterval interval : bins) {
            binTrees.putIfAbsent(interval.getContig(), new IntervalTree<>());
            binTrees.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), new Tuple2<>(0, 0));
        }

        final Iterator<SimpleCount> countsIterator = counts.getRecords().iterator();
        while (countsIterator.hasNext()) {
            final SimpleCount count = countsIterator.next();
            final IntervalTree<Tuple2<Integer,Integer>> tree = binTrees.get(count.getContig());
            if (tree == null) continue;
            final Iterator<IntervalTree.Node<Tuple2<Integer,Integer>>> overlappers = tree.overlappers(count.getStart(), count.getEnd());
            if (overlappers.hasNext()) {
                final IntervalTree.Node<Tuple2<Integer,Integer>> node = overlappers.next();
                if (overlappers.hasNext()) {
                    throw new GATKException("Counts bin " + count.getInterval().toString() + " had multiple overlapping target bins");
                }
                if (count.getStart() < node.getStart() || count.getEnd() > node.getEnd()) {
                    throw new GATKException("Counts bin " + count.getInterval().toString() + " partially overlaps target bin " + node.toString());
                }
                final Tuple2<Integer,Integer> value = node.getValue();
                node.setValue(new Tuple2<>(value._1 + count.getLengthOnReference(), value._2 + count.getCount()));
            }
        }
        final List<SimpleCount> rebinnedCounts = binTrees.entrySet().stream()
                .flatMap(entry -> Utils.stream(entry.getValue().iterator())
                            .filter(node -> node.getValue()._1 == node.getLength())
                            .map(node -> new SimpleCount(new SimpleInterval(entry.getKey(), node.getStart(), node.getEnd()), node.getValue()._2)))
                .collect(Collectors.toList());
        final int numRebinnedCounts = rebinnedCounts.size();
        if (numRebinnedCounts < bins.size()) {
            logger.warn((bins.size() - numRebinnedCounts) + " of " + bins.size() + " bins were not fully covered by counts and will be dropped.");
        }
        return rebinnedCounts;
    }

    private void validateCountCollection(final SimpleCountCollection counts) {
        if (counts.size() == 0) return;
        if (!dictionary.isSameDictionary(counts.getMetadata().getSequenceDictionary())) {
            throw new UserException("Counts dictionary did not match the provided sequence dictionary");
        }
        final List<Integer> binSizes = counts.getIntervalsStream().filter(interval -> interval.getEnd() < dictionary.getSequence(interval.getContig()).getSequenceLength()).mapToInt(SimpleInterval::getLengthOnReference).distinct().boxed().collect(Collectors.toList());
        if (binSizes.size() > 1) {
            final String sizes = String.join(", ", binSizes.stream().map(String::valueOf).collect(Collectors.toList()));
            throw new UserException.BadInput("Counts must have uniform bin size but found sizes " + sizes);
        }
        final int countsBinSize = binSizes.get(0);
        for (final Integer binSize : binSizes) {
            if (binSize % countsBinSize != 0) {
                throw new UserException.BadInput("Counts must have bin size that divides evenly into the target bin sizes");
            }
        }
        final boolean binsNotAligned = counts.getIntervalsStream().anyMatch(interval -> interval.getStart() % countsBinSize != 1);
        if (binsNotAligned) {
            throw new UserException.BadInput("Counts must start at multiples of their bin size, beginning at position 1");
        }
    }

    private List<SimpleInterval> stratifyAndBinIntervals(final List<SimpleInterval> intervals) {
        final int binSize = arguments.binSize;
        final int minSize = binSize * arguments.minBins;
        final int maxSize = arguments.maxBins < 0 ? Integer.MAX_VALUE : binSize * arguments.maxBins;
        final List<GenomeLoc> stratifiedIntervals = intervals.stream()
                .filter(interval -> interval.getLengthOnReference() >= minSize && interval.getLengthOnReference() < maxSize)
                .map(interval -> SVIntervalUtils.convertToSVInterval(interval, dictionary))
                .map(interval -> SVIntervalUtils.getPaddedInterval(interval, arguments.breakpointPadding * binSize, dictionary))
                .map(interval -> SVIntervalUtils.convertToGenomeLoc(interval, dictionary))
                .collect(Collectors.toList());
        final List<SimpleInterval> mergedIntervals = IntervalUtils.mergeIntervalLocations(stratifiedIntervals, IntervalMergingRule.ALL).stream()
                .map(SVIntervalUtils::convertToSimpleInterval)
                .collect(Collectors.toList());
        final List<SimpleInterval> binnedIntervals = mergedIntervals.stream().flatMap(interval -> {
            final int start = interval.getStart() - (interval.getStart() % binSize);
            final int end = Math.min(interval.getEnd() - (interval.getStart() % binSize) + binSize, dictionary.getSequence(interval.getContig()).getSequenceLength());
            final int numBins = (end - start) / binSize;
            return IntStream.range(0, numBins).mapToObj(binIndex -> new SimpleInterval(interval.getContig(), start + binIndex * binSize + 1, start + (binIndex + 1) * binSize));
        }).collect(Collectors.toList());
        return binnedIntervals;
    }


}
