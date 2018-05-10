package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTree;
import org.broadinstitute.barclay.argparser.Argument;
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
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.collections.Lists;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
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

    @Argument(doc = "Sample counts files", fullName = "counts")
    private List<String> countsPathList;

    @Argument(doc = "Output directory", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPath;

    @Argument(doc = "Sample names", fullName = "sample-name")
    private List<String> sampleNameList;

    private final StructuralVariationDiscoveryArgumentCollection.PreprocessSVCohortIntervals arguments = new StructuralVariationDiscoveryArgumentCollection.PreprocessSVCohortIntervals();
    private SAMSequenceDictionary dictionary;

    public static final List<Integer> BIN_SIZES = Lists.newArrayList(100, 1000, 10000);
    public static final int MIN_BINS_PER_EVENT = 10;

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    private String getBinDirectory(final int binSize) {
        return Paths.get(outputPath, "bin-" + binSize + "bp").toAbsolutePath().toString();
    }

    private String getSampleDirectory(final int binSize, final int sampleIndex) {
        return Paths.get(getBinDirectory(binSize), "SAMPLE_" + sampleIndex).toAbsolutePath().toString();
    }

    @Override
    public void traverse() {
        dictionary = getBestAvailableSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(dictionary); // intervalPathList.stream().map(path -> getIntervals(path, dictionary)).collect(Collectors.toList());
        final List<Tuple2<Integer,List<SimpleInterval>>> binsSizesAndIntervals = stratifyAndBinIntervals(intervals);
        for (final Tuple2<Integer,List<SimpleInterval>> pair : binsSizesAndIntervals) {
            logger.info("Generated " + pair._2.size() + " intervals at bin size " + pair._1);
        }

        try {
            for (final Integer binSize : BIN_SIZES) {
                final String binDirName = getBinDirectory(binSize);
                IOUtils.createDirectory(binDirName);
                for (int i = 0; i < countsPathList.size(); i++) {
                    final String sampleDirName = getSampleDirectory(binSize, i);
                    IOUtils.createDirectory(sampleDirName);
                }
            }
        } catch (final IOException e) {
            throw new GATKException("Could not create output directory: " + e.getMessage());
        }

        /*
        for (final Tuple2<Integer,List<SimpleInterval>> pair : binsSizesAndIntervals) {
            final int binSize = pair._1;
            final Map<String,List<SimpleInterval>> intervalListMap = new HashMap<>(SVUtils.hashMapCapacity(dictionary.size()));
            for (final SimpleInterval interval : pair._2) {
                intervalListMap.putIfAbsent(interval.getContig(), new ArrayList<>());
                intervalListMap.get(interval.getContig()).add(interval);
            }
            final String binDirName = getBinDirectory(binSize);
            for (final String key : intervalListMap.keySet()) {
                final SimpleIntervalCollection rebinnedIntervals = new SimpleIntervalCollection(locatableMetadata, intervalListMap.get(key));
                final Path sampleIntervalOutput = Paths.get(binDirName, "SV-intervals-" + key + ".interval_list");
                rebinnedIntervals.write(sampleIntervalOutput.toFile());
            }
        }
        */

        for (int i = 0; i < countsPathList.size(); i++) {
            final String countsPath = countsPathList.get(i);
            final String sampleName = sampleNameList.get(i);
            logger.info("Loading " + countsPath);
            final SimpleCountCollection counts = SimpleCountCollection.read(new File(countsPath));
            validateCountCollection(counts);
            for (final Tuple2<Integer,List<SimpleInterval>> pair : binsSizesAndIntervals) {
                final int binSize = pair._1;
                logger.info("Rebinning at bin size " + binSize + "...");
                final String binDirName = getBinDirectory(binSize);
                final String sampleDirName = getSampleDirectory(binSize, i);
                final Map<String,SimpleCountCollection> rebinnedCountsMap = rebinCounts(counts, pair._2);
                for (final String key : rebinnedCountsMap.keySet()) {
                    final SimpleCountCollection rebinnedCounts = rebinnedCountsMap.get(key);
                    if (rebinnedCounts.size() > 0) {
                        final Path sampleCountsOutput = Paths.get(sampleDirName, sampleName + "-counts-" + key + ".tsv");
                        rebinnedCounts.write(sampleCountsOutput.toFile());
                    }
                    if (i == 0) { //TODO assuming for now that counts have exact same bins for all samples
                        final IntervalList rebinnedIntervals = new IntervalList(dictionary);
                        rebinnedIntervals.addall(rebinnedCounts.getIntervalsStream().map(interval -> new Interval(interval.getContig(), interval.getStart(), interval.getEnd())).collect(Collectors.toList()));
                        final Path sampleIntervalOutput = Paths.get(binDirName, "SV-intervals-" + key + ".interval_list");
                        rebinnedIntervals.write(sampleIntervalOutput.toFile());
                    }
                }
            }
        }
    }

    private Map<String,SimpleCountCollection> rebinCounts(final SimpleCountCollection counts, final List<SimpleInterval> bins) {
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
        final Map<String,SimpleCountCollection> rebinnedCounts = binTrees.entrySet().stream()
                .map(entry -> {
                    final List<SimpleCount> entryCounts = Utils.stream(entry.getValue().iterator())
                            .filter(node -> node.getValue()._1 == node.getLength())
                            .map(node -> new SimpleCount(new SimpleInterval(entry.getKey(), node.getStart(), node.getEnd()), node.getValue()._2))
                            .collect(Collectors.toList());
                    return new HashMap.SimpleEntry<>(entry.getKey(), entryCounts);
                }).map(entry -> new HashMap.SimpleEntry<>(entry.getKey(), new SimpleCountCollection(counts.getMetadata(), entry.getValue())))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        final int numRebinnedCounts = rebinnedCounts.values().stream().mapToInt(SimpleCountCollection::size).sum();
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

    private List<Tuple2<Integer,List<SimpleInterval>>> stratifyAndBinIntervals(final List<SimpleInterval> intervals) {
        final int numBinSizes = BIN_SIZES.size();
        final List<Tuple2<Integer,List<SimpleInterval>>> binsSizesAndIntervals = new ArrayList<>(numBinSizes);
        for (int i = 0; i < numBinSizes; i++) {
            final int binSize = BIN_SIZES.get(i);
            final int minSize = binSize * MIN_BINS_PER_EVENT;
            final int maxSize;
            if (i == BIN_SIZES.size() - 1) {
                maxSize = Integer.MAX_VALUE;
            } else {
                maxSize = minSize * BIN_SIZES.get(i + 1);
            }
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
            binsSizesAndIntervals.add(new Tuple2<>(binSize, binnedIntervals));
        }
        return binsSizesAndIntervals;
    }


}
