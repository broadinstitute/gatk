package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.WeightedInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;
import picard.util.IntervalList.IntervalListScatterMode;
import picard.util.IntervalList.IntervalListScatterer;
import spire.math.interval.Overlap;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
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
 *  * <pre>
 *  * gatk SplitIntervals \
 *  *   -R ref_fasta.fa \
 *  *   -L adjacent_intervals.list \
 *  *   --scatter-count 10 \
 *  *   --interval-merging-rule OVERLAPPING_ONLY \
 *  *   --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
 *  *   -O interval-files-folder
 *  * </pre>
 *
 * <p>
 *     Note that adjacent intervals will be merged by default.  In cases where the desired behavior is to apportion a set
 *     of small adjacent intervals with nearly uniform runtime among X new interval lists, the argument
 *     `--interval-merging-rule OVERLAPPING_ONLY` should be included.
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

    public static final String DONT_MIX_CONTIGS_LONG_NAME = "dont-mix-contigs";

    public static final String MIN_CONTIG_SIZE_LONG_NAME = "min-contig-size";

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

    @Argument(fullName = MIN_CONTIG_SIZE_LONG_NAME, doc = "Minimum contig size to keep if getting intervals from the reference", optional = true)
    private int minContigSize = 0;

    @Argument(fullName = SUBDIVISION_MODE_lONG_NAME, shortName = SUBDIVISION_MODE_SHORT_NAME, doc = "How to divide intervals.")
    private IntervalListScatterMode subdivisionMode = IntervalListScatterMode.INTERVAL_SUBDIVISION;

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

    @Argument(doc = "Scattered interval files do not contain intervals from multiple contigs.  This is applied after the initial scatter, so that the requested scatter count is a lower bound on the number of actual scattered files.", fullName = DONT_MIX_CONTIGS_LONG_NAME, optional = true)
    public boolean dontMixContigs = false;

    @Argument(doc = "BED file of genomic weights, represented by the score field", fullName = WEIGHTS_BED_FILE_FULL_NAME, optional = true)
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

        final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? intervalArgumentCollection.getIntervals(sequenceDictionary)
                : IntervalUtils.getAllIntervalsForReference(sequenceDictionary).stream()
                .filter(contig -> contig.getLengthOnReference() >= minContigSize)
                .collect(Collectors.toList());

        final IntervalList intervalList = new IntervalList(sequenceDictionary);

        // the non-weighted mode
        // TODO: replace with test against enum of scatter mode
        if (true) { // non-weighted mode
            intervals.stream().map(si -> new Interval(si.getContig(), si.getStart(), si.getEnd())).forEach(intervalList::add);
        } else { // weighted mode

            // first read the BED
            final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(weightsBedFile.getAbsolutePath(), new BEDCodec(), false);

            // and construct an OverlapDetector with it
            final List<WeightedInterval> weights = new ArrayList<>();

            try {
                bedReader.iterator().stream().map(f -> new WeightedInterval(f.getContig(), f.getStart(), f.getEnd(), f.getScore())).forEach(weights::add);
            } catch (IOException e) {
                throw new GATKException("Error reading BED file", e);
            }
            OverlapDetector<WeightedInterval> od = OverlapDetector.create(weights);


                // weights should be entirely disjoint sets of intervals
            IntervalUtils.validateNoOverlappingIntervals(weights);

            for (SimpleInterval si : intervals) {
                List<WeightedInterval> l = applyWeightsToInterval(intervalList.getHeader(), new Interval(si.getContig(), si.getStart(), si.getEnd()), od, defaultWeightPerBase);
                intervalList.addall(Collections.unmodifiableList(l));
            }
        }

        final IntervalListScatterer scatterer = subdivisionMode.make();
        final List<IntervalList> scattered = scatterer.scatter(intervalList, scatterCount);

        // optionally split interval lists that contain intervals from multiple contigs
        final List<IntervalList> scatteredFinal = !dontMixContigs ? scattered :
                scattered.stream().flatMap(il -> il.getIntervals().stream()
                        .collect(Collectors.groupingBy(Interval::getContig)).entrySet().stream()    // group each interval list into sublists
                        .sorted(Comparator.comparingInt(entry -> sequenceDictionary.getSequenceIndex(entry.getKey())))  // sort entries by contig
                        .map(entry -> entry.getValue()) // discard the keys and just keep the lists of intervals
                        .map(list -> {
                            final IntervalList singleContigList = new IntervalList(sequenceDictionary);
                            singleContigList.addall(list);
                            return singleContigList;
                        })  // turn the intervals back into an IntervalList
                ).collect(Collectors.toList());

        final int maxNumberOfPlaces = Math.max((int)Math.floor(Math.log10(scatterCount-1))+1, numDigits);
        final String formatString = "%0" + maxNumberOfPlaces + "d";
        IntStream.range(0, scatteredFinal.size()).forEach(n -> scatteredFinal.get(n).write(new File(outputDir, prefix + String.format(formatString, n) + extension)));
    }

    @Override
    public void traverse() { }  // no traversal for this tool!

    // TODO: sort out real plumbing for SAMFileHeader
    public static List<WeightedInterval> applyWeightsToInterval(SAMFileHeader header, Interval interval, OverlapDetector<WeightedInterval> weightsOverlapDetector, long defaultWeightPerBase) {
        List<WeightedInterval> outputIntervals = new ArrayList<>();

        for( WeightedInterval w : weightsOverlapDetector.getOverlaps(interval) ) {
            float weightPerBase = w.getWeight() / (float) w.length();
            // System.out.println("Weight per base of " + weightPerBase + " from length " + w.length() + " from interval " + w );

            // get the weighted piece and add it to the output
            Interval piece = interval.intersect(w);
            outputIntervals.add(new WeightedInterval(piece, (float) piece.length() * weightPerBase ));
        }

        // now we need to emit any pieces of the interval that were not covered by a weight
        IntervalList original = new IntervalList(header);
        original.add(interval);

        IntervalList weighted = new IntervalList(header);
        weighted.addall(Collections.unmodifiableList(outputIntervals));

        IntervalList uncovered = IntervalList.subtract(original, weighted);
        for( Interval piece : uncovered ) {
            outputIntervals.add(new WeightedInterval(piece, (float) piece.length() * defaultWeightPerBase ));
        }

        // TODO: sanity check that all bases in the input interval are covered by the output Intervals

        IntervalUtils.validateNoOverlappingIntervals(outputIntervals);

        return outputIntervals;
    }
}
