package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;
import org.broadinstitute.hellbender.metrics.PerUnitMetricCollector;
import org.broadinstitute.hellbender.metrics.SAMRecordMultiLevelCollector;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.*;

/**
 * TargetMetrics, are metrics to measure how well we hit specific targets (or baits) when using a targeted sequencing process like hybrid selection
 * or Targeted PCR Techniques (TSCA).  TargetMetrics at the moment are the metrics that are shared by both HybridSelection and TargetedPcrMetrics.
 *
 * TargetMetricsCollector collects for a run these common metrics and can be sub-classed to provide metrics more specific to a targeted sequencing
 * run.
 *
 * Note: Probe is the name I've used to indicate the bait set or amplicon set (e.g. the individual technological units used to target specific
 * sites).
 *
 * @author Jonathan Burke
 */
public abstract class TargetMetricsCollector<METRIC_TYPE extends MultiLevelMetrics> extends SAMRecordMultiLevelCollector<METRIC_TYPE, Integer> {

    // What is considered "near" to the bait
    private static final int NEAR_PROBE_DISTANCE = 250;

    //If perTargetCoverage != null then coverage is computed for each specified target and output to this file
    private final File perTargetCoverage;

    //The name of the set of probes used
    private final String probeSetName;

    private static final Log log = Log.getInstance(TargetMetricsCollector.class);

    //The interval list indicating the regions targeted by all probes
    private final IntervalList allProbes;

    //The interval list of the the regions we intend to cover
    private final IntervalList allTargets;

    // Overlap detector for finding overlaps between reads and the experimental targets
    private final OverlapDetector<Interval> targetDetector;

    // Overlap detector for finding overlaps between the reads and the baits (and the near bait space)
    private final OverlapDetector<Interval> probeDetector;

    private Map<Interval,Double> intervalToGc = null;

    //The number of bases within all unique intervals in allProbes
    private final long probeTerritory;

    //The number of bases within all unique intervals found in allTargets
    private final long targetTerritory;

    private final long genomeSize;

    //A map of coverage by target in which coverage is reset every read, this is done
    //so that we can calculate overlap for a read once and the resulting coverage is
    //than added to the cumulative coverage of every collector that collects
    //information on that read
    private Map<Interval, Coverage> coverageByTargetForRead;
    private Coverage [] cov;

    //Converts a targetMetric into a more specific metric of METRIC_TYPE
    public abstract METRIC_TYPE convertMetric(final TargetMetrics targetMetrics);

    /**
     * Since the targeted metrics (HsMetrics, TargetedPcrMetrics,...) share many of the same values as TargetMetrics, this copy will copy all public attributes in targetMetrics
     * to the outputMetrics' attributes of the same name.  If no matching attribute exists in the outputMetrics or the attribute of the target metrics class also is found
     * in targetKeys then it's value is not copied.  Further more, targetKeys and outputKeys are attribute name arrays synchronized by the index.
     * For each target key, targetMetrics.<targetKeys[i]> is assigned to outputMetrics.<outputKeys[i]>
     *
     * @param targetMetrics A metric with values to be copied
     * @param outputMetrics A metrics intended to receive values from targetMetrics
     * @param targetKeys Specific names of attributes of targetMetrics to copy to outputMetrics, each key has a corresponding one in outputKeys
     * @param outputKeys Specific names of the destination attributes of outputMetrics that will be filled with values of outputMetrics, each key has a corresponding one in targetKeys
     * @param <MT> The type of metric of outputMetrics
     */
    protected static <MT extends MetricBase> void reflectiveCopy(final TargetMetrics targetMetrics, final MT outputMetrics, final String[] targetKeys, final String[] outputKeys) {

        if(targetKeys == null || outputKeys == null) {
            if(outputKeys != null) {
                throw new GATKException("Target keys is null but output keys == " + StringUtil.join(",", outputKeys));
            }

            if(targetKeys != null) {
                throw new GATKException("Output keys is null but target keys == " + StringUtil.join(",", targetKeys));
            }
        } else {
            if(targetKeys.length != outputKeys.length) {
                throw new GATKException("Target keys and output keys do not have the same length: " +
                        "targetKeys == (" + StringUtil.join(",", targetKeys) + ") " +
                        "outputKeys == (" + StringUtil.join(",", outputKeys) + ")");
            }
        }

        final Class<?> mtClass = outputMetrics.getClass();
        final Set<Field> targetSet = CollectionUtil.makeSet(TargetMetrics.class.getFields());

        for(final String targetKey : targetKeys) {
            if(targetSet.contains(targetKey)) {
                targetSet.remove(targetKey);
            }
        }

        final Set<String> outputSet = new HashSet<>();
        for(final Field field : outputMetrics.getClass().getFields()) {
            outputSet.add(field.getName());
        }

        for(final Field field : targetSet) {
            if(outputSet.contains(field.getName())) {
                try {
                    final Field outputField = mtClass.getField(field.getName());
                    outputField.set(outputMetrics, field.get(targetMetrics));
                } catch (Exception e) {
                    throw new GATKException("Exception while copying targetMetrics to " + outputMetrics.getClass().getName(), e);
                }
            }
        }

        for(int i = 0; i < targetKeys.length; i++) {
            try {
                Field targetMetricField = TargetMetrics.class.getField(targetKeys[i]);
                Field outputMetricField = mtClass.getField(outputKeys[i]);
                outputMetricField.set(outputMetrics, targetMetricField.get(targetMetrics));
            } catch(final Exception exc) {
                throw new GATKException("Exception while copying TargetMetrics." + targetKeys[i] + " to " + mtClass.getName() + "." + outputKeys[i], exc);
            }
        }
    }

    public TargetMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords, final ReferenceSequenceFile refFile,
                                  final File perTargetCoverage, final IntervalList targetIntervals, final IntervalList probeIntervals, final String probeSetName) {
        this.perTargetCoverage = perTargetCoverage;
        this.probeSetName = probeSetName;

        this.allProbes  = probeIntervals;
        this.allTargets = targetIntervals;

        final List<Interval> uniqueBaits = this.allProbes.uniqued().getIntervals();
        this.probeDetector = new OverlapDetector<>(-NEAR_PROBE_DISTANCE, 0);
        this.probeDetector.addAll(uniqueBaits, uniqueBaits);
        this.probeTerritory = Interval.countBases(uniqueBaits);

        final List<Interval> uniqueTargets = this.allTargets.uniqued().getIntervals();
        targetDetector = new OverlapDetector<>(0,0);
        this.targetDetector.addAll(uniqueTargets, uniqueTargets);
        this.targetTerritory = Interval.countBases(uniqueTargets);

        // Populate the coverage by target map
        int i = 0;
        cov = new Coverage[uniqueTargets.size()];
        this.coverageByTargetForRead = new LinkedHashMap<>(uniqueTargets.size() * 2, 0.5f);
        for (final Interval target : uniqueTargets) {
            final Coverage coverage = new Coverage(target, 0);
            this.coverageByTargetForRead.put(target, coverage);
            cov[i++] = coverage;
        }

        long genomeSizeAccumulator = 0;
        for (final SAMSequenceRecord seq : this.allProbes.getHeader().getSequenceDictionary().getSequences()) {
            genomeSizeAccumulator += seq.getSequenceLength();
        }
        this.genomeSize = genomeSizeAccumulator;


        if (refFile != null) {
            intervalToGc = new HashMap<>();
            for (final Interval target : uniqueTargets) {
                final ReferenceSequence rs = refFile.getSubsequenceAt(target.getContig(), target.getStart(), target.getEnd());
                intervalToGc.put(target, SequenceUtil.calculateGc(rs.getBases()));
            }
        }

        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
        final PerUnitTargetMetricCollector collector =  new PerUnitTargetMetricCollector(probeSetName, coverageByTargetForRead.keySet(),
                                                                                         sample, library, readGroup, probeTerritory, targetTerritory, genomeSize,
                                                                                         intervalToGc);
        if (this.probeSetName != null) {
            collector.setBaitSetName(probeSetName);
        }

        return collector;
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeAllReadCollector() {
        final PerUnitTargetMetricCollector collector = (PerUnitTargetMetricCollector) makeChildCollector(null, null, null);
        if (perTargetCoverage != null) {
            collector.setPerTargetOutput(perTargetCoverage);
        }

        return collector;
    }

    /**
     * Collect the Target Metrics for one unit of "accumulation" (i.e. for one sample, or for one library ...)
     */
    public final class PerUnitTargetMetricCollector implements PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> {

        private final Map<Interval,Double> intervalToGc;
        private File perTargetOutput;

        // A Map to accumulate per-bait-region (i.e. merge of overlapping targets) coverage. */
        private final Map<Interval, Coverage> coverageByTarget;

        private final TargetMetrics metrics = new TargetMetrics();

        /**
         * Constructor that parses the squashed reference to genome reference file and stores the
         * information in a map for later use.
         */
        public PerUnitTargetMetricCollector(final String probeSetName, final Set<Interval> coverageTargets,
                                            final String sample, final String library, final String readGroup,
                                            final long probeTerritory, final long targetTerritory, final long genomeSize,
                                            final Map<Interval, Double> intervalToGc) {
            this.metrics.SAMPLE           = sample;
            this.metrics.LIBRARY          = library;
            this.metrics.READ_GROUP       = readGroup;
            this.metrics.PROBE_SET        = probeSetName;

            metrics.PROBE_TERRITORY  = probeTerritory;
            metrics.TARGET_TERRITORY = targetTerritory;
            metrics.GENOME_SIZE      = genomeSize;

            this.coverageByTarget = new LinkedHashMap<>(coverageTargets.size() * 2, 0.5f);
            for (Interval target : coverageTargets) {
                this.coverageByTarget.put(target, new Coverage(target,0));
            }

            this.intervalToGc = intervalToGc;
        }

        /** If set, the metrics collector will output per target coverage information to this file. */
        public void setPerTargetOutput(final File perTargetOutput) {
            this.perTargetOutput = perTargetOutput;
        }

        /** Sets the name of the bait set explicitly instead of inferring it from the bait file. */
        public void setBaitSetName(final String name) {
            this.metrics.PROBE_SET = name;
        }

        /** Adds information about an individual SAMRecord to the statistics. */
        public void acceptRecord(final SAMRecord rec) {
            // Just plain avoid records that are marked as not-primary
            if (rec.isSecondaryOrSupplementary()) return;

            this.metrics.TOTAL_READS += 1;

            // Check for PF reads
            if (rec.getReadFailsVendorQualityCheckFlag()) {
                return;
            }

            // Prefetch the list of target and bait overlaps here as they're needed multiple times.
            final Collection<Interval> targets;
            final Collection<Interval> probes;

            if (!rec.getReadUnmappedFlag()) {
                final Interval read = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
                targets = targetDetector.getOverlaps(read);
                probes   = probeDetector.getOverlaps(read);
            }
            else {
                targets = null;
                probes = null;
            }

            ++this.metrics.PF_READS;
            this.metrics.PF_BASES += rec.getReadLength();

            // And now calculate the values we need for HS_LIBRARY_SIZE
            if (rec.getReadPairedFlag() && rec.getFirstOfPairFlag() && !rec.getReadUnmappedFlag() && !rec.getMateUnmappedFlag()) {
                if (probes != null && !probes.isEmpty()) {
                    ++this.metrics.PF_SELECTED_PAIRS;
                    if (!rec.getDuplicateReadFlag()) ++this.metrics.PF_SELECTED_UNIQUE_PAIRS;
                }
            }

            // Check for reads that are marked as duplicates
            if (rec.getDuplicateReadFlag()) {
                return;
            }
            else {
                ++this.metrics.PF_UNIQUE_READS;
            }

            // Don't bother with reads that didn't align uniquely
            if (rec.getReadUnmappedFlag() || rec.getMappingQuality() == 0) {
                return;
            }

            this.metrics.PF_UQ_READS_ALIGNED += 1;
            for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                this.metrics.PF_UQ_BASES_ALIGNED += block.getLength();
            }

            final boolean mappedInPair = rec.getReadPairedFlag() && !rec.getMateUnmappedFlag();

            // Find the target overlaps
            if (targets != null && !targets.isEmpty()) {
                for (final Interval target : targets) {
                    final Coverage coverage = this.coverageByTarget.get(target);

                    for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                        final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                        for (int pos=block.getReferenceStart(); pos<=end; ++ pos) {
                            if (pos >= target.getStart() && pos <= target.getEnd()) {
                                ++this.metrics.ON_TARGET_BASES;
                                if (mappedInPair) ++this.metrics.ON_TARGET_FROM_PAIR_BASES;
                                coverage.addBase(pos - target.getStart());
                            }
                        }
                    }
                }
            }

            // Now do the bait overlaps
            int mappedBases = 0;
            for (final AlignmentBlock block : rec.getAlignmentBlocks()) mappedBases += block.getLength();
            int onBaitBases = 0;

            if (probes != null && !probes.isEmpty()) {
                for (final Interval bait : probes) {
                    for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                        final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());

                        for (int pos=block.getReferenceStart(); pos<=end; ++pos) {
                            if (pos >= bait.getStart() && pos <= bait.getEnd()) ++onBaitBases;
                        }
                    }
                }

                this.metrics.ON_PROBE_BASES   += onBaitBases;
                this.metrics.NEAR_PROBE_BASES += (mappedBases - onBaitBases);
            }
            else {
                this.metrics.OFF_PROBE_BASES += mappedBases;
            }

        }

        @Override
        public void finish() {
            metrics.PCT_PF_READS         = metrics.PF_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS      = metrics.PF_UNIQUE_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS_ALIGNED = metrics.PF_UQ_READS_ALIGNED / (double) metrics.PF_UNIQUE_READS;

            final double denominator   = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES + metrics.OFF_PROBE_BASES);

            metrics.PCT_SELECTED_BASES = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES) / denominator;
            metrics.PCT_OFF_PROBE         = metrics.OFF_PROBE_BASES / denominator;
            metrics.ON_PROBE_VS_SELECTED = metrics.ON_PROBE_BASES / (double) (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES);
            metrics.MEAN_PROBE_COVERAGE   = metrics.ON_PROBE_BASES / (double) metrics.PROBE_TERRITORY;
            metrics.FOLD_ENRICHMENT       = (metrics.ON_PROBE_BASES/ denominator) / ((double) metrics.PROBE_TERRITORY / metrics.GENOME_SIZE);

            calculateTargetCoverageMetrics();
            calculateGcMetrics();
        }

        /** Calculates how much additional sequencing is needed to raise 80% of bases to the mean for the lane. */
        private void calculateTargetCoverageMetrics() {
            final short[] depths = new short[(int) this.metrics.TARGET_TERRITORY];  // may not use entire array
            int zeroCoverageTargets = 0;
            int depthIndex = 0;
            double totalCoverage = 0;
            int basesConsidered = 0;

            for (final Coverage c : this.coverageByTarget.values()) {
                if (!c.hasCoverage()) {
                    ++zeroCoverageTargets;
                    continue;
                }

                final short[] targetDepths = c.getDepths();
                basesConsidered += targetDepths.length;

                for (final short depth : targetDepths) {
                    depths[depthIndex++] = depth;
                    totalCoverage += depth;
                }
            }

            this.metrics.MEAN_TARGET_COVERAGE = totalCoverage / basesConsidered;

            // Sort the array (ASCENDING) and then find the base the coverage value that lies at the 80%
            // line, which is actually at 20% into the array now
            Arrays.sort(depths);
            // Note.  basesConsidered can be between 0 and depths.length inclusive.  indexOf80thPercentile will be -1 in the latter case
            final int indexOf80thPercentile = Math.max((depths.length - 1 - basesConsidered) + (int) (basesConsidered * 0.2), 0);
            final int coverageAt80thPercentile = depths[indexOf80thPercentile];
            this.metrics.FOLD_80_BASE_PENALTY = this.metrics.MEAN_TARGET_COVERAGE / coverageAt80thPercentile;
            this.metrics.ZERO_CVG_TARGETS_PCT = zeroCoverageTargets / (double) allTargets.getIntervals().size();

            // Now do the "how many bases at X" calculations.
            int totalTargetBases = 0;
            int targetBases2x  = 0;
            int targetBases10x = 0;
            int targetBases20x = 0;
	        int targetBases30x = 0;
	        int targetBases40x = 0;
	        int targetBases50x = 0;
	        int targetBases100x = 0;

            for (final Coverage c : this.coverageByTarget.values()) {
                for (final short depth : c.getDepths()) {
                    ++totalTargetBases;

                    if (depth >= 2) {
                        ++targetBases2x;
                        if (depth >=10) {
                            ++targetBases10x;
                            if (depth >= 20) {
                                ++targetBases20x;
                                if (depth >=30) {
                                    ++targetBases30x;
	                                if (depth >=40) {
		                                ++targetBases40x;
		                                if (depth >=50) {
			                                ++targetBases50x;
			                                if (depth >=100) {
				                                ++targetBases100x;
			                                }
		                                }
	                                }
                                }
                            }
                        }
                    }
                }
            }

            this.metrics.PCT_TARGET_BASES_2X  = (double) targetBases2x  / (double) totalTargetBases;
            this.metrics.PCT_TARGET_BASES_10X = (double) targetBases10x / (double) totalTargetBases;
            this.metrics.PCT_TARGET_BASES_20X = (double) targetBases20x / (double) totalTargetBases;
	        this.metrics.PCT_TARGET_BASES_30X = (double) targetBases30x / (double) totalTargetBases;
	        this.metrics.PCT_TARGET_BASES_40X = (double) targetBases40x / (double) totalTargetBases;
	        this.metrics.PCT_TARGET_BASES_50X = (double) targetBases50x / (double) totalTargetBases;
	        this.metrics.PCT_TARGET_BASES_100X = (double) targetBases100x / (double) totalTargetBases;
        }

        private void calculateGcMetrics() {
            if (this.intervalToGc != null) {
                log.info("Calculating GC metrics");

                // Setup the output file if we're outputting per-target coverage
                FormatUtil fmt = new FormatUtil();
                final PrintWriter out;
                try {
                    if (perTargetOutput != null) {
                        out = new PrintWriter(perTargetOutput);
                        out.println("chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage\tnormalized_coverage");
                    }
                    else {
                        out = null;
                    }
                }
                catch (IOException ioe) { throw new RuntimeIOException(ioe); }

                final int bins = 101;
                final long[] targetBasesByGc  = new long[bins];
                final long[] alignedBasesByGc = new long[bins];

                for (final Map.Entry<Interval,Coverage> entry : this.coverageByTarget.entrySet()) {
                    final Interval interval = entry.getKey();
                    final Coverage cov = entry.getValue();

                    final double gcDouble = this.intervalToGc.get(interval);
                    final int gc = (int) Math.round(gcDouble * 100);

                    targetBasesByGc[gc]  += interval.length();
                    alignedBasesByGc[gc] += cov.getTotal();

                    if (out != null) {
                        final double coverage = cov.getTotal() / (double) interval.length();

                        out.println(interval.getContig() + "\t" +
                                    interval.getStart() + "\t" +
                                    interval.getEnd() + "\t" +
                                    interval.length() + "\t" +
                                    interval.getName() + "\t" +
                                    fmt.format(gcDouble) + "\t" +
                                    fmt.format(coverage) + "\t" +
                                    fmt.format(coverage / this.metrics.MEAN_TARGET_COVERAGE)
                        );
                    }
                }

                if (out != null) out.close();

                // Total things up
                long totalTarget = 0;
                long totalBases  = 0;
                for (int i=0; i<targetBasesByGc.length; ++i) {
                    totalTarget += targetBasesByGc[i];
                    totalBases  += alignedBasesByGc[i];
                }

                // Re-express things as % of the totals and calculate dropout metrics
                for (int i=0; i<targetBasesByGc.length; ++i) {
                    final double targetPct  = targetBasesByGc[i]  / (double) totalTarget;
                    final double alignedPct = alignedBasesByGc[i] / (double) totalBases;

                    double dropout = (alignedPct - targetPct) * 100d;
                    if (dropout < 0) {
                        dropout = Math.abs(dropout);

                        if (i <=50) this.metrics.AT_DROPOUT += dropout;
                        if (i >=50) this.metrics.GC_DROPOUT += dropout;
                    }
                }
            }
        }


        @Override
        public void addMetricsToFile(MetricsFile<METRIC_TYPE, Integer> hsMetricsComparableMetricsFile) {
            hsMetricsComparableMetricsFile.addMetric(convertMetric(this.metrics));
        }
    }

    /**
     * A simple class that is used to store the coverage information about an interval.
     *
     * @author Tim Fennell
     */
    public static class Coverage {
        private final Interval interval;
        private final short[] depths;

        /** Constructs a new coverage object for the provided mapping with the desired padding either side. */
        public Coverage(final Interval i, final int padding) {
            this.interval = i;
            this.depths = new short[interval.length() + 2*padding];
        }

        /** Adds a single point of depth at the desired offset into the coverage array. */
        public void addBase(final int offset) {
            if (offset >= 0 && offset < this.depths.length) {
                // Prevent overflow if depth is too great, while avoiding doubling memory requirement.
                if (this.depths[offset] < Short.MAX_VALUE) {
                    this.depths[offset] += 1;
                }
            }
        }

        /** Returns true if any base in the range has coverage of > 1 */
        public boolean hasCoverage() {
            for (final short s : depths) {
                if (s > 1) return true;
            }

            return false;
        }

        /** Gets the coverage depths as an array of shorts. */
        public short[] getDepths() { return this.depths; }

        public int getTotal() {
            int total = 0;
            for (int i=0; i<depths.length; ++i) total += depths[i];
            return total;
        }

        @Override
        public String toString() {
            return "TargetedMetricCollector(interval=" + interval + ", depths = [" + StringUtil.intValuesToString(this.depths) + "])";
        }
    }
}
