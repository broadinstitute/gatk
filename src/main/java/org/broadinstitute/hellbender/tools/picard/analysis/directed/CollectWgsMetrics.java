package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 *
 * @author tfennell
 */
@CommandLineProgramProperties(
        summary = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "whole genome sequencing experiments.",
        oneLineSummary = "Produces metrics related to whole genome sequencing for a SAM/BAM file",
        programGroup = QCProgramGroup.class
)
public final class CollectWgsMetrics extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input SAM/BAM file.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output metrics file.")
    public File OUTPUT;

    @Argument(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.")
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Argument(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.")
    public int MINIMUM_BASE_QUALITY = 20;

    @Argument(shortName = "CAP", doc = "Treat bases with coverage exceeding this value as if they had coverage at this value.")
    public int COVERAGE_CAP = 250;

    @Argument(doc = "For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;

    @Argument(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetrics extends MetricBase {
        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
        public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
        public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
        public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
        public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
        public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
        public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
        public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
        public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
        public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
        public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
        public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
        public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
        public double PCT_100X;
    }

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        // Setup all the inputs
        final ProgressLogger progress = new ProgressLogger(logger, 10000000, "Processed", "loci");
        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).validationStringency(VALIDATION_STRINGENCY).open(INPUT);

        // Load up the reference sequence and double check sequence dictionaries
        if (!in.getFileHeader().getSequenceDictionary().isEmpty()) {
            SequenceUtil.assertSequenceDictionariesEqual(in.getFileHeader().getSequenceDictionary(), refWalker.getSequenceDictionary());
        }

        final SamLocusIterator iterator = new SamLocusIterator(in);
        final List<SamRecordFilter> filters = new ArrayList<>();
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingFilter mapqFilter = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
        final CountingPairedFilter pairFilter = new CountingPairedFilter();
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        filters.add(pairFilter);
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        iterator.setSamFilters(filters);
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setQualityScoreCutoff(0);        // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);

        final int max = COVERAGE_CAP;
        final long[] HistogramArray = new long[max + 1];
        final long[] baseQHistogramArray = new long[Byte.MAX_VALUE];
        final boolean usingStopAfter = STOP_AFTER > 0;
        final long stopAfter = STOP_AFTER - 1;
        long counter = 0;

        long basesExcludedByBaseq = 0;
        long basesExcludedByOverlap = 0;
        long basesExcludedByCapping = 0;

        // Loop through all the loci
        while (iterator.hasNext()) {
            final SamLocusIterator.LocusInfo info = iterator.next();

            // Check that the reference is not N
            final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (base == 'N') continue;

            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<>(info.getRecordAndPositions().size());
            int pileupSize = 0;
            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {

                if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY)                   { ++basesExcludedByBaseq;   continue; }
                if (!readNames.add(recs.getRecord().getReadName()))                 { ++basesExcludedByOverlap; continue; }
                pileupSize++;
                if (pileupSize <= max) {
                    baseQHistogramArray[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                }
            }

            final int depth = Math.min(readNames.size(), max);
            if (depth < readNames.size()) basesExcludedByCapping += readNames.size() - max;
            HistogramArray[depth]++;

            // Record progress and perhaps stop
            progress.record(info.getSequenceName(), info.getPosition());
            if (usingStopAfter && ++counter > stopAfter) break;
        }

        // Construct and write the outputs
        final Histogram<Integer> histo = new Histogram<>("coverage", "count");
        for (int i = 0; i < HistogramArray.length; ++i) {
            histo.increment(i, HistogramArray[i]);
        }

        // Construct and write the outputs
        final Histogram<Integer> baseQHisto = new Histogram<>("value", "baseq_count");
        for (int i=0; i<baseQHistogramArray.length; ++i) {
            baseQHisto.increment(i, baseQHistogramArray[i]);
        }

        final WgsMetrics metrics = generateWgsMetrics();
        metrics.GENOME_TERRITORY = (long) histo.getSumOfValues();
        metrics.MEAN_COVERAGE = histo.getMean();
        metrics.SD_COVERAGE = histo.getStandardDeviation();
        metrics.MEDIAN_COVERAGE = histo.getMedian();
        metrics.MAD_COVERAGE = histo.getMedianAbsoluteDeviation();

        final long basesExcludedByDupes = dupeFilter.getFilteredBases();
        final long basesExcludedByMapq = mapqFilter.getFilteredBases();
        final long basesExcludedByPairing = pairFilter.getFilteredBases();
        final double total = histo.getSum();
        final double totalWithExcludes = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;
        metrics.PCT_EXC_DUPE = basesExcludedByDupes / totalWithExcludes;
        metrics.PCT_EXC_MAPQ = basesExcludedByMapq / totalWithExcludes;
        metrics.PCT_EXC_UNPAIRED = basesExcludedByPairing / totalWithExcludes;
        metrics.PCT_EXC_BASEQ = basesExcludedByBaseq / totalWithExcludes;
        metrics.PCT_EXC_OVERLAP = basesExcludedByOverlap / totalWithExcludes;
        metrics.PCT_EXC_CAPPED = basesExcludedByCapping / totalWithExcludes;
        metrics.PCT_EXC_TOTAL = (totalWithExcludes - total) / totalWithExcludes;

        metrics.PCT_5X = MathUtils.sum(HistogramArray, 5, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_10X = MathUtils.sum(HistogramArray, 10, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_15X = MathUtils.sum(HistogramArray, 15, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_20X = MathUtils.sum(HistogramArray, 20, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_25X = MathUtils.sum(HistogramArray, 25, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_30X = MathUtils.sum(HistogramArray, 30, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_40X = MathUtils.sum(HistogramArray, 40, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_50X = MathUtils.sum(HistogramArray, 50, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_60X = MathUtils.sum(HistogramArray, 60, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_70X = MathUtils.sum(HistogramArray, 70, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_80X = MathUtils.sum(HistogramArray, 80, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_90X = MathUtils.sum(HistogramArray, 90, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_100X = MathUtils.sum(HistogramArray, 100, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;

        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        out.addMetric(metrics);
        out.addHistogram(histo);
        if (INCLUDE_BQ_HISTOGRAM) {
            out.addHistogram(baseQHisto);
        }
        out.write(OUTPUT);

        return null;
    }

    protected WgsMetrics generateWgsMetrics() {
        return new WgsMetrics();
    }
}

/**
 * A SamRecordFilter that counts the number of aligned bases in the reads which it filters out. Abstract and designed
 * to be subclassed to implement the desired filter.
 */
abstract class CountingFilter implements SamRecordFilter {
    private long filteredRecords = 0;
    private long filteredBases = 0;

    /** Gets the number of records that have been filtered out thus far. */
    public long getFilteredRecords() { return this.filteredRecords; }

    /** Gets the number of bases that have been filtered out thus far. */
    public long getFilteredBases() { return this.filteredBases; }

    @Override
    public final boolean filterOut(final SAMRecord record) {
        final boolean filteredOut = reallyFilterOut(record);
        if (filteredOut) {
            ++filteredRecords;
            for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                this.filteredBases += block.getLength();
            }
        }
        return filteredOut;
    }

    public abstract boolean reallyFilterOut(final SAMRecord record);

    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new UnsupportedOperationException();
    }
}

/** Counting filter that discards reads that have been marked as duplicates. */
final class CountingDuplicateFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getDuplicateReadFlag(); }
}

/** Counting filter that discards reads below a configurable mapping quality threshold. */
final class CountingMapQFilter extends CountingFilter {
    private final int minMapq;

    CountingMapQFilter(final int minMapq) { this.minMapq = minMapq; }

    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getMappingQuality() < minMapq; }
}

/** Counting filter that discards reads that are unpaired in sequencing and paired reads who's mates are not mapped. */
final class CountingPairedFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return !record.getReadPairedFlag() || record.getMateUnmappedFlag(); }
}

