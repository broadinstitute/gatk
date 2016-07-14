package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.PerUnitMetricCollector;
import org.broadinstitute.hellbender.metrics.SAMRecordMultiLevelCollector;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.gene.Gene;
import org.broadinstitute.hellbender.utils.gene.LocusFunction;

import java.io.File;
import java.util.*;

public final class RnaSeqMetricsCollector extends SAMRecordMultiLevelCollector<RnaSeqMetrics, Integer> {
    public enum StrandSpecificity {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}

    private final int minimumLength;
    private final StrandSpecificity strandSpecificity;
    private final double rrnaFragmentPercentage;
    private final Long ribosomalInitialValue;

    private final Set<Integer> ignoredSequenceIndices;

    private final OverlapDetector<Gene> geneOverlapDetector;
    private final OverlapDetector<Interval> ribosomalSequenceOverlapDetector;
    private final boolean collectCoverageStatistics;
    
    public RnaSeqMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords,
                                  final Long ribosomalBasesInitialValue, OverlapDetector<Gene> geneOverlapDetector, OverlapDetector<Interval> ribosomalSequenceOverlapDetector,
                                  final Set<Integer> ignoredSequenceIndices, final int minimumLength, final StrandSpecificity strandSpecificity,
                                  final double rrnaFragmentPercentage, boolean collectCoverageStatistics) {
        this.ribosomalInitialValue  = ribosomalBasesInitialValue;
        this.ignoredSequenceIndices = ignoredSequenceIndices;
        this.geneOverlapDetector    = geneOverlapDetector;
        this.ribosomalSequenceOverlapDetector = ribosomalSequenceOverlapDetector;
        this.minimumLength          = minimumLength;
        this.strandSpecificity      = strandSpecificity;
        this.rrnaFragmentPercentage = rrnaFragmentPercentage;
        this.collectCoverageStatistics = collectCoverageStatistics;
        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitRnaSeqMetricsCollector(sample, library, readGroup, ribosomalInitialValue);
    }

    public static OverlapDetector<Interval> makeOverlapDetector(final File samFile, final SAMFileHeader header, final File ribosomalIntervalsFile) {

        OverlapDetector<Interval> ribosomalSequenceOverlapDetector = new OverlapDetector<>(0, 0);
        if (ribosomalIntervalsFile != null) {

            final IntervalList ribosomalIntervals = IntervalList.fromFile(ribosomalIntervalsFile);
            try {
                SequenceUtil.assertSequenceDictionariesEqual(header.getSequenceDictionary(), ribosomalIntervals.getHeader().getSequenceDictionary());
            } catch (SequenceUtil.SequenceListsDifferException e) {
                throw new UserException("Sequence dictionaries differ in " + samFile.getAbsolutePath() + " and " + ribosomalIntervalsFile.getAbsolutePath(), e);
            }
            final IntervalList uniquedRibosomalIntervals = ribosomalIntervals.uniqued();
            final List<Interval> intervals = uniquedRibosomalIntervals.getIntervals();
            ribosomalSequenceOverlapDetector.addAll(intervals, intervals);
        }
        return ribosomalSequenceOverlapDetector;
    }

    public static Set<Integer> makeIgnoredSequenceIndicesSet(final SAMFileHeader header, final Set<String> ignoredSequence) {
        final Set<Integer> ignoredSequenceIndices = new HashSet<>();
        for (final String sequenceName: ignoredSequence) {
            final SAMSequenceRecord sequenceRecord = header.getSequence(sequenceName);
            if (sequenceRecord == null) {
                throw new UserException("Unrecognized sequence " + sequenceName + " passed as argument to IGNORE_SEQUENCE");
            }
            ignoredSequenceIndices.add(sequenceRecord.getSequenceIndex());
        }
        return ignoredSequenceIndices;
    }

    private class PerUnitRnaSeqMetricsCollector implements PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> {

        final RnaSeqMetrics metrics = new RnaSeqMetrics();
        
        private final Map<Gene.Transcript, int[]> coverageByTranscript = new HashMap<>();

        public PerUnitRnaSeqMetricsCollector(final String sample,
                                             final String library,
                                             final String readGroup,
                                             final Long ribosomalBasesInitialValue) {
            this.metrics.SAMPLE = sample;
            this.metrics.LIBRARY = library;
            this.metrics.READ_GROUP = readGroup;
            this.metrics.RIBOSOMAL_BASES = ribosomalBasesInitialValue;
            
        }

        public void acceptRecord(SAMRecord rec) {
            // Filter out some reads, and collect the total number of PF bases
            if (rec.getReadFailsVendorQualityCheckFlag() || rec.isSecondaryOrSupplementary()) return;

            this.metrics.PF_BASES += rec.getReadLength();
            if (rec.getReadUnmappedFlag()) return;

            if (ignoredSequenceIndices.contains(rec.getReferenceIndex())) {
                ++this.metrics.IGNORED_READS;
                return;
            }

            // Grab information about the alignment and overlapping genes etc.
            final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());

            // Attempt to get an interval for the entire fragment (if paired read) else just use the read itself.
            // If paired read is chimeric or has one end unmapped, don't create an interval.
            final Interval fragmentInterval;
            if (!rec.getReadPairedFlag()) {
                fragmentInterval = readInterval;
            } else if (rec.getMateUnmappedFlag() || rec.getReferenceIndex() != rec.getMateReferenceIndex()) {
                fragmentInterval = null;
            } else {
                final int fragmentStart = Math.min(rec.getAlignmentStart(), rec.getMateAlignmentStart());
                final int fragmentEnd = CoordMath.getEnd(fragmentStart, Math.abs(rec.getInferredInsertSize()));
                fragmentInterval = new Interval(rec.getReferenceName(), fragmentStart, fragmentEnd);
            }
            if (fragmentInterval != null) {
                final Collection<Interval> overlappingRibosomalIntervals = ribosomalSequenceOverlapDetector.getOverlaps(fragmentInterval);
                int intersectionLength = 0;
                for (final Interval overlappingInterval : overlappingRibosomalIntervals) {
                    final int thisIntersectionLength = overlappingInterval.getIntersectionLength(fragmentInterval);
                    intersectionLength = Math.max(intersectionLength, thisIntersectionLength);
                }
                if (intersectionLength/(double)fragmentInterval.length() >= rrnaFragmentPercentage) {
                    // Assume entire read is ribosomal.
                    // TODO: Should count reads, not bases?
                    metrics.RIBOSOMAL_BASES += rec.getReadLength();
                    int numAlignedBases = 0;
                    for (final AlignmentBlock alignmentBlock : rec.getAlignmentBlocks()) {
                        numAlignedBases += alignmentBlock.getLength();
                    }
                    metrics.PF_ALIGNED_BASES += numAlignedBases;
                    return;
                }
            }

            final Collection<Gene> overlappingGenes                  = geneOverlapDetector.getOverlaps(readInterval);
            final List<AlignmentBlock> alignmentBlocks               = rec.getAlignmentBlocks();
            boolean overlapsExon = false;

            for (final AlignmentBlock alignmentBlock : alignmentBlocks) {
                // Get functional class for each position in the alignment block.
                final LocusFunction[] locusFunctions = new LocusFunction[alignmentBlock.getLength()];

                // By default, if base does not overlap with rRNA or gene, it is intergenic.
                Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

                for (final Gene gene : overlappingGenes) {
                    for (final Gene.Transcript transcript : gene) {
                        transcript.assignLocusFunctionForRange(alignmentBlock.getReferenceStart(), locusFunctions);
                        // if you want to gather coverage statistics, this variable should be true.
                        // added for cases with many units [samples/read groups] which overwhelm memory.                
                        // Add coverage to our coverage counter for this transcript
                        if (collectCoverageStatistics) {
                            int[] coverage = this.coverageByTranscript.get(transcript);
                            if (coverage == null) {
                                coverage = new int[transcript.length()];
                                this.coverageByTranscript.put(transcript, coverage);
                            }
                            transcript.addCoverageCounts(alignmentBlock.getReferenceStart(),
                                    CoordMath.getEnd(alignmentBlock.getReferenceStart(), alignmentBlock.getLength()),
                                    coverage);
                        }

                    }
                }

                // Tally the function of each base in the alignment block.
                for (final LocusFunction locusFunction : locusFunctions) {
                    ++metrics.PF_ALIGNED_BASES;
                    switch (locusFunction) {
                        case INTERGENIC:
                            ++metrics.INTERGENIC_BASES;
                            break;
                        case INTRONIC:
                            ++metrics.INTRONIC_BASES;
                            break;
                        case UTR:
                            ++metrics.UTR_BASES;
                            overlapsExon = true;
                            break;
                        case CODING:
                            ++metrics.CODING_BASES;
                            overlapsExon = true;
                            break;
                        case RIBOSOMAL:
                            ++metrics.RIBOSOMAL_BASES;
                            break;
                    }
                }
            }

            // Strand-specificity is tallied on read basis rather than base at a time.  A read that aligns to more than one
            // gene is not counted.
            if (overlapsExon && strandSpecificity != StrandSpecificity.NONE && overlappingGenes.size() == 1) {
                final boolean negativeTranscriptionStrand = overlappingGenes.iterator().next().isNegativeStrand();
                final boolean negativeReadStrand = rec.getReadNegativeStrandFlag();
                final boolean readAndTranscriptStrandsAgree = negativeReadStrand == negativeTranscriptionStrand;
                final boolean readOneOrUnpaired = !rec.getReadPairedFlag() || rec.getFirstOfPairFlag();
                final boolean firstReadExpectedToAgree = strandSpecificity == StrandSpecificity.FIRST_READ_TRANSCRIPTION_STRAND;
                final boolean thisReadExpectedToAgree = readOneOrUnpaired == firstReadExpectedToAgree;
                // If the read strand is the same as the strand of the transcript, and the end is the one that is supposed to agree,
                // then the strand specificity for this read is correct.
                // -- OR --
                // If the read strand is not the same as the strand of the transcript, and the end is not the one that is supposed
                // to agree, then the strand specificity for this read is correct.
                if (readAndTranscriptStrandsAgree == thisReadExpectedToAgree) {
                    ++metrics.CORRECT_STRAND_READS;
                } else {
                    ++metrics.INCORRECT_STRAND_READS;
                }
            }

        }

        public void finish() {
            if (metrics.PF_ALIGNED_BASES > 0) {
                if (metrics.RIBOSOMAL_BASES != null) {
                    metrics.PCT_RIBOSOMAL_BASES =  metrics.RIBOSOMAL_BASES  / (double) metrics.PF_ALIGNED_BASES;
                }
                metrics.PCT_CODING_BASES =     metrics.CODING_BASES     / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_UTR_BASES =        metrics.UTR_BASES        / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_INTRONIC_BASES =   metrics.INTRONIC_BASES   / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_INTERGENIC_BASES = metrics.INTERGENIC_BASES / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_MRNA_BASES =       metrics.PCT_CODING_BASES + metrics.PCT_UTR_BASES;
                metrics.PCT_USABLE_BASES =     (metrics.CODING_BASES + metrics.UTR_BASES) / (double) metrics.PF_BASES;
            }

            if (metrics.CORRECT_STRAND_READS > 0 || metrics.INCORRECT_STRAND_READS > 0) {
                metrics.PCT_CORRECT_STRAND_READS = metrics.CORRECT_STRAND_READS/(double)(metrics.CORRECT_STRAND_READS + metrics.INCORRECT_STRAND_READS);
            }
        }

        @Override
        public void addMetricsToFile(final MetricsFile<RnaSeqMetrics, Integer> file) {
            // Compute metrics based on coverage of top 1000 genes
            final Histogram<Integer> normalizedCovByPos = computeCoverageMetrics();
            file.addMetric(metrics);
            file.addHistogram(normalizedCovByPos);
        }

        /**
         * Computes a set of coverage based metrics on the mostly highly expressed genes' most highly
         * expressed transcripts.
         */
        private Histogram<Integer> computeCoverageMetrics() {
            final Histogram<Double> cvs = new Histogram<>();
            final Histogram<Double> fivePrimeSkews = new Histogram<>();
            final Histogram<Double> threePrimeSkews = new Histogram<>();
            final Histogram<Double> gapBasesPerKb = new Histogram<>();
            final Histogram<Double> fiveToThreeSkews = new Histogram<>();
            String prefix = null;
            if (this.metrics.READ_GROUP != null) {
                prefix = this.metrics.READ_GROUP + ".";
            }
            else if (this.metrics.LIBRARY != null) {
                prefix = this.metrics.LIBRARY + ".";
            }
            else if (this.metrics.SAMPLE != null) {
                prefix = this.metrics.SAMPLE + ".";
            }
            else {
                prefix = "All_Reads.";
            }

            final Histogram<Integer> normalizedCoverageByNormalizedPosition = new Histogram<>("normalized_position", prefix + "normalized_coverage");

            final Map<Gene.Transcript,int[]> transcripts = pickTranscripts(coverageByTranscript);
            final double transcriptCount = transcripts.size();

            for (final Map.Entry<Gene.Transcript,int[]> entry : transcripts.entrySet()) {
                final Gene.Transcript tx = entry.getKey();
                final double[] coverage;
                {
                    final double[] tmp = MathUtils.promote(entry.getValue());
                    if (tx.getGene().isPositiveStrand())  coverage = tmp;
                    else coverage = copyAndReverse(tmp);
                }
                final double mean = MathUtils.mean(coverage, 0, coverage.length);

                // Calculate the CV of coverage for this tx
                final double stdev = new StandardDeviation().evaluate(coverage, mean);
                final double cv    = stdev / mean;
                cvs.increment(cv);

                // Calculate the 5' and 3' biases
                {
                    final int PRIME_BASES = 100;
                    final double fivePrimeCoverage = MathUtils.mean(coverage, 0, PRIME_BASES);
                    final double threePrimeCoverage = MathUtils.mean(coverage, coverage.length - PRIME_BASES, coverage.length);

                    fivePrimeSkews.increment(fivePrimeCoverage / mean);
                    threePrimeSkews.increment(threePrimeCoverage / mean);
                    fiveToThreeSkews.increment(fivePrimeCoverage / threePrimeCoverage);
                }

                // Calculate normalized coverage vs. normalized position
                {
                    final int lastIndex = coverage.length - 1;

                    for (int percent=0; percent<=100; ++percent) {
                        final double p = percent / 100d;
                        final int start  = (int) Math.max(0, lastIndex * (p - 0.005));
                        final int end    = (int) Math.min(lastIndex, lastIndex * (p + 0.005));
                        final int length = end - start + 1;

                        double sum = 0;
                        for (int i=start; i<=end; ++i) sum += coverage[i];
                        final double normalized = (sum / length) / mean;
                        normalizedCoverageByNormalizedPosition.increment(percent, normalized / transcriptCount);
                    }
                }

                // Calculate gap bases per kilobase
                //            {
                //                int gapBases = 0;
                //                final double minCoverage = mean * 0.1;
                //                for (int i=0; i<coverage.length; ++i) {
                //                    if (coverage[i] < minCoverage) ++gapBases;
                //                }
                //                gapBasesPerKb.increment(gapBases / (coverage.length / 1000d));
                //            }
            }

            this.metrics.MEDIAN_CV_COVERAGE = cvs.getMedian();
            this.metrics.MEDIAN_5PRIME_BIAS = fivePrimeSkews.getMedian();
            this.metrics.MEDIAN_3PRIME_BIAS = threePrimeSkews.getMedian();
            this.metrics.MEDIAN_5PRIME_TO_3PRIME_BIAS = fiveToThreeSkews.getMedian();

            return normalizedCoverageByNormalizedPosition;
        }

        /** Little method to copy an array and reverse it at the same time. */
        private double[] copyAndReverse(final double[] in) {
            final double[] out = new double[in.length];
            for (int i=0, j=in.length-1; i<in.length; ++i, --j) out[j] = in[i];
            return out;
        }

        /** Picks the set of transcripts on which the coverage metrics are to be calculated. */
        public Map<Gene.Transcript, int[]> pickTranscripts(final Map<Gene.Transcript, int[]> transcriptCoverage) {
            final Map<Gene.Transcript, Double> bestPerGene = new HashMap<>();

            // Make a map of the best transcript per gene to it's mean coverage
            for (final Gene gene : geneOverlapDetector.getAll()) {
                Gene.Transcript best = null;
                double bestMean = 0;

                for (final Gene.Transcript tx : gene) {
                    final int[] cov = transcriptCoverage.get(tx);

                    if (tx.length() < Math.max(minimumLength, 100)) continue;
                    if (cov == null) continue;

                    final double mean = MathUtils.mean(MathUtils.promote(cov), 0, cov.length);
                    if (mean < 1d) continue;
                    if (best == null || mean > bestMean) {
                        best = tx;
                        bestMean = mean;
                    }
                }

                if (best != null) bestPerGene.put(best, bestMean);
            }

            // Find the 1000th best coverage value
            final double[] coverages = new double[bestPerGene.size()];
            int i=0;
            for (final double d : bestPerGene.values()) coverages[i++] = d;
            Arrays.sort(coverages);
            final double min = coverages.length == 0 ? 0 : coverages[Math.max(0, coverages.length - 1001)];

            // And finally build the output map
            final Map<Gene.Transcript, int[]> retval = new HashMap<>();
            for (final Map.Entry<Gene.Transcript,Double> entry : bestPerGene.entrySet()) {
                final Gene.Transcript tx = entry.getKey();
                final double coverage = entry.getValue();

                if (coverage >= min) {
                    retval.put(tx, transcriptCoverage.get(tx));
                }
            }

            return retval;
        }

    }
}
