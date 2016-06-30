package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.PerUnitMetricCollector;
import org.broadinstitute.hellbender.metrics.SAMRecordAndReference;
import org.broadinstitute.hellbender.metrics.SAMRecordAndReferenceMultiLevelCollector;

import java.util.*;

public final class AlignmentSummaryMetricsCollector extends SAMRecordAndReferenceMultiLevelCollector<AlignmentSummaryMetrics,  Long> {
    // If we have a reference sequence, collect metrics on how well we aligned to it
    private final boolean doRefMetrics;

    //the adapter sequences converted to byte arrays
    private final byte[][] adapterKmers;

    //A list of Strings representing the sequence of bases in an adapter
    private final List<String> adapterSequence;

    //Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.
    private final int maxInsertSize;


    //Whether the SAM/BAM file consists of bisulfite sequenced reads.
    private final boolean isBisulfiteSequenced;

    //The minimum mapping quality a base has to meet in order to be considered high quality
    private final int MAPPING_QUALITY_THRESOLD = 20;

    //The minimum quality a base has to meet in order to be consider hq_20
    private static final int BASE_QUALITY_THRESHOLD = 20;

    //The number of bases to check in order to map a read to an adapter
    private static final int ADAPTER_MATCH_LENGTH = 16;

    // The maximum number of mismatches a read can have and still be considered as matching an adapter
    private static final int MAX_ADAPTER_ERRORS = 1;

    public AlignmentSummaryMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords,
                                            final boolean doRefMetrics, final List<String> adapterSequence, final int maxInsertSize, boolean isBisulfiteSequenced) {
        this.doRefMetrics       = doRefMetrics;
        this.adapterSequence    = adapterSequence;
        this.adapterKmers = prepareAdapterSequences();
        this.maxInsertSize = maxInsertSize;
        this.isBisulfiteSequenced = isBisulfiteSequenced;
        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<AlignmentSummaryMetrics, Long, SAMRecordAndReference> makeChildCollector(String sample, String library, String readGroup) {
        return new GroupAlignmentSummaryMetricsPerUnitMetricCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {
        if (!rec.isSecondaryOrSupplementary()) {
            super.acceptRecord(rec, ref);
        }
    }

    /** Converts the supplied adapter sequences to byte arrays in both fwd and rc. */
    private byte [][] prepareAdapterSequences() {
        final Set<String> kmers = new HashSet<>();

        // Make a set of all kmers of adapterMatchLength
        for (final String seq : adapterSequence) {
            for (int i=0; i<=seq.length() - ADAPTER_MATCH_LENGTH; ++i) {
                final String kmer = seq.substring(i, i+ADAPTER_MATCH_LENGTH).toUpperCase();

                int ns = 0;
                for (final char ch : kmer.toCharArray()) if (ch == 'N') ++ns;
                if (ns <= MAX_ADAPTER_ERRORS) {
                    kmers.add(kmer);
                    kmers.add(SequenceUtil.reverseComplement(kmer));
                }
            }
        }

        // Make an array of byte[] for the kmers
        final byte [][] adapterKmers = new byte[kmers.size()][];
        int i=0;
        for (final String kmer : kmers) {
            adapterKmers[i++] = StringUtil.stringToBytes(kmer);
        }
        return adapterKmers;
    }

    /**
     * Checks the first ADAPTER_MATCH_LENGTH bases of the read against known adapter sequences and returns
     * true if the read matches an adapter sequence with MAX_ADAPTER_ERRORS mismsatches or fewer.
     *
     * @param read the basecalls for the read in the order and orientation the machine read them
     * @return true if the read matches an adapter and false otherwise
     */
    private boolean isAdapterSequence(final byte[] read) {
        if (read.length < ADAPTER_MATCH_LENGTH) return false;

        for (final byte[] adapter : adapterKmers) {
            int errors = 0;

            for (int i=0; i<adapter.length; ++i) {
                if (read[i] != adapter[i]) {
                    if (++errors > MAX_ADAPTER_ERRORS) break;
                }
            }

            if (errors <= MAX_ADAPTER_ERRORS) return true;
        }

        return false;
    }

    private class GroupAlignmentSummaryMetricsPerUnitMetricCollector implements PerUnitMetricCollector<AlignmentSummaryMetrics, Long, SAMRecordAndReference> {
        final IndividualAlignmentSummaryMetricsCollector unpairedCollector;
        final IndividualAlignmentSummaryMetricsCollector firstOfPairCollector;
        final IndividualAlignmentSummaryMetricsCollector secondOfPairCollector;
        final IndividualAlignmentSummaryMetricsCollector pairCollector;
        final String sample;
        final String library;
        final String readGroup;

        public GroupAlignmentSummaryMetricsPerUnitMetricCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            unpairedCollector     = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.UNPAIRED, sample, library, readGroup);
            firstOfPairCollector  = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.FIRST_OF_PAIR, sample, library, readGroup);
            secondOfPairCollector = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.SECOND_OF_PAIR, sample, library, readGroup);
            pairCollector         = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.PAIR, sample, library, readGroup);
        }

        public void acceptRecord(final SAMRecordAndReference args) {
            final SAMRecord rec         = args.getSamRecord();
            final ReferenceSequence ref = args.getReferenceSequence();

            if (rec.getReadPairedFlag()) {
                if (rec.getFirstOfPairFlag()) {
                    firstOfPairCollector.addRecord(rec, ref);
                }
                else {
                    secondOfPairCollector.addRecord(rec, ref);
                }

                pairCollector.addRecord(rec, ref);
            }
            else {
                unpairedCollector.addRecord(rec, ref);
            }
        }

        @Override
        public void finish() {
            // Let the collectors do any summary computations etc.
            unpairedCollector.onComplete();
            firstOfPairCollector.onComplete();
            secondOfPairCollector.onComplete();
            pairCollector.onComplete();
        }

        @Override
        public void addMetricsToFile(final MetricsFile<AlignmentSummaryMetrics, Long> file) {
            if (firstOfPairCollector.getMetrics().TOTAL_READS > 0) {
                // override how bad cycle is determined for paired reads, it should be
                // the sum of first and second reads
                pairCollector.getMetrics().BAD_CYCLES = firstOfPairCollector.getMetrics().BAD_CYCLES +
                        secondOfPairCollector.getMetrics().BAD_CYCLES;

                file.addMetric(firstOfPairCollector.getMetrics());
                file.addMetric(secondOfPairCollector.getMetrics());
                file.addMetric(pairCollector.getMetrics());
            }

            //if there are no reads in any category then we will returned an unpaired alignment summary metric with all zero values
            if (unpairedCollector.getMetrics().TOTAL_READS > 0 || firstOfPairCollector.getMetrics().TOTAL_READS == 0) {
                file.addMetric(unpairedCollector.getMetrics());
            }
        }

        /**
         * Class that counts reads that match various conditions
         */
        private class IndividualAlignmentSummaryMetricsCollector {
            private long numPositiveStrand = 0;
            private final Histogram<Integer> readLengthHistogram = new Histogram<>();
            private final AlignmentSummaryMetrics metrics;
            private long chimeras;
            private long chimerasDenominator;
            private long adapterReads;
            private long indels;

            private long nonBisulfiteAlignedBases = 0;
            private long hqNonBisulfiteAlignedBases = 0;
            private final Histogram<Long> mismatchHistogram = new Histogram<>();
            private final Histogram<Long> hqMismatchHistogram = new Histogram<>();
            private final Histogram<Integer> badCycleHistogram = new Histogram<>();

            public IndividualAlignmentSummaryMetricsCollector(final AlignmentSummaryMetrics.Category pairingCategory,
                                                              final String sample,
                                                              final String library,
                                                              final String readGroup) {
                metrics = new AlignmentSummaryMetrics();
                metrics.CATEGORY = pairingCategory;
                metrics.SAMPLE = sample;
                metrics.LIBRARY = library;
                metrics.READ_GROUP = readGroup;
            }

            public void addRecord(final SAMRecord record, final ReferenceSequence ref) {
                if (record.isSecondaryOrSupplementary()) {
                    // only want 1 count per read so skip non primary alignments
                    return;
                }

                collectReadData(record, ref);
                collectQualityData(record, ref);
            }
            @SuppressWarnings("unchecked")
            public void onComplete() {
                //summarize read data
                if (metrics.TOTAL_READS > 0)
                {
                    metrics.PCT_PF_READS = (double) metrics.PF_READS / (double) metrics.TOTAL_READS;
                    metrics.PCT_ADAPTER = this.adapterReads / (double) metrics.PF_READS;
                    metrics.MEAN_READ_LENGTH = readLengthHistogram.getMean();

                    //Calculate BAD_CYCLES
                    metrics.BAD_CYCLES = 0;

                    for (final Histogram.Bin<Integer> cycleBin : badCycleHistogram.values()) {
                        final double badCyclePercentage = cycleBin.getValue() / metrics.TOTAL_READS;
                        if (badCyclePercentage >= .8) {
                            metrics.BAD_CYCLES++;
                        }
                    }

                    if(doRefMetrics) {
                        if (metrics.PF_READS > 0)         metrics.PCT_PF_READS_ALIGNED = (double) metrics.PF_READS_ALIGNED / (double) metrics.PF_READS;
                        if (metrics.PF_READS_ALIGNED > 0) metrics.PCT_READS_ALIGNED_IN_PAIRS = (double) metrics.READS_ALIGNED_IN_PAIRS/ (double) metrics.PF_READS_ALIGNED;
                        if (metrics.PF_READS_ALIGNED > 0) metrics.STRAND_BALANCE = numPositiveStrand / (double) metrics.PF_READS_ALIGNED;
                        if (this.chimerasDenominator > 0) metrics.PCT_CHIMERAS = this.chimeras / (double) this.chimerasDenominator;

                        if (nonBisulfiteAlignedBases > 0) metrics.PF_MISMATCH_RATE = mismatchHistogram.getSum() / (double) nonBisulfiteAlignedBases;
                        metrics.PF_HQ_MEDIAN_MISMATCHES = hqMismatchHistogram.getMedian();
                        if (hqNonBisulfiteAlignedBases > 0) metrics.PF_HQ_ERROR_RATE = hqMismatchHistogram.getSum() / (double) hqNonBisulfiteAlignedBases;
                        if (metrics.PF_ALIGNED_BASES > 0) metrics.PF_INDEL_RATE = this.indels / (double) metrics.PF_ALIGNED_BASES;
                    }
                }
            }

            private void collectReadData(final SAMRecord record, final ReferenceSequence ref) {
                metrics.TOTAL_READS++;
                readLengthHistogram.increment(record.getReadBases().length);

                if (!record.getReadFailsVendorQualityCheckFlag()) {
                    metrics.PF_READS++;
                    if (isNoiseRead(record)) metrics.PF_NOISE_READS++;

                    if (record.getReadUnmappedFlag()) {
                        // If the read is unmapped see if it's adapter sequence
                        final byte[] readBases = record.getReadBases();
                        if (!(record instanceof BAMRecord)) StringUtil.toUpperCase(readBases);

                        if (isAdapterSequence(readBases)) {
                            this.adapterReads++;
                        }
                    }
                    else if(doRefMetrics) {
                        metrics.PF_READS_ALIGNED++;
                        if (!record.getReadNegativeStrandFlag()) numPositiveStrand++;

                        if (record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
                            metrics.READS_ALIGNED_IN_PAIRS++;

                            // Check that both ends have mapq > minimum
                            final Integer mateMq = record.getIntegerAttribute("MQ");
                            if (mateMq == null || mateMq >= MAPPING_QUALITY_THRESOLD && record.getMappingQuality() >= MAPPING_QUALITY_THRESOLD) {
                                ++this.chimerasDenominator;

                                // With both reads mapped we can see if this pair is chimeric
                                if (Math.abs(record.getInferredInsertSize()) > maxInsertSize ||
                                        !record.getReferenceIndex().equals(record.getMateReferenceIndex())) {
                                    ++this.chimeras;
                                }
                            }
                        }
                    }
                }
            }

            private void collectQualityData(final SAMRecord record, final ReferenceSequence reference) {
                // If the read isnt an aligned PF read then look at the read for no-calls
                if (record.getReadUnmappedFlag() || record.getReadFailsVendorQualityCheckFlag() || !doRefMetrics) {
                    final byte[] readBases = record.getReadBases();
                    for (int i = 0; i < readBases.length; i++) {
                        if (SequenceUtil.isNoCall(readBases[i])) {
                            badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                        }
                    }
                }
                else if (!record.getReadFailsVendorQualityCheckFlag()) {
                    final boolean highQualityMapping = isHighQualityMapping(record);
                    if (highQualityMapping) metrics.PF_HQ_ALIGNED_READS++;

                    final byte[] readBases = record.getReadBases();
                    final byte[] refBases = reference.getBases();
                    final byte[] qualities  = record.getBaseQualities();
                    final int refLength = refBases.length;
                    long mismatchCount   = 0;
                    long hqMismatchCount = 0;

                    for (final AlignmentBlock alignmentBlock : record.getAlignmentBlocks()) {
                        final int readIndex = alignmentBlock.getReadStart() - 1;
                        final int refIndex  = alignmentBlock.getReferenceStart() - 1;
                        final int length    = alignmentBlock.getLength();

                        for (int i=0; i<length && refIndex+i<refLength; ++i) {
                            final int readBaseIndex = readIndex + i;
                            boolean mismatch = !SequenceUtil.basesEqual(readBases[readBaseIndex], refBases[refIndex + i]);
                            boolean bisulfiteBase = false;
                            if (mismatch && isBisulfiteSequenced) {
                                if ( (record.getReadNegativeStrandFlag() &&
                                        (refBases[refIndex+i] == 'G' || refBases[refIndex+i] =='g') &&
                                        (readBases[readBaseIndex] == 'A' || readBases[readBaseIndex] == 'a'))
                                        || ((!record.getReadNegativeStrandFlag()) &&
                                        (refBases[refIndex+i] == 'C' || refBases[refIndex+i] == 'c') &&
                                        (readBases[readBaseIndex] == 'T') || readBases[readBaseIndex] == 't') ) {

                                    bisulfiteBase = true;
                                    mismatch = false;
                                }
                            }

                            if(mismatch) mismatchCount++;

                            metrics.PF_ALIGNED_BASES++;
                            if(!bisulfiteBase) nonBisulfiteAlignedBases++;

                            if (highQualityMapping) {
                                metrics.PF_HQ_ALIGNED_BASES++;
                                if (!bisulfiteBase) hqNonBisulfiteAlignedBases++;
                                if (qualities[readBaseIndex] >= BASE_QUALITY_THRESHOLD) metrics.PF_HQ_ALIGNED_Q20_BASES++;
                                if (mismatch) hqMismatchCount++;
                            }

                            if (mismatch || SequenceUtil.isNoCall(readBases[readBaseIndex])) {
                                badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                            }
                        }
                    }

                    mismatchHistogram.increment(mismatchCount);
                    hqMismatchHistogram.increment(hqMismatchCount);

                    // Add any insertions and/or deletions to the global count
                    for (final CigarElement elem : record.getCigar().getCigarElements()) {
                        final CigarOperator op = elem.getOperator();
                        if (op == CigarOperator.INSERTION || op == CigarOperator.DELETION) ++ this.indels;
                    }
                }
            }

            private boolean isNoiseRead(final SAMRecord record) {
                final Object noiseAttribute = record.getAttribute(ReservedTagConstants.XN);
                return (noiseAttribute != null && noiseAttribute.equals(1));
            }

            private boolean isHighQualityMapping(final SAMRecord record) {
                return !record.getReadFailsVendorQualityCheckFlag() &&
                        record.getMappingQuality() >= MAPPING_QUALITY_THRESOLD;
            }

            public AlignmentSummaryMetrics getMetrics() {
                return this.metrics;
            }
        }
    }
}
