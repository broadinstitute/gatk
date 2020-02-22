package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.IntervalPileup;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class DragstrCasesSampler {

    private SAMSequenceDictionary dictionary;
    private ReadsDataSource readsSource;
    private ReferenceDataSource referenceSource;
    private static final Logger logger = LogManager.getLogger(EstimateDragstrModelParameters.class);
    private DragstrCasesSamplerArgumentCollection config;

    public DragstrCasesSampler(final DragstrCasesSamplerArgumentCollection dragstrCasesSamplerArgumentCollection,
                               final ReferenceDataSource referenceSource,
                               final ReadsDataSource readsSource) {
        this.config = dragstrCasesSamplerArgumentCollection;
        this.readsSource = readsSource;
        this.referenceSource = referenceSource;
        this.dictionary = referenceSource.getSequenceDictionary();
    }

    void sample(final DragstrModelEstimator.RepeatCases dest, final List<DragstrLocus> loci) {
        final int period = dest.getPeriod();
        final int repeats = dest.getRepeats();
        logger.info("Sampling period = " + period + " and repeat count = " + repeats);
        final Random rdn = new Random(((config.randomSeed * 31) + period * 31) + repeats * 31);
        Collections.shuffle(loci, rdn);
        final int size = loci.size();
        int nonAllRefCases = 0; // number of cases that contain at least some non-ref read.
        int noReads = 0;
        int lowMinMQ = 0;
        for (int i = 0; i < size; i++) {
            final DragstrLocus locus = loci.get(i);
            final SimpleInterval startLocation = locus.getStartInterval(dictionary, 0);
            final SimpleInterval repeatInterval = locus.getRepeatInterval(dictionary, config.pileupPadding, config.pileupPadding + period);
            final List<GATKRead> reads = Utils.stream(readsSource
                    .query(startLocation))
                    .filter(m -> (m.getFlags() & 0x0f04) == 0)
                    .filter(m -> m.getMappingQuality() != 0)
                    .filter(m -> !m.isUnmapped())
                    .collect(Collectors.toList());
            if (reads.isEmpty()) {
                if (logger.isDebugEnabled()) {
                    logger.debug(String.format("P=%d, L=%d, j=%d, %s: no alignments found", period, locus.getRepeats(),
                            i + 1, startLocation));
                }
                noReads++;
            } else if (reads.stream().mapToInt(GATKRead::getMappingQuality).min().getAsInt() < config.samplingMinMQ) {
                if (logger.isDebugEnabled()) {
                    logger.debug(String.format("P=%d, L=%d, j=%d, %s: MAPQ check failed", period, locus.getRepeats(),
                            i + 1, startLocation));
                }
                lowMinMQ++;
            } else {
                byte[] referenceSequence = referenceSource.queryAndPrefetch(repeatInterval).getBases();
                final ReferenceBases referenceBases = new ReferenceBases(referenceSequence, repeatInterval);
                final IntervalPileup pileup = IntervalPileup.of(reads, referenceBases);
                final int startColumn = (int) locus.getStart() - repeatInterval.getStart();
                int j, k, l;
                for (j = startColumn, k = j + locus.getRepeats() * period, l = 0; l < period && k < referenceSequence.length; l++, j++, k++) {
                    if (referenceSequence[j] != referenceSequence[k]) {
                        break;
                    }
                }
                final int endOfRepeatWithExtraBases = l + startColumn + locus.getRepeats() * period - 1;
                final int firstRelevantColumn = 0;
                final int lastRelevantColumn = Math.min(endOfRepeatWithExtraBases + config.pileupPadding, pileup.width() - 1);
                int total = 0;
                int nonRefTotal = 0;
                int readsWithNoBases = 0;
                int readsWithTooManyBadBQs = 0;
                int readsWithNoQuals = 0;
                row_for:
                for (int row = 0; row < pileup.height(); row++) {
                    if (!pileup.reads().get(row).hasBaseQualities()) {
                        readsWithNoBases++;
                        continue;
                    }
                    int disqualifyingBaseCalls = 0;
                    int indelLength = 0;
                    final IntervalPileup.Element element = pileup.element(row);
                    for (int column = firstRelevantColumn; column <= lastRelevantColumn; column++) {
                        final byte qual = pileup.qualAt(row, column);
                        final byte base = pileup.baseAt(row, column);
                        if (base == IntervalPileup.NO_BASE) {
                            //if (++disqualifyingBaseCalls > config.baseQualExceptionsAllowed) {
                            //readsWithTooManyBadBQs++;
                            readsWithNoBases++;
                            continue row_for;
                            //}
                        } else if (base == IntervalPileup.GAP) {
                            if (column >= startColumn && column <= endOfRepeatWithExtraBases) {
                                indelLength--;
                            }
                        } else if (qual != IntervalPileup.NO_BQ) {
                            if (qual < config.baseQualThreshold
                                    && ++disqualifyingBaseCalls > config.baseQualExceptionsAllowed) {
                                readsWithTooManyBadBQs++;
                                continue row_for;
                            }
                        } else { //if (qual == IntervalPileup.NO_BQ) {
                            readsWithNoQuals++;
                            continue row_for;
                        }
                    }
                    for (final IntervalPileup.Insert insert : element.inserts(startColumn - 1, endOfRepeatWithExtraBases)) {
                        indelLength += insert.length();
                    }
                    total++;
                    if (indelLength != 0) {
                        nonRefTotal++;
                    }
                }
                if (logger.isDebugEnabled()) {
                    logger.debug(String.format("P=%d, L=%d, j=%d, %s: taken with non-ref/total (%d/%d) rows disqualified (%d,%d,%d)", period, locus.getRepeats(),
                            i + 1, startLocation, nonRefTotal, total, readsWithNoBases, readsWithNoQuals, readsWithTooManyBadBQs));
                }
                dest.add(total, nonRefTotal);
                if (dest.size() >= config.maximumNumberOfCases) {
                    logger.info("Finished sampling period = " + period + " and repeat count = " + repeats + " as maximum number of case was reached: " + dest.size());
                    logger.info("Number of cases containing non-ref reads is " + nonAllRefCases + " over " + dest.size());
                    logger.info("Number of dismissed cases are " + noReads + " due to lack of mapped reads and " + lowMinMQ + " due to low minimum MQ");
                    return;
                } else if (nonRefTotal > 0) {
                    if (++nonAllRefCases >= config.targetMinimumNonRefCases && dest.size() >= config.minimumNumberOfCases) {
                        logger.info("Finished sampling period = " + period + " and repeat count = " + repeats + " as minimum non-ref containing cases with a minimum number of cases was reached: " + nonAllRefCases + "/" + dest.size());
                        logger.info("Number of cases containing non-ref reads is " + nonAllRefCases + " over " + dest.size());
                        logger.info("Number of dismissed cases are " + noReads + " due to lack of mapped reads and " + lowMinMQ + " due to low minimum MQ");
                        return;
                    }
                }
            }
        }
        logger.info("Finished sampling period = " + period + " and repeat count = " + repeats + " exhausted all possible cases with non-ref containing " + nonAllRefCases + " over a total of " + dest.size());
        logger.info("Number of cases containing non-ref reads is " + nonAllRefCases + " over " + dest.size());
        logger.info("Number of dismissed cases are " + noReads + " due to lack of mapped reads and " + lowMinMQ + " due to low minimum MQ");
    }
}
