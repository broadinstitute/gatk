package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents the DRAGstr model fitting relevant stats at a given locus on the genome for the target sample.
 */
public final class DragstrLocusCase {

    /**
     * The target locus.
     * @return never {@code null}.
     */
    public DragstrLocus getLocus() {
        return locus;
    }

    /**
     * The period length of this STR.
     * @return never {@code null}.
     */
    public int getPeriod() {
        return locus.getPeriod();
    }

    /**
     * The repeat lenght in integer repeated units.
     * @return never {@code null}.
     */
    public int getRepeatLength() {
        return locus.getRepeats();
    }

    /**
     * Returns the total number of qualifying reads that overlap the STR that starts at this locus for the target sample.
     * @return 0 or greater.
     */
    public int getDepth() {
        return depth;
    }

    /**
     * Returns the total number of indel events on reads that overlap the STR a this locus.
     * @return 0 or greater.
     */
    public int getIndels() {
        return indels;
    }

    /**
     * The minimum MQ amongst the reads qualifying reads that overlap this position.
     * @return 0 or greater.
     */
    public int getMinMQ() {
        return minMQ;
    }

    /**
     * Number of qualifying reads records that overlap the STR that are supplementary alignments.
     * @return 0 or greater.
     */
    public int getNSup() {
        return nSup;
    }

    private final DragstrLocus locus;
    private final int depth;
    private final int indels;
    private final int minMQ;
    private final int nSup;

    public static DragstrLocusCase create(final DragstrLocus locus, final int depth, final int indels, final int minMQ, final int nSup) {
        Utils.nonNull(locus);
        Utils.validateArg(depth >= 0, "depth must be 0 or positive");
        Utils.validateArg(indels >= 0, "the non-ref indel count must be 0 or positive");
        Utils.validateArg( minMQ >= 0, "the min-MQ must be 0 or positive");
        Utils.validateArg(nSup >= 0, "the number of supplementaries must 0 or positive");
        return new DragstrLocusCase(locus, depth, indels, minMQ, nSup);
    }

    public DragstrLocusCase(final DragstrLocus locus, final int total, final int nonRef, final int minMQ, final int nSup) {
        this.locus = locus;
        this.depth = total;
        this.indels = nonRef;
        this.minMQ = minMQ;
        this.nSup = nSup;
    }

    /**
     * Tests whether this site can be used for parameter estimation given the limits for relevant stats.
     * @param minDepth minimum depth to qualify.
     * @param samplingMinMQ minimum min-MQ.
     * @param maxSup maximum number of supplementary alignments.
     * @return {@code true} iff it qualifies.
     */
    public boolean qualifies(final int minDepth, final int samplingMinMQ, final int maxSup) {
        return depth >= minDepth && nSup <= maxSup && minMQ >= samplingMinMQ;
    }

    /**
     * Returns a locatable for this caze's locus given the sequence dictionary.
     * @param dictionary
     * @return never {@code null}.
     */
    public Locatable getLocation(final SAMSequenceDictionary dictionary) {
        return new SimpleInterval(dictionary.getSequence(locus.getChromosomeIndex()).getSequenceName(),
                (int) locus.getStart(), (int) locus.getEnd());
    }
}
