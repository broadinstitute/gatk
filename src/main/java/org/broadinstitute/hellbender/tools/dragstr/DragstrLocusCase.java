package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public final class DragstrLocusCase {
    public DragstrLocus getLocus() {
        return locus;
    }

    public int getPeriod() {
        return locus.getPeriod();
    }

    public int getRepeatLength() {
        return locus.getRepeats();
    }

    public int getN() {
        return n;
    }

    public int getK() {
        return k;
    }

    public int getMinMQ() {
        return minMQ;
    }

    public int getNSup() {
        return nSup;
    }

    private final DragstrLocus locus;
    private final int n;
    private final int k;
    private final int minMQ;
    private final int nSup;

    public static DragstrLocusCase create(final DragstrLocus locus, final int n, final int k, final int minMQ, final int nSup) {
        Utils.nonNull(locus);
        Utils.validateArg(n >= 0, "depth must be 0 or positive");
        Utils.validateArg(k >= 0, "the non-ref count must be 0 or positive");
        Utils.validateArg( minMQ >= 0, "the min-MQ must be 0 or positive");
        Utils.validateArg(nSup >= 0, "the number of supplementaries must 0 or positive");
        return new DragstrLocusCase(locus, n, k, minMQ, nSup);
    }

    public DragstrLocusCase(final DragstrLocus locus, final int total, final int nonRef, final int minMQ, final int nSup) {
        this.locus = locus;
        this.n = total;
        this.k = nonRef;
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
        return n >= minDepth && nSup <= maxSup && minMQ >= samplingMinMQ;
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
