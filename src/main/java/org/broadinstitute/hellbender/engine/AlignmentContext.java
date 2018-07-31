package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Bundles together a pileup and a location.
 */
public final class AlignmentContext implements Locatable, HasGenomeLocation {
    // Definitions:
    //   COMPLETE = full alignment context
    //   FORWARD  = reads on forward strand
    //   REVERSE  = reads on reverse strand
    //
    public enum ReadOrientation { COMPLETE, FORWARD, REVERSE }

    private final Locatable loc;
    private final ReadPileup basePileup;

    private boolean hasPileupBeenDownsampled;

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup) {
        this(loc, basePileup, false);
    }

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup, final boolean hasPileupBeenDownsampled ) {
        Utils.nonNull(loc, "BUG: GenomeLoc in Alignment context is null");
        Utils.nonNull(basePileup, "BUG: ReadBackedPileup in Alignment context is null");

        this.loc = loc;
        this.basePileup = basePileup;
        this.hasPileupBeenDownsampled = hasPileupBeenDownsampled;
    }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int size() {
        return basePileup.size();
    }

    @Override
    public String getContig() {
        return getLocation().getContig();
    }

    @Override
    public int getStart() {
        return getLocation().getStart();
    }

    @Override
    public int getEnd() {
        return getLocation().getEnd();
    }

    public long getPosition() { return getStart(); }

    @Override
    public Locatable getLocation() { return loc; }

    /**
     * Returns true if any reads have been filtered out of the pileup due to excess DoC.
     * @return True if reads have been filtered out.  False otherwise.
     */
    public boolean hasPileupBeenDownsampled() { return hasPileupBeenDownsampled; }

    public ReadPileup getBasePileup() {
        return basePileup;
    }

    /**
     * Returns a potentially derived subcontext containing only forward, reverse, or in fact all reads
     * in alignment context context.
     */
    public AlignmentContext stratify(final ReadOrientation type) {
        switch(type) {
            case COMPLETE:
                return this;
            case FORWARD:
                return new AlignmentContext(loc, basePileup.makeFilteredPileup(pe -> !pe.getRead().isReverseStrand()));
            case REVERSE:
                return new AlignmentContext(loc, basePileup.makeFilteredPileup(pe -> pe.getRead().isReverseStrand()));
            default:
                throw new IllegalArgumentException("Unable to get alignment context for type = " + type);
        }
    }

    public Map<String, AlignmentContext> splitContextBySampleName(final SAMFileHeader header) {
        return this.splitContextBySampleName((String) null, header);
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample, but referenced by sample name instead
     * of sample object.
     *
     * @param assumedSingleSample If non-null, assume this is the only sample in our pileup and return immediately.
     *                            If null, get the list of samples from the provided header and do the work of splitting by sample.
     * @return a Map of sample name to StratifiedAlignmentContext; samples without coverage are not included
     **/
    public Map<String, AlignmentContext> splitContextBySampleName(final String assumedSingleSample, final SAMFileHeader header) {
        if (assumedSingleSample != null){
            return Collections.singletonMap(assumedSingleSample, this);
        }
        final Locatable loc = this.getLocation();
        // this will throw an user error if there are samples without RG/sampleName
        final Map<String, ReadPileup> pileups = this.getBasePileup().splitBySample(header, assumedSingleSample);
        final Map<String, AlignmentContext> contexts = new LinkedHashMap<>(pileups.size());
        for (final Map.Entry<String, ReadPileup> entry : pileups.entrySet()) {
            // Don't add empty pileups to the split context.
            if (entry.getValue().isEmpty()) {
                continue;
            }
            contexts.put(entry.getKey(), new AlignmentContext(loc, entry.getValue()));
        }
        return contexts;
    }

    public static Map<String, AlignmentContext> splitContextBySampleName(final ReadPileup pileup, final SAMFileHeader header) {
        return new AlignmentContext(pileup.getLocation(), pileup).splitContextBySampleName(header);
    }

    @Override
    public String toString() {
        return "AlignmentContext{" +
                "loc=" + loc +
                ", basePileup=" + basePileup +
                ", hasPileupBeenDownsampled=" + hasPileupBeenDownsampled +
                '}';
    }
}
