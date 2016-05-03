package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * Bundles together a pileup and a location.
 */
public final class AlignmentContext implements Locatable, HasGenomeLocation {
    // Definitions:
    //   COMPLETE = full alignment context
    //   FORWARD  = reads on forward strand
    //   REVERSE  = reads on forward strand
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
        return this.splitContextBySampleName((String)null, header);
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample, but referencd by sample name instead
     * of sample object.
     *
     * @param assumedSingleSample If non-null, assume this is the only sample in our pileup and return immediately.
     *                            If null, get the list of samples from the provided header and do the work of splitting by sample.
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public Map<String, AlignmentContext> splitContextBySampleName(final String assumedSingleSample, final SAMFileHeader header) {
        if (assumedSingleSample != null){
            return Collections.singletonMap(assumedSingleSample, this);
        }

        final Map<String, AlignmentContext> contexts = new HashMap<>();

        for(final String sample: basePileup.getSamples(header)) {
            final ReadPileup pileupForSample = basePileup.makeFilteredPileup(pe -> Objects.equals(sample, ReadUtils.getSampleName(pe.getRead(), header)));

            if (sample == null){
                throw new UserException.ReadMissingReadGroup(pileupForSample.iterator().next().getRead());
            }

            // Don't add empty pileups to the split context.
            if(! pileupForSample.isEmpty()) {
                contexts.put(sample, new AlignmentContext(loc, pileupForSample));
            }
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
