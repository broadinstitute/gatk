package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Bundles together a pileup and a location.
 */
public final class AlignmentContext implements HasGenomeLocation {
    // Definitions:
    //   COMPLETE = full alignment context
    //   FORWARD  = reads on forward strand
    //   REVERSE  = reads on forward strand
    //
    public enum ReadOrientation { COMPLETE, FORWARD, REVERSE }

    private final Locatable loc;
    private final ReadPileup basePileup;

    protected boolean hasPileupBeenDownsampled;

    /**
     * The number of bases we've skipped over in the reference since the last map invocation.
     * By default, nothing is being skipped, so skippedBases == 0.
     */
    private long skippedBases = 0;

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup) {
        this(loc, basePileup, 0, false);
    }

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup, final boolean hasPileupBeenDownsampled) {
        this(loc, basePileup, 0, hasPileupBeenDownsampled);
    }

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup, final long skippedBases) {
        this(loc, basePileup, skippedBases, false);
    }

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup, final long skippedBases, final boolean hasPileupBeenDownsampled ) {
        Utils.nonNull(loc, "BUG: GenomeLoc in Alignment context is null");
        Utils.nonNull(basePileup, "BUG: ReadBackedPileup in Alignment context is null");
        if ( skippedBases < 0 ) {
            throw new IllegalArgumentException("BUG: skippedBases is -1 in Alignment context");
        }

        this.loc = loc;
        this.basePileup = basePileup;
        this.skippedBases = skippedBases;
        this.hasPileupBeenDownsampled = hasPileupBeenDownsampled;
    }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int size() {
        return basePileup.size();
    }

    public String getContig() { return getLocation().getContig(); }

    public long getPosition() { return getLocation().getStart(); }

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
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public Map<String, AlignmentContext> splitContextBySampleName(final String assumedSingleSample, final SAMFileHeader header) {
        final Locatable loc = this.getLocation();
        final Map<String, AlignmentContext> contexts = new HashMap<>();

        for(final String sample: this.getBasePileup().getSamples(header)) {
            final ReadPileup pileupBySample = this.getBasePileup().makeFilteredPileup(pe -> sample.equals(ReadUtils.getSampleName(pe.getRead(), header)));

            // Don't add empty pileups to the split context.
            if(pileupBySample.isEmpty()) {
                continue;
            }

            if(sample != null) {
                contexts.put(sample, new AlignmentContext(loc, pileupBySample));
            } else {
                if(assumedSingleSample == null) {
                    throw new UserException.ReadMissingReadGroup(pileupBySample.iterator().next().getRead());
                }
                contexts.put(assumedSingleSample,new AlignmentContext(loc, pileupBySample));
            }
        }

        return contexts;
    }

    /**
     * Splits the AlignmentContext into one context per read group
     *
     * @return a Map of ReadGroup to AlignmentContext, or an empty map if context has no base pileup
     *
     **/
    public Map<SAMReadGroupRecord, AlignmentContext> splitContextByReadGroup(final Collection<SAMReadGroupRecord> readGroups) {
        final Map<SAMReadGroupRecord, AlignmentContext> contexts = new HashMap<>();

        for (final SAMReadGroupRecord rg : readGroups) {
            final ReadPileup rgPileup = this.getBasePileup().makeFilteredPileup(pe -> rg.getReadGroupId().equals(pe.getRead().getReadGroup()));
            if ( rgPileup != null ) // there we some reads for RG
            {
                contexts.put(rg, new AlignmentContext(this.getLocation(), rgPileup));
            }
        }

        return contexts;
    }

    public static Map<String, AlignmentContext> splitContextBySampleName(final ReadPileup pileup, final SAMFileHeader header) {
        return new AlignmentContext(pileup.getLocation(), pileup).splitContextBySampleName(header);
    }

    public static AlignmentContext joinContexts(final Collection<AlignmentContext> contexts) {
        // validation
        final Locatable loc = contexts.iterator().next().getLocation();
        for(final AlignmentContext context: contexts) {
            if(!loc.equals(context.getLocation())) {
                throw new GATKException("Illegal attempt to join contexts from different genomic locations");
            }
        }

        final List<PileupElement> pe = new ArrayList<>();
        for(final AlignmentContext context: contexts) {
            for(final PileupElement pileupElement: context.basePileup) {
                pe.add(pileupElement);
            }
        }
        return new AlignmentContext(loc, new ReadPileup(loc,pe));
    }
}
