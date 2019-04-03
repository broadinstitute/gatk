package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Extract simple VariantContext events from a single haplotype
 */
public final class EventMap extends TreeMap<Integer, VariantContext> {
    private static final long serialVersionUID = 1L;

    private static final Logger logger = LogManager.getLogger(EventMap.class);
    protected static final int MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION = 3;
    private static final int MAX_EVENTS_PER_HAPLOTYPE = 3;
    private static final int MAX_INDELS_PER_HAPLOTYPE = 2;
    public static final Allele SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele.create("<UNASSEMBLED_EVENT>", false);

    private final Haplotype haplotype;
    private final byte[] ref;
    private final Locatable refLoc;
    private final String sourceNameToAdd;

    public EventMap(final Haplotype haplotype, final byte[] ref, final Locatable refLoc, final String sourceNameToAdd, final int maxMnpDistance) {
        super();
        this.haplotype = haplotype;
        this.ref = ref;
        this.refLoc = refLoc;
        this.sourceNameToAdd = sourceNameToAdd;
        processCigarForInitialEvents(maxMnpDistance);
    }

    /**
     * For testing.  Let's you set up a explicit configuration without having to process a haplotype and reference
     * @param stateForTesting
     */
    public EventMap(final Collection<VariantContext> stateForTesting) {
        haplotype = null;
        ref = null;
        refLoc = null;
        sourceNameToAdd = null;
        for ( final VariantContext vc : stateForTesting )
            addVC(vc);
    }

    /**
     *
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occurring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     */
    protected void processCigarForInitialEvents(final int maxMnpDistance) {
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
        final Cigar cigar = haplotype.getCigar();
        final byte[] alignment = haplotype.getBases();

        int refPos = haplotype.getAlignmentStartHapwrtRef();
        if( refPos < 0 ) {
            return;
        } // Protection against SW failures

        final List<VariantContext> proposedEvents = new ArrayList<>();

        int alignmentPos = 0;

        for( int cigarIndex = 0; cigarIndex < cigar.numCigarElements(); cigarIndex++ ) {
            final CigarElement ce = cigar.getCigarElement(cigarIndex);
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    if( refPos > 0 ) { // protect against trying to create insertions/deletions at the beginning of a contig
                        final List<Allele> insertionAlleles = new ArrayList<>();
                        final int insertionStart = refLoc.getStart() + refPos - 1;
                        final byte refByte = ref[refPos-1];
                        if( BaseUtils.isRegularBase(refByte) ) {
                            insertionAlleles.add( Allele.create(refByte, true) );
                        }
                        if( cigarIndex == 0 || cigarIndex == cigar.numCigarElements() - 1 ) {
                            // if the insertion isn't completely resolved in the haplotype, skip it
                            // note this used to emit SYMBOLIC_UNASSEMBLED_EVENT_ALLELE but that seems dangerous
                        } else {
                            byte[] insertionBases = {};
                            insertionBases = ArrayUtils.add(insertionBases, ref[refPos - 1]); // add the padding base
                            insertionBases = ArrayUtils.addAll(insertionBases, Arrays.copyOfRange(alignment, alignmentPos, alignmentPos + elementLength));
                            if( BaseUtils.isAllRegularBases(insertionBases) ) {
                                insertionAlleles.add( Allele.create(insertionBases, false) );
                            }
                        }
                        if( insertionAlleles.size() == 2 ) { // found a proper ref and alt allele
                            proposedEvents.add(new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), insertionStart, insertionStart, insertionAlleles).make());
                        }
                    }
                    alignmentPos += elementLength;
                    break;
                }
                case S:
                {
                    alignmentPos += elementLength;
                    break;
                }
                case D:
                {
                    if( refPos > 0 ) { // protect against trying to create insertions/deletions at the beginning of a contig
                        final byte[] deletionBases = Arrays.copyOfRange( ref, refPos - 1, refPos + elementLength );  // add padding base
                        final List<Allele> deletionAlleles = new ArrayList<>();
                        final int deletionStart = refLoc.getStart() + refPos - 1;
                        final byte refByte = ref[refPos-1];
                        if( BaseUtils.isRegularBase(refByte) && BaseUtils.isAllRegularBases(deletionBases) ) {
                            deletionAlleles.add( Allele.create(deletionBases, true) );
                            deletionAlleles.add( Allele.create(refByte, false) );
                            proposedEvents.add(new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), deletionStart, deletionStart + elementLength, deletionAlleles).make());
                        }
                    }
                    refPos += elementLength;
                    break;
                }
                case M:
                case EQ:
                case X:
                {
                    final Queue<Integer> mismatchOffsets = new ArrayDeque<>();    // a FIFO queue of substitutions
                    for( int offset = 0; offset < elementLength; offset++ ) {
                        final byte refByte = ref[refPos + offset ];
                        final byte altByte = alignment[alignmentPos + offset];
                        final boolean mismatch = refByte != altByte && BaseUtils.isRegularBase(refByte) && BaseUtils.isRegularBase(altByte);
                        if (mismatch) {
                            mismatchOffsets.add(offset);
                        }
                    }

                    while (!mismatchOffsets.isEmpty()) {
                        final int start = mismatchOffsets.poll();
                        int end = start;
                        while (!mismatchOffsets.isEmpty() && mismatchOffsets.peek() - end <= maxMnpDistance) {
                            end = mismatchOffsets.poll();
                        }
                        final Allele refAllele = Allele.create(Arrays.copyOfRange(ref, refPos + start, refPos + end + 1), true);
                        final Allele altAllele = Allele.create(Arrays.copyOfRange(alignment, alignmentPos + start, alignmentPos + end + 1), false);
                        proposedEvents.add(new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), refLoc.getStart() + refPos + start, refLoc.getStart() + refPos + end, Arrays.asList(refAllele, altAllele)).make());
                    }

                    // move refPos and alignmentPos forward to the end of this cigar element
                    refPos += elementLength;
                    alignmentPos += elementLength;
                    break;
                }
                case N:
                case H:
                case P:
                default:
                    throw new GATKException( "Unsupported cigar operator created during SW alignment: " + ce.getOperator() );
            }
        }

        for ( final VariantContext proposedEvent : proposedEvents )
            addVC(proposedEvent, true);
    }

    /**
     * Add VariantContext vc to this map, merging events with the same start sites if necessary
     * @param vc the variant context to add
     */
    public void addVC(final VariantContext vc) {
        addVC(vc, true);
    }

    /**
     * Add VariantContext vc to this map
     * @param vc the variant context to add
     * @param merge should we attempt to merge it with an already existing element, or should we throw an error in that case?
     */
    public void addVC(final VariantContext vc, final boolean merge) {
        Utils.nonNull(vc);
        if ( containsKey(vc.getStart()) ) {
            Utils.validate(merge, () -> "Will not merge previously bound variant contexts as merge is false at " + vc);
            final VariantContext prev = get(vc.getStart());
            put(vc.getStart(), makeBlock(prev, vc));
        } else
            put(vc.getStart(), vc);
    }

    /**
     * Create a block substitution out of two variant contexts that start at the same position
     *
     * vc1 can be SNP, and vc2 can then be either a insertion or deletion.
     * If vc1 is an indel, then vc2 must be the opposite type (vc1 deletion => vc2 must be an insertion)
     *
     * @param vc1 the first variant context we want to merge
     * @param vc2 the second
     * @return a block substitution that represents the composite substitution implied by vc1 and vc2
     */
    protected VariantContext makeBlock(final VariantContext vc1, final VariantContext vc2) {
        Utils.validateArg( vc1.getStart() == vc2.getStart(), () -> "vc1 and 2 must have the same start but got " + vc1 + " and " + vc2);
        Utils.validateArg( vc1.isBiallelic(), "vc1 must be biallelic");
        if ( ! vc1.isSNP() ) {
            Utils.validateArg ( (vc1.isSimpleDeletion() && vc2.isSimpleInsertion()) || (vc1.isSimpleInsertion() && vc2.isSimpleDeletion()),
                    () -> "Can only merge single insertion with deletion (or vice versa) but got " + vc1 + " merging with " + vc2);
        } else {
            Utils.validateArg(!vc2.isSNP(), () -> "vc1 is " + vc1 + " but vc2 is a SNP, which implies there's been some terrible bug in the cigar " + vc2);
        }

        final Allele ref, alt;
        final VariantContextBuilder b = new VariantContextBuilder(vc1);
        if ( vc1.isSNP() ) {
            // we have to repair the first base, so SNP case is special cased
            if ( vc1.getReference().equals(vc2.getReference()) ) {
                // we've got an insertion, so we just update the alt to have the prev alt
                ref = vc1.getReference();
                alt = Allele.create(vc1.getAlternateAllele(0).getDisplayString() + vc2.getAlternateAllele(0).getDisplayString().substring(1), false);
            } else {
                // we're dealing with a deletion, so we patch the ref
                ref = vc2.getReference();
                alt = vc1.getAlternateAllele(0);
                b.stop(vc2.getEnd());
            }
        } else {
            final VariantContext insertion = vc1.isSimpleInsertion() ? vc1 : vc2;
            final VariantContext deletion  = vc1.isSimpleInsertion() ? vc2 : vc1;
            ref = deletion.getReference();
            alt = insertion.getAlternateAllele(0);
            b.stop(deletion.getEnd());
        }

        return b.alleles(Arrays.asList(ref, alt)).make();
    }

    // TODO -- warning this is an O(N^3) algorithm because I'm just lazy.  If it's valuable we need to reengineer it
    protected void replaceClumpedEventsWithBlockSubstitutions() {
        if ( getNumberOfEvents() >= MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION) {
            int lastStart = -1;
            for ( boolean foundOne = true; foundOne; ) {
                foundOne = false;
                for ( final VariantContext vc : getVariantContexts() ) {
                    if ( vc.getStart() > lastStart ) {
                        lastStart = vc.getStart();
                        final List<VariantContext> neighborhood = getNeighborhood(vc, 10);
                        if ( updateToBlockSubstitutionIfBetter(neighborhood) ) {
                            foundOne = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    protected boolean updateToBlockSubstitutionIfBetter(final List<VariantContext> neighbors) {
        if (neighbors.size() < MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION)
            return false;
        // TODO -- need more tests to decide if this is really so good

        final VariantContext first = neighbors.get(0);
        final int refStartOffset = first.getStart() - refLoc.getStart();
        final int refEndOffset = neighbors.get(neighbors.size() - 1).getEnd() - refLoc.getStart();

        final byte[] refBases = Arrays.copyOfRange(ref, refStartOffset, refEndOffset + 1);
        final byte[] hapBases = AlignmentUtils.getBasesCoveringRefInterval(refStartOffset, refEndOffset, haplotype.getBases(), haplotype.getAlignmentStartHapwrtRef(), haplotype.getCigar());

        final VariantContextBuilder builder = new VariantContextBuilder(first);
        builder.stop(first.getStart() + refBases.length - 1);
        builder.alleles(Arrays.asList(Allele.create(refBases, true), Allele.create(hapBases)));
        final VariantContext block = builder.make();

        // remove all merged events
        for ( final VariantContext merged : neighbors ) {
            if ( remove(merged.getStart()) == null )
                throw new IllegalArgumentException("Expected to remove variant context from the event map but remove said there wasn't any element there: " + merged);
        }

        // note must be after we remove the previous events as the treeset only allows one key per start
        logger.info("Transforming into block substitution at " + block);
        addVC(block, false);

        return true;
    }

    /**
     * Get all of the variant contexts starting at leftMost that are within maxBP of each other
     *
     * @param leftMost the left most (smallest position) variant context that will start the neighborhood
     * @param maxBPBetweenEvents the maximum distance in BP between the end of one event the start of the next
     *                           to be included the the resulting list
     * @return a list that contains at least one element (leftMost)
     */
    protected List<VariantContext> getNeighborhood(final VariantContext leftMost, final int maxBPBetweenEvents) {
        final List<VariantContext> neighbors = new LinkedList<>();

        VariantContext left = leftMost;
        for ( final VariantContext vc : getVariantContexts() ) {
            if ( vc.getStart() < leftMost.getStart() )
                continue;

            if ( vc.getStart() - left.getEnd() < maxBPBetweenEvents ) {
                // this vc is within max distance to the end of the left event, so accumulate it
                neighbors.add(vc);
                left = vc;
            }
        }

        return neighbors;
    }

    /**
     * Get the starting positions of events in this event map
     * @return
     */
    public Set<Integer> getStartPositions() {
        return keySet();
    }

    /**
     * Get the variant contexts in order of start position in this event map
     * @return
     */
    public Collection<VariantContext> getVariantContexts() {
        return values();
    }

    /**
     * How many events do we have?
     * @return
     */
    public int getNumberOfEvents() {
        return size();
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("EventMap{");
        for ( final VariantContext vc : getVariantContexts() )
            b.append(String.format("%s:%d-%d %s,", vc.getContig(), vc.getStart(), vc.getEnd(), vc.getAlleles()));
        b.append("}");
        return b.toString();
    }

    /**
     * Build event maps for each haplotype, returning the sorted set of all of the starting positions of all
     * events across all haplotypes
     *
     * @param haplotypes a list of haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     * @param debug if true, we'll emit debugging information during this operation
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occuring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     * @return a sorted set of start positions of all events among all haplotypes
     */
    public static TreeSet<Integer> buildEventMapsForHaplotypes( final List<Haplotype> haplotypes,
                                                                final byte[] ref,
                                                                final Locatable refLoc,
                                                                final boolean debug,
                                                                final int maxMnpDistance) {
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<>();
        int hapNumber = 0;

        if( debug ) logger.info("=== Best Haplotypes ===");
        for( final Haplotype h : haplotypes ) {
            // Walk along the alignment and turn any difference from the reference into an event
            h.setEventMap(new EventMap(h, ref, refLoc, "HC" + hapNumber++, maxMnpDistance));
            startPosKeySet.addAll(h.getEventMap().getStartPositions());

            if( debug ) {
                logger.info(h.toString());
                logger.info("> Cigar = " + h.getCigar());
                logger.info(">> Events = " + h.getEventMap());
            }
        }

        return startPosKeySet;
    }

    /**
     * Returns any events in the map that overlap loc, including spanning deletions and events that start at loc.
     */
    public List<VariantContext> getOverlappingEvents(final int loc) {
        return headMap(loc, true).values().stream().filter(v -> v.getEnd() >= loc).collect(Collectors.toList());
    }

    private static class VariantContextComparator implements Comparator<VariantContext>, Serializable {
        private static final long serialVersionUID = -2549166273822365485L;

        @Override
        public int compare(VariantContext vc1, VariantContext vc2) {
            return vc1.getStart() - vc2.getStart();
        }
    }

    /**
     * Get all of the VariantContexts in the event maps for all haplotypes, sorted by their start position
     * @param haplotypes the set of haplotypes to grab the VCs from
     * @return a sorted set of variant contexts
     */
    public static SortedSet<VariantContext> getAllVariantContexts( final List<Haplotype> haplotypes ) {
        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<VariantContext> vcs = new TreeSet<>(new VariantContextComparator());

        for( final Haplotype h : haplotypes ) {
            vcs.addAll(h.getEventMap().getVariantContexts());
        }

        return vcs;
    }
}
