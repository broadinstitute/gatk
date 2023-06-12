package org.broadinstitute.hellbender.utils.haplotype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Extract simple VariantContext events from a single haplotype
 */
public final class EventMap extends TreeMap<Integer, Event> {
    private static final long serialVersionUID = 1L;

    private static final Logger logger = LogManager.getLogger(EventMap.class);

    public EventMap(final Collection<Event> events) {
        super();
        events.forEach(this::addEvent);
    }

    public static EventMap fromHaplotype(final Haplotype haplotype, final byte[] ref, final Locatable refLoc, final int maxMnpDistance) {
        return new EventMap(getEvents(haplotype, ref, refLoc, maxMnpDistance));
    }

    public static EventMap fromHaplotype(final Haplotype haplotype, final byte[] ref, final int maxMnpDistance) {
        return new EventMap(getEvents(haplotype, ref, haplotype.getLocation(), maxMnpDistance));
    }

    // this is really just a convenient way to make EventMap objects in unit tests
    @VisibleForTesting
    public static EventMap of(final Event ... events) {
        return new EventMap(Arrays.asList(events));
    }

    /**
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occurring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     */
    private static List<Event> getEvents(final Haplotype haplotype, final byte[] ref, final Locatable refLoc, final int maxMnpDistance) {
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
        final Cigar cigar = haplotype.getCigar();
        final byte[] alignment = haplotype.getBases();

        int refPos = haplotype.getAlignmentStartHapwrtRef();
        if( refPos < 0 ) {
            return Collections.emptyList();
        } // Protection against SW failures

        final List<Event> proposedEvents = new ArrayList<>();

        int alignmentPos = 0;

        for( int cigarIndex = 0; cigarIndex < cigar.numCigarElements(); cigarIndex++ ) {
            final CigarElement ce = cigar.getCigarElement(cigarIndex);
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
                case I:
                {
                    if( refPos > 0 && cigarIndex > 0 && cigarIndex < cigar.numCigarElements() - 1) { // forbid insertions at start of contig or not resolved within the haplotype
                        final int insertionStart = refLoc.getStart() + refPos - 1;
                        final byte refByte = ref[refPos-1];

                        byte[] insertionBases = {refByte};  // add the padding base
                        insertionBases = ArrayUtils.addAll(insertionBases, Arrays.copyOfRange(alignment, alignmentPos, alignmentPos + elementLength));
                        if( BaseUtils.isAllRegularBases(insertionBases) ) {
                            proposedEvents.add(new Event(refLoc.getContig(), insertionStart, Allele.create(refByte, true), Allele.create(insertionBases)));
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
                    if( refPos > 0 ) { // forbid deletions at the beginning of a contig
                        final byte[] deletionBases = Arrays.copyOfRange( ref, refPos - 1, refPos + elementLength );  // add padding base
                        final int deletionStart = refLoc.getStart() + refPos - 1;
                        final byte refByte = ref[refPos-1];
                        if( BaseUtils.isRegularBase(refByte) && BaseUtils.isAllRegularBases(deletionBases) ) {
                            proposedEvents.add(new Event(refLoc.getContig(), deletionStart, Allele.create(deletionBases, true), Allele.create(refByte)));
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
                        final Event event = new Event(refLoc.getContig(), refLoc.getStart() + refPos + start, refAllele, altAllele);

                        if ( haplotype.isCollapsed() ) {
                            event.setVariantAttribute(AssemblyBasedCallerUtils.EXT_COLLAPSED_TAG, "1");
                        }
                        proposedEvents.add(event);
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
        return proposedEvents;
    }

    // add an event, combining it into a compound event if an event already exists at the same position
    private void addEvent(final Event newEvent) {
        Utils.nonNull(newEvent);
        computeIfPresent(newEvent.getStart(), (pos, oldEvent) -> makeCompoundEvents(oldEvent, newEvent));
        putIfAbsent(newEvent.getStart(), newEvent);
    }

    /**
     * Merge two events with the same start into a single compound event.  The resulting event will not be a SNP
     * or a simple indel.  This should not in any way be conflated with making a multiallelic VariantContext.  Here we
     * apply the substitution encoded by one event and then apply the substitution encoded by a second event to obtain
     * a single complex event.
     *
     * Examples: A -> C + A -> AT = A -> CT; C -> CGGG + CAAA -> C = CAAA -> CGGG
     *
     * e1 can be SNP, and e2 can then be either a insertion or deletion.
     * If e1 is an indel, then e2 must be the opposite type (e1 deletion => e2 must be an insertion)
     */
    protected static Event makeCompoundEvents(final Event e1, final Event e2) {
        Utils.validateArg( e1.getStart() == e2.getStart(), "e1 and e2 must have the same start");

        if ( e1.isSNP() || e2.isSNP()) {
            Utils.validateArg(!(e1.isSNP() && e2.isSNP()), "Trying to put two overlapping SNPs in one EventMap.  This could be a CIGAR bug.");
            final Event snp = e1.isSNP() ? e1 : e2;
            final Event indel = e1.isSNP() ? e2 : e1;

            // SNP + insertion. Example: A -> G (e1) + A -> CT (e2) = A -> GT
            if ( snp.refAllele().equals(indel.refAllele()) ) {
                return new Event(snp.getContig(), snp.getStart(), snp.refAllele(), Allele.create(snp.altAllele().getDisplayString() + indel.altAllele().getDisplayString().substring(1)));
            } else {
                // SNP + deletion. Example:  A -> T + AC -> A = AC -> T
                return new Event(snp.getContig(), snp.getStart(), indel.refAllele(), snp.altAllele());
            }
        } else {    // insertion + deletion.  Example: AC -> A + A -> AGT = AC -> AGT
            Utils.validateArg ( (e1.isSimpleDeletion() && e2.isSimpleInsertion()) || (e1.isSimpleInsertion() && e2.isSimpleDeletion()),
                    () -> "Can only merge single insertion with deletion (or vice versa)");
            final Event insertion = e1.isSimpleInsertion() ? e1 : e2;
            final Event deletion  = e1.isSimpleInsertion() ? e2 : e1;
            return new Event(e1.getContig(), e1.getStart(), deletion.refAllele(), insertion.altAllele());
        }
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
    public Collection<Event> getEvents() {
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
        for ( final Event event : getEvents() )
            b.append(String.format("%s:%d-%d %s,", event.getContig(), event.getStart(), event.getEnd(), Arrays.asList(event.refAllele(), event.altAllele())));
        b.append("}");
        return b.toString();
    }

    /**
     * Build event maps for each haplotype
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
     */
    public static void buildEventMapsForHaplotypes( final Collection<Haplotype> haplotypes,
                                                                final byte[] ref,
                                                                final Locatable refLoc,
                                                                final boolean debug,
                                                                final int maxMnpDistance) {
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
        int hapNumber = 0;

        if( debug ) logger.info("=== Best Haplotypes ===");
        for( final Haplotype h : haplotypes ) {
            //TODO h.recomputeAndSetEventMap() with overrides for the two haplotype classes should replace this casting+checking code here which is prone to error.

            // Since PD haplotypes Know what alleles are variants, simply ask it and generate the map that way.
            final EventMap events = h.isPartiallyDetermined() ? new EventMap(((PartiallyDeterminedHaplotype) h).getDeterminedAlleles()) :
                    fromHaplotype(h, ref, refLoc, maxMnpDistance);
            h.setEventMap(events);

            if (debug) {
                logger.info(h.toString());
                logger.info("> Cigar = " + h.getCigar());
                logger.info(">> Events = " + h.getEventMap());
            }
        }
    }

    /**
     * Return the sorted set of all of the starting positions of all events across all haplotypes
     */
    public static TreeSet<Integer> getEventStartPositions(final Collection<Haplotype> haplotypes) {
        final TreeSet<Integer> result = new TreeSet<>();
        for( final Haplotype h : haplotypes ) {
            Utils.nonNull(h.getEventMap(), "Haplotype event map has not been set");
            result.addAll(h.getEventMap().getStartPositions());
        }
        return result;
    }

    /**
     * Returns any events in the map that overlap loc, including spanning deletions and events that start at loc.
     */
    public List<Event> getOverlappingEvents(final int loc) {
        final List<Event> overlappingEvents = headMap(loc, true).values().stream().filter(v -> v.getEnd() >= loc).collect(Collectors.toList());
        // if we're at the start of an insertion, exclude deletions that end here; otherwise keep everything
        final Predicate<Event> filter = overlappingEvents.stream().anyMatch(Event::isSimpleInsertion) ?
                v -> !(v.isSimpleDeletion() && v.getEnd() == loc) : v -> true;
        return overlappingEvents.stream().filter(filter).collect(Collectors.toList());
    }
}
