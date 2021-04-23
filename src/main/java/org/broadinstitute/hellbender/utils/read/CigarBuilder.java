package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * This class allows code that manipulates cigars to do so naively by handling complications such as merging consecutive
 * identical operators within the builder.  A CigarBuilder takes care of the following:
 *
 * 1)  Merging consecutive identical operators, eg 10M5M -> 15M
 * 2)  Eliminating leading and trailing deletions, eg 10D10M -> 10M and 10M10D -> 10M
 * 3)  Shifting deletions to the left of adjacent insertions, eg 10M1ID10D -> 10M10D10I
 * 4)  Validating the overall structure of [hard clip] [soft clip] non-clip [soft clip] [hard clip]
 *
 * Edge cases, such as removing a deletion that immediately follows a leading insertion, *are* handled correctly.  See the unit tests.
 *
 * Leading and trailing deletions may be kept by using the non-default CigarBuilder(false) constructor.
 *
 * All of this is achieved simply by invoking add() repeatedly, followed by make().
 */
public class CigarBuilder {

    private final List<CigarElement> cigarElements = new ArrayList<>();

    // track the last operator so we can merge consecutive elements with the same operator
    // for example, adding 3M and 4M is equivalent to adding 7M
    // also we ignore leading deletions so for example 10S + 5D = 10S
    private CigarOperator lastOperator = null;

    private Section section = Section.LEFT_HARD_CLIP;

    private final boolean removeDeletionsAtEnds;

    private int leadingDeletionBasesRemoved = 0;
    private int trailingDeletionBasesRemoved = 0;
    private int trailingDeletionBasesRemovedInMake = 0;

    public CigarBuilder(final boolean removeDeletionsAtEnds) {
        this.removeDeletionsAtEnds = removeDeletionsAtEnds;
    }

    public CigarBuilder() {
        this(true);
    }

    public CigarBuilder add(final CigarElement element) {
        if (element.getLength() == 0) {
            return this;
        }

        final CigarOperator operator = element.getOperator();

        // skip a deletion after clipping ie at the beginning of the read
        // note the edge case of a deletion following a leading insertion, which we also skip
        if (removeDeletionsAtEnds && operator == CigarOperator.DELETION) {
            if(lastOperator == null || lastOperator.isClipping()) {
                leadingDeletionBasesRemoved += element.getLength();
                return this;
            } else if (lastOperator == CigarOperator.INSERTION && (cigarElements.size() == 1 || cigarElements.get(cigarElements.size() - 2).getOperator().isClipping())) {
                leadingDeletionBasesRemoved += element.getLength();
                return this;
            }
        }

        advanceSectionAndValidateCigarOrder(operator);

        // merge consecutive elements with the same operator
        if (operator == lastOperator) {
            final int n = cigarElements.size() - 1;
            cigarElements.set(n, new CigarElement(cigarElements.get(n).getLength() + element.getLength(), operator));
        } else {
            if (lastOperator == null) {
                cigarElements.add(element);
                lastOperator = operator;
            } else if (operator.isClipping()) {
                // if we have just started clipping on the right and realize the last operator was a deletion, remove it
                // if we have just started clipping on the right and the last two operators were a deletion and insertion, remove the deletion
                if (removeDeletionsAtEnds && !lastOperator.consumesReadBases() && !lastOperator.isClipping()) {
                    trailingDeletionBasesRemoved += cigarElements.get(cigarElements.size() - 1).getLength();
                    cigarElements.set(cigarElements.size() - 1, element);
                    lastOperator = operator;
                } else if (removeDeletionsAtEnds && lastTwoElementsWereDeletionAndInsertion()) {
                    trailingDeletionBasesRemoved += cigarElements.get(cigarElements.size() - 2).getLength();
                    cigarElements.set(cigarElements.size() - 2, cigarElements.get(cigarElements.size() - 1));
                    cigarElements.set(cigarElements.size() - 1, element);
                } else {
                    cigarElements.add(element);
                    lastOperator = operator;
                }
            } else if (operator == CigarOperator.DELETION && lastOperator == CigarOperator.INSERTION) {
                // The order of deletion and insertion elements is arbitrary, so to standardize we shift deletions to the left
                // that is, we place the deletion before the insertion and shift the insertion right
                // if the element before the insertion is another deletion, we merge in the new deletion
                // note that the last operator remains an insertion
                final int size = cigarElements.size();
                if (size > 1 && cigarElements.get(size - 2).getOperator() == CigarOperator.DELETION) {
                    cigarElements.set(size - 2, new CigarElement(cigarElements.get(size-2).getLength() + element.getLength(), CigarOperator.DELETION));
                } else {
                    cigarElements.add(cigarElements.size() - 1, element);
                }
            } else {
                cigarElements.add(element);
                lastOperator = operator;
            }
        }

        return this;
    }

    private boolean lastTwoElementsWereDeletionAndInsertion() {
        return lastOperator == CigarOperator.INSERTION && cigarElements.size() > 1 && cigarElements.get(cigarElements.size() - 2).getOperator() == CigarOperator.DELETION;
    }

    public CigarBuilder addAll(final Iterable<CigarElement> elements) {
        for (final CigarElement element : elements) {
            add(element);
        }
        return this;
    }

    public Cigar make(final boolean allowEmpty) {
        Utils.validate(!(section == Section.LEFT_SOFT_CLIP && cigarElements.get(0).getOperator() == CigarOperator.SOFT_CLIP), "cigar is completely soft-clipped");
        trailingDeletionBasesRemovedInMake = 0;
        if (removeDeletionsAtEnds && lastOperator == CigarOperator.DELETION) {
            trailingDeletionBasesRemovedInMake = cigarElements.get(cigarElements.size() - 1).getLength();
            cigarElements.remove(cigarElements.size() - 1);
        } else if (removeDeletionsAtEnds && lastTwoElementsWereDeletionAndInsertion()) {
            trailingDeletionBasesRemovedInMake = cigarElements.get(cigarElements.size() - 2).getLength();
            cigarElements.remove(cigarElements.size() - 2);
        }
        Utils.validate(allowEmpty || !cigarElements.isEmpty(), "No cigar elements left after removing leading and trailing deletions.");
        return new Cigar(cigarElements);    // removing flanking deletions may cause an empty cigar to be output.  We do not throw an error or return null.
    }

    public Cigar make() {
        return make(false);
    }

    private enum Section {LEFT_HARD_CLIP, LEFT_SOFT_CLIP, MIDDLE, RIGHT_SOFT_CLIP, RIGHT_HARD_CLIP}

    // validate that cigar structure is hard clip, soft clip, unclipped, soft clip, hard clip
    private void advanceSectionAndValidateCigarOrder(CigarOperator operator) {
        if (operator == CigarOperator.HARD_CLIP) {
            if (section == Section.LEFT_SOFT_CLIP || section == Section.MIDDLE || section == Section.RIGHT_SOFT_CLIP) {
                section = Section.RIGHT_HARD_CLIP;
            }
        } else if (operator == CigarOperator.SOFT_CLIP) {
            Utils.validate(section != Section.RIGHT_HARD_CLIP, "cigar has already reached its right hard clip");
            if (section == Section.LEFT_HARD_CLIP) {
                section = Section.LEFT_SOFT_CLIP;
            } else if(section == Section.MIDDLE) {
                section = Section.RIGHT_SOFT_CLIP;
            }
        } else {
            Utils.validate(section != Section.RIGHT_SOFT_CLIP && section != Section.RIGHT_HARD_CLIP, "cigar has already reached right clip");
            if (section == Section.LEFT_HARD_CLIP || section == Section.LEFT_SOFT_CLIP) {
                section = Section.MIDDLE;
            }
        }
    }

    /**
     * Count the number of leading deletion bases that have been removed by this builder and that will not show up in any call to make().
     * Note that all leading deletions are removed prior to calling make().  For example, successively adding 3S2D10I7D10M would result in
     * the 2D and 7D elements being discarded, for a total of 9 removed deletion bases.
     */
    public int getLeadingDeletionBasesRemoved() {
        return leadingDeletionBasesRemoved;
    }

    /**
     * Counts the number of trailing deletion bases that were removed in the last call to make().  These may be removed
     * before or during make().  For example, adding 3M and 3D does not removed the 3D because the builder does not know that 3D
     * is a terminal element.  If make() is then called, the builder will record the discarded 3D and this method will return 3.
     * Subsequently adding 3M, calling make(), and then calling this method will result in 0.
     */
    public int getTrailingDeletionBasesRemoved() {
        return trailingDeletionBasesRemoved + trailingDeletionBasesRemovedInMake;
    }

    /**
     * Return a Result object containing the output of make() as well as the number of leading and trailing deletion bases
     * removed relative to the cigar elements that were add()ed.  This is very useful when in addition to transforming a cigar we must also
     * keep track of an alignment start or end.
     */
    public Result makeAndRecordDeletionsRemovedResult() {
        final Cigar cigar = make();
        return new Result(cigar, getLeadingDeletionBasesRemoved(), getTrailingDeletionBasesRemoved());
    }

    public static final class Result {
        private Cigar cigar;
        private final int leadingDeletionBasesRemoved;
        private final int trailingDeletionBasesRemoved;

        public Result(final Cigar cigar, final int leadingDeletionBasesRemoved, final int trailingDeletionBasesRemoved) {
            this.cigar = cigar;
            this.leadingDeletionBasesRemoved = leadingDeletionBasesRemoved;
            this.trailingDeletionBasesRemoved = trailingDeletionBasesRemoved;
        }

        public Cigar getCigar() {
            return cigar;
        }

        public int getLeadingDeletionBasesRemoved() {
            return leadingDeletionBasesRemoved;
        }

        public int getTrailingDeletionBasesRemoved() {
            return trailingDeletionBasesRemoved;
        }
    }
}