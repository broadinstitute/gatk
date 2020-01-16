package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class CigarBuilder {

    private final List<CigarElement> cigarElements = new ArrayList<>();

    // track the last operator so we can merge consecutive elements with the same operator
    // for example, adding 3M and 4M is equivalent to adding 7M
    // also we ignore leading deletions so for example 10S + 5D = 10S
    private CigarOperator lastOperator = null;

    private Section section = Section.LEFT_HARD_CLIP;

    public CigarBuilder() { }

    public CigarBuilder add(final CigarElement element) {
        final CigarOperator operator = element.getOperator();

        // skip a deletion after clipping ie at the beginning of the read
        // note the edge case of a deletion following a leading insertion, which we also skip
        if (operator == CigarOperator.DELETION) {
            if(lastOperator == null || lastOperator.isClipping()) {
                return this;
            } else if (lastOperator == CigarOperator.INSERTION && (cigarElements.size() == 1 || cigarElements.get(cigarElements.size() - 2).getOperator().isClipping())) {
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
            } else if (operator.isClipping() && !lastOperator.consumesReadBases() && !lastOperator.isClipping()) {
                // if we have just start clipping on the right and realize the last operator was a deletion, remove it
                cigarElements.set(cigarElements.size() - 1, element);
                lastOperator = operator;
            } else if (operator == CigarOperator.INSERTION && lastOperator == CigarOperator.DELETION) {
                // The order of deletion and insertion elements is arbitrary, so to standardize we shift insertions to the left
                // that is, we place the insertion before the deletion and shift the deletion right
                // if the element before the deletion is another insertion, we merge in the new insertion
                // note that the last operator remains a deletion
                final int size = cigarElements.size();
                if (size > 1 && cigarElements.get(size - 2).getOperator() == CigarOperator.INSERTION) {
                    cigarElements.set(size - 2, new CigarElement(cigarElements.get(size-2).getLength() + element.getLength(), CigarOperator.INSERTION));
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

    public CigarBuilder addAll(final Collection<CigarElement> elements) {
        for (final CigarElement element : elements) {
            add(element);
        }
        return this;
    }

    public Cigar make() {
        Utils.validate(!(section == Section.LEFT_SOFT_CLIP && cigarElements.get(0).getOperator() == CigarOperator.SOFT_CLIP), "cigar is completely soft-clipped");
        if (lastOperator == CigarOperator.DELETION) {
            cigarElements.remove(cigarElements.size() - 1);
        }
        return new Cigar(cigarElements);
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
}