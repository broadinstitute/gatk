package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Utilities for unit-testing cigars.
 */
public final class CigarTestUtils {

    @Test
    public static void testRandomValidCigarsAreValid() {
        final List<Cigar> randomCigars = randomValidCigars(new Random(113), 1_000, 5, 50);
        for (final Cigar cigar : randomCigars) {
            Assert.assertNotNull(cigar);
            Assert.assertNull(cigar.isValid("read-name", 0));
        }
    }

    /**
     * Returns a list of randomly generted valid CIGAR instances.
     * @param rdn random number generator to use to generate random elements in the result cigars.
     * @param count number of random cigars to produce.
     * @param maximumNumberOfCoreElements maximum number of core (anything except clipping) operations in any given cigar.
     * @param maximumElementLength maximum length for any cigar element.
     * @param prepend additional cigar to be added a the beginning of the result list.
     * @return never {@code null}.
     */
    public static List<Cigar> randomValidCigars(final Random rdn, final int count, final int maximumNumberOfCoreElements,
                                                final int maximumElementLength, final Cigar ... prepend) {
        Utils.nonNull(rdn);
        ParamUtils.isPositiveOrZero(count, "number of cigars");
        Utils.nonNull(prepend);
        final List<Cigar> result = new ArrayList<>(prepend.length + count);
        for (final Cigar cigar : prepend) {
            result.add(cigar);
        }
        for (int i = 0; i < count; i++) {
            final boolean leftClipping = rdn.nextBoolean();
            final boolean rightClipping = rdn.nextBoolean();
            final boolean leftHardClipping = leftClipping && rdn.nextBoolean();
            final boolean rightHardClipping = rightClipping && rdn.nextBoolean();
            final int leftClippingLength = leftClipping ? rdn.nextInt(maximumElementLength) + 1 : 0;
            final int rightClippingLength = rightClipping ? rdn.nextInt(maximumElementLength) + 1: 0;
            final int leftHardClippingLength = leftHardClipping ? (rdn.nextBoolean() ? leftClippingLength : rdn.nextInt(leftClippingLength) + 1) : 0;
            final int rightHardClippingLength = rightHardClipping ? (rdn.nextBoolean() ? rightClippingLength : rdn.nextInt(rightClippingLength + 1) ): 0;
            final int leftSoftClippingLength = leftClippingLength - leftHardClippingLength;
            final int rightSoftClippingLength = rightClippingLength - rightHardClippingLength;
            final List<CigarElement> coreElements = new ArrayList<>();
            final int coreElementCount = rdn.nextInt(maximumNumberOfCoreElements + 1);
            coreElements.add(new CigarElement(rdn.nextInt(maximumElementLength) + 1, CigarOperator.M));
            for (int j = 0; j < coreElementCount; j++) {
                final CigarOperator op = randomCoreOperator(rdn);
                coreElements.add(new CigarElement(rdn.nextInt(maximumElementLength) + 1, op));
            }
            Collections.shuffle(coreElements, rdn);
            if (!coreElements.get(0).getOperator().isAlignment()) {
                coreElements.add(0, new CigarElement(rdn.nextInt(maximumElementLength) + 1, CigarOperator.M));
            }
            if (!coreElements.get(coreElements.size() - 1).getOperator().isAlignment()) {
                coreElements.add(new CigarElement(rdn.nextInt(maximumElementLength) + 1, CigarOperator.M));
            }
            final List<CigarElement> elements = new ArrayList<>();
            if (leftHardClippingLength > 0) elements.add(new CigarElement(leftHardClippingLength, CigarOperator.H));
            if (leftSoftClippingLength > 0) elements.add(new CigarElement(leftSoftClippingLength, CigarOperator.S));
            elements.addAll(coreElements);
            if (rightSoftClippingLength > 0) elements.add(new CigarElement(rightSoftClippingLength, CigarOperator.S));
            if (rightHardClippingLength > 0) elements.add(new CigarElement(rightHardClippingLength, CigarOperator.H));
            final Cigar cigar = CigarUtils.combineAdjacentCigarElements(new Cigar(elements));
            result.add(cigar);
        }
        return result;
    }

    private static CigarOperator randomCoreOperator(final Random rdn) {
        while (true) {
            final int idx = rdn.nextInt(CigarOperator.values().length);
            final CigarOperator op = CigarOperator.values()[idx];
            if (!op.isClipping()) {
                return op;
            }
        }
    }
}
