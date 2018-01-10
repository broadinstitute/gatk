package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.locusiterator.LIBSTest;
import org.broadinstitute.hellbender.utils.locusiterator.LIBS_position;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByStateBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static htsjdk.samtools.CigarOperator.*;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public final class PileupElementUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "PileupElementTest")
    public Object[][] makePileupElementTest() {
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "PileupElementTest")
    public void testPileupElementTest(final LIBSTest params) {
        final GATKRead read = params.makeRead();
        final AlignmentStateMachine state = new AlignmentStateMachine(read);
        final LIBS_position tester = new LIBS_position(read);

        while (state.stepForwardOnGenome() != null) {
            tester.stepForwardOnGenome();
            final PileupElement pe = state.makePileupElement();

            Assert.assertEquals(pe.getRead(), read);
            Assert.assertEquals(pe.getMappingQual(), read.getMappingQuality());
            Assert.assertEquals(pe.getOffset(), state.getReadOffset());

            Assert.assertEquals(pe.isDeletion(), state.getCigarOperator() == D);
            Assert.assertEquals(pe.isAfterInsertion(), tester.isAfterInsertion);
            Assert.assertEquals(pe.isBeforeInsertion(), tester.isBeforeInsertion);
            Assert.assertEquals(pe.isNextToSoftClip(), tester.isNextToSoftClip);

            if (!hasNeighboringPaddedOps(params.getElements(), pe.getCurrentCigarOffset())) {
                Assert.assertEquals(pe.isAfterDeletionEnd(), tester.isAfterDeletionEnd);
                Assert.assertEquals(pe.isBeforeDeletionStart(), tester.isBeforeDeletionStart);
            }

            Assert.assertEquals(pe.getOffsetInCurrentCigar(), tester.getCurrentPositionOnOperatorBase0(), "CigarElement index failure");
            Assert.assertTrue(pe.getOffsetInCurrentCigar() >= 0, "Offset into current cigar too small");
            Assert.assertTrue(pe.getOffsetInCurrentCigar() < pe.getCurrentCigarElement().getLength(), "Offset into current cigar too big");

            Assert.assertNotNull(pe.toString());

            Assert.assertEquals(pe.atEndOfCurrentCigar(), state.getOffsetIntoCurrentCigarElement() == state.getCurrentCigarElement().getLength() - 1, "atEndOfCurrentCigar failed");
            Assert.assertEquals(pe.atStartOfCurrentCigar(), state.getOffsetIntoCurrentCigarElement() == 0, "atStartOfCurrentCigar failed");

            Assert.assertEquals(pe.getBase(), pe.isDeletion() ? PileupElement.DELETION_BASE : read.getBases()[state.getReadOffset()]);
            Assert.assertEquals(pe.getQual(), pe.isDeletion() ? PileupElement.DELETION_QUAL : read.getBaseQualities()[state.getReadOffset()]);

            Assert.assertEquals(pe.getCurrentCigarElement(), state.getCurrentCigarElement());
            Assert.assertEquals(pe.getCurrentCigarOffset(), state.getCurrentCigarElementOffset());

            final int lengthOfImmediatelyFollowingIndel = pe.getLengthOfImmediatelyFollowingIndel();
            final String basesOfImmediatelyFollowingInsertion = pe.getBasesOfImmediatelyFollowingInsertion();
            Assert.assertTrue(lengthOfImmediatelyFollowingIndel != 0 || basesOfImmediatelyFollowingInsertion == null);
            Assert.assertTrue(basesOfImmediatelyFollowingInsertion == null || basesOfImmediatelyFollowingInsertion.length() == lengthOfImmediatelyFollowingIndel);

            // Don't test -- pe.getBaseIndex();
            if ( pe.atEndOfCurrentCigar() && state.getCurrentCigarElementOffset() < read.numCigarElements() - 1 ) {
                final CigarElement nextElement = read.getCigar().getCigarElement(state.getCurrentCigarElementOffset() + 1);
                if (nextElement.getOperator() == CigarOperator.I) {
                    Assert.assertTrue(pe.getBetweenNextPosition().size() >= 1);
                    Assert.assertEquals(pe.getBetweenNextPosition().get(0), nextElement);
                }
                if (nextElement.getOperator() == M) {
                    Assert.assertTrue(pe.getBetweenNextPosition().isEmpty());
                }
            } else {
                Assert.assertTrue(pe.getBetweenNextPosition().isEmpty());
            }

            if (pe.atStartOfCurrentCigar() && state.getCurrentCigarElementOffset() > 0) {
                final CigarElement prevElement = read.getCigar().getCigarElement(state.getCurrentCigarElementOffset() - 1);
                if (prevElement.getOperator() == CigarOperator.I) {
                    Assert.assertTrue(pe.getBetweenPrevPosition().size() >= 1);
                    Assert.assertEquals(pe.getBetweenPrevPosition().get(pe.getBetweenPrevPosition().size() - 1), prevElement);
                }
                if (prevElement.getOperator() == M) {
                    Assert.assertTrue(pe.getBetweenPrevPosition().isEmpty());
                }
            } else {
                Assert.assertTrue(pe.getBetweenPrevPosition().isEmpty());
            }

            // TODO -- add meaningful tests
            pe.getBaseInsertionQual();
            pe.getBaseDeletionQual();
        }
    }


    @DataProvider(name = "PrevAndNextTest")
    public Object[][] makePrevAndNextTest() {
        final List<Object[]> tests = new LinkedList<>();

        final List<CigarOperator> operators = Arrays.asList(CigarOperator.I, CigarOperator.P, CigarOperator.S);

        for (final CigarOperator firstOp : Arrays.asList(M)) {
            for (final CigarOperator lastOp : Arrays.asList(M, D)) {
                for (final int nIntermediate : Arrays.asList(1, 2, 3)) {
                    for (final List<CigarOperator> combination : Utils.makePermutations(operators, nIntermediate, false)) {
                        final int readLength = 2 + combination.size();
                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 1, readLength);
                        read.setBases(Utils.dupBytes((byte) 'A', readLength));
                        read.setBaseQualities(Utils.dupBytes((byte) 30, readLength));

                        String cigar = "1" + firstOp;
                        for (final CigarOperator op : combination) {
                            cigar += "1" + op;
                        }
                        cigar += "1" + lastOp;
                        read.setCigar(cigar);

                        tests.add(new Object[]{read, firstOp, lastOp, combination});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PrevAndNextTest")
    public void testPrevAndNextTest(final GATKRead read, final CigarOperator firstOp, final CigarOperator lastOp, final List<CigarOperator> ops) {
        final AlignmentStateMachine state = new AlignmentStateMachine(read);

        //before the first 'op' from the list
        state.stepForwardOnGenome();
        final PileupElement pe = state.makePileupElement();
        Assert.assertEquals(pe.getBetweenNextPosition().size(), ops.size());
        Assert.assertEquals(pe.getBetweenPrevPosition().size(), 0);
        assertEqualsOperators(pe.getBetweenNextPosition(), ops);
        Assert.assertEquals(pe.getPreviousOnGenomeCigarElement(), null);
        Assert.assertNotNull(pe.getNextOnGenomeCigarElement());
        Assert.assertEquals(pe.getNextOnGenomeCigarElement().getOperator(), lastOp);

        //after the first 'op' from the list
        state.stepForwardOnGenome();
        final PileupElement pe2 = state.makePileupElement();
        Assert.assertEquals(pe2.getBetweenPrevPosition().size(), ops.size());
        Assert.assertEquals(pe2.getBetweenNextPosition().size(), 0);
        assertEqualsOperators(pe2.getBetweenPrevPosition(), ops);
        Assert.assertNotNull(pe2.getPreviousOnGenomeCigarElement());
        Assert.assertEquals(pe2.getPreviousOnGenomeCigarElement().getOperator(), firstOp);
        Assert.assertEquals(pe2.getNextOnGenomeCigarElement(), null);
    }

    @Test(dataProvider = "PrevAndNextTest")
    public void testImmediateBeforeAndAfterTest(final GATKRead read, final CigarOperator firstOp, final CigarOperator lastOp, final List<CigarOperator> ops) {
        final AlignmentStateMachine state = new AlignmentStateMachine(read);

        //before the first 'op' from the list
        state.stepForwardOnGenome();
        final PileupElement pe1 = state.makePileupElement();
        Assert.assertEquals(pe1.getAdjacentOperator(PileupElement.Direction.PREV), null);
        Assert.assertEquals(pe1.getAdjacentOperator(PileupElement.Direction.NEXT), ops.get(0));
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertEquals(pe1.isImmediatelyBefore(op), ops.get(0) == op, op.toString());
            Assert.assertFalse(pe1.isImmediatelyAfter(op), op.toString());
        }

        //after the first 'op' from the list
        state.stepForwardOnGenome();
        final PileupElement pe2 = state.makePileupElement();
        Assert.assertEquals(pe2.getAdjacentOperator(PileupElement.Direction.PREV), ops.get(ops.size() - 1));
        Assert.assertEquals(pe2.getAdjacentOperator(PileupElement.Direction.NEXT), null);
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe2.isImmediatelyBefore(op), "immediately before " + op);
            Assert.assertEquals(pe2.isImmediatelyAfter(op), op == ops.get(ops.size() - 1), "immediately after " + op);
        }
    }

    @DataProvider(name = "PrevAndNextTest_simple")
    public Object[][] makePrevAndNextTest_simple() {
        final List<Object[]> tests = new LinkedList<>();

        for (final CigarOperator firstOp : Arrays.asList(M)) {
            for (final CigarOperator lastOp : Arrays.asList(M)) {
                for (final CigarOperator middleOp : Arrays.asList(CigarOperator.I, CigarOperator.P, CigarOperator.S)) {
                    final int readLength = 3;
                    final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 1, readLength);
                    read.setBases(Utils.dupBytes((byte) 'A', readLength));
                    read.setBaseQualities(Utils.dupBytes((byte) 30, readLength));

                    String cigar = "1" + firstOp;
                    cigar += "1" + middleOp;
                    cigar += "1" + lastOp;
                    read.setCigar(cigar);

                    tests.add(new Object[]{read, middleOp});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PrevAndNextTest_simple")
    public void testImmediateBeforeAndAfterTest_simple(final GATKRead read, final CigarOperator middleOp) {
        final AlignmentStateMachine state = new AlignmentStateMachine(read);

        state.stepForwardOnGenome();

        //before the 'middleOp'
        final PileupElement pe1 = state.makePileupElement();
        Assert.assertEquals(pe1.getAdjacentOperator(PileupElement.Direction.PREV), null, "PREV");
        Assert.assertEquals(pe1.getAdjacentOperator(PileupElement.Direction.NEXT), middleOp, "NEXT");
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertEquals(pe1.isImmediatelyBefore(op), middleOp == op, op.toString());
            Assert.assertFalse(pe1.isImmediatelyAfter(op), op.toString());
        }

        state.stepForwardOnGenome();

        //after the 'middleOp'
        final PileupElement pe2 = state.makePileupElement();
        Assert.assertEquals(pe2.getAdjacentOperator(PileupElement.Direction.PREV), middleOp, "PREV");
        Assert.assertEquals(pe2.getAdjacentOperator(PileupElement.Direction.NEXT), null, "NEXT");
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe2.isImmediatelyBefore(op), op.toString());
            Assert.assertEquals(pe2.isImmediatelyAfter(op), op == middleOp, op.toString());
        }
    }

    private void assertEqualsOperators(final List<CigarElement> elements, final List<CigarOperator> ops) {
        for (int i = 0; i < elements.size(); i++) {
            Assert.assertEquals(elements.get(i).getOperator(), ops.get(i), "elements doesn't have expected operator at position " + i);
        }
    }

    @Test
    public void testCreatePileupForReadAndOffset1() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
        final PileupElement pe0 = PileupElement.createPileupForReadAndOffset(read, 0);
        final List<CigarElement> betweenNextPosition = pe0.getBetweenNextPosition();
        Assert.assertTrue(betweenNextPosition.isEmpty());
        Assert.assertEquals(pe0.getOffset(), 0);

        final PileupElement pe10 = PileupElement.createPileupForReadAndOffset(read, 10);
        Assert.assertEquals(pe10.getOffset(), 10);
    }

    @Test
    public void testCreatePileupForReadAndOffset2() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M10D10M");
        final PileupElement pe10 = PileupElement.createPileupForReadAndOffset(read, 10);
        Assert.assertEquals(pe10.getOffset(), 10);

        final PileupElement pe15 = PileupElement.createPileupForReadAndOffset(read, 15);
        Assert.assertEquals(pe15.getOffset(), 15);
    }

    @Test
    public void testIsImmediatelyAfter_insideMs() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M10D10M");
        final PileupElement pe = PileupElement.createPileupForReadAndOffset(read, 5);
        Assert.assertFalse(pe.atStartOfCurrentCigar(), "atStartOfCurrentCigar");
        Assert.assertFalse(pe.atEndOfCurrentCigar(), "atEndOfCurrentCigar");
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe.isImmediatelyAfter(op), "isImmediatelyAfter " + op);
        }
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe.isImmediatelyBefore(op), "isImmediatelyBefore " + op);
        }
    }

    @Test
    public void testIsImmediatelyAfter_afterLastM() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M10D10M");
        final PileupElement pe = PileupElement.createPileupForReadAndOffset(read, 9);
        Assert.assertFalse(pe.atStartOfCurrentCigar(), "atStartOfCurrentCigar");
        Assert.assertTrue(pe.atEndOfCurrentCigar(), "atEndOfCurrentCigar");
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe.isImmediatelyAfter(op), "isImmediatelyAfter " + op);
        }

        Assert.assertTrue(pe.isImmediatelyBefore(D));
        for (final CigarOperator op : CigarOperator.values()) {
            if (op != D) {
                Assert.assertFalse(pe.isImmediatelyBefore(op), "isImmediatelyBefore " + op);
            }
        }
    }

    @Test
    public void testIsImmediatelyAfter_afterFirstD() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M10D10M");
        final PileupElement pe = PileupElement.createPileupForReadAndOffset(read, 10);
        Assert.assertFalse(pe.atEndOfCurrentCigar(), "atEndOfCurrentCigar");
        Assert.assertTrue(pe.atStartOfCurrentCigar(), "atStartOfCurrentCigar");
        Assert.assertTrue(pe.isImmediatelyAfter(D));
        for (final CigarOperator op : CigarOperator.values()) {
            if (op != D) {
                Assert.assertFalse(pe.isImmediatelyAfter(op), "isImmediatelyAfter " + op);
            }
        }

        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe.isImmediatelyBefore(op));
        }
    }

    @Test
    public void testIsImmediatelyAfter_insideFinalMs() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M10D10M");
        final PileupElement pe = PileupElement.createPileupForReadAndOffset(read, 15);
        Assert.assertFalse(pe.atStartOfCurrentCigar(), "atStartOfCurrentCigar");
        Assert.assertFalse(pe.atEndOfCurrentCigar(), "atEndOfCurrentCigar");
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe.isImmediatelyAfter(op), "isImmediatelyAfter " + op);
        }
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertFalse(pe.isImmediatelyBefore(op), "isImmediatelyBefore " + op);
        }
    }

    @DataProvider(name = "adjacentElementTestData")
    private Object[][] adjacentElementTestData() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M" + "10D" + "10N" + "10P" + "10=" + "10X" + "10M" + "10S" + "10H");
        return new Object[][]{
                {read, 0, true, false, null, D},
                {read, 1, false, false, null, D},
                {read, 2, false, false, null, D},
                {read, 3, false, false, null, D},
                {read, 4, false, false, null, D},
                {read, 5, false, false, null, D},
                {read, 6, false, false, null, D},
                {read, 7, false, false, null, D},
                {read, 8, false, false, null, D},
                {read, 9, false, true, null, D},
                {read, 10, true, false, P, X},
                {read, 11, false, false, P, X},
                {read, 12, false, false, P, X},
                {read, 13, false, false, P, X},
                {read, 14, false, false, P, X},
                {read, 15, false, false, P, X},
                {read, 16, false, false, P, X},
                {read, 17, false, false, P, X},
                {read, 18, false, false, P, X},
                {read, 19, false, true, P, X},
                {read, 20, true, false, EQ, M},
                {read, 21, false, false, EQ, M},
                {read, 22, false, false, EQ, M},
                {read, 23, false, false, EQ, M},
                {read, 24, false, false, EQ, M},
                {read, 25, false, false, EQ, M},
                {read, 26, false, false, EQ, M},
                {read, 27, false, false, EQ, M},
                {read, 28, false, false, EQ, M},
                {read, 29, false, true, EQ, M},
                {read, 30, true, false, X, S},
                {read, 31, false, false, X, S},
                {read, 32, false, false, X, S},
                {read, 33, false, false, X, S},
                {read, 34, false, false, X, S},
                {read, 35, false, false, X, S},
                {read, 36, false, false, X, S},
                {read, 37, false, false, X, S},
                {read, 38, false, false, X, S},
                {read, 39, false, true, X, S},
        };
    }

    @Test(dataProvider = "adjacentElementTestData")
    public void testAdjacentElements(GATKRead read, int goodOffset, boolean atStartOfCurrentCigar, boolean atEndOfCurrentCigar, CigarOperator prev, CigarOperator next) throws Exception {
        final PileupElement pe = PileupElement.createPileupForReadAndOffset(read, goodOffset);
        Assert.assertEquals(pe.atStartOfCurrentCigar(), atStartOfCurrentCigar, "atStartOfCurrentCigar");
        Assert.assertEquals(pe.atEndOfCurrentCigar(), atEndOfCurrentCigar, "atEndOfCurrentCigar");
        Assert.assertEquals(pe.getAdjacentOperator(PileupElement.Direction.PREV), prev, "prev");
        Assert.assertEquals(pe.getAdjacentOperator(PileupElement.Direction.NEXT), next, "next");
        for (final CigarOperator op : CigarOperator.values()) {
            Assert.assertEquals(pe.isImmediatelyAfter(op), atStartOfCurrentCigar && op == prev);
            Assert.assertEquals(pe.isImmediatelyBefore(op), atEndOfCurrentCigar && op == next);
        }
    }

    @DataProvider(name = "elementData_badOffsetWithinBounds")
    private Object[][] elementData_badOffsetWithinBounds() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M" + "10D" + "10N" + "10P" + "10=" + "10X" + "10M" + "10S" + "10H");
        return new Object[][]{
                {read, 40},
                {read, 41},
                {read, 42},
                {read, 43},
                {read, 44},
                {read, 45},
                {read, 46},
                {read, 47},
                {read, 48},
                {read, 49},
        };
    }

    @Test(dataProvider = "elementData_badOffsetWithinBounds", expectedExceptions = IllegalStateException.class)
    public void testBadOffsetWithinBounds(GATKRead read, int badOffset) throws Exception {
        PileupElement.createPileupForReadAndOffset(read, badOffset);
    }

    @DataProvider(name = "elementData_badOffsetOutOfBounds")
    private Object[][] elementData_badOffsetOutOfBounds() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M" + "10D" + "10N" + "10P" + "10=" + "10X" + "10M" + "10S" + "10H");
        return new Object[][]{
                {read, 50},     //these are not even within the span of the read
                {read, 1000},
        };
    }

    @Test(dataProvider = "elementData_badOffsetOutOfBounds", expectedExceptions = IllegalStateException.class)
    public void testBadOffsetOutOfBounds(GATKRead read, int badOffset) throws Exception {
        PileupElement.createPileupForReadAndOffset(read, badOffset);
    }
}
