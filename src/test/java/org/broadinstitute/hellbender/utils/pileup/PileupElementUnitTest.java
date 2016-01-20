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

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

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

        while ( state.stepForwardOnGenome() != null ) {
            tester.stepForwardOnGenome();
            final PileupElement pe = state.makePileupElement();

            Assert.assertEquals(pe.getRead(), read);
            Assert.assertEquals(pe.getMappingQual(), read.getMappingQuality());
            Assert.assertEquals(pe.getOffset(), state.getReadOffset());

            Assert.assertEquals(pe.isDeletion(), state.getCigarOperator() == CigarOperator.D);
            Assert.assertEquals(pe.isAfterInsertion(), tester.isAfterInsertion);
            Assert.assertEquals(pe.isBeforeInsertion(), tester.isBeforeInsertion);
            Assert.assertEquals(pe.isNextToSoftClip(), tester.isNextToSoftClip);

            if ( ! hasNeighboringPaddedOps(params.getElements(), pe.getCurrentCigarOffset()) ) {
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
            if ( pe.atEndOfCurrentCigar() && state.getCurrentCigarElementOffset() < read.getCigar().numCigarElements() - 1 ) {
                final CigarElement nextElement = read.getCigar().getCigarElement(state.getCurrentCigarElementOffset() + 1);
                if ( nextElement.getOperator() == CigarOperator.I ) {
                    Assert.assertTrue(pe.getBetweenNextPosition().size() >= 1);
                    Assert.assertEquals(pe.getBetweenNextPosition().get(0), nextElement);
                }
                if ( nextElement.getOperator() == CigarOperator.M ) {
                    Assert.assertTrue(pe.getBetweenNextPosition().isEmpty());
                }
            } else {
                Assert.assertTrue(pe.getBetweenNextPosition().isEmpty());
            }

            if ( pe.atStartOfCurrentCigar() && state.getCurrentCigarElementOffset() > 0 ) {
                final CigarElement prevElement = read.getCigar().getCigarElement(state.getCurrentCigarElementOffset() - 1);
                if ( prevElement.getOperator() == CigarOperator.I ) {
                    Assert.assertTrue(pe.getBetweenPrevPosition().size() >= 1);
                    Assert.assertEquals(pe.getBetweenPrevPosition().getLast(), prevElement);
                }
                if ( prevElement.getOperator() == CigarOperator.M ) {
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

        for ( final CigarOperator firstOp : Arrays.asList(CigarOperator.M) ) {
            for ( final CigarOperator lastOp : Arrays.asList(CigarOperator.M, CigarOperator.D) ) {
                for ( final int nIntermediate : Arrays.asList(1, 2, 3) ) {
                    for ( final List<CigarOperator> combination : Utils.makePermutations(operators, nIntermediate, false) ) {
                        final int readLength = 2 + combination.size();
                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 1, readLength);
                        read.setBases(Utils.dupBytes((byte) 'A', readLength));
                        read.setBaseQualities(Utils.dupBytes((byte) 30, readLength));

                        String cigar = "1" + firstOp;
                        for ( final CigarOperator op : combination ) {
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

        state.stepForwardOnGenome();
        final PileupElement pe = state.makePileupElement();
        Assert.assertEquals(pe.getBetweenNextPosition().size(), ops.size());
        Assert.assertEquals(pe.getBetweenPrevPosition().size(), 0);
        assertEqualsOperators(pe.getBetweenNextPosition(), ops);
        Assert.assertEquals(pe.getPreviousOnGenomeCigarElement(), null);
        Assert.assertNotNull(pe.getNextOnGenomeCigarElement());
        Assert.assertEquals(pe.getNextOnGenomeCigarElement().getOperator(), lastOp);

        state.stepForwardOnGenome();
        final PileupElement pe2 = state.makePileupElement();
        Assert.assertEquals(pe2.getBetweenPrevPosition().size(), ops.size());
        Assert.assertEquals(pe2.getBetweenNextPosition().size(), 0);
        assertEqualsOperators(pe2.getBetweenPrevPosition(), ops);
        Assert.assertNotNull(pe2.getPreviousOnGenomeCigarElement());
        Assert.assertEquals(pe2.getPreviousOnGenomeCigarElement().getOperator(), firstOp);
        Assert.assertEquals(pe2.getNextOnGenomeCigarElement(), null);
    }

    private void assertEqualsOperators(final List<CigarElement> elements, final List<CigarOperator> ops) {
        for ( int i = 0; i < elements.size(); i++ ) {
            Assert.assertEquals(elements.get(i).getOperator(), ops.get(i), "elements doesn't have expected operator at position " + i);
        }
    }

    @Test
    public void testCreatePileupForReadAndOffset1(){
        final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
        final PileupElement pe0 = PileupElement.createPileupForReadAndOffset(read, 0);
        final LinkedList<CigarElement> betweenNextPosition = pe0.getBetweenNextPosition();
        Assert.assertTrue(betweenNextPosition.isEmpty());
        Assert.assertEquals(pe0.getOffset(), 0);

        final PileupElement pe10 = PileupElement.createPileupForReadAndOffset(read, 10);
        Assert.assertEquals(pe10.getOffset(), 10);
    }

    @Test
    public void testCreatePileupForReadAndOffset2(){
        final GATKRead read = ArtificialReadUtils.createArtificialRead("10M10D10M");
        final PileupElement pe10 = PileupElement.createPileupForReadAndOffset(read, 10);
        Assert.assertEquals(pe10.getOffset(), 10);

        final PileupElement pe15 = PileupElement.createPileupForReadAndOffset(read, 15);
        Assert.assertEquals(pe15.getOffset(), 15);

    }
}
