package org.broadinstitute.hellbender.utils.locusiterator;

import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class AlignmentStateMachineUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "AlignmentStateMachineTest")
    public Object[][] makeAlignmentStateMachineTest() {
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "AlignmentStateMachineTest")
    public void testAlignmentStateMachineTest(LIBSTest params) {
        final GATKRead read = params.makeRead();
        final AlignmentStateMachine state = new AlignmentStateMachine(read);
        final LIBS_position tester = new LIBS_position(read);

        // min is one because always visit something, even for 10I reads
        final int expectedBpToVisit = read.getEnd() - read.getStart() + 1;

        Assert.assertSame(state.getRead(), read);
        Assert.assertNotNull(state.toString());

        int bpVisited = 0;
        int lastOffset = -1;

        // TODO -- more tests about test state machine state before first step?
        Assert.assertTrue(state.isLeftEdge());
        Assert.assertNull(state.getCigarOperator());
        Assert.assertNotNull(state.toString());
        Assert.assertEquals(state.getReadOffset(), -1);
        Assert.assertEquals(state.getGenomeOffset(), -1);
        Assert.assertEquals(state.getCurrentCigarElementOffset(), -1);
        Assert.assertEquals(state.getCurrentCigarElement(), null);

        while ( state.stepForwardOnGenome() != null ) {
            Assert.assertNotNull(state.toString());

            tester.stepForwardOnGenome();

            Assert.assertTrue(state.getReadOffset() >= lastOffset, "Somehow read offsets are decreasing: lastOffset " + lastOffset + " current " + state.getReadOffset());
            Assert.assertEquals(state.getReadOffset(), tester.getCurrentReadOffset(), "Read offsets are wrong at " + bpVisited);

            Assert.assertFalse(state.isLeftEdge());

            Assert.assertEquals(state.getCurrentCigarElement(), read.getCigar().getCigarElement(tester.currentOperatorIndex), "CigarElement index failure");
            Assert.assertEquals(state.getOffsetIntoCurrentCigarElement(), tester.getCurrentPositionOnOperatorBase0(), "CigarElement index failure");

            Assert.assertEquals(read.getCigar().getCigarElement(state.getCurrentCigarElementOffset()), state.getCurrentCigarElement(), "Current cigar element isn't what we'd get from the read itself");

            Assert.assertTrue(state.getOffsetIntoCurrentCigarElement() >= 0, "Offset into current cigar too small");
            Assert.assertTrue(state.getOffsetIntoCurrentCigarElement() < state.getCurrentCigarElement().getLength(), "Offset into current cigar too big");

            Assert.assertEquals(state.getGenomeOffset(), tester.getCurrentGenomeOffsetBase0(), "Offset from alignment start is bad");
            Assert.assertEquals(state.getGenomePosition(), tester.getCurrentGenomeOffsetBase0() + read.getStart(), "GenomePosition start is bad");
            Assert.assertEquals(state.getLocation().size(), 1, "GenomeLoc position should have size == 1");
            Assert.assertEquals(state.getLocation().getStart(), state.getGenomePosition(), "GenomeLoc position is bad");
            // most tests of this functionality are in LIBS
            Assert.assertNotNull(state.makePileupElement());

            lastOffset = state.getReadOffset();
            bpVisited++;
        }

        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
        Assert.assertEquals(state.getReadOffset(), read.getLength());
        Assert.assertEquals(state.getCurrentCigarElementOffset(), read.numCigarElements());
        Assert.assertEquals(state.getCurrentCigarElement(), null);
        Assert.assertNotNull(state.toString());
    }
}
