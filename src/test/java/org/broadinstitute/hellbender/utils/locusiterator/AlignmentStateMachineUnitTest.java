/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class AlignmentStateMachineUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "AlignmentStateMachineTest")
    public Object[][] makeAlignmentStateMachineTest() {
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "AlignmentStateMachineTest")
    public void testAlignmentStateMachineTest(LIBSTest params) {
        final SAMRecord read = params.makeRead();
        final AlignmentStateMachine state = new AlignmentStateMachine(read);
        final LIBS_position tester = new LIBS_position(read);

        // min is one because always visit something, even for 10I reads
        final int expectedBpToVisit = read.getAlignmentEnd() - read.getAlignmentStart() + 1;

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
            Assert.assertEquals(state.getGenomePosition(), tester.getCurrentGenomeOffsetBase0() + read.getAlignmentStart(), "GenomePosition start is bad");
            Assert.assertEquals(state.getLocation(genomeLocParser).size(), 1, "GenomeLoc position should have size == 1");
            Assert.assertEquals(state.getLocation(genomeLocParser).getStart(), state.getGenomePosition(), "GenomeLoc position is bad");
            // most tests of this functionality are in LIBS
            Assert.assertNotNull(state.makePileupElement());

            lastOffset = state.getReadOffset();
            bpVisited++;
        }

        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
        Assert.assertEquals(state.getReadOffset(), read.getReadLength());
        Assert.assertEquals(state.getCurrentCigarElementOffset(), read.getCigarLength());
        Assert.assertEquals(state.getCurrentCigarElement(), null);
        Assert.assertNotNull(state.toString());
    }
}
