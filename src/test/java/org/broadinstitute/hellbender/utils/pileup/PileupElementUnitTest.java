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

package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.locusiterator.LIBSTest;
import org.broadinstitute.hellbender.utils.locusiterator.LIBS_position;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByStateBaseTest;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class PileupElementUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "PileupElementTest")
    public Object[][] makePileupElementTest() {
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "PileupElementTest")
    public void testPileupElementTest(LIBSTest params) {
        final SAMRecord read = params.makeRead();
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



            Assert.assertEquals(pe.atEndOfCurrentCigar(), state.getOffsetIntoCurrentCigarElement() == state.getCurrentCigarElement().getLength() - 1, "atEndOfCurrentCigar failed");
            Assert.assertEquals(pe.atStartOfCurrentCigar(), state.getOffsetIntoCurrentCigarElement() == 0, "atStartOfCurrentCigar failed");

            Assert.assertEquals(pe.getBase(), pe.isDeletion() ? PileupElement.DELETION_BASE : read.getReadBases()[state.getReadOffset()]);
            Assert.assertEquals(pe.getQual(), pe.isDeletion() ? PileupElement.DELETION_QUAL : read.getBaseQualities()[state.getReadOffset()]);

            Assert.assertEquals(pe.getCurrentCigarElement(), state.getCurrentCigarElement());
            Assert.assertEquals(pe.getCurrentCigarOffset(), state.getCurrentCigarElementOffset());

            // tested in libs
            //pe.getLengthOfImmediatelyFollowingIndel();
            //pe.getBasesOfImmediatelyFollowingInsertion();

            // Don't test -- pe.getBaseIndex();
            if ( pe.atEndOfCurrentCigar() && state.getCurrentCigarElementOffset() < read.getCigarLength() - 1 ) {
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
        final List<Object[]> tests = new LinkedList<Object[]>();

        final List<CigarOperator> operators = Arrays.asList(CigarOperator.I, CigarOperator.P, CigarOperator.S);

        for ( final CigarOperator firstOp : Arrays.asList(CigarOperator.M) ) {
            for ( final CigarOperator lastOp : Arrays.asList(CigarOperator.M, CigarOperator.D) ) {
                for ( final int nIntermediate : Arrays.asList(1, 2, 3) ) {
                    for ( final List<CigarOperator> combination : Utils.makePermutations(operators, nIntermediate, false) ) {
                        final int readLength = 2 + combination.size();
                        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, 1, readLength);
                        read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
                        read.setBaseQualities(Utils.dupBytes((byte) 30, readLength));

                        String cigar = "1" + firstOp;
                        for ( final CigarOperator op : combination ) cigar += "1" + op;
                        cigar += "1" + lastOp;
                        read.setCigarString(cigar);

                        tests.add(new Object[]{read, firstOp, lastOp, combination});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PrevAndNextTest")
    public void testPrevAndNextTest(final SAMRecord read, final CigarOperator firstOp, final CigarOperator lastOp, final List<CigarOperator> ops) {
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
}
