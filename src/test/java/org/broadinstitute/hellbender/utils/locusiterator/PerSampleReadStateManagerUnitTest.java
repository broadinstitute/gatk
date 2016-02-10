/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.locusiterator;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class PerSampleReadStateManagerUnitTest extends LocusIteratorByStateBaseTest {
    private class PerSampleReadStateManagerTest extends TestDataProvider {
        private List<Integer> readCountsPerAlignmentStart;
        private List<SAMRecord> reads;
        private List<ArrayList<AlignmentStateMachine>> recordStatesByAlignmentStart;
        private int removalInterval;

        public PerSampleReadStateManagerTest( List<Integer> readCountsPerAlignmentStart, int removalInterval ) {
            super(PerSampleReadStateManagerTest.class);

            this.readCountsPerAlignmentStart = readCountsPerAlignmentStart;
            this.removalInterval = removalInterval;

            reads = new ArrayList<SAMRecord>();
            recordStatesByAlignmentStart = new ArrayList<ArrayList<AlignmentStateMachine>>();

            setName(String.format("%s: readCountsPerAlignmentStart: %s  removalInterval: %d",
                    getClass().getSimpleName(), readCountsPerAlignmentStart, removalInterval));
        }

        public void run() {
            PerSampleReadStateManager perSampleReadStateManager = new PerSampleReadStateManager(LocusIteratorByState.NO_DOWNSAMPLING);

            makeReads();

            for ( ArrayList<AlignmentStateMachine> stackRecordStates : recordStatesByAlignmentStart ) {
                perSampleReadStateManager.addStatesAtNextAlignmentStart(new LinkedList<AlignmentStateMachine>(stackRecordStates));
            }

            // read state manager should have the right number of reads
            Assert.assertEquals(reads.size(), perSampleReadStateManager.size());

            Iterator<SAMRecord> originalReadsIterator = reads.iterator();
            Iterator<AlignmentStateMachine> recordStateIterator = perSampleReadStateManager.iterator();
            int recordStateCount = 0;
            int numReadStatesRemoved = 0;

            // Do a first-pass validation of the record state iteration by making sure we get back everything we
            // put in, in the same order, doing any requested removals of read states along the way
            while ( recordStateIterator.hasNext() ) {
                AlignmentStateMachine readState = recordStateIterator.next();
                recordStateCount++;
                SAMRecord readFromPerSampleReadStateManager = readState.getRead();

                Assert.assertTrue(originalReadsIterator.hasNext());
                SAMRecord originalRead = originalReadsIterator.next();

                // The read we get back should be literally the same read in memory as we put in
                Assert.assertTrue(originalRead == readFromPerSampleReadStateManager);

                // If requested, remove a read state every removalInterval states
                if ( removalInterval > 0 && recordStateCount % removalInterval == 0 ) {
                    recordStateIterator.remove();
                    numReadStatesRemoved++;
                }
            }

            Assert.assertFalse(originalReadsIterator.hasNext());

            // If we removed any read states, do a second pass through the read states to make sure the right
            // states were removed
            if ( numReadStatesRemoved > 0 ) {
                Assert.assertEquals(perSampleReadStateManager.size(), reads.size() - numReadStatesRemoved);

                originalReadsIterator = reads.iterator();
                recordStateIterator = perSampleReadStateManager.iterator();
                int readCount = 0;
                int readStateCount = 0;

                // Match record states with the reads that should remain after removal
                while ( recordStateIterator.hasNext() ) {
                    AlignmentStateMachine readState = recordStateIterator.next();
                    readStateCount++;
                    SAMRecord readFromPerSampleReadStateManager = readState.getRead();

                    Assert.assertTrue(originalReadsIterator.hasNext());

                    SAMRecord originalRead = originalReadsIterator.next();
                    readCount++;

                    if ( readCount % removalInterval == 0 ) {
                        originalRead = originalReadsIterator.next(); // advance to next read, since the previous one should have been discarded
                        readCount++;
                    }

                    // The read we get back should be literally the same read in memory as we put in (after accounting for removals)
                    Assert.assertTrue(originalRead == readFromPerSampleReadStateManager);
                }

                Assert.assertEquals(readStateCount, reads.size() - numReadStatesRemoved);
            }

            // Allow memory used by this test to be reclaimed
            readCountsPerAlignmentStart = null;
            reads = null;
            recordStatesByAlignmentStart = null;
        }

        private void makeReads() {
            int alignmentStart = 1;

            for ( int readsThisStack : readCountsPerAlignmentStart ) {
                ArrayList<GATKSAMRecord> stackReads = new ArrayList<GATKSAMRecord>(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(readsThisStack, header, "foo", 0, alignmentStart, MathUtils.randomIntegerInRange(50, 100)));
                ArrayList<AlignmentStateMachine> stackRecordStates = new ArrayList<AlignmentStateMachine>();

                for ( GATKSAMRecord read : stackReads ) {
                    stackRecordStates.add(new AlignmentStateMachine(read));
                }

                reads.addAll(stackReads);
                recordStatesByAlignmentStart.add(stackRecordStates);
            }
        }
    }

    @DataProvider(name = "PerSampleReadStateManagerTestDataProvider")
    public Object[][] createPerSampleReadStateManagerTests() {
        for ( List<Integer> thisTestReadStateCounts : Arrays.asList( Arrays.asList(1),
                Arrays.asList(2),
                Arrays.asList(10),
                Arrays.asList(1, 1),
                Arrays.asList(2, 2),
                Arrays.asList(10, 10),
                Arrays.asList(1, 10),
                Arrays.asList(10, 1),
                Arrays.asList(1, 1, 1),
                Arrays.asList(2, 2, 2),
                Arrays.asList(10, 10, 10),
                Arrays.asList(1, 1, 1, 1, 1, 1),
                Arrays.asList(10, 10, 10, 10, 10, 10),
                Arrays.asList(1, 2, 10, 1, 2, 10)
        ) ) {

            for ( int removalInterval : Arrays.asList(0, 2, 3) ) {
                new PerSampleReadStateManagerTest(thisTestReadStateCounts, removalInterval);
            }
        }

        return PerSampleReadStateManagerTest.getTests(PerSampleReadStateManagerTest.class);
    }

    @Test(dataProvider = "PerSampleReadStateManagerTestDataProvider")
    public void runPerSampleReadStateManagerTest( PerSampleReadStateManagerTest test ) {
        logger.warn("Running test: " + test);

        test.run();
    }
}
