package org.broadinstitute.hellbender.utils.locusiterator;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public final class PerSampleReadStateManagerUnitTest extends LocusIteratorByStateBaseTest {

    private static final class PerSampleReadStateManagerTester extends TestDataProvider {
        private final List<Integer> readCountsPerAlignmentStart;
        private final int removalInterval;

        public PerSampleReadStateManagerTester(final List<Integer> readCountsPerAlignmentStart, final int removalInterval ) {
            super(PerSampleReadStateManagerTester.class);

            this.readCountsPerAlignmentStart = readCountsPerAlignmentStart;
            this.removalInterval = removalInterval;

            setName(String.format("%s: readCountsPerAlignmentStart: %s  removalInterval: %d",
                    getClass().getSimpleName(), readCountsPerAlignmentStart, removalInterval));
        }

        public void run() {
            final List<GATKRead> reads = new ArrayList<>();
            final List<ArrayList<AlignmentStateMachine>> recordStatesByAlignmentStart = new ArrayList<>();
            makeReads(reads, recordStatesByAlignmentStart);

            final PerSampleReadStateManager perSampleReadStateManager = new PerSampleReadStateManager(LocusIteratorByState.NO_DOWNSAMPLING);
            recordStatesByAlignmentStart.stream().map(LinkedList<AlignmentStateMachine>::new).forEach(perSampleReadStateManager::addStatesAtNextAlignmentStart);

            Assert.assertEquals(reads.size(), perSampleReadStateManager.size());

            final Iterator<GATKRead> originalReadsIterator = reads.iterator();
            final Iterator<AlignmentStateMachine> stateIterator = perSampleReadStateManager.iterator();
            int recordStateIndex = 0;
            int numReadStatesRemoved = 0;

            // Do a first-pass validation of the record state iteration by making sure we get back everything we
            // put in, in the same order, doing any requested removals of read states along the way
            while ( stateIterator.hasNext() ) {
                // The read we get back should be literally the same read in memory as we put in
                Assert.assertTrue(originalReadsIterator.hasNext());
                Assert.assertTrue(originalReadsIterator.next() == stateIterator.next().getRead());

                // If requested, remove a read state every removalInterval states
                if ( removalInterval > 0 && recordStateIndex % removalInterval == 0 ) {
                    stateIterator.remove();
                    numReadStatesRemoved++;
                }
                recordStateIndex++;
            }

            Assert.assertFalse(originalReadsIterator.hasNext());
            Assert.assertEquals(perSampleReadStateManager.size(), reads.size() - numReadStatesRemoved);

            // second pass through the read states to make sure the right states were removed
            if ( numReadStatesRemoved > 0 ) {

                Iterator<AlignmentStateMachine> secondPassStateIterator = perSampleReadStateManager.iterator();
                for ( int readIndex = 0; readIndex < reads.size(); readIndex++ ) {
                    // skip reads that were removed in the first pass
                    if ( readIndex % removalInterval == 0 ) {
                        continue;
                    }
                    Assert.assertTrue(secondPassStateIterator.hasNext());
                    Assert.assertTrue(reads.get(readIndex) == secondPassStateIterator.next().getRead());
                }
                Assert.assertFalse(secondPassStateIterator.hasNext());
            }
        }

        private void makeReads(List<GATKRead> reads, List<ArrayList<AlignmentStateMachine>> recordStatesByAlignmentStart) {
            final int alignmentStart = 1;
            for ( int readsThisStack : readCountsPerAlignmentStart ) {
                final int readLength = Utils.getRandomGenerator().nextInt(51) + 50;
                ArrayList<GATKRead> stackReads = new ArrayList<>(ArtificialReadUtils.createIdenticalArtificialReads(readsThisStack, header, "foo", 0, alignmentStart, readLength));
                ArrayList<AlignmentStateMachine> stackRecordStates = new ArrayList<>();

                for ( GATKRead read : stackReads ) {
                    stackRecordStates.add(new AlignmentStateMachine(read));
                }

                reads.addAll(stackReads);
                recordStatesByAlignmentStart.add(stackRecordStates);
            }
        }
    }

    @DataProvider(name = "PerSampleReadStateManagerTestDataProvider")
    public Object[][] createPerSampleReadStateManagerTests() {
        final List<List<Integer>> testReadStateCounts = Arrays.asList(Arrays.asList(1),
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
        );

        final List<Integer> removeIntervals = Arrays.asList(0, 2, 3);

        final List<Object[]> tests = testReadStateCounts.stream()
                .flatMap(counts -> removeIntervals.stream()
                        .map(interval -> new Object[] {new PerSampleReadStateManagerTester(counts, interval)}))
                .collect(Collectors.toList());

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PerSampleReadStateManagerTestDataProvider")
    public void runPerSampleReadStateManagerTest( PerSampleReadStateManagerTester test ) {
        logger.warn("Running test: " + test);
        Utils.resetRandomGenerator();
        test.run();
    }
}
