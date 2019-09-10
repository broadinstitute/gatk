package org.broadinstitute.hellbender.utils.locusiterator;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public final class PerSampleReadStateManagerUnitTest extends LocusIteratorByStateBaseTest {

    private static final class PerSampleReadStateManagerTester extends TestDataProvider {
        private final List<Integer> readCountsPerAlignmentStart;
        private final List<Integer> readLengths;
        private final int removalInterval;

        public PerSampleReadStateManagerTester(final List<Integer> readCountsPerAlignmentStart, final List<Integer> readLengths, final int removalInterval ) {
            super(PerSampleReadStateManagerTester.class);

            this.readCountsPerAlignmentStart = readCountsPerAlignmentStart;
            this.readLengths = readLengths;
            this.removalInterval = removalInterval;

            setName(String.format("%s: readCountsPerAlignmentStart: %s  removalInterval: %d",
                    getClass().getSimpleName(), readCountsPerAlignmentStart, removalInterval));
        }

        public void run() {
            final List<Collection<GATKRead>> readStacks = makeReadStacks();
            final List<GATKRead> reads = readStacks.stream().flatMap(Collection::stream).collect(Collectors.toList());
            final List<List<AlignmentStateMachine>> recordStatesByAlignmentStart = readStacks.stream()
                    .map(stackReads -> stackReads.stream().map(AlignmentStateMachine::new).collect(Collectors.toList()))
                    .collect(Collectors.toList());

            final PerSampleReadStateManager perSampleReadStateManager = new PerSampleReadStateManager(LocusIteratorByState.NO_DOWNSAMPLING);
            recordStatesByAlignmentStart.stream().map(LinkedList<AlignmentStateMachine>::new).forEach(perSampleReadStateManager::addStatesAtNextAlignmentStart);

            Assert.assertEquals(reads.size(), perSampleReadStateManager.size());

            final Iterator<AlignmentStateMachine> stateIterator = perSampleReadStateManager.iterator();
            int numReadStatesRemoved = 0;

            // Do a first-pass validation of the record state iteration by making sure we get back everything we
            // put in, in the same order, doing any requested removals of read states along the way
            for ( int readIndex = 0; readIndex < reads.size(); readIndex++ ) {
                // The read we get back should be literally the same read in memory as we put in
                Assert.assertTrue(stateIterator.hasNext());
                final GATKRead stateIteratorRead = stateIterator.next().getRead();
                Assert.assertTrue(reads.get(readIndex) == stateIteratorRead);

                // If requested, remove a read state every removalInterval states
                if ( removalInterval > 0 && readIndex % removalInterval == 0 ) {
                    stateIterator.remove();
                    numReadStatesRemoved++;
                }
            }

            Assert.assertFalse(stateIterator.hasNext());
            Assert.assertEquals(perSampleReadStateManager.size(), reads.size() - numReadStatesRemoved);

            // second pass through the read states to make sure the right states were removed
            if ( numReadStatesRemoved > 0 ) {

                final Iterator<AlignmentStateMachine> secondPassStateIterator = perSampleReadStateManager.iterator();
                for ( int readIndex = 0; readIndex < reads.size(); readIndex++ ) {
                    // skip reads that were removed in the first pass
                    if ( readIndex % removalInterval == 0 ) {
                        continue;
                    }
                    Assert.assertTrue(secondPassStateIterator.hasNext());
                    final GATKRead stateIteratorRead = secondPassStateIterator.next().getRead();
                    Assert.assertTrue(reads.get(readIndex) == stateIteratorRead);
                }
                Assert.assertFalse(secondPassStateIterator.hasNext());
            }
        }

        private List<Collection<GATKRead>> makeReadStacks() {
            final int alignmentStart = 1;
            return IntStream.range(0, readCountsPerAlignmentStart.size())
                    .mapToObj(n -> ArtificialReadUtils.createIdenticalArtificialReads(readCountsPerAlignmentStart.get(n), header, "foo", 0, alignmentStart, readLengths.get(n)))
                    .collect(Collectors.toList());
        }
    }

    @DataProvider(name = "PerSampleReadStateManagerTestDataProvider")
    public Object[][] createPerSampleReadStateManagerTests() {
        Utils.resetRandomGenerator();
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

        // each alignment start has reads of all the same length.  Thus each instance of readLengths must be a list of integers
        // with the same size as the list of read state counts.
        final List<List<Integer>> readLengths = testReadStateCounts.stream()
                .map(counts -> counts.stream().map(c -> Utils.getRandomGenerator().nextInt(51) + 50).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final List<Integer> removalIntervals = Arrays.asList(0, 2, 3);

        final List<Object[]> tests = IntStream.range(0, testReadStateCounts.size()).boxed()
                .flatMap(n -> removalIntervals.stream()
                        .map(interval -> new Object[] {new PerSampleReadStateManagerTester(testReadStateCounts.get(n), readLengths.get(n), interval)}))
                .collect(Collectors.toList());

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PerSampleReadStateManagerTestDataProvider")
    public void runPerSampleReadStateManagerTest( PerSampleReadStateManagerTester test ) {
        test.run();
    }
}
