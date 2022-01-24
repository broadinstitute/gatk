package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Manages and updates mapping from sample -> Iterable<AlignmentStateMachine>
 */
final class ReadStateManager implements Iterable<Map.Entry<String, PerSampleReadStateManager>> {
    private final List<String> samples;
    private final PeekableIterator<GATKRead> iterator;
    private final SamplePartitioner samplePartitioner;

    /**
     * A mapping from sample name -> the per sample read state manager that manages
     *
     * IT IS CRITICAL THAT THIS BE A LINKED HASH MAP, SO THAT THE ITERATION OF THE MAP OCCURS IN THE SAME
     * ORDER AS THE ORIGINAL SAMPLES
     */
    private final Map<String, PerSampleReadStateManager> readStatesBySample = new LinkedHashMap<>();

    private int totalReadStates = 0;

    public ReadStateManager(final Iterator<GATKRead> source,
                            final List<String> samples,
                            final LIBSDownsamplingInfo info,
                            final SAMFileHeader header) {
        Utils.nonNull(source, "source");
        Utils.nonNull(samples, "samples");
        Utils.nonNull(info, "downsampling info");
        Utils.nonNull(header, "header");
        this.samples = samples;
        this.iterator = new PeekableIterator<>(source);

        for (final String sample : samples) {
            // because this is a linked hash map the order of iteration will be in sample order
            readStatesBySample.put(sample, new PerSampleReadStateManager(info));
        }

        samplePartitioner = new SamplePartitioner(info, samples, header);
    }

    /**
     * Returns a iterator over all the sample -> per-sample read state managers with each sample in this read state manager.
     *
     * The order of iteration is the same as the order of the samples provided upon construction to this
     * ReadStateManager.
     *
     * @return Iterator over sample + per sample read state manager pairs for this read state manager.
     */
    @Override
    public Iterator<Map.Entry<String, PerSampleReadStateManager>> iterator() {
        return readStatesBySample.entrySet().iterator();
    }

    public boolean isEmpty() {
        return totalReadStates == 0;
    }

    /**
     * Retrieves the total number of reads in the manager across all samples.
     *
     * @return Total number of reads over all samples.
     */
    public int size() {
        return totalReadStates;
    }

    /**
     * Retrieves the total number of reads in the manager in the given sample.
     *
     * @param sample The sample.
     * @return Total number of reads in the given sample.
     */
    public int size(final String sample) {
        Utils.nonNull(sample);
        return readStatesBySample.get(sample).size();
    }

    public AlignmentStateMachine getFirst() {
        for ( final PerSampleReadStateManager manager : readStatesBySample.values() ) {
            if ( ! manager.isEmpty() ) {
                return manager.getFirst();
            }
        }
        return null;
    }

    public boolean hasNext() {
        return totalReadStates > 0 || iterator.hasNext();
    }

    /**
     * Advances all of the read states by one bp.  After this call the read states are reflective
     * of the next pileup.
     */
    public void updateReadStates() {
        for (final PerSampleReadStateManager perSampleReadStateManager : readStatesBySample.values() ) {
            totalReadStates -= perSampleReadStateManager.updateReadStates();
        }
    }

    /**
     * Does read start at the same position as described by currentContextIndex and currentAlignmentStart?
     *
     * @param read the read we want to test
     * @param currentContig the contig of the reads in this state manager
     * @param currentAlignmentStart the alignment start of the of the left-most position on the
     *                           genome of the reads in this read state manager
     * @return true if read has contig index and start equal to the current ones
     */
    private boolean readStartsAtCurrentPosition(final GATKRead read, final String currentContig, final int currentAlignmentStart) {
        return read.getStart() == currentAlignmentStart && read.getContig().equals(currentContig);
    }

    /**
     * Pull all of the reads off the iterator that overlap the left-most position among all
     * reads this ReadStateManager
     */
    public void collectPendingReads() {
        if (!iterator.hasNext()) {
            return;
        }

        // determine the left-most boundary that determines which reads to keep in this new pileup
        final String firstContig;
        final int firstAlignmentStart;
        if ( isEmpty() ) {
            // there are no reads here, so our next state is the next read in the stream
            firstContig = iterator.peek().getContig();
            firstAlignmentStart = iterator.peek().getStart();
        } else {
            // there's a read in the system, so it's our targeted first read
            final AlignmentStateMachine firstState = getFirst();
            firstContig = firstState.getContig();
            // note this isn't the alignment start of the read, but rather the alignment start position
            firstAlignmentStart = firstState.getGenomePosition();
        }

        while ( iterator.hasNext() && readStartsAtCurrentPosition(iterator.peek(), firstContig, firstAlignmentStart) ) {
            submitRead(iterator.next());
        }

        samplePartitioner.doneSubmittingReads();

        for (final String sample : samples) {
            final Collection<GATKRead> newReads = samplePartitioner.getReadsForSample(sample);

            final PerSampleReadStateManager statesBySample = readStatesBySample.get(sample);
            addReadsToSample(statesBySample, newReads);
        }

        samplePartitioner.reset();
    }

    /**
     * Add a read to the sample partitioner, potentially adding it to all submitted reads, if appropriate
     * @param read a non-null read
     */
    void submitRead(final GATKRead read) {
        samplePartitioner.submitRead(read);
    }

    /**
     * Add reads with the given sample name to the given hanger entry.
     *
     * @param readStates The list of read states to add this collection of reads.
     * @param reads      Reads to add.  Selected reads will be pulled from this source.
     */
    private void addReadsToSample(final PerSampleReadStateManager readStates, final Collection<GATKRead> reads) {
        if (reads.isEmpty()) {
            return;
        }

        final LinkedList<AlignmentStateMachine> newReadStates = new LinkedList<>();

        for (final GATKRead read : reads) {
            final AlignmentStateMachine state = new AlignmentStateMachine(read);
            if ( state.stepForwardOnGenome() != null ){ // todo -- should be an assertion not a skip
                // explicitly filter out reads that are all insertions / soft clips
                newReadStates.add(state);
            }
        }

        totalReadStates += readStates.addStatesAtNextAlignmentStart(newReadStates);
    }
}
