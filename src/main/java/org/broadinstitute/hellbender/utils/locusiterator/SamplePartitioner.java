package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.Downsampler;
import org.broadinstitute.hellbender.utils.downsampling.PassThroughDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReservoirDownsampler;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Divides reads by sample and (if requested) does a preliminary downsampling pass
 * with a ReservoirDownsampler.
 *
 * Note: stores reads by sample ID string, not by sample object
 */
final class SamplePartitioner {
    /**
     * Map from sample name (as a string) to a downsampler of reads for that sample
     */
    private final Map<String, Downsampler<GATKRead>> readsBySample;

    /**
     * Are we in a state where we're done submitting reads and have semi-finalized the
     * underlying per sample downsampler?
     */
    boolean doneSubmittingReads = false;

    private final SAMFileHeader header;

    /**
     * Create a new SamplePartitioner capable of splitting reads up into buckets of reads for
     * each sample in samples, and perform a preliminary downsampling of these reads
     * (separately for each sample) if downsampling is requested in LIBSDownsamplingInfo
     *
     * Note that samples must be comprehensive, in that all reads every submitted to this
     * partitioner must come from one of the samples provided here.  If not, submitRead
     * will throw an exception.  Duplicates in the list of samples will be ignored
     *
     * @param LIBSDownsamplingInfo do we want to downsample, and if so to what coverage?
     * @param samples the complete list of samples we're going to partition reads into. Can be
     *                empty, but in that case this code cannot function properly if you
     *                attempt to add data to it.
     */
    public SamplePartitioner(final LIBSDownsamplingInfo LIBSDownsamplingInfo, final List<String> samples, final SAMFileHeader header) {
        Utils.nonNull(LIBSDownsamplingInfo, "LIBSDownsamplingInfo cannot be null");
        Utils.nonNull(samples, "samples must be a non-null list");
        Utils.nonNull(header, "header must be not null");

        readsBySample = new LinkedHashMap<>(samples.size());
        for ( final String sample : samples ) {
            readsBySample.put(sample, createDownsampler(LIBSDownsamplingInfo));
        }
        this.header = header;
    }

    /**
     * Create a new, ready to use downsampler based on the parameters in LIBSDownsamplingInfo
     * @param LIBSDownsamplingInfo the parameters to use in creating the downsampler
     * @return a downsampler appropriate for LIBSDownsamplingInfo.  If no downsampling is requested,
     *   uses the PassThroughDownsampler, which does nothing at all.
     */
    private Downsampler<GATKRead> createDownsampler(final LIBSDownsamplingInfo LIBSDownsamplingInfo) {
        return LIBSDownsamplingInfo.isPerformDownsampling()
                ? new ReservoirDownsampler(LIBSDownsamplingInfo.getToCoverage(), true)
                : new PassThroughDownsampler();
    }

    /**
     * Offer this read to the partitioner, putting it into the bucket of reads for the sample
     * of read (obtained via the read's read group).
     *
     * If the read group is missing, uses the special "null" read group
     *
     * @throws IllegalStateException if the sample of read wasn't present in the original
     *   set of samples provided to this SamplePartitioner at construction
     *
     * @param read the read to add to the sample's list of reads
     */
    public void submitRead(final GATKRead read) {
        Utils.nonNull(read);
        final String sampleName = read.getReadGroup() != null ? ReadUtils.getSampleName(read, header) : null;
        final Downsampler<GATKRead> downsampler = readsBySample.get(sampleName);
        Utils.validate(downsampler != null, () -> "Offered read with sample name " + sampleName + " to SamplePartitioner " +
                    "but this sample wasn't provided as one of possible samples at construction");

        downsampler.submit(read);
        doneSubmittingReads = false;
    }

    /**
     * Tell this partitioner that all reads in this cycle have been submitted, so that we
     * can finalize whatever downsampling is required by each sample.
     *
     * Note that we *must* call this function before getReadsForSample, or else that
     * function will exception out.
     */
    public void doneSubmittingReads() {
        for ( final Downsampler<GATKRead> downsampler : readsBySample.values() ) {
            downsampler.signalEndOfInput();
        }
        doneSubmittingReads = true;
    }

    /**
     * Get the final collection of reads for this sample for this cycle
     *
     * The cycle is defined as all of the reads that occur between
     * the first call to submitRead until doneSubmittingReads is called.  At that
     * point additional downsampling may occur (depending on construction arguments)
     * and that set of reads is returned here.
     *
     * Note that this function can only be called once per cycle, as underlying
     * collection of reads is cleared.
     *
     * @param sampleName the sample we want reads for, must be present in the original samples. Can be null.
     * @return a non-null collection of reads for sample in this cycle
     */
    public Collection<GATKRead> getReadsForSample(final String sampleName) {
        Utils.validate(doneSubmittingReads, "getReadsForSample called before doneSubmittingReads was called");

        final Downsampler<GATKRead> downsampler = readsBySample.get(sampleName);
        if ( downsampler == null ) throw new NoSuchElementException("Sample name not found");

        return downsampler.consumeFinalizedItems();
    }

    /**
     * Resets this SamplePartitioner, indicating that we're starting a new
     * cycle of adding reads to each underlying downsampler.
     */
    public void reset() {
        for ( final Downsampler<GATKRead> downsampler : readsBySample.values() ) {
            downsampler.clearItems();
            downsampler.resetStats();
        }
        doneSubmittingReads = false;
    }
}
