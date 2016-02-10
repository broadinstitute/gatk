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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.downsampling.Downsampler;
import org.broadinstitute.gatk.utils.downsampling.PassThroughDownsampler;
import org.broadinstitute.gatk.utils.downsampling.ReservoirDownsampler;

import java.util.*;

/**
 * Divides reads by sample and (if requested) does a preliminary downsampling pass
 * with a ReservoirDownsampler.
 *
 * Note: stores reads by sample ID string, not by sample object
 */
class SamplePartitioner<T extends SAMRecord> {
    /**
     * Map from sample name (as a string) to a downsampler of reads for that sample
     */
    final private Map<String, Downsampler<T>> readsBySample;

    /**
     * Are we in a state where we're done submitting reads and have semi-finalized the
     * underlying per sample downsampler?
     */
    boolean doneSubmittingReads = false;

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
    @Ensures({
            "readsBySample != null",
            "readsBySample.size() == new HashSet(samples).size()"
    })
    public SamplePartitioner(final LIBSDownsamplingInfo LIBSDownsamplingInfo, final List<String> samples) {
        if ( LIBSDownsamplingInfo == null ) throw new IllegalArgumentException("LIBSDownsamplingInfo cannot be null");
        if ( samples == null ) throw new IllegalArgumentException("samples must be a non-null list");

        readsBySample = new LinkedHashMap<String, Downsampler<T>>(samples.size());
        for ( final String sample : samples ) {
            readsBySample.put(sample, createDownsampler(LIBSDownsamplingInfo));
        }
    }

    /**
     * Create a new, ready to use downsampler based on the parameters in LIBSDownsamplingInfo
     * @param LIBSDownsamplingInfo the parameters to use in creating the downsampler
     * @return a downsampler appropriate for LIBSDownsamplingInfo.  If no downsampling is requested,
     *   uses the PassThroughDownsampler, which does nothing at all.
     */
    @Requires("LIBSDownsamplingInfo != null")
    @Ensures("result != null")
    private Downsampler<T> createDownsampler(final LIBSDownsamplingInfo LIBSDownsamplingInfo) {
        return LIBSDownsamplingInfo.isPerformDownsampling()
                ? new ReservoirDownsampler<T>(LIBSDownsamplingInfo.getToCoverage(), true)
                : new PassThroughDownsampler<T>();
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
    @Requires("read != null")
    @Ensures("doneSubmittingReads == false")
    public void submitRead(final T read) {
        final String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;
        final Downsampler<T> downsampler = readsBySample.get(sampleName);
        if ( downsampler == null )
            throw new IllegalStateException("Offered read with sample name " + sampleName + " to SamplePartitioner " +
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
    @Ensures("doneSubmittingReads == true")
    public void doneSubmittingReads() {
        for ( final Downsampler<T> downsampler : readsBySample.values() ) {
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
     * @param sampleName the sample we want reads for, must be present in the original samples
     * @return a non-null collection of reads for sample in this cycle
     */
    @Ensures("result != null")
    public Collection<T> getReadsForSample(final String sampleName) {
        if ( ! doneSubmittingReads ) throw new IllegalStateException("getReadsForSample called before doneSubmittingReads was called");

        final Downsampler<T> downsampler = readsBySample.get(sampleName);
        if ( downsampler == null ) throw new NoSuchElementException("Sample name not found");

        return downsampler.consumeFinalizedItems();
    }

    /**
     * Resets this SamplePartitioner, indicating that we're starting a new
     * cycle of adding reads to each underlying downsampler.
     */
    @Ensures("doneSubmittingReads == false")
    public void reset() {
        for ( final Downsampler<T> downsampler : readsBySample.values() ) {
            downsampler.clearItems();
            downsampler.resetStats();
        }
        doneSubmittingReads = false;
    }
}
