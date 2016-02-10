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

package org.broadinstitute.gatk.utils.downsampling;

import htsjdk.samtools.SAMRecord;

/**
 * An extension of the basic downsampler API with reads-specific operations
 *
 * @author David Roazen
 */
public abstract class ReadsDownsampler<T extends SAMRecord> extends Downsampler<T> {

    /**
     * Does this downsampler require that reads be fed to it in coordinate order?
     *
     * @return true if reads must be submitted to this downsampler in coordinate order, otherwise false
     */
    public abstract boolean requiresCoordinateSortOrder();

    /**
     * Tell this downsampler that no more reads located before the provided read (according to
     * the sort order of the read stream) will be fed to it.
     *
     * Allows position-aware downsamplers to finalize pending reads earlier than they would
     * otherwise be able to, particularly when doing per-sample downsampling and reads for
     * certain samples are sparser than average.
     *
     * @param read the downsampler will assume that no reads located before this read will ever
     *             be submitted to it in the future
     */
    public abstract void signalNoMoreReadsBefore( final T read );
}
