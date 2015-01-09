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

package org.broadinstitute.hellbender.utils;

/**
 * How should we represent a clipped bases in a read?
 */
public enum ClippingRepresentation {
    /** Clipped bases are changed to Ns */
    WRITE_NS,

    /** Clipped bases are changed to have Q0 quality score */
    WRITE_Q0S,

    /** Clipped bases are change to have both an N base and a Q0 quality score */
    WRITE_NS_Q0S,

    /**
     * Change the read's cigar string to soft clip (S, see sam-spec) away the bases.
     * Note that this can only be applied to cases where the clipped bases occur
     * at the start or end of a read.
     */
    SOFTCLIP_BASES,

    /**
     * WARNING: THIS OPTION IS STILL UNDER DEVELOPMENT AND IS NOT SUPPORTED.
     *
     * Change the read's cigar string to hard clip (H, see sam-spec) away the bases.
     * Hard clipping, unlike soft clipping, actually removes bases from the read,
     * reducing the resulting file's size but introducing an irrevesible (i.e.,
     * lossy) operation.  Note that this can only be applied to cases where the clipped
     * bases occur at the start or end of a read.
     */
    HARDCLIP_BASES,

    /**
     * Turn all soft-clipped bases into matches
     */
    REVERT_SOFTCLIPPED_BASES,
}
