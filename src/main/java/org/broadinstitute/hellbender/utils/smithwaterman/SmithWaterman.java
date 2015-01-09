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

package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;

/**
 * Generic interface for SmithWaterman calculations
 *
 * This interface allows clients to use a generic SmithWaterman variable, without propagating the specific
 * implementation of SmithWaterman throughout their code:
 *
 * SmithWaterman sw = new SpecificSmithWatermanImplementation(ref, read, params)
 * sw.getCigar()
 * sw.getAlignmentStart2wrt1()
 */
public interface SmithWaterman {

    /**
     * Get the cigar string for the alignment of this SmithWaterman class
     * @return a non-null cigar
     */
    public Cigar getCigar();

    /**
     * Get the starting position of the read sequence in the reference sequence
     * @return a positive integer >= 0
     */
    public int getAlignmentStart2wrt1();
}
