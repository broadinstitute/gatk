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
 * GenomeLocs are very useful objects to keep track of genomic locations and perform set operations
 * with them.
 *
 * However, GenomeLocs are bound to strict validation through the GenomeLocParser and cannot
 * be created easily for small tasks that do not require the rigors of the GenomeLocParser validation
 *
 * UnvalidatingGenomeLoc is a simple utility to create GenomeLocs without going through the parser.
 *
 * WARNING: SHOULD BE USED ONLY BY EXPERT USERS WHO KNOW WHAT THEY ARE DOING!
 */
public class UnvalidatingGenomeLoc extends GenomeLoc {

    public UnvalidatingGenomeLoc(String contigName, int contigIndex, int start, int stop) {
        super(contigName, contigIndex, start, stop);
    }
}
