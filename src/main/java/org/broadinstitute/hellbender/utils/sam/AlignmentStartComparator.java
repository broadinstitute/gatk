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

package org.broadinstitute.hellbender.utils.sam;

import htsjdk.samtools.SAMRecord;

import java.util.Comparator;

/**
 * Compares two SAMRecords only the basis on alignment start.  Note that
 * comparisons are performed ONLY on the basis of alignment start; any
 * two SAM records with the same alignment start will be considered equal.
 *
 * Unmapped alignments will all be considered equal.
 */
public class AlignmentStartComparator implements Comparator<SAMRecord> {
    public int compare(SAMRecord lhs, SAMRecord rhs) {
        if(!lhs.getReferenceIndex().equals(rhs.getReferenceIndex()))
            return lhs.getReferenceIndex() - rhs.getReferenceIndex();

        // Note: no integer overflow here because alignment starts are >= 0.
        return lhs.getAlignmentStart() - rhs.getAlignmentStart();
    }
}
