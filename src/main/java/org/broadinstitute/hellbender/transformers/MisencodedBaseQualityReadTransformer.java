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

package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Checks for and errors out (or fixes if requested) when it detects reads with base qualities that are not encoded with
 * phred-scaled quality scores.  Q0 == ASCII 33 according to the SAM specification, whereas Illumina encoding starts at
 * Q64.  The idea here is simple: if we are asked to fix the scores then we just subtract 31 from every quality score.
 */
public final class MisencodedBaseQualityReadTransformer implements ReadTransformer {

    private static final int ILLUMINA_ENCODING_FIX_VALUE = 31;  // Illumina_64 - PHRED_33

    @Override
    public SAMRecord apply(final SAMRecord read) {
        final byte[] quals = read.getBaseQualities();
        for ( int i = 0; i < quals.length; i++ ) {
            quals[i] -= ILLUMINA_ENCODING_FIX_VALUE;
            if ( quals[i] < 0 )
                throw new UserException.BadInput("while fixing mis-encoded base qualities we encountered a read that was correctly encoded; we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some other tool");
        }
        read.setBaseQualities(quals);
        return read;
    }
}
