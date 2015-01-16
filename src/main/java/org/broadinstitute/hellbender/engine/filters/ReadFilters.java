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

package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.function.Predicate;

public final class ReadFilters  {

    private ReadFilters(){}
    /*
     * Below is a list of standard pre-packaged read filters.
     */

    public static final Predicate<SAMRecord> UNMAPPED =  read -> read.getReadUnmappedFlag() || read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;

    public static final Predicate<SAMRecord> NOT_PRIMARY_ALIGNMENT = SAMRecord::getNotPrimaryAlignmentFlag;

    public static final Predicate<SAMRecord> DUPLICATE = SAMRecord::getDuplicateReadFlag;

    public static final Predicate<SAMRecord> FAILS_VENDOR_QUALITY_CHECK = SAMRecord::getReadFailsVendorQualityCheckFlag;

    public static final Predicate<SAMRecord> MAPPING_QUALITY_UNAVAIALBLE = read -> read.getMappingQuality() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE;

    public static final Predicate<SAMRecord> MAPPING_QUALITY_ZERO = read -> read.getMappingQuality() == 0;
}
