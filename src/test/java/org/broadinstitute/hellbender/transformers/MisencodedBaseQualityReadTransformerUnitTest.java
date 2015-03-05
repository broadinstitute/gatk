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


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * Basic unit test for misencoded quals
 */
public class MisencodedBaseQualityReadTransformerUnitTest extends BaseTest {

    private static final byte[] badQuals = { 59, 60, 62, 63, 64, 61, 62, 58, 57, 56 };
    private static final byte[] goodQuals = { 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 };
    private static final byte[] fixedQuals = { 28, 29, 31, 32, 33, 30, 31, 27, 26, 25 };
    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
    }

    private SAMRecord createRead(final byte[] quals) {
        final String readBases = "AAAAAAAAAA";
        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 10, readBases.getBytes(), quals);
        read.setCigarString("10M");
        return read;
    }

    @Test(enabled = true)
    public void testGoodQuals() {
        final ReadTransformer tr = new MisencodedBaseQualityReadTransformer();
        SAMRecord read = createRead(goodQuals);
        SAMRecord newRead = tr.apply(read);
        Assert.assertEquals(read, newRead);
    }

    @Test(enabled = true)
    public void testFixBadQuals() {
        final ReadTransformer tr = new MisencodedBaseQualityReadTransformer();
        final SAMRecord read = createRead(badQuals);
        final SAMRecord fixedRead = tr.apply(read);
        Assert.assertEquals(fixedQuals, fixedRead.getBaseQualities());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFixGoodQualsBlowUp() {
        final ReadTransformer tr = new MisencodedBaseQualityReadTransformer();
        final SAMRecord read = createRead(fixedQuals);
        tr.apply(read);
    }
}