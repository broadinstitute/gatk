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

package org.broadinstitute.hellbender.utils.text;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

public class TextFormattingUtilsUnitTest extends BaseTest{
    @Test(expectedExceptions = GATKException.class)
    public void testSplitWhiteSpaceNullLine() {
        TextFormattingUtils.splitWhiteSpace(null);
    }

    @Test
    public void testSplitWhiteSpace() {
        Assert.assertEquals(TextFormattingUtils.splitWhiteSpace("foo bar baz"), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitWhiteSpace("foo  bar  baz"), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitWhiteSpace(" foo bar baz"), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitWhiteSpace(" foo bar baz "), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitWhiteSpace("foo bar baz "), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitWhiteSpace("\tfoo\tbar\tbaz\t"), new String[]{"foo", "bar", "baz"});
    }

    @Test(expectedExceptions = GATKException.class)
    public void testGetWordStartsNullLine() {
        TextFormattingUtils.getWordStarts(null);
    }

    @Test
    public void testGetWordStarts() {
        Assert.assertEquals(TextFormattingUtils.getWordStarts("foo bar baz"), Arrays.asList(4, 8));
        Assert.assertEquals(TextFormattingUtils.getWordStarts("foo  bar  baz"), Arrays.asList(5, 10));
        Assert.assertEquals(TextFormattingUtils.getWordStarts(" foo bar baz"), Arrays.asList(1, 5, 9));
        Assert.assertEquals(TextFormattingUtils.getWordStarts(" foo bar baz "), Arrays.asList(1, 5, 9));
        Assert.assertEquals(TextFormattingUtils.getWordStarts("foo bar baz "), Arrays.asList(4, 8));
        Assert.assertEquals(TextFormattingUtils.getWordStarts("\tfoo\tbar\tbaz\t"), Arrays.asList(1, 5, 9));
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSplitFixedWidthNullLine() {
        TextFormattingUtils.splitFixedWidth(null, Collections.emptyList());
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSplitFixedWidthNullColumnStarts() {
        TextFormattingUtils.splitFixedWidth("foo bar baz", null);
    }

    @Test
    public void testSplitFixedWidth() {
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth("foo bar baz", Arrays.asList(4, 8)), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth("foo  bar  baz", Arrays.asList(5, 10)), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth(" foo bar baz", Arrays.asList(5, 9)), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth(" foo bar baz ", Arrays.asList(5, 9)), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth("foo bar baz ", Arrays.asList(4, 8)), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth("\tfoo\tbar\tbaz\t", Arrays.asList(5, 9)), new String[]{"foo", "bar", "baz"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth("f o b r b z", Arrays.asList(4, 8)), new String[]{"f o", "b r", "b z"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth(" f o b r b z", Arrays.asList(4, 8)), new String[]{"f o", "b r", "b z"});
        Assert.assertEquals(TextFormattingUtils.splitFixedWidth("  f o b r b z", Arrays.asList(4, 8)), new String[]{"f", "o b", "r b z"});
    }
}
