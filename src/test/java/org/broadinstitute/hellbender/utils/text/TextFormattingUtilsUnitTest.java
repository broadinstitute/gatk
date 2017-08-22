package org.broadinstitute.hellbender.utils.text;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

public final class TextFormattingUtilsUnitTest extends GATKBaseTest {
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
