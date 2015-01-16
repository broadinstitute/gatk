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

package org.broadinstitute.hellbender.utils.R;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class RUtilsUnitTest {
    @DataProvider(name = "stringLists")
    public Object[][] getStringLists() {
        return new Object[][] {
                new Object[] { null, "NA" },
                new Object[] { Collections.EMPTY_LIST, "c()" },
                new Object[] { Arrays.asList("1", "2", "3"), "c('1','2','3')" }
        };
    }

    @Test(dataProvider = "stringLists")
    public void testToStringList(List<? extends CharSequence> actual, String expected) {
        Assert.assertEquals(RUtils.toStringList(actual), expected);
    }

    @DataProvider(name = "numberLists")
    public Object[][] getNumberLists() {
        return new Object[][] {
                new Object[] { null, "NA" },
                new Object[] { Collections.EMPTY_LIST, "c()" },
                new Object[] { Arrays.asList(1, 2, 3), "c(1,2,3)" },
                new Object[] { Arrays.asList(1D, 2D, 3D), "c(1.0,2.0,3.0)" }
        };
    }

    @Test(dataProvider = "numberLists")
    public void testToNumberList(List<? extends Number> actual, String expected) {
        Assert.assertEquals(RUtils.toNumberList(actual), expected);
    }
}
