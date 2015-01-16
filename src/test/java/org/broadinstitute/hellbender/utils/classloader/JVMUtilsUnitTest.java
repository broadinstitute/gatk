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

package org.broadinstitute.hellbender.utils.classloader;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

@SuppressWarnings("rawtypes")
public class JVMUtilsUnitTest {

    // Test classes used by the tests for JVMUtils.getCallingClass():
    private static class DummyTestClass1 {
        public static Class getCaller( final Class callee ) {
            return DummyTestClass2.getCaller(callee);
        }
    }

    private static class DummyTestClass2 {
        public static Class getCaller( final Class callee ) {
            return DummyTestClass3.getCaller(callee);
        }
    }

    private static class DummyTestClass3 {
        public static Class getCaller( final Class callee ) {
            return JVMUtils.getCallingClass(callee);
        }
    }

    @DataProvider( name = "TestGetCallingClassDataProvider" )
    public Object[][] getTestCallingClassTestData() {
        return new Object[][] {
            { DummyTestClass1.class, JVMUtilsUnitTest.class },
            { DummyTestClass2.class, DummyTestClass1.class },
            { DummyTestClass3.class, DummyTestClass2.class }
        };
    }

    @Test( dataProvider = "TestGetCallingClassDataProvider" )
    public void testGetCallingClass( final Class callee, final Class expectedCaller ) {
        final Class reportedCaller = DummyTestClass1.getCaller(callee);

        Assert.assertEquals(reportedCaller, expectedCaller,
                String.format("Wrong calling class returned from DummyTestClass1.getCaller(%s)", callee.getSimpleName()));
    }

    @Test( expectedExceptions = IllegalArgumentException.class )
    public void testGetCallingClassCalleeNotFound() {
        // Trying to get the calling class of a class not on the runtime stack should produce an exception.
        JVMUtils.getCallingClass(DummyTestClass1.class);
    }
}
