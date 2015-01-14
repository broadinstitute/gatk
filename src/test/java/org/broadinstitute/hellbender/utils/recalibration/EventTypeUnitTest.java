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

package org.broadinstitute.hellbender.utils.recalibration;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashSet;
import java.util.Set;

public final class EventTypeUnitTest extends BaseTest {
    @Test
    public void testEventTypes() {
        for ( final EventType et : EventType.values() ) {
            Assert.assertNotNull(et.toString());
            Assert.assertNotNull(et.prettyPrint());
            Assert.assertFalse("".equals(et.toString()));
            Assert.assertFalse("".equals(et.prettyPrint()));
            Assert.assertEquals(EventType.eventFrom(et.ordinal()), et);
            Assert.assertEquals(EventType.eventFrom(et.toString()), et);
        }
    }

    @Test
    public void testEventTypesEnumItself() {
        final Set<String> shortReps = new HashSet<String>();
        for ( final EventType et : EventType.values() ) {
            Assert.assertFalse(shortReps.contains(et.toString()), "Short representative for EventType has duplicates for " + et);
            shortReps.add(et.toString());
        }
        Assert.assertEquals(shortReps.size(), EventType.values().length, "Short representatives for EventType aren't unique");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadString() {
        EventType.eventFrom("asdfhalsdjfalkjsdf");
    }
}
