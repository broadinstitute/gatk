package org.broadinstitute.hellbender.utils.recalibration;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.LinkedHashSet;
import java.util.Set;

public final class EventTypeUnitTest extends GATKBaseTest {
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
        final Set<String> shortReps = new LinkedHashSet<>();
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
