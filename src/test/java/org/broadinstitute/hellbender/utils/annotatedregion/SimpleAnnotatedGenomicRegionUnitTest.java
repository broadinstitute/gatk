package org.broadinstitute.hellbender.utils.annotatedregion;

import com.google.common.collect.ImmutableMap;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SimpleAnnotatedGenomicRegionUnitTest extends BaseTest {

    @Test(dataProvider = "equalsTests")
    public void testEquals(SimpleAnnotatedGenomicRegion test1, SimpleAnnotatedGenomicRegion test2, boolean gt) {
        Assert.assertEquals(test1.equals(test2), gt);
    }

    @DataProvider(name = "equalsTests")
    public Object [][] createEqualsTests() {

        return new Object[][] {

                {
                    new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                            ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        true
                }, {
                    new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                            ImmutableMap.of("Foo", "bar", "Foo1", "bar2")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        false
                }, {
                    new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                            ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableMap.of("Foo", "bar", "Foobar1", "bar1")),
                        false
                }, {
                    new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 1000, 2000),
                            ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        false
                }, {
                    new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                            ImmutableMap.of("Foo", "bar", "Foo1", "bar1")),
                        new SimpleAnnotatedGenomicRegion( new SimpleInterval("1", 100, 200),
                                ImmutableMap.of("Foo", "bar")),
                        false
                }
        };
    }
}
