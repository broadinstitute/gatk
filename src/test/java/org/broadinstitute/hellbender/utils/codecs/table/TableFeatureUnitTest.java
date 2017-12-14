package org.broadinstitute.hellbender.utils.codecs.table;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;


/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class TableFeatureUnitTest extends GATKBaseTest {

    @Test
    public void testTableFeatureGetters() {
        final TableFeature feature = new TableFeature(
                new SimpleInterval("1", 10, 100),
                Arrays.asList("1", "2", "3"),
                Arrays.asList("a", "b", "c"));

        // test all values and columns
        Assert.assertEquals(feature.getAllValues(), Arrays.asList("1", "2", "3"));
        Assert.assertEquals(feature.getHeader(), Arrays.asList("a", "b", "c"));

        // test retrieval of all the elements one by one and its mapping
        for (int i = 0; i < 3; i++) {
            final String colName = feature.getHeader().get(i);
            Assert.assertEquals(feature.getValue(i), feature.get(colName));
        }

        // test getValuesTo
        Assert.assertEquals(feature.getValuesTo(1), Arrays.asList("1"));
        Assert.assertEquals(feature.getValuesTo(2), Arrays.asList("1", "2"));

        // test that invalid columns throw an Illegal argument exception
        Assert.assertThrows(IllegalArgumentException.class, () -> feature.getValue(3));
        Assert.assertThrows(IllegalArgumentException.class, () -> feature.getValue(-1));
        Assert.assertThrows(IllegalArgumentException.class, () -> feature.getValuesTo(4));
        Assert.assertThrows(IllegalArgumentException.class, () -> feature.getValuesTo(-1));
    }

}