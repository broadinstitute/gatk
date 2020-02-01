package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.*;

/**
 * Checks on Caller argument collection cloning and getting/setting of values.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class StandardCallerArgumentCollectionUnitTest extends GATKBaseTest {

    @Test
    public void testParameterLessInitialization() {
        try {
            new StandardCallerArgumentCollection();
        } catch (final Exception ex) {
            Assert.fail();
        }
    }

    @Test
    public void testGetSampleContaminationUninitializedArgs() {
        final StandardCallerArgumentCollection args = new StandardCallerArgumentCollection();

        Assert.assertFalse(args.isSampleContaminationPresent());
        Map<String, Double> returnedContaminationMap = args.getSampleContamination();

        // The returned map should be officially empty, but return 0.0 on query for any sample (it's a DefaultedMap)
        Assert.assertTrue(returnedContaminationMap.isEmpty());
        Assert.assertEquals((double) returnedContaminationMap.get("MySample"), 0.0);
    }

    @Test
    public void testGetSampleContaminationInitializedFraction() {
        final StandardCallerArgumentCollection args = new StandardCallerArgumentCollection();

        args.CONTAMINATION_FRACTION = 0.1;
        Assert.assertTrue(args.isSampleContaminationPresent());
        Map<String, Double> returnedContaminationMap = args.getSampleContamination();

        // The returned map should be officially empty, but return 0.1 on query for any sample (it's a DefaultedMap)
        Assert.assertTrue(returnedContaminationMap.isEmpty());
        Assert.assertEquals((double) returnedContaminationMap.get("MySample"), 0.1);
    }

    @Test
    public void testGetSampleContaminationInitializedMap() {
        final StandardCallerArgumentCollection args = new StandardCallerArgumentCollection();

        Map<String,Double> contaminationMap = new HashMap<>();
        contaminationMap.put("Sample1", 0.1);
        contaminationMap.put("Sample2", 0.2);
        args.setSampleContamination(contaminationMap);

        Assert.assertTrue(args.isSampleContaminationPresent());
        Map<String, Double> returnedContaminationMap = args.getSampleContamination();
        
        // The returned map should be of size 2, and return 0.0 on query for unknown samples (it's a DefaultedMap)
        Assert.assertEquals(returnedContaminationMap.size(), 2);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample1"), 0.1);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample2"), 0.2);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample3"), 0.0);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample4"), 0.0);
    }

    @Test
    public void testGetSampleContaminationInitializedMapAndFraction() {
        final StandardCallerArgumentCollection args = new StandardCallerArgumentCollection();

        args.CONTAMINATION_FRACTION = 0.05;
        Map<String,Double> contaminationMap = new HashMap<>();
        contaminationMap.put("Sample1", 0.1);
        contaminationMap.put("Sample2", 0.2);
        args.setSampleContamination(contaminationMap);

        Assert.assertTrue(args.isSampleContaminationPresent());
        Map<String, Double> returnedContaminationMap = args.getSampleContamination();

        // The returned map should be of size 2, and return 0.05 on query for unknown samples (it's a DefaultedMap)
        Assert.assertEquals(returnedContaminationMap.size(), 2);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample1"), 0.1);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample2"), 0.2);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample3"), 0.05);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample4"), 0.05);
    }

    @Test
    public void testGetSampleContaminationMapWithNoContamination() {
        final StandardCallerArgumentCollection args = new StandardCallerArgumentCollection();

        // Create a map that doesn't actually have any contamination set
        Map<String,Double> contaminationMap = new HashMap<>();
        contaminationMap.put("Sample1", 0.0);
        contaminationMap.put("Sample2", 0.0);
        args.setSampleContamination(contaminationMap);

        Assert.assertFalse(args.isSampleContaminationPresent());
        Map<String, Double> returnedContaminationMap = args.getSampleContamination();

        // The returned map should be of size 2, and return 0.0 on queries for any sample (it's a DefaultedMap)
        Assert.assertEquals(returnedContaminationMap.size(), 2);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample1"), 0.0);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample2"), 0.0);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample3"), 0.0);
        Assert.assertEquals((double) returnedContaminationMap.get("Sample4"), 0.0);
    }
}
