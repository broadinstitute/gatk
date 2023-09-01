package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.OccurrenceMatrix;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class OccurrenceMatrixUnitTest<C> {

    @Test(dataProvider = "createMatrixTestProvider")
    public void testCreateMatrix(Map<Integer, Collection<Integer>> testInput, List<Pair<Integer, Integer>> testOutput, List<Set<Integer>> testComponents){
        OccurrenceMatrix<Integer, Integer> occm = new OccurrenceMatrix<>(testInput);
        List<Pair<Integer, Integer>> result = occm.nonCoOcurringColumns();
        List<Set<Integer>> components = occm.getIndependentSets(result);

        Assert.assertEquals(result.equals(testOutput),true);
        Assert.assertEquals(components.equals(testComponents),true);

    }


    @DataProvider(name = "createMatrixTestProvider")
    public Object[][] createMatrixTest() {
        List<Object[]> tests = new ArrayList<>();

        Map<Integer, Collection<Integer>> testInput = new HashMap<>();
        Integer[] vals = new Integer[] {1,2};
        testInput.put(1, Arrays.asList(vals));
        vals = new Integer[] {2,3};
        testInput.put(2, Arrays.asList(vals));
        vals = new Integer[] {1,3};
        testInput.put(3, Arrays.asList(vals));

        List<Pair<Integer, Integer>> testOutput = new ArrayList<>();
        List<Set<Integer>> testCC = new ArrayList<>();
        testCC.add(new HashSet<>(Arrays.asList(1)));
        testCC.add(new HashSet<>(Arrays.asList(2)));
        testCC.add(new HashSet<>(Arrays.asList(3)));


        tests.add(new Object[]{testInput, testOutput, testCC});

        testInput = new HashMap<>();
        vals = new Integer[] {1,2};
        testInput.put(1, Arrays.asList(vals));
        vals = new Integer[] {3};
        testInput.put(2, Arrays.asList(vals));

        testOutput = new ArrayList<>();
        testOutput.add(new ImmutablePair<>(1,3));
        testOutput.add(new ImmutablePair<>(2,3));

        testCC = new ArrayList<>();
        testCC.add(new HashSet<>(Arrays.asList(1,2,3)));
        tests.add(new Object[]{testInput, testOutput, testCC});

        testInput = new HashMap<>();
        vals = new Integer[] {1,2};
        testInput.put(1, Arrays.asList(vals));
        vals = new Integer[] {3,4};
        testInput.put(2, Arrays.asList(vals));

        testOutput = new ArrayList<>();
        testOutput.add(new ImmutablePair<>(1,3));
        testOutput.add(new ImmutablePair<>(1,4));
        testOutput.add(new ImmutablePair<>(2,3));
        testOutput.add(new ImmutablePair<>(2,4));

        testCC = new ArrayList<>();
        testCC.add(new HashSet<>(Arrays.asList(1,2,3,4)));
        tests.add(new Object[]{testInput, testOutput, testCC});


        testInput = new HashMap<>();
        vals = new Integer[] {1,2};
        testInput.put(1, Arrays.asList(vals));
        vals = new Integer[] {4};
        testInput.put(2, Arrays.asList(vals));

        testOutput = new ArrayList<>();
        testOutput.add(new ImmutablePair<>(1,4));
        testOutput.add(new ImmutablePair<>(2,4));
        testCC = new ArrayList<>();
        testCC.add(new HashSet<>(Arrays.asList(1,2,4)));
        tests.add(new Object[]{testInput, testOutput, testCC});

        return tests.toArray(new Object[][]{});

    }


}
