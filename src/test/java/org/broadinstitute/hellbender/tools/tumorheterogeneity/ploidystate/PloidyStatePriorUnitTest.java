package org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

public class PloidyStatePriorUnitTest {
    private static final double EPSILON = 1E-10;

    @Test
    public void testPloidyStatePrior() throws Exception {
        final List<PloidyState> ploidyStates = Arrays.asList(
                new PloidyState(0, 0),
                new PloidyState(0, 1),
                new PloidyState(1, 0),
                new PloidyState(1, 1));
        final List<Double> unnormalizedLogProbabilities = Arrays.asList(
                Math.log(5.),
                Math.log(10.),
                Math.log(10.),
                Math.log(75.));
        final Map<PloidyState, Double> ploidyStateToUnnormalizedLogProbabilityMap = new LinkedHashMap<>();
        IntStream.range(0, ploidyStates.size()).forEach(i -> ploidyStateToUnnormalizedLogProbabilityMap.put(ploidyStates.get(i), unnormalizedLogProbabilities.get(i)));
        final PloidyStatePrior ploidyStatePrior = new PloidyStatePrior(ploidyStateToUnnormalizedLogProbabilityMap);

        Assert.assertEquals(ploidyStatePrior.ploidyStates(), ploidyStates);
        Assert.assertEquals(ploidyStatePrior.maxCopyNumber(), 2);
        Assert.assertEquals(ploidyStatePrior.logProbability(new PloidyState(0, 0)), Math.log(0.05), EPSILON);
        Assert.assertEquals(ploidyStatePrior.logProbability(new PloidyState(0, 1)), Math.log(0.1), EPSILON);
        Assert.assertEquals(ploidyStatePrior.logProbability(new PloidyState(1, 0)), Math.log(0.1), EPSILON);
        Assert.assertEquals(ploidyStatePrior.logProbability(new PloidyState(1, 1)), Math.log(0.75), EPSILON);
        Assert.assertEquals(ploidyStatePrior.logProbability(0), Math.log(0.05), EPSILON);
        Assert.assertEquals(ploidyStatePrior.logProbability(1), Math.log(0.2), EPSILON);
        Assert.assertEquals(ploidyStatePrior.logProbability(2), Math.log(0.75), EPSILON);
    }
}