package org.broadinstitute.hellbender.utils.dragstr;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link DoubleSequence}
 */
public class DoubleSequenceTest {

    @Test(dataProvider = "correctSequenceData")
    public void testCorrectSequence(final String spec, final double[] expectedValues, final double epsilon) {
        final DoubleSequence subject = new DoubleSequence(spec);
        Assert.assertEquals(subject.size(), expectedValues.length);
        // check single value access:
        for (int i = 0; i < expectedValues.length; i++) {
            Assert.assertEquals(subject.get(i), expectedValues[i], epsilon);
        }
        // check array access:
        final double[] actualValues = subject.toDoubleArray();
        Assert.assertNotNull(actualValues);
        Assert.assertEquals(actualValues, expectedValues, epsilon);
    }

    @Test(dataProvider = "badSpecsData", expectedExceptions = RuntimeException.class)
    public void testBadSpecs(final String spec) {
        new DoubleSequence(spec);
    }

    @DataProvider(name = "badSpecsData")
    public Iterator<Object[]> badSpecsData() {
        return Arrays.stream(new String[] {null, "", "::", "1:2", "2:2:2:2:2", "2,2", "4;5;100", "10:-1:100"})
                .map(s -> new Object[] { s }).iterator();
    }

    @DataProvider(name = "correctSequenceData")
    public Iterator<Object[]> correctSequenceData() {
        final String STD_EPSILON = "0.001";
        final String[] testCases = {
                "0:1:10",      "0,1,2,3,4,5,6,7,8,9,10", STD_EPSILON,
                "-1:2.5:5",    "-1,1.5,4", STD_EPSILON,
                "1:1e-1:2.3",  "1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3", STD_EPSILON,
                "5:-1:1",      "5,4,3,2,1", STD_EPSILON,
                "1:0.00001:1", "1", STD_EPSILON,
                "5:-123:5",    "5", STD_EPSILON,
                "8:1:10.00001", "8,9,10.00001", "0.0000000001",
        };
        return IntStream.range(0, testCases.length / 3)
                 .map(i -> 3 * i)
                 .mapToObj(i -> {
                     final String spec = testCases[i++];
                     final String strValues = testCases[i++];
                     final double epsilon = Double.parseDouble(testCases[i]);
                     final double[] dblValues = Arrays.stream(strValues.split("\\s*,\\s*")).mapToDouble(Double::parseDouble).toArray();
                     return new Object[]{spec, dblValues, epsilon};
                 }).iterator();
    }
}
