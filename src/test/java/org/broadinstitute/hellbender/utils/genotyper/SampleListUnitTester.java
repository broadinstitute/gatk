package org.broadinstitute.hellbender.utils.genotyper;

import org.testng.Assert;

import java.util.*;

/**
 * Helper class for those unit-test classes that test on implementations of SampleList.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleListUnitTester {

    /**
     * Test that the contents of a sample-list are the ones expected.
     *
     * <p>
     *     This method perform various consistency check involving all the {@link org.broadinstitute.gatk.utils.genotyper.SampleList} interface methods.
     *     Therefore calling this method is equivalent to a thorough check of the {@link org.broadinstitute.gatk.utils.genotyper.SampleList} aspect of
     *     the {@code actual} argument.
     * </p>
     *
     * @param actual the sample-list to assess.
     * @param expected the expected sample-list.
     *
     * @throws IllegalArgumentException if {@code expected} is {@code null} or contains
     *   {@code null}s which is an indication of an bug in the testing code.
     *
     * @throws java.lang.RuntimeException if there is some testing assertion exception which
     *   is an indication of an actual bug the code that is been tested.
     */
    public static void assertSampleList(final SampleList actual, final List<String> expected) {
        if (expected == null)
            throw new IllegalArgumentException("the expected list cannot be null");
        final Set<String> expectedNames = new LinkedHashSet<>(expected.size());
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.numberOfSamples(), expected.size());
        for (int i = 0; i < expected.size(); i++) {
            final String expectedSample = expected.get(i);
            if (expectedSample == null)
                throw new IllegalArgumentException("the expected sample cannot be null");
            if (expectedSample.equals(NEVER_USE_SAMPLE_NAME))
                throw new IllegalArgumentException("you cannot use the forbidden sample name");
            if (expectedNames.contains(expected.get(i)))
                throw new IllegalArgumentException("repeated names in the expected list, this is a test bug");
            final String actualSample = actual.getSample(i);
            Assert.assertNotNull(actualSample, "sample name cannot be null");
            Assert.assertFalse(expectedNames.contains(actualSample), "repeated sample name: " + actualSample);
            Assert.assertEquals(actualSample, expectedSample, "wrong sample name order; index = " + i);
            Assert.assertEquals(actual.indexOfSample(actualSample), i, "sample index mismatch");
            expectedNames.add(actualSample);
        }

        Assert.assertEquals(actual.indexOfSample(NEVER_USE_SAMPLE_NAME), -1);
    }

    /**
     * Creates a sample list for testing given the number of samples in it.
     * @param sampleCount the required sample count.
     * @return never {@code null}.
     */
    public static SampleList sampleList(final int sampleCount) {
        if (sampleCount < 0)
            throw new IllegalArgumentException("the number of sample cannot be negative");
        final List<String> result = new ArrayList<>(sampleCount);
        for (int i =0; i < sampleCount; i++)
            result.add("SAMPLE_" + i);
        return new IndexedSampleList(result);
    }

    /**
     * Save to assume that this sample name will never be used.
     */
    private static final String NEVER_USE_SAMPLE_NAME = "WHY_WOULD_YOU_CALL_A_SAMPLE_LIKE_THIS?  ArE yOu Crazzzzy? " + new Date().toString();
}
