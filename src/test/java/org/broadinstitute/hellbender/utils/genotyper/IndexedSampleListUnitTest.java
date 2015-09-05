package org.broadinstitute.hellbender.utils.genotyper;


import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests {@link org.broadinstitute.gatk.utils.genotyper.IndexedSampleList}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IndexedSampleListUnitTest {

    @Test
    public void testEmptyConstructor() {
        final IndexedSampleList subject = new IndexedSampleList();
        SampleListUnitTester.assertSampleList(subject, Collections.emptyList());
    }

    @Test(dataProvider="sampleCountMaxSampleIndexData")
    public void testArrayConstructor(final int sampleCount, final int maxSampleIndex) {
        final String[] sampleNames = generateSampleNames(sampleCount,maxSampleIndex);

        final LinkedHashSet<String> nonRepeatedNames = new LinkedHashSet<>(Arrays.asList(sampleNames));
        final IndexedSampleList subject = new IndexedSampleList(sampleNames);
        SampleListUnitTester.assertSampleList(subject, Arrays.asList(nonRepeatedNames.toArray(new String[nonRepeatedNames.size()])));
    }

    @Test(dataProvider="sampleCountMaxSampleIndexData")
    public void testCollectionConstructor(final int sampleCount, final int maxSampleIndex) {
        final String[] sampleNames = generateSampleNames(sampleCount,maxSampleIndex);

        final List<String> sampleNameList = Arrays.asList(sampleNames);
        final LinkedHashSet<String> nonRepeatedNames = new LinkedHashSet<>(Arrays.asList(sampleNames));
        final IndexedSampleList subject = new IndexedSampleList(sampleNameList);
        SampleListUnitTester.assertSampleList(subject, Arrays.asList(nonRepeatedNames.toArray(new String[nonRepeatedNames.size()])));
    }

    /**
     * Generate testing sample names.
     *
     * <p>
     *     Basically all have a common prefix "SAMPLE_" followed by a numeric index.
     * </p>
     *
     * <p>
     *     With {@code maxSampleIndex} you can force to have some repeated sample names;
     *     (if {@code sampleCount < maxSampleIndex}.
     * </p>
     *
     * @param sampleCount number of sample names to generate.
     * @param maxSampleIndex the maximum sample numeric index.
     *
     * @throws RuntimeException if {@code sampleCount} or {@code maxSampleIndex} are negative.
     * @return never {@code null}.
     */
    private String[] generateSampleNames(final int sampleCount, final int maxSampleIndex) {
        final String[] result = new String[sampleCount];
        for (int i = 0; i < sampleCount; i++)
            result[i] = "SAMPLE_" + rnd.nextInt(maxSampleIndex + 1);
        return result;
    }

    private static final int[] SAMPLE_COUNT = { 0, 1, 5, 10, 20};

    private static final int[] MAX_SAMPLE_INDEX = { 0, 1, 4, 9, 10000};

    private static final Random rnd = Utils.getRandomGenerator();


    @DataProvider(name="sampleCountMaxSampleIndexData")
    public Object[][] sampleCountMaxSampleIndexData() {
        final Object[][] result = new Object[SAMPLE_COUNT.length * MAX_SAMPLE_INDEX.length][];
        int nextIndex = 0;
        for (int i = 0; i < SAMPLE_COUNT.length; i++)
            for (int j = 0; j < MAX_SAMPLE_INDEX.length; j++)
                result[nextIndex++] = new Object[] { SAMPLE_COUNT[i], MAX_SAMPLE_INDEX[j]};
        return result;
    }



}
