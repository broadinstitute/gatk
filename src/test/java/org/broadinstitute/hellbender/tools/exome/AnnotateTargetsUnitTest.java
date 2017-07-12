package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link AnnotateTargets}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AnnotateTargetsUnitTest extends BaseTest {

    @Test(dataProvider = "calculateBaitCountsTestData")
    public void testCalculateBaitCounts(final List<Target> targets,
                                        final List<Set<SimpleInterval>> overlappingBaitSetPerTarget,
                                        final double[] expected) {
        final double[] result = AnnotateTargets.calculateBaitCounts(targets, overlappingBaitSetPerTarget);
        ArrayAsserts.assertArrayEquals(expected, result, 1e-12);
    }

    @DataProvider(name = "calculateBaitCountsTestData")
    public Object[][] getCalculateBaitCountsTestData() {
        final List<Target> targets = Arrays.stream(new Target[] {
                new Target(new SimpleInterval("1", 10, 50)),
                new Target(new SimpleInterval("1", 60, 100)),
                new Target(new SimpleInterval("2", 20, 30))
        }).collect(Collectors.toList());

        final SimpleInterval bait_1 = new SimpleInterval("1", 15, 40);
        final SimpleInterval bait_2 = new SimpleInterval("1", 40, 70);
        final SimpleInterval bait_3 = new SimpleInterval("1", 10, 60);
        final SimpleInterval bait_4 = new SimpleInterval("2", 30, 40);

        return new Object[][] {
                {targets, Arrays.asList(
                        new HashSet<>(Arrays.asList(bait_1)),
                        new HashSet<>(),
                        new HashSet<>()),
                        new double[] {1.0, 0.0, 0.0}},
                {targets, Arrays.asList(
                        new HashSet<>(Arrays.asList(bait_1, bait_2)),
                        new HashSet<>(Arrays.asList(bait_2)),
                        new HashSet<>()),
                        new double[] {1.5, 0.5, 0.0}},
                {targets, Arrays.asList(
                        new HashSet<>(Arrays.asList(bait_1, bait_2)),
                        new HashSet<>(Arrays.asList(bait_2)),
                        new HashSet<>(Arrays.asList(bait_4))),
                        new double[] {1.5, 0.5, 1.0}},
                {targets, Arrays.asList(
                        new HashSet<>(Arrays.asList(bait_1, bait_2)),
                        new HashSet<>(Arrays.asList(bait_2)),
                        new HashSet<>(Arrays.asList(bait_4))),
                        new double[] {1.5, 0.5, 1.0}},
                {targets, Arrays.asList(
                        new HashSet<>(Arrays.asList(bait_1, bait_2)),
                        new HashSet<>(Arrays.asList(bait_2, bait_3)),
                        new HashSet<>(Arrays.asList(bait_4))),
                        new double[] {1.5, 1.5, 1.0}}
        };
    }
}
