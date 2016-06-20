package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TargetPadderUnitTest extends BaseTest {

    @Test(dataProvider = "simple")
    public void testSimplePaddingWithoutOverlap(List<Target> targets){
        final List<Target> tc = TargetPadder.padTargets(targets, 100);
        Assert.assertTrue(tc.get(0).equals(new Target("target1", new SimpleInterval("1", 100, 275))));
        Assert.assertTrue(tc.get(1).equals(new Target("target2", new SimpleInterval("1", 276, 450))));
        Assert.assertTrue(tc.get(2).equals(new Target("target3", new SimpleInterval("1", 800, 1050))));
    }

    @Test(dataProvider = "simpleNearStartOfChromosome")
    public void testSimplePaddingWithZero(List<Target> targets){
        final List<Target> tc = TargetPadder.padTargets(targets, 100);
        Assert.assertTrue(tc.get(0).equals(new Target("target1", new SimpleInterval("1", 1, 275))));
        Assert.assertTrue(tc.get(1).equals(new Target("target2", new SimpleInterval("1", 276, 450))));
        Assert.assertTrue(tc.get(2).equals(new Target("target3", new SimpleInterval("1", 800, 1050))));
    }

    @DataProvider(name="simple")
    public Object [][] simpleData() {
        Target t1 = new Target("target1", new SimpleInterval("1", 200, 250));
        Target t2 = new Target("target2", new SimpleInterval("1", 300, 350));
        Target t3 = new Target("target3", new SimpleInterval("1", 900, 950));
        final List<Target> targetList = Arrays.asList(t1,t2,t3);
        final List<Object[]> result = new ArrayList<>();

        result.add(new Object[]{targetList});
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="simpleNearStartOfChromosome")
    public Object [][] simpleDataNearStartOfContig() {
        Target t1 = new Target("target1", new SimpleInterval("1", 20, 250));
        Target t2 = new Target("target2", new SimpleInterval("1", 300, 350));
        Target t3 = new Target("target3", new SimpleInterval("1", 900, 950));
        final List<Target> targetList = Arrays.asList(t1,t2,t3);
        final List<Object[]> result = new ArrayList<>();

        result.add(new Object []{targetList});
        return result.toArray(new Object[result.size()][]);
    }
}
