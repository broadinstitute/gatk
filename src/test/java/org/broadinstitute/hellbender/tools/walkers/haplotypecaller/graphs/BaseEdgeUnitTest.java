package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class BaseEdgeUnitTest extends BaseTest {
    @DataProvider(name = "EdgeCreationData")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int multiplicity : Arrays.asList(1, 2, 3) ) {
            for ( final boolean isRef : Arrays.asList(true, false) ) {
                tests.add(new Object[]{isRef, multiplicity});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "EdgeCreationData")
    public void testBasic(final boolean isRef, final int mult) {
        final BaseEdge e = new BaseEdge(isRef, mult);
        Assert.assertEquals(e.isRef(), isRef);
        Assert.assertEquals(e.getMultiplicity(), mult);
        Assert.assertEquals(e.getPruningMultiplicity(), mult);
        Assert.assertEquals(e.getDotLabel(), Integer.toString(mult));

        e.toString();//just check not blowing up

        e.setIsRef(!isRef);
        Assert.assertEquals(e.isRef(), !isRef);

        e.toString();//just check not blowing up

        e.setMultiplicity(mult + 1);
        Assert.assertEquals(e.getMultiplicity(), mult + 1);
        Assert.assertEquals(e.getPruningMultiplicity(), mult + 1);
        Assert.assertEquals(e.getDotLabel(), Integer.toString(mult + 1));

        e.toString();//just check not blowing up

        e.incMultiplicity(2);
        Assert.assertEquals(e.getMultiplicity(), mult + 3);
        Assert.assertEquals(e.getPruningMultiplicity(), mult + 3);
        Assert.assertEquals(e.getDotLabel(), Integer.toString(mult + 3));

        e.toString();//just check not blowing up

        final BaseEdge copy = e.copy();
        Assert.assertEquals(copy.isRef(), e.isRef());
        Assert.assertEquals(copy.getMultiplicity(), e.getMultiplicity());
        Assert.assertEquals(copy.getPruningMultiplicity(), e.getPruningMultiplicity());
        Assert.assertEquals(copy.getDotLabel(), e.getDotLabel());

        e.toString();//just check not blowing up
    }

    @Test(dataProvider = "EdgeCreationData")
    public void testAdd(final boolean isRef, final int mult) {
        final BaseEdge e1 = new BaseEdge(isRef, mult);
        final BaseEdge e2 = new BaseEdge(isRef, mult);
        final BaseEdge e3 = e1.add(e2);
        Assert.assertTrue(e1 == e3);//identity
        Assert.assertEquals(e1.isRef(), isRef);
        Assert.assertEquals(e1.getMultiplicity(), mult*2);
        Assert.assertEquals(e1.getPruningMultiplicity(), mult*2);
        Assert.assertEquals(e1.getDotLabel(), Integer.toString(mult*2));

        final BaseEdge e4 = new BaseEdge(!isRef, mult);
        e1.add(e4);
        Assert.assertEquals(e1.isRef(), true); //one or the other was ref
    }

    @Test
    public void testAddOr() {
        final BaseEdge e1f = new BaseEdge(false, 1);
        final BaseEdge e2f = new BaseEdge(false, 2);
        final BaseEdge e1e2 = BaseEdge.makeOREdge(Arrays.asList(e1f, e2f), 4);
        Assert.assertEquals(e1e2.getMultiplicity(), 4);
        Assert.assertEquals(e1e2.isRef(), false);

        final BaseEdge e3t = new BaseEdge(true, 3);
        final BaseEdge e1e2e3 = BaseEdge.makeOREdge(Arrays.asList(e1f, e2f, e3t), 4);
        Assert.assertEquals(e1e2e3.getMultiplicity(), 4);
        Assert.assertEquals(e1e2e3.isRef(), true);
    }

        @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddNull() {
        final BaseEdge e = new BaseEdge(false, 1);
        e.add(null);
    }

    @Test
    public void testEdgeWeightComparator() {
        final BaseEdge e10 = new BaseEdge(false, 10);
        final BaseEdge e5 = new BaseEdge(true, 5);
        final BaseEdge e2 = new BaseEdge(false, 2);
        final BaseEdge e1 = new BaseEdge(false, 1);

        final List<BaseEdge> edges = new ArrayList<>(Arrays.asList(e1, e2, e5, e10));
        Collections.sort(edges, BaseEdge.EDGE_MULTIPLICITY_ORDER);
        Assert.assertEquals(edges.get(0), e10);
        Assert.assertEquals(edges.get(1), e5);
        Assert.assertEquals(edges.get(2), e2);
        Assert.assertEquals(edges.get(3), e1);
    }
}
