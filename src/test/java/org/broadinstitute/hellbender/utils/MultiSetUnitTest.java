package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Comparator;

public class MultiSetUnitTest {

    @Test
    public void testAddSingly(){
        MultiSet<Integer> intList = new MultiSet<>();
        intList.add(2);
        intList.add(5);
        intList.add(5);
        intList.add(2);
        intList.add(2);
        intList.add(3);
        intList.add(2);
        intList.add(2);
        intList.add(4);
        intList.add(4);
        intList.add(4);
        intList.add(2);

        Assert.assertEquals(intList.isEmpty(), false);
        Assert.assertEquals(intList.toString(Comparator.naturalOrder(), ",", ":"), "2:6,3:1,4:3,5:2");
    }

    @Test
    public void testAddValueCounts(){
        MultiSet<Integer> intList = new MultiSet<>();
        intList.add(5,2);
        intList.add(2,6);
        intList.add(3,1);
        intList.add(4,3);

        Assert.assertEquals(intList.isEmpty(), false);
        Assert.assertEquals(intList.toString(Comparator.naturalOrder(), ",", ":"), "2:6,3:1,4:3,5:2");
    }

    @Test
    public void testAddBothWays(){
        MultiSet<Integer> intList = new MultiSet<>();
        intList.add(2);
        intList.add(5,2);
        intList.add(2);
        intList.add(2);
        intList.add(3);
        intList.add(2);
        intList.add(2);
        intList.add(4,2);
        intList.add(2);
        intList.add(4,1);

        Assert.assertEquals(intList.toString(Comparator.naturalOrder(), ",", ":"), "2:6,3:1,4:3,5:2");
    }

    @Test
    public void testCombineLists(){
        MultiSet<Integer> intList1 = new MultiSet<>();
        intList1.add(5,2);
        intList1.add(2,6);
        intList1.add(3,1);
        intList1.add(4,3);

        MultiSet<Integer> intList2 = new MultiSet<>();
        intList2.add(2,5);
        intList2.add(6,2);
        intList2.add(1,3);
        intList2.add(3,4);

        intList1.addAll(intList2);

        Assert.assertEquals(intList1.toString(Comparator.naturalOrder(), ",", ":"), "1:3,2:11,3:5,4:3,5:2,6:2");

    }

    @Test
    public void testIterator(){
        MultiSet<Integer> intList1 = new MultiSet<>();
        intList1.add(5,2);
        intList1.add(2,6);
        intList1.add(3,1);
        intList1.add(4,3);

        MultiSet<Integer> intList2 = new MultiSet<>();
        for(Integer i : intList1) {
            intList2.add(i);
        }

        Assert.assertEquals(intList1.toString(Comparator.naturalOrder()), intList2.toString(Comparator.naturalOrder()));
    }

}