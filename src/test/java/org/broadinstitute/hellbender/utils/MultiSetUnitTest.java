package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

public class MultiSetUnitTest {

    @Test
    public void testAddSingly(){
        MultiSet<Integer> intList = new MultiSet<>();
        intList.addAll(2);
        intList.addAll(5);
        intList.addAll(5);
        intList.addAll(2);
        intList.addAll(2);
        intList.addAll(3);
        intList.addAll(2);
        intList.addAll(2);
        intList.addAll(4);
        intList.addAll(4);
        intList.addAll(4);
        intList.addAll(2);

        Assert.assertEquals(intList.isEmpty(), false);
        Assert.assertEquals(intList.toString(), "2,6,3,1,4,3,5,2");
    }

    @Test
    public void testAddValueCounts(){
        MultiSet<Integer> intList = new MultiSet<>();
        intList.addAll(5,2);
        intList.addAll(2,6);
        intList.addAll(3,1);
        intList.addAll(4,3);

        Assert.assertEquals(intList.isEmpty(), false);
        Assert.assertEquals(intList.toString(), "2,6,3,1,4,3,5,2");
    }

    @Test
    public void testAddBothWays(){
        MultiSet<Integer> intList = new MultiSet<>();
        intList.addAll(2);
        intList.addAll(5,2);
        intList.addAll(2);
        intList.addAll(2);
        intList.addAll(3);
        intList.addAll(2);
        intList.addAll(2);
        intList.addAll(4,2);
        intList.addAll(2);
        intList.addAll(4,1);

        Assert.assertEquals(intList.toString(), "2,6,3,1,4,3,5,2");
    }

    @Test
    public void testCombineLists(){
        MultiSet<Integer> intList1 = new MultiSet<>();
        intList1.addAll(5,2);
        intList1.addAll(2,6);
        intList1.addAll(3,1);
        intList1.addAll(4,3);

        MultiSet<Integer> intList2 = new MultiSet<>();
        intList2.addAll(2,5);
        intList2.addAll(6,2);
        intList2.addAll(1,3);
        intList2.addAll(3,4);

        intList1.addAll(intList2);

        Assert.assertEquals(intList1.toString(), "1,3,2,11,3,5,4,3,5,2,6,2");

    }

    @Test
    public void testIterator(){
        MultiSet<Integer> intList1 = new MultiSet<>();
        intList1.addAll(5,2);
        intList1.addAll(2,6);
        intList1.addAll(3,1);
        intList1.addAll(4,3);

        MultiSet<Integer> intList2 = new MultiSet<>();
        for(Integer i : intList1) {
            intList2.addAll(i);
        }

        Assert.assertEquals(intList1.toString(),intList2.toString());
    }

}