package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.CompressedDataList;
import org.testng.Assert;
import org.testng.annotations.Test;

public class CompressedDataListUnitTest {

    @Test
    public void testAddSingly(){
        CompressedDataList<Integer> intList = new CompressedDataList<>();
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
        Assert.assertEquals(intList.toString(), "2,6,3,1,4,3,5,2");
    }

    @Test
    public void testAddValueCounts(){
        CompressedDataList<Integer> intList = new CompressedDataList<>();
        intList.add(5,2);
        intList.add(2,6);
        intList.add(3,1);
        intList.add(4,3);

        Assert.assertEquals(intList.isEmpty(), false);
        Assert.assertEquals(intList.toString(), "2,6,3,1,4,3,5,2");
    }

    @Test
    public void testAddBothWays(){
        CompressedDataList<Integer> intList = new CompressedDataList<>();
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

        Assert.assertEquals(intList.toString(), "2,6,3,1,4,3,5,2");
    }

    @Test
    public void testCombineLists(){
        CompressedDataList<Integer> intList1 = new CompressedDataList<>();
        intList1.add(5,2);
        intList1.add(2,6);
        intList1.add(3,1);
        intList1.add(4,3);

        CompressedDataList<Integer> intList2 = new CompressedDataList<>();
        intList2.add(2,5);
        intList2.add(6,2);
        intList2.add(1,3);
        intList2.add(3,4);

        intList1.add(intList2);

        Assert.assertEquals(intList1.toString(), "1,3,2,11,3,5,4,3,5,2,6,2");

    }

    @Test
    public void testIterator(){
        CompressedDataList<Integer> intList1 = new CompressedDataList<>();
        intList1.add(5,2);
        intList1.add(2,6);
        intList1.add(3,1);
        intList1.add(4,3);

        CompressedDataList<Integer> intList2 = new CompressedDataList<>();
        for(Integer i : intList1) {
            intList2.add(i);
        }

        Assert.assertEquals(intList1.toString(),intList2.toString());
    }

}