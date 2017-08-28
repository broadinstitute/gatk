package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;

import static org.testng.Assert.*;

public class BarcodeSetByIntervalIteratorTest extends BaseTest {

    @DataProvider(name="iteratorTests")
    public Object[][] getIteratorTestData() {
        final List<Object[]> tests = new ArrayList<>();

        final SVIntervalTree<List<String>> intervalTree1 = new SVIntervalTree<>();
        intervalTree1.put(new SVInterval(1, 100, 200), Arrays.asList("ACTG"));
        final List<Tuple2<SVInterval, Set<String>>> expectedList1 = new ArrayList<>();
        expectedList1.add(new Tuple2<>(new SVInterval(1, 100, 200), new HashSet<>(Arrays.asList("ACTG"))));
        tests.add(new Object[] { 1, intervalTree1, expectedList1 });

        final SVIntervalTree<List<String>> intervalTree2 = new SVIntervalTree<>();
        intervalTree2.put(new SVInterval(1, 100, 200), Arrays.asList("ACTG"));
        intervalTree2.put(new SVInterval(1, 150, 250), Arrays.asList("GCTA"));
        final List<Tuple2<SVInterval, Set<String>>> expectedList2 = new ArrayList<>();
        expectedList2.add(new Tuple2<>(new SVInterval(1, 100, 150), new HashSet<>(Arrays.asList("ACTG"))));
        expectedList2.add(new Tuple2<>(new SVInterval(1, 150, 200), new HashSet<>(Arrays.asList("ACTG", "GCTA"))));
        expectedList2.add(new Tuple2<>(new SVInterval(1, 200, 250), new HashSet<>(Arrays.asList("GCTA"))));
        tests.add(new Object[] { 1, intervalTree2, expectedList2 });

        final SVIntervalTree<List<String>> intervalTree3 = new SVIntervalTree<>();
        intervalTree3.put(new SVInterval(1, 100, 200), Arrays.asList("ACTG"));
        intervalTree3.put(new SVInterval(1, 100, 175), Arrays.asList("CCCC"));
        intervalTree3.put(new SVInterval(1, 150, 250), Arrays.asList("GCTA"));
        final List<Tuple2<SVInterval, Set<String>>> expectedList3 = new ArrayList<>();
        expectedList3.add(new Tuple2<>(new SVInterval(1, 100, 150), new HashSet<>(Arrays.asList("ACTG", "CCCC"))));
        expectedList3.add(new Tuple2<>(new SVInterval(1, 150, 175), new HashSet<>(Arrays.asList("ACTG", "GCTA", "CCCC"))));
        expectedList3.add(new Tuple2<>(new SVInterval(1, 175, 200), new HashSet<>(Arrays.asList("ACTG", "GCTA"))));
        expectedList3.add(new Tuple2<>(new SVInterval(1, 200, 250), new HashSet<>(Arrays.asList("GCTA"))));
        tests.add(new Object[] { 1, intervalTree3, expectedList3 });

        final SVIntervalTree<List<String>> intervalTree4 = new SVIntervalTree<>();
        intervalTree4.put(new SVInterval(1, 100, 200), Arrays.asList("ACTG"));
        intervalTree4.put(new SVInterval(1, 150, 250), Arrays.asList("GCTA"));
        intervalTree4.put(new SVInterval(1, 225, 250), Arrays.asList("AAAA"));
        intervalTree4.put(new SVInterval(1, 225, 275), Arrays.asList("GGGG"));
        final List<Tuple2<SVInterval, Set<String>>> expectedList4 = new ArrayList<>();
        expectedList4.add(new Tuple2<>(new SVInterval(1, 100, 150), new HashSet<>(Arrays.asList("ACTG"))));
        expectedList4.add(new Tuple2<>(new SVInterval(1, 150, 200), new HashSet<>(Arrays.asList("ACTG", "GCTA"))));
        expectedList4.add(new Tuple2<>(new SVInterval(1, 200, 225), new HashSet<>(Arrays.asList("GCTA"))));
        expectedList4.add(new Tuple2<>(new SVInterval(1, 225, 250), new HashSet<>(Arrays.asList("GCTA", "AAAA", "GGGG"))));
        expectedList4.add(new Tuple2<>(new SVInterval(1, 250, 275), new HashSet<>(Arrays.asList("GGGG"))));
        tests.add(new Object[] { 1, intervalTree4, expectedList4 });

        final SVIntervalTree<List<String>> intervalTree5 = new SVIntervalTree<>();
        intervalTree5.put(new SVInterval(1, 100, 200), Arrays.asList("ACTG", "GCTA"));
        final List<Tuple2<SVInterval, Set<String>>> expectedList5 = new ArrayList<>();
        expectedList5.add(new Tuple2<>(new SVInterval(1, 100, 200), new HashSet<>(Arrays.asList("ACTG", "GCTA"))));
        tests.add(new Object[] { 1, intervalTree5, expectedList5 });

        final SVIntervalTree<List<String>> intervalTree6 = new SVIntervalTree<>();
        intervalTree6.put(new SVInterval(1, 100, 200), Arrays.asList("ACTG"));
        intervalTree6.put(new SVInterval(1, 200, 300), Arrays.asList("GCTA"));
        final List<Tuple2<SVInterval, Set<String>>> expectedList6 = new ArrayList<>();
        expectedList6.add(new Tuple2<>(new SVInterval(1, 100, 200), new HashSet<>(Arrays.asList("ACTG"))));
        expectedList6.add(new Tuple2<>(new SVInterval(1, 200, 300), new HashSet<>(Arrays.asList("GCTA"))));
        tests.add(new Object[] { 1, intervalTree6, expectedList6 });

        return tests.toArray(new Object[][]{});

    }

    @Test(dataProvider = "iteratorTests")
    public void testBarcodeIterator(final int contig, final SVIntervalTree<List<String>> intervalTree,
                                    final List<Tuple2<SVInterval, Set<String>>> expectedResult) {
        BarcodeSetByIntervalIterator iterator = new BarcodeSetByIntervalIterator(contig, intervalTree);
        int expectedIdx = 0;
        while(iterator.hasNext()) {
            Assert.assertEquals(iterator.next(), expectedResult.get(expectedIdx));
            expectedIdx++;
        }
        Assert.assertFalse(iterator.hasNext());
    }

}