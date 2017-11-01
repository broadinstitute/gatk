package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.collections.ListUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class IntervalOverlappingIteratorUnitTest extends GATKBaseTest {

    @DataProvider(name="data")
    public Object[][] getData() {
        // the sequence dictionary
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("1", 1000000));
        dictionary.addSequence(new SAMSequenceRecord("2", 1000000));
        // the set of intervals
        final List<SimpleInterval> intervals_1 = Arrays.asList(new SimpleInterval("1:500-600"),	new SimpleInterval("1:700-800"));
        final List<SimpleInterval> intervals_2 = Arrays.asList(new SimpleInterval("2:100-200"), new SimpleInterval("2:400-1000"));
        // some records
        final SimpleInterval record_1_1_100 = new SimpleInterval("1:1-100");
        final SimpleInterval record_1_1_800 = new SimpleInterval("1:1-800");
        final SimpleInterval record_1_500_600 = new SimpleInterval("1:500-600");
        final SimpleInterval record_1_700_750 = new SimpleInterval("1:700-750");
        final SimpleInterval record_2_100_150 = new SimpleInterval("2:100-150");
        final SimpleInterval record_2_900_999 = new SimpleInterval("2:900-999");
        // test cases
        return new Object[][] {
            // first record starts before the first interval, second record overlaps the first interval
            {intervals_1, dictionary, new SimpleInterval[]{record_1_1_100, record_1_500_600, record_2_900_999}, new SimpleInterval[]{record_1_500_600}},
            // first record starts after the first interval, second interval overlaps the first record
            {intervals_1, dictionary, new SimpleInterval[]{record_1_700_750, record_2_900_999}, new SimpleInterval[]{record_1_700_750}},
            // first interval is on a later contig than the first record, but overlaps later records
            {intervals_2, dictionary, new SimpleInterval[]{record_1_1_100, record_2_900_999}, new SimpleInterval[]{record_2_900_999}},
            // first interval is on an earlier contig than the first record, but later records overlap later intervals
            {ListUtils.union(intervals_1, intervals_2), dictionary, new SimpleInterval[]{record_2_100_150, record_2_900_999}, new SimpleInterval[]{record_2_100_150, record_2_900_999}},
            // no records overlap any intervals
            {intervals_1, dictionary, new SimpleInterval[]{record_2_900_999}, new SimpleInterval[0]},
            // an interval overlaps multiple records
            {intervals_1, dictionary, new SimpleInterval[]{record_1_1_800, record_1_500_600, record_2_900_999}, new SimpleInterval[]{record_1_1_800, record_1_500_600}},
            // a record overlaps multiple intervals
            {intervals_1, dictionary, new SimpleInterval[]{record_1_1_800, record_2_100_150}, new SimpleInterval[]{record_1_1_800}}
        };
    }

    @Test(dataProvider = "data")
    public void testIterator(List<SimpleInterval> intervals, SAMSequenceDictionary dictionary, Locatable[] records,	Locatable[] expected) {
        IntervalOverlappingIterator<Locatable> iterator = new IntervalOverlappingIterator<>(Arrays.asList(records).iterator(), intervals, dictionary);

        final Iterator<Locatable> expectedIterator = Arrays.asList(expected).iterator();
        final List<Locatable> iteratorAsList = new ArrayList<>();
        iterator.forEachRemaining(loc -> {
            iteratorAsList.add(loc);
            // assert that the locations are the expected
            Assert.assertEquals(loc, expectedIterator.next());
        });
        // assert that all the expected values are present
        Assert.assertEquals(iteratorAsList.size(), expected.length);
    }
}