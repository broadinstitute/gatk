package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public class IntervalArgumentCollectionTest extends BaseTest{

    @Test(expectedExceptions = GATKException.class)
    public void emptyIntervalsTest(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        Assert.assertFalse(iac.intervalsSpecified());
        iac.getIntervals(hg19GenomeLocParser);  //should throw an exception
    }

    @Test
    public void testExcludeWithNoIncludes(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.excludeIntervalStrings.addAll(Arrays.asList("1", "2", "3"));
        Assert.assertTrue(iac.intervalsSpecified());
        GenomeLoc chr4GenomeLoc = hg19GenomeLocParser.createOverEntireContig("4");
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval(chr4GenomeLoc)));
    }

    @Test
    public void testIncludeWithExclude(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.intervalStrings.add("1:1-100");
        iac.excludeIntervalStrings.add("1:90-200");
        Assert.assertTrue(iac.intervalsSpecified());
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 1, 89)));
    }

    @Test
    public void testIntervalSetRule(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.intervalStrings.add("1:1-100");
        iac.intervalStrings.add("1:90-200");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 90, 100)));
        iac.intervalSetRule = IntervalSetRule.UNION;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 1, 200)));
    }

    @Test
    public void testPadding(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.intervalStrings.add("1:20-30");
        iac.intervalPadding = 10;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 10, 40)));
    }

    @Test( expectedExceptions = UserException.BadArgumentValue.class)
    public void testAllExcluded(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.intervalStrings.add("1:10-20");
        iac.excludeIntervalStrings.add("1:1-200");
        iac.getIntervals(hg19GenomeLocParser);
    }

    @Test( expectedExceptions= UserException.BadArgumentValue.class)
    public void testNoIntersection(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.intervalStrings.add("1:10-20");
        iac.intervalStrings.add("1:50-200");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        iac.getIntervals(hg19GenomeLocParser);
    }

    @Test( expectedExceptions = UserException.BadArgumentValue.class)
    public void testUmapped(){
        IntervalArgumentCollection iac = new IntervalArgumentCollection();
        iac.intervalStrings.add("unmapped");
        iac.getIntervals(hg19GenomeLocParser);
    }

}
