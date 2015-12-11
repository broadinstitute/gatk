package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.parser.CommandLineParser;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class IntervalArgumentCollectionTest extends BaseTest{

    @DataProvider(name = "optionalOrNot")
    public Object[][] optionalOrNot(){
        return new Object[][]{
                { new OptionalIntervalArgumentCollection()},
                { new RequiredIntervalArgumentCollection()}
        };
    }

    private static class WithOptionalIntervals{
        @ArgumentCollection
        IntervalArgumentCollection iac = new OptionalIntervalArgumentCollection();
    }

    private static class WithRequiredIntervals{
        @ArgumentCollection
        IntervalArgumentCollection iac = new RequiredIntervalArgumentCollection();
    }

    @Test
    public void testOptionalIsOptional(){
        WithOptionalIntervals opt = new WithOptionalIntervals();
        CommandLineParser clp = new CommandLineParser(opt);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test(expectedExceptions = UserException.class)
    public void testRequiredIsRequired(){
        WithRequiredIntervals opt = new WithRequiredIntervals();
        CommandLineParser clp = new CommandLineParser(opt);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }


    @Test(dataProvider = "optionalOrNot",expectedExceptions = GATKException.class)
    public void emptyIntervalsTest(IntervalArgumentCollection iac){
        Assert.assertFalse(iac.intervalsSpecified());
        iac.getIntervals(hg19GenomeLocParser);  //should throw an exception
    }

    @Test(dataProvider = "optionalOrNot")
    public void testExcludeWithNoIncludes(IntervalArgumentCollection iac){
        iac.excludeIntervalStrings.addAll(Arrays.asList("1", "2", "3"));
        Assert.assertTrue(iac.intervalsSpecified());
        GenomeLoc chr4GenomeLoc = hg19GenomeLocParser.createOverEntireContig("4");
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval(chr4GenomeLoc)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIncludeWithExclude(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.excludeIntervalStrings.add("1:90-200");
        Assert.assertTrue(iac.intervalsSpecified());
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 1, 89)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIntervalSetRule(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.addToIntervalStrings("1:90-200");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 90, 100)));
        iac.intervalSetRule = IntervalSetRule.UNION;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 1, 200)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testPadding(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:20-30");
        iac.intervalPadding = 10;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser), Arrays.asList(new SimpleInterval("1", 10, 40)));
    }

    @Test(dataProvider = "optionalOrNot", expectedExceptions = UserException.BadArgumentValue.class)
    public void testAllExcluded(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:10-20");
        iac.excludeIntervalStrings.add("1:1-200");
        iac.getIntervals(hg19GenomeLocParser);
    }

    @Test(dataProvider = "optionalOrNot", expectedExceptions= UserException.BadArgumentValue.class)
    public void testNoIntersection(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:10-20");
        iac.addToIntervalStrings("1:50-200");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        iac.getIntervals(hg19GenomeLocParser);
    }

    @Test(dataProvider = "optionalOrNot", expectedExceptions = UserException.BadArgumentValue.class)
    public void testUmapped(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("unmapped");
        iac.getIntervals(hg19GenomeLocParser);
    }

}
