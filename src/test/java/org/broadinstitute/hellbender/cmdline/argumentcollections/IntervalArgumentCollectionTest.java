package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class IntervalArgumentCollectionTest extends GATKBaseTest {

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
        CommandLineParser clp = new CommandLineArgumentParser(opt);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testRequiredIsRequired(){
        WithRequiredIntervals opt = new WithRequiredIntervals();
        CommandLineParser clp = new CommandLineArgumentParser(opt);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test(dataProvider = "optionalOrNot",expectedExceptions = GATKException.class)
    public void emptyIntervalsTest(IntervalArgumentCollection iac){
        Assert.assertFalse(iac.intervalsSpecified());
        iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary());  //should throw an exception
    }

    @Test(dataProvider = "optionalOrNot")
    public void testExcludeWithNoIncludes(IntervalArgumentCollection iac){
        iac.excludeIntervalStrings.addAll(Arrays.asList("1", "2", "3"));
        Assert.assertTrue(iac.intervalsSpecified());
        GenomeLoc chr4GenomeLoc = hg19GenomeLocParser.createOverEntireContig("4");
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval(chr4GenomeLoc)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testExcludeWithPadding(IntervalArgumentCollection iac){
        iac.intervalExclusionPadding = 10;
        iac.addToIntervalStrings("1:1-100");
        iac.excludeIntervalStrings.add("1:90-100");
        Assert.assertTrue(iac.intervalsSpecified());
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval("1", 1, 79)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIncludeWithExclude(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.excludeIntervalStrings.add("1:90-200");
        Assert.assertTrue(iac.intervalsSpecified());
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval("1", 1, 89)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIncludeWithPadding(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:20-30");
        iac.intervalPadding = 10;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval("1", 10, 40)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIntervalMergingRuleAdjacentMerge(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.addToIntervalStrings("1:101-200");
        iac.intervalMergingRule = IntervalMergingRule.ALL;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval("1", 1, 200)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIntervalMergingRuleAdjacentNoMerge(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.addToIntervalStrings("1:101-200");
        iac.intervalMergingRule = IntervalMergingRule.OVERLAPPING_ONLY;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()).size(), 2);
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()).get(0), new SimpleInterval("1", 1, 100));
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()).get(1), new SimpleInterval("1", 101, 200));
    }

    /**
     * Asserts that the interval set rule is applied first, then the interval ordering rule. This should give an error because the overlap is empty.
     * @param iac
     */
    @Test(dataProvider = "optionalOrNot", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testIntervalSetAndMergingOverlap(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.addToIntervalStrings("1:101-200");
        iac.addToIntervalStrings("1:90-110");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        iac.intervalMergingRule = IntervalMergingRule.ALL;
        iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary());
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIntervalSetRuleIntersection(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.addToIntervalStrings("1:90-200");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval("1", 90, 100)));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testIntervalSetRuleUnion(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:1-100");
        iac.addToIntervalStrings("1:90-200");
        iac.intervalSetRule = IntervalSetRule.UNION;
        Assert.assertEquals(iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary()), Arrays.asList(new SimpleInterval("1", 1, 200)));
    }


    @Test(dataProvider = "optionalOrNot", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testAllExcluded(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:10-20");
        iac.excludeIntervalStrings.add("1:1-200");
        iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary());
    }

    @Test(dataProvider = "optionalOrNot", expectedExceptions= CommandLineException.BadArgumentValue.class)
    public void testNoIntersection(IntervalArgumentCollection iac){
        iac.addToIntervalStrings("1:10-20");
        iac.addToIntervalStrings("1:50-200");
        iac.intervalSetRule = IntervalSetRule.INTERSECTION;
        iac.getIntervals(hg19GenomeLocParser.getSequenceDictionary());
    }

    @Test(dataProvider = "optionalOrNot")
    public void testUnmappedInclusion(IntervalArgumentCollection iac) {
        iac.addToIntervalStrings("unmapped");
        final TraversalParameters traversalParameters = iac.getTraversalParameters(hg19GenomeLocParser.getSequenceDictionary());
        Assert.assertTrue(traversalParameters.traverseUnmappedReads());
        Assert.assertTrue(traversalParameters.getIntervalsForTraversal().isEmpty());
    }

    @Test(dataProvider = "optionalOrNot")
    public void testUnmappedAndMappedInclusion(IntervalArgumentCollection iac) {
        iac.addToIntervalStrings("1:10-20");
        iac.addToIntervalStrings("2:1-5");
        iac.addToIntervalStrings("unmapped");
        final TraversalParameters traversalParameters = iac.getTraversalParameters(hg19GenomeLocParser.getSequenceDictionary());
        Assert.assertTrue(traversalParameters.traverseUnmappedReads());
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().size(), 2);
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().get(0), new SimpleInterval("1", 10, 20));
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().get(1), new SimpleInterval("2", 1, 5));
    }

    @Test(dataProvider = "optionalOrNot")
    public void testUnmappedAndMappedInclusionPlusMappedExclusion(IntervalArgumentCollection iac) {
        iac.addToIntervalStrings("1:10-20");
        iac.addToIntervalStrings("2:1-5");
        iac.addToIntervalStrings("unmapped");
        iac.excludeIntervalStrings.addAll(Arrays.asList("1"));
        final TraversalParameters traversalParameters = iac.getTraversalParameters(hg19GenomeLocParser.getSequenceDictionary());
        Assert.assertTrue(traversalParameters.traverseUnmappedReads());
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().size(), 1);
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().get(0), new SimpleInterval("2", 1, 5));
    }

    @Test(dataProvider = "optionalOrNot", expectedExceptions = UserException.class)
    public void testThrowOnUnmappedExclusion(IntervalArgumentCollection iac) {
        iac.excludeIntervalStrings.addAll(Arrays.asList("unmapped"));
        iac.getTraversalParameters(hg19GenomeLocParser.getSequenceDictionary());
    }

    @Test(dataProvider = "optionalOrNot")
    public void testMultipleUnmappedInclusion(IntervalArgumentCollection iac) {
        iac.addToIntervalStrings("unmapped");
        iac.addToIntervalStrings("1:10-20");
        iac.addToIntervalStrings("unmapped");
        iac.addToIntervalStrings("unmapped");
        final TraversalParameters traversalParameters = iac.getTraversalParameters(hg19GenomeLocParser.getSequenceDictionary());
        Assert.assertTrue(traversalParameters.traverseUnmappedReads());
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().size(), 1);
        Assert.assertEquals(traversalParameters.getIntervalsForTraversal().get(0), new SimpleInterval("1", 10, 20));
    }
}
