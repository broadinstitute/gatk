package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.util.List;

@CommandLineProgramProperties(summary="Compare two interval lists to see if they are equal",
        oneLineSummary = "Compare two interval lists for equality",
        programGroup = IntervalsManipulationProgramGroup.class)
public class CompareIntervalLists extends CommandLineProgram {

    @Argument(fullName ="L")
    public String firstIntervalFile;

    @Argument(fullName ="L2")
    public String secondIntervalFile;

    @ArgumentCollection
    final ReferenceInputArgumentCollection reference = new RequiredReferenceInputArgumentCollection();

    @Override
    public Object doWork() {
        final SAMSequenceDictionary samSequenceDictionary;
        try(final ReferenceDataSource ref = ReferenceDataSource.of(reference.getReferencePath())) {
            samSequenceDictionary = ref.getSequenceDictionary();
        }

        final GenomeLocParser genomeLocParser = new GenomeLocParser(samSequenceDictionary);
        final GenomeLocSortedSet firstIntervals = getGenomeLocs(genomeLocParser, firstIntervalFile);
        final GenomeLocSortedSet secondIntervals = getGenomeLocs(genomeLocParser, secondIntervalFile);
        final String result = IntervalUtils.equateIntervals(firstIntervals.toList(), secondIntervals.toList());
        if(result == null) {
            System.out.println("Intervals are equal");
            return 0;
        } else {
            throw new UserException("Intervals are not equal: \n" + result);
        }
    }

    private static GenomeLocSortedSet getGenomeLocs(final GenomeLocParser genomeLocParser, final String intervalFile) {

        final List<GenomeLoc> genomeLocs = IntervalUtils.parseIntervalArguments(genomeLocParser, intervalFile);
        return  IntervalUtils.sortAndMergeIntervals(genomeLocParser, genomeLocs, IntervalMergingRule.ALL);
    }
}
