package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.io.File;

@CommandLineProgramProperties(
        summary = "asdf",
        oneLineSummary = "asdf",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class Bam2IntervalFile extends LocusWalker {
    private static final int DEFAULT_MIN_COVERAGE = 10;
    int currentStart = 0;
    int currentLength = 0;
    boolean inDesert = true;
    int minIntervalLength = 10;


    @Argument(doc = "", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputIntervalFile;

    @Argument(doc = "", fullName = "min-coverage")
    private int minCoverage = DEFAULT_MIN_COVERAGE;

    private SAMSequenceDictionary sequenceDictionary;
    private IntervalList intervalList;

    @Override
    public void onTraversalStart(){
        sequenceDictionary = getBestAvailableSequenceDictionary();
        intervalList = new IntervalList(sequenceDictionary);
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final long coverage = MathUtils.sum(alignmentContext.getBasePileup().getBaseCounts());
        if (coverage < minCoverage){
            if (currentLength >= minIntervalLength){
                intervalList.add(new Interval(alignmentContext.getContig(), currentStart, alignmentContext.getStart()));
            }

            inDesert = true;
            currentLength = 0;
        }

        // If we get here, coverage is greater than the minCoverage
        if (inDesert){
            // Start a new interval
            currentStart = alignmentContext.getStart();
            inDesert = false;
            currentLength++;
            return;
        } else {
            // Extend
            currentLength++;
            return;
        }
    }

    @Override
    public Object onTraversalSuccess(){
        intervalList.write(outputIntervalFile);
        return "SUCCESS";
    }
}
