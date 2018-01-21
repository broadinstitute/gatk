package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyReferenceLocations;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFReader;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

public class SvDiscoveryUtils {

    //==================================================================================================================

    public static void evaluateIntervalsAndNarls(final List<SVInterval> assembledIntervals,
                                                 final List<NovelAdjacencyReferenceLocations> narls,
                                                 final SAMSequenceDictionary referenceSequenceDictionary,
                                                 final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters,
                                                 final Logger toolLogger) {
        if ( parameters.truthVCF != null ) {
            final SVIntervalTree<String> trueBreakpoints =
                    SVVCFReader.readBreakpointsFromTruthVCF(parameters.truthVCF, referenceSequenceDictionary, parameters.truthIntervalPadding);

            if ( assembledIntervals != null ) {
                evaluateIntervalsAgainstTruth(assembledIntervals, trueBreakpoints, toolLogger);
            }

            final SVIntervalTree<String> narlyBreakpoints =
                    readBreakpointsFromNarls(narls, referenceSequenceDictionary, parameters.truthIntervalPadding);

            evaluateNarlsAgainstTruth(narlyBreakpoints, trueBreakpoints, toolLogger);
        }
    }

    private static SVIntervalTree<String> readBreakpointsFromNarls(final List<NovelAdjacencyReferenceLocations> narls,
                                                                   final SAMSequenceDictionary dictionary,
                                                                   final int breakpointPadding) {
        final SVIntervalTree<String> breakpoints = new SVIntervalTree<>();
        for ( final NovelAdjacencyReferenceLocations narl : narls ) {
            final int padding = breakpointPadding + narl.complication.getLength();

            final SimpleInterval si1 = narl.leftJustifiedLeftRefLoc;
            breakpoints.put(
                    new SVInterval(dictionary.getSequenceIndex(si1.getContig()), si1.getStart()-padding, si1.getStart()+padding), null);

            final SimpleInterval si2 = narl.leftJustifiedRightRefLoc;
            breakpoints.put(
                    new SVInterval(dictionary.getSequenceIndex(si2.getContig()), si2.getStart()-padding, si2.getStart()+padding), null);
        }
        return breakpoints;
    }

    private static void evaluateNarlsAgainstTruth(final SVIntervalTree<String> narlyBreakpoints,
                                                  final SVIntervalTree<String> trueBreakpoints,
                                                  final Logger localLogger) {
        final float falsePos = 1.f - narlyBreakpoints.overlapFraction(trueBreakpoints);
        final int nNarly = narlyBreakpoints.size();
        localLogger.info("Breakpoint false positive rate = " + falsePos + " (" + Math.round(falsePos*nNarly) + "/" + nNarly + ")");
        final float falseNeg = 1.f - trueBreakpoints.overlapFraction(narlyBreakpoints);
        final int nTrue = trueBreakpoints.size();
        localLogger.info("Breakpoint false negative rate = " + falseNeg + " (" + Math.round(falseNeg*nTrue) + "/" + nTrue + ")");
    }

    private static void evaluateIntervalsAgainstTruth(final List<SVInterval> assembledIntervals,
                                                      final SVIntervalTree<String> trueBreakpoints,
                                                      final Logger localLogger) {
        final SVIntervalTree<Integer> intervals = new SVIntervalTree<>();
        final int nIntervals = assembledIntervals.size();
        for ( int idx = 0; idx != nIntervals; ++idx ) {
            intervals.put(assembledIntervals.get(idx), idx);
        }
        final float falsePos = 1.f - intervals.overlapFraction(trueBreakpoints);
        localLogger.info("Interval false positive rate = " + falsePos + " (" + Math.round(falsePos*nIntervals) + "/" + nIntervals + ")");
        final float falseNeg = 1.f - trueBreakpoints.overlapFraction(intervals);
        final int nTrue = trueBreakpoints.size();
        localLogger.info("Interval false negative rate = " + falseNeg + " (" + Math.round(falseNeg*nTrue) + "/" + nTrue + ")");
    }
}
