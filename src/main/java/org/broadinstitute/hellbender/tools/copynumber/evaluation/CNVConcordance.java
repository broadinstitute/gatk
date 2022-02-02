package org.broadinstitute.hellbender.tools.copynumber.evaluation;

import htsjdk.samtools.util.CoordMath;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Supported Data Types
 * .called.seg (ModelSegments)
 * ??? (ichorCNA)
 */
@CommandLineProgramProperties(
        summary = "CNV Concordance",
        oneLineSummary = "CNV Concordance",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class CNVConcordance extends GATK {

    @Argument(doc = "Eval", fullName = "eval")
    private File evalFile;

    @Argument(doc = "Truth", fullName = "truth")
    private File truthFile;

    // Truth calls annotated with eval calls
    @Argument(doc = "output1", fullName = "output1")
    private File output1;

    // Eval calls annotated with truth calls
    @Argument(doc = "output2", fullName = "output2")
    private File output2;

    // Eval calls annotated with truth calls
    @Argument(doc = "output3", fullName = "output3")
    private File output3;

    @Argument(doc = "minimum event size", fullName = "min-event-size", optional = true)
    private int minEventSize = 50;

    @Argument(doc = "An eval segment that overlap a truth segment more than this fraction is considered concordant",
            fullName = "min-overlapping-bases", optional = true)
    private double minOverlappingFractionThreshold = 0.3;

    @Override
    public void traverse() {
        final CalledCopyRatioSegmentCollection evalSegments = new CalledCopyRatioSegmentCollection(evalFile);
        final CalledCopyRatioSegmentCollection truthSegments = new CalledCopyRatioSegmentCollection(truthFile);

        final CalledCopyRatioSegmentCollection filteredTruthSegments = applySegmentFilters(truthSegments, minEventSize);
        // TODO: inside createOverlapMap, change to a TreeMap, which entails implementing Comparable within CalledCopyRatioSegment
        // This is to be implemented. But the per-eval segment map should be useful.
        final Map<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> mapEvalAsKey =
                IntervalUtils.createOverlapMap(evalSegments.getRecords(), truthSegments.getRecords(),
                        getBestAvailableSequenceDictionary());
        Utils.validate(evalSegments.size() == mapEvalAsKey.size(), "");

        final Map<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> perTruthOverlaps =
                IntervalUtils.createOverlapMap(filteredTruthSegments.getRecords(), evalSegments.getRecords(),
                        getBestAvailableSequenceDictionary());
        Utils.validate(filteredTruthSegments.size() == perTruthOverlaps.size(), "");

        // Use AnnotatedIntervalCollection to output the result (for now)
        final AnnotationKey<String> TRUTH_CALL = new AnnotationKey<>("truth_call", String.class, s -> s.length() > 0);
        final AnnotationKey<Integer> NUM_EVAL_SEGMENTS_KEY = new AnnotationKey<>("num_eval_segments", Integer.class, n -> n >= 0);
        final AnnotationKey<Integer> NUM_OVERLAPPING_EVAL_BASES_KEY = new AnnotationKey<>("num_overlapping_eval_bases", Integer.class, n -> n >= 0);
        final AnnotationKey<Integer> EVENT_WIDTH_KEY = new AnnotationKey<>("event_width", Integer.class, n -> n > 0);
        final AnnotationKey<Double> FRACTION_OVERLAP_KEY = new AnnotationKey<>("fraction_overlap", Double.class, p -> p >= 0);
        final AnnotationKey<String> STATUS_KEY = new AnnotationKey<>("status", String.class, s -> s.length() > 0);

        final List<AnnotatedInterval> result = new ArrayList<>();
        for (final Map.Entry<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> entry : perTruthOverlaps.entrySet()) {
            final CalledCopyRatioSegment truthSegment = entry.getKey();
            final List<CalledCopyRatioSegment> overlappingEvalSegments = entry.getValue();
            final CalledCopyRatioSegment.Call truthCall = truthSegment.getCall();

            // How about return { num segments, overlapping bases, width }
            int numOverlappingEvalSegments = overlappingEvalSegments.size();
            int numOverlappingBases = 0;
            for (final CalledCopyRatioSegment overlappingEvalSegment : overlappingEvalSegments){
                final CalledCopyRatioSegment.Call evalCall = overlappingEvalSegment.getCall();
                // Is this the right way?
                if (evalCall != truthCall){
                    continue;
                }

                // Must check that the call is the same between tr and ev
                final int overlappingBases = CoordMath.getOverlap(truthSegment.getStart(), truthSegment.getEnd(),
                        overlappingEvalSegment.getStart(), overlappingEvalSegment.getEnd());
                numOverlappingBases += overlappingBases;
            }
            // +1 because the interval includes both ends. e.g. [10, 20] contains 11 bases
            final int eventWidth = truthSegment.getEnd() - truthSegment.getStart() + 1;
            Utils.validate(eventWidth > 0, "event width must be positive but got " + eventWidth);
            final double fractionOverlap = numOverlappingBases/ Math.max(eventWidth, 1.0);
            final String status = getStatus(truthCall, fractionOverlap);

            final List<Pair<AnnotationKey<?>, Object>> annotationValues = new ArrayList<>();
            annotationValues.add(new ImmutablePair<>(TRUTH_CALL, truthCall.toString()));
            annotationValues.add(new ImmutablePair<>(NUM_EVAL_SEGMENTS_KEY, numOverlappingEvalSegments));
            annotationValues.add(new ImmutablePair<>(NUM_OVERLAPPING_EVAL_BASES_KEY, numOverlappingBases));
            annotationValues.add(new ImmutablePair<>(EVENT_WIDTH_KEY, eventWidth));
            annotationValues.add(new ImmutablePair<>(FRACTION_OVERLAP_KEY, fractionOverlap));
            annotationValues.add(new ImmutablePair<>(STATUS_KEY, status));

            // Which AnnotatedInterval should I use?
            final AnnotatedInterval annotatedInterval = new AnnotatedInterval(truthSegment.getInterval(),
                    new AnnotationMap(annotationValues));
            result.add(annotatedInterval);
        }

        final AnnotatedIntervalCollection aic = new AnnotatedIntervalCollection(truthSegments.getMetadata(), result);

        // These are things that R would be much better at computing; should I just not return this information?
        final long numAmplificationsInTruth = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(TRUTH_CALL).equals(CalledCopyRatioSegment.Call.AMPLIFICATION.toString()))
                .count();
        final long numAmplificationsTPs = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(TRUTH_CALL).equals(CalledCopyRatioSegment.Call.AMPLIFICATION.toString()))
                .filter(r -> r.getAnnotationMap().getValue(STATUS_KEY).equals("TP"))
                .count();

        final long numDeletionsInTruth = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(TRUTH_CALL).equals(CalledCopyRatioSegment.Call.DELETION.toString()))
                .count();
        final long numDeletionsTPs = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(TRUTH_CALL).equals(CalledCopyRatioSegment.Call.DELETION.toString()))
                .filter(r -> r.getAnnotationMap().getValue(STATUS_KEY).equals("TP"))
                .count();

        final double amplificationSensitivity = numAmplificationsInTruth == 0 ? 0 : (double) numAmplificationsTPs / numAmplificationsInTruth;
        final double deletionSensitivity = numDeletionsInTruth == 0 ? 0 : (double) numDeletionsTPs / numDeletionsInTruth;

        aic.write(output1);

        // Write a class of these metrics---eventually
        try (final PrintWriter pw = new PrintWriter(output3)){
            pw.println("amplification_sensitivity\tdeletion_sensitivity");
            pw.println(String.format("%f\t%f", amplificationSensitivity, deletionSensitivity));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }

    private String getStatus(CalledCopyRatioSegment.Call truthCall, double fractionOverlap) {
        if (truthCall == CalledCopyRatioSegment.Call.NEUTRAL){
            return fractionOverlap > minOverlappingFractionThreshold ? "TN" : "FP"; // This is probably right but still feels funny.
        } else {
            return fractionOverlap > minOverlappingFractionThreshold ? "TP" : "FN";
        }

    }

    private CalledCopyRatioSegmentCollection applySegmentFilters(final CalledCopyRatioSegmentCollection truthSegments,
                                                                 final int minimumEventSize) {
        List<CalledCopyRatioSegment> filteredSegments = truthSegments.getRecords().stream()
                .filter(seg -> seg.getEnd() - seg.getStart() > minimumEventSize)
                .collect(Collectors.toList());
        return new CalledCopyRatioSegmentCollection(truthSegments.getMetadata(), filteredSegments);


    }
}
