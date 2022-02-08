package org.broadinstitute.hellbender.tools.copynumber.evaluation;

import htsjdk.samtools.util.CoordMath;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
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
public class CNVConcordance extends GATKTool {

    @Argument(doc = "Eval", fullName = "eval")
    private File evalFile;

    @Argument(doc = "Truth", fullName = "truth")
    private File truthFile;

    public static final String TRUTH_ANNOTATED_WITH_EVAL_ARG_NAME = "truth-annotated-with-eval";
    public static final String TRUTH_ANNOTATED_WITH_EVAL_SHORT_NAME = "te";
    @Argument(doc = "Truth calls annotated with eval calls", fullName = TRUTH_ANNOTATED_WITH_EVAL_ARG_NAME, shortName = "te")
    private File truthAnnotatedWithEvalFile;

    // To be implemented. Optional for now.
    public static final String EVAL_ANNOTATED_WITH_TRUTH_ARG_NAME = "eval-annotated-with-truth";
    public static final String EVAL_ANNOTATED_WITH_TRUTH_SHORT_NAME = "et";
    @Argument(doc = "Eval calls annotated with truth calls", fullName = EVAL_ANNOTATED_WITH_TRUTH_ARG_NAME, optional = true)
    private File evalAnnotatedWithTruthFile;

    public static final String SUMMARY_ARG_NAME = "summary";
    @Argument(doc = "Sumamry statsitics", fullName = SUMMARY_ARG_NAME)
    private File summaryFile;

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

        // A dictionary with eval segments as key, and overlapping truth segments as value
        final Map<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> perEvalOverlapMap =
                IntervalUtils.createOverlapMap(evalSegments.getRecords(), truthSegments.getRecords(),
                        getBestAvailableSequenceDictionary());
        Utils.validate(evalSegments.size() == perEvalOverlapMap.size(),
                "The number of eval segments must equal the size of perEvalOverlapMap");

        final Map<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> perTruthOverlapMap =
                IntervalUtils.createOverlapMap(filteredTruthSegments.getRecords(), evalSegments.getRecords(),
                        getBestAvailableSequenceDictionary());
        Utils.validate(filteredTruthSegments.size() == perTruthOverlapMap.size(),
                "The number of truth segments must equal the size of perEvalOverlapMap");

        // Use AnnotatedIntervalCollection to output the result (for now)
        final AnnotationKey<String> truthCallKey = new AnnotationKey<>("truth_call", String.class, s -> s.length() > 0);
        final AnnotationKey<Integer> numEvalSegmentsKey = new AnnotationKey<>("num_eval_segments", Integer.class, n -> n >= 0);
        final AnnotationKey<Integer> numOverlappingEvalBasesKey = new AnnotationKey<>("num_overlapping_eval_bases", Integer.class, n -> n >= 0);
        final AnnotationKey<Integer> eventWidthKey = new AnnotationKey<>("event_width", Integer.class, n -> n > 0);
        final AnnotationKey<Double> fractionOverlapKey = new AnnotationKey<>("fraction_overlap", Double.class, p -> p >= 0);
        final AnnotationKey<String> statusKey = new AnnotationKey<>("status", String.class, s -> s.length() > 0);

        // Traverse the truth calls
        // TODO: Extract Method
        final List<AnnotatedInterval> result = new ArrayList<>();
        for (final Map.Entry<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> entry : perTruthOverlapMap.entrySet()) {
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
            final int eventSize = truthSegment.getEnd() - truthSegment.getStart() + 1;
            Utils.validate(eventSize > 0, "event size must be positive but got " + eventSize);
            final double fractionOverlap = numOverlappingBases/ Math.max(eventSize, 1.0);
            final String status = getStatus(truthCall, fractionOverlap);

            final List<Pair<AnnotationKey<?>, Object>> annotationValues = new ArrayList<>();
            annotationValues.add(new ImmutablePair<>(truthCallKey, truthCall.toString()));
            annotationValues.add(new ImmutablePair<>(numEvalSegmentsKey, numOverlappingEvalSegments));
            annotationValues.add(new ImmutablePair<>(numOverlappingEvalBasesKey, numOverlappingBases));
            annotationValues.add(new ImmutablePair<>(eventWidthKey, eventSize));
            annotationValues.add(new ImmutablePair<>(fractionOverlapKey, fractionOverlap));
            annotationValues.add(new ImmutablePair<>(statusKey, status));

            // Which AnnotatedInterval should I use?
            final AnnotatedInterval annotatedInterval = new AnnotatedInterval(truthSegment.getInterval(),
                    new AnnotationMap(annotationValues));
            result.add(annotatedInterval);
        }

        final AnnotatedIntervalCollection aic = new AnnotatedIntervalCollection(truthSegments.getMetadata(), result);

        // TODO: do the same with eval as base. This will allow us to compute FPs.

        // These are things that R would be much better at computing; should I just not return this information?
        final long numAmplificationsInTruth = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(truthCallKey).equals(CalledCopyRatioSegment.Call.AMPLIFICATION.toString()))
                .count();
        final long numAmplificationsTPs = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(truthCallKey).equals(CalledCopyRatioSegment.Call.AMPLIFICATION.toString()))
                .filter(r -> r.getAnnotationMap().getValue(statusKey).equals("TP"))
                .count();

        final long numDeletionsInTruth = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(truthCallKey).equals(CalledCopyRatioSegment.Call.DELETION.toString()))
                .count();
        final long numDeletionsTPs = aic.getRecords().stream()
                .filter(r -> r.getAnnotationMap().getValue(truthCallKey).equals(CalledCopyRatioSegment.Call.DELETION.toString()))
                .filter(r -> r.getAnnotationMap().getValue(statusKey).equals("TP"))
                .count();

        final double amplificationSensitivity = numAmplificationsInTruth == 0 ? 0 : (double) numAmplificationsTPs / numAmplificationsInTruth;
        final double deletionSensitivity = numDeletionsInTruth == 0 ? 0 : (double) numDeletionsTPs / numDeletionsInTruth;
        final long numTruthCalls = numAmplificationsInTruth + numDeletionsInTruth;
        final long numTP = numAmplificationsTPs + numDeletionsTPs;

        aic.write(truthAnnotatedWithEvalFile);

        // Write a class of these metrics---eventually
        try (final PrintWriter pw = new PrintWriter(summaryFile)){
            pw.println("num_truth_calls\ttp\tamplification_sensitivity\tdeletion_sensitivity\t");
            pw.println(String.format("%d\t%d\t%f\t%f",
                    numTruthCalls,
                    numTP,
                    amplificationSensitivity,
                    deletionSensitivity));
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

    /**
     * Subset the truth segments by various criteria:
     * - Event size should be larger than a threshold (minimumEventSize)
     * @param truthSegments
     * @param minimumEventSize
     * @return
     */
    private CalledCopyRatioSegmentCollection applySegmentFilters(final CalledCopyRatioSegmentCollection truthSegments,
                                                                 final int minimumEventSize) {
        List<CalledCopyRatioSegment> filteredSegments = truthSegments.getRecords().stream()
                .filter(seg -> seg.getEnd() - seg.getStart() > minimumEventSize)
                .collect(Collectors.toList());
        return new CalledCopyRatioSegmentCollection(truthSegments.getMetadata(), filteredSegments);


    }
}
