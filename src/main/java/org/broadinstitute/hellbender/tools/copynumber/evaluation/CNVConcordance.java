package org.broadinstitute.hellbender.tools.copynumber.evaluation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CoordMath;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Supported Data Types
 * .called.seg (ModelSegments)
 * ??? (ichorCNA)
 *
 *
 * To be implemented
 * - Arm-level. Yes, preprocess cytoband, make an arm file.
 * - Make per-arm truth.
 * - The output should be per-arm calls.
 * - Make a heatmap.
 *
 * ARM mode is better for comparison e.g. ULP-WGS vs TP
 *
 * Output:
 * 1. Arm level: (eval vs truth calls, w/ tp etc)
 * 2. Arm level: (truth)
 * 3. Arm Level: (eval)
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

    public static final String FILE_TYPE_NAME = "file-type";
    @Argument(doc = "The CNV caller used to generate the seg file",
            fullName = FILE_TYPE_NAME)
    private String fileType;

    // Make it an argument in case different reference is used
    public static final String CHR_ARMS_NAME = "chr-arms";
    @Argument(doc = "A tsv file containing chm arms",
            fullName = CHR_ARMS_NAME, optional = true)
    private String chrArmFile = "/Users/tsato/workspace/cfCNV/resources/chr_arms_hg19.tsv";

    // Do this later
    // ichor takes in the *seg.txt file.
    enum FileType {
        MODEL_SEGMENTS, ICHOR_CNA
    }

    // Annotation keys
    final AnnotationKey<String> chrArmKey = new AnnotationKey<>("chr_arm", String.class, s -> s.length() > 0);

    // chr arm only
    final AnnotationKey<String> truthCallKey = new AnnotationKey<>("truth_call", String.class, s -> s.length() > 0);
    final AnnotationKey<String> evalCallKey = new AnnotationKey<>("eval_call", String.class, s -> s.length() > 0);
    final AnnotationKey<String> truthFractionOverlapKey = new AnnotationKey<>("truth_fraction_overlap", String.class, s -> s.length() > 0);
    final AnnotationKey<String> evalFractionOverlapKey = new AnnotationKey<>("eval_fraction_overlap", String.class, s -> s.length() > 0);



    final AnnotationKey<Integer> numEvalSegmentsKey = new AnnotationKey<>("num_eval_segments", Integer.class, n -> n >= 0);
    final AnnotationKey<Integer> numOverlappingEvalBasesKey = new AnnotationKey<>("num_overlapping_eval_bases", Integer.class, n -> n >= 0);
    final AnnotationKey<Integer> eventWidthKey = new AnnotationKey<>("event_width", Integer.class, n -> n > 0);
    final AnnotationKey<Double> fractionOverlapKey = new AnnotationKey<>("fraction_overlap", Double.class, p -> p >= 0);

    final AnnotationKey<String> statusKey = new AnnotationKey<>("status", String.class, s -> s.length() > 0);


    // Need to require reference for ichorCNA, since its output does not have the sequencing dictionary header
    // Or should we loosen this restriction?
    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void traverse() {
        CalledCopyRatioSegmentCollection evalSegments = null;
        CalledCopyRatioSegmentCollection truthSegments = null;

        // Read chrarm files (what data structure?)
        final List<AnnotatedInterval> chrArms = readChmArms(chrArmFile);

        if (fileType.equalsIgnoreCase("MODELSEGMENTS")) {
            evalSegments = new CalledCopyRatioSegmentCollection(evalFile);
            truthSegments = new CalledCopyRatioSegmentCollection(truthFile);
        } else if (fileType.equalsIgnoreCase("ICHORCNA")) {
            evalSegments = parseICHORCalls(evalFile);
            truthSegments = parseICHORCalls(truthFile);
        }

        // BEGIN CHR ARM MODE
        final Map<AnnotatedInterval, List<CalledCopyRatioSegment>> perChrArmTruthOverlapMap =
                IntervalUtils.createOverlapMap(chrArms, truthSegments.getRecords(),
                        getBestAvailableSequenceDictionary());
        final List<AnnotatedInterval> perChrArmTruthCalls = getArmLevelCalls(perChrArmTruthOverlapMap);
        // TODO: Output per-chromosome truth file

        // Get per CHR ARM eval calls
        final Map<AnnotatedInterval, List<CalledCopyRatioSegment>> perChrArmEvalOverlapMap =
                IntervalUtils.createOverlapMap(chrArms, evalSegments.getRecords(),
                        getBestAvailableSequenceDictionary());
        final List<AnnotatedInterval> perChrArmEvalCalls = getArmLevelCalls(perChrArmEvalOverlapMap);
        final AnnotatedIntervalComparator comparator = new AnnotatedIntervalComparator(new SAMFileHeader(getBestAvailableSequenceDictionary()));

        Collections.sort(chrArms, comparator);
        Collections.sort(perChrArmEvalCalls, comparator);
        Collections.sort(perChrArmTruthCalls, comparator);

        // Combine truth and eval
        int d = 3;
        // Can we have the list be sorted?
        // Perhaps better to use tree set...https://stackoverflow.com/questions/8725387/why-is-there-no-sortedlist-in-java

        final List<AnnotatedInterval> result2 = new ArrayList<>();
        int i = 0;
        for (final AnnotatedInterval chrArm : chrArms){
            final List<Pair<AnnotationKey<?>, Object>> annotationValues = new ArrayList<>();
            final AnnotatedInterval truth = perChrArmTruthCalls.get(i); // Be more describe with var name
            final AnnotatedInterval eval = perChrArmEvalCalls.get(i); // Be more describe with var name

            // Best to check here that the things match up etc. But ok.

            final String truthCall = truth.getAnnotationMap().getValue(statusKey);
            final String evalCall = eval.getAnnotationMap().getValue(statusKey);

            // Should be a way to "inherit" annotation values of an exising annotation map but w.e.
            annotationValues.add(new ImmutablePair<>(chrArmKey, chrArm.getAnnotationMap().getValue(chrArmKey)));
            annotationValues.add(new ImmutablePair<>(truthCallKey, truthCall));
            annotationValues.add(new ImmutablePair<>(evalCallKey, evalCall));
            annotationValues.add(new ImmutablePair<>(truthFractionOverlapKey, chrArm.getAnnotationMap().getValue(fractionOverlapKey)));
            annotationValues.add(new ImmutablePair<>(evalFractionOverlapKey, chrArm.getAnnotationMap().getValue(fractionOverlapKey)));

            final String callStatus;
            final boolean truthHasEvent = truthCall.equals(CalledCopyRatioSegment.Call.AMPLIFICATION.toString()) ||
                    truthCall.equals(CalledCopyRatioSegment.Call.DELETION.toString());
            if (truthCall.equals(evalCall)){
                if (truthHasEvent){
                    callStatus = "TP";
                } else {
                    callStatus = "TN";
                }
            } else {
                if (truthHasEvent){
                    callStatus = "FN";
                } else {
                    callStatus = "FP";
                }
            }

            // This is weird, two different things are both called status, fix it.
            annotationValues.add(new ImmutablePair<>(statusKey, callStatus));

            // annotationValues.add(new ImmutablePair<>(statusKey, eval.getAnnotationMap().getValue(statusKey)));
            // END MORE ANNOTATIONS

            final AnnotatedInterval annotatedInterval = new AnnotatedInterval(chrArm.getInterval(),
                    new AnnotationMap(annotationValues));
            result2.add(annotatedInterval);
            i++;
        }

        final String dir = "/Users/tsato/workspace/cfCNV/analysis/ichor_output";
        final AnnotatedIntervalCollection aic2 = new AnnotatedIntervalCollection(truthSegments.getMetadata(), result2);
        aic2.write(new File(dir, "per_chr_test.tsv"));

        final AnnotatedIntervalCollection aic3 = new AnnotatedIntervalCollection(truthSegments.getMetadata(), perChrArmTruthCalls);
        aic3.write(new File(dir, "per_chr_truth_test.tsv"));

        final AnnotatedIntervalCollection aic4 = new AnnotatedIntervalCollection(truthSegments.getMetadata(), perChrArmEvalCalls);
        aic4.write(new File(dir, "per_chr_eval_test.tsv"));
        // END CHR ARM MODE

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

        // Traverse the truth calls. Use AnnotatedIntervalCollection to output the result (for now)
        final List<AnnotatedInterval> result = getCNVConcordance(perTruthOverlapMap);
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

    private List<AnnotatedInterval> getArmLevelCalls(final Map<AnnotatedInterval, List<CalledCopyRatioSegment>> perArmOverlapMap) {
        final List<AnnotatedInterval> result = new ArrayList<>();
        for (final Map.Entry<AnnotatedInterval, List<CalledCopyRatioSegment>> entry : perArmOverlapMap.entrySet()) {
            final AnnotatedInterval chrArm = entry.getKey();
            final List<CalledCopyRatioSegment> overlappingSegments = entry.getValue();

            int[] segmentSizeByEvent = new int[3];
            int NEUTRAL_INDEX = 0;
            int AMPLIFICATION_INDEX = 1;
            int DELETION_INDEX = 2;
            String[] eventNames = new String[]{ "NEUTRAL", "AMPLIFICATION", "DELETION"};
            for (final CalledCopyRatioSegment overlappingTruthSegment : overlappingSegments){
                final CalledCopyRatioSegment.Call call = overlappingTruthSegment.getCall();
                // Must check that the call is the same between tr and ev
                final int overlappingSegmentSize = CoordMath.getOverlap(chrArm.getStart(), chrArm.getEnd(),
                        overlappingTruthSegment.getStart(), overlappingTruthSegment.getEnd());
                if (call == CalledCopyRatioSegment.Call.NEUTRAL){
                    segmentSizeByEvent[NEUTRAL_INDEX] += overlappingSegmentSize;
                } else if (call == CalledCopyRatioSegment.Call.AMPLIFICATION){
                    segmentSizeByEvent[AMPLIFICATION_INDEX] += overlappingSegmentSize;
                } else {
                    segmentSizeByEvent[DELETION_INDEX] += overlappingSegmentSize;
                }
            }

            final int eventIndex = MathUtils.maxElementIndex(segmentSizeByEvent);
            final String eventName = eventNames[eventIndex];
            final int eventSize = segmentSizeByEvent[eventIndex];

            // +1 because the interval includes both ends. e.g. [10, 20] contains 11 bases
            final int chmArmSize = chrArm.getEnd() - chrArm.getStart() + 1;
            final double fractionOverlap = eventSize/ Math.max(chmArmSize, 1.0);

            final List<Pair<AnnotationKey<?>, Object>> annotationValues = new ArrayList<>();
            annotationValues.add(new ImmutablePair<>(chrArmKey, chrArm.getAnnotationMap().getValue(chrArmKey)));
            annotationValues.add(new ImmutablePair<>(numOverlappingEvalBasesKey, eventSize)); // TODO: rename numOver...to EventSize
            annotationValues.add(new ImmutablePair<>(eventWidthKey, chmArmSize));
            annotationValues.add(new ImmutablePair<>(fractionOverlapKey, fractionOverlap));
            annotationValues.add(new ImmutablePair<>(statusKey, eventName));

            final AnnotatedInterval annotatedInterval = new AnnotatedInterval(chrArm.getInterval(),
                    new AnnotationMap(annotationValues));
            result.add(annotatedInterval);
        }

        return result;
    }

    private List<AnnotatedInterval> getCNVConcordance(final Map<CalledCopyRatioSegment, List<CalledCopyRatioSegment>> perTruthOverlapMap) {
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

        return result;
    }

    private List<AnnotatedInterval> readChmArms(String chrArmFile) {
        final List<AnnotatedInterval> result = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(chrArmFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                final String[] row = line.split("\t");
                if (row[0].equalsIgnoreCase("chr")){
                    continue;
                }

                final String chrom = row[0];
                final String arm = row[1];
                final int start = Math.max(Integer.parseInt(row[2]), 1);
                final int end = Integer.parseInt(row[3]);

                final List<Pair<AnnotationKey<?>, Object>> annotationValues = new ArrayList<>();
                annotationValues.add(new ImmutablePair<>(chrArmKey, arm));

                final AnnotatedInterval interval = new AnnotatedInterval(new SimpleInterval(chrom, start, end), new AnnotationMap(annotationValues));
                result.add(interval);
            }
        } catch (IOException e){
            throw new UserException("Hello");
        }

        return result;
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

    // Eventually use TableReader.... //
    private CalledCopyRatioSegmentCollection parseICHORCalls(final File file){
        final List<CalledCopyRatioSegment> rows = new ArrayList<>();
        String sampleName = null;
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                final String[] row = line.split("\t");
                if (row[0].equalsIgnoreCase("ID")){
                    continue;
                }

                if (sampleName == null){
                    sampleName = row[0];
                }
                final String chrom = row[1];
                final int start = Integer.parseInt(row[2]);
                final int end = Integer.parseInt(row[3]);
                final int numMarkers = Integer.parseInt(row[4]);
                final double segmentMedianLogRatio = Double.parseDouble(row[5]);
                final int copyNumber = Integer.parseInt(row[6]);
                final String ichorCall = row[7];
                final String subcloneStatus = row[8];
                // What's logR copy number?
                // Corrected call?
                final SimpleInterval interval = new SimpleInterval(chrom, start, end);
                final CopyRatioSegment segment = new CopyRatioSegment(interval, numMarkers, segmentMedianLogRatio);
                // This is where we lose some resolution...which should be ok
                // But if not return
                final CalledCopyRatioSegment.Call msCall = ichor2MSCallConversion(ichorCall);
                final CalledCopyRatioSegment calledSegment = new CalledCopyRatioSegment(segment, msCall);
                rows.add(calledSegment);
            }
        } catch (IOException e){
            throw new UserException("Hello");
        }

        final SimpleSampleLocatableMetadata metaData = new SimpleSampleLocatableMetadata(sampleName, getBestAvailableSequenceDictionary());
        return new CalledCopyRatioSegmentCollection(metaData, rows);
    }

    private CalledCopyRatioSegment.Call ichor2MSCallConversion(final String ichorCNACall){
        if (ichorCNACall.equals("NEUT")){
            return CalledCopyRatioSegment.Call.NEUTRAL;
        } else if (ichorCNACall.equals("HETD")){
            return CalledCopyRatioSegment.Call.DELETION;
        } else {
            return CalledCopyRatioSegment.Call.AMPLIFICATION;
        }

    }
}
