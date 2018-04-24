package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

final class ReadDepthCluster {

    private final List<ReadDepthEvent> eventsList;
    private final SVIntervalTree<ReadDepthEvent> eventsTree;
    private final List<Tuple2<List<ReadDepthModel.OverlapInfo>, Double>> copyNumberInfo;
    private final List<Tuple2<Integer, Integer>> nearestCallDistances;
    private final SimpleSVType.TYPES type;

    public ReadDepthCluster(final List<ReadDepthEvent> events, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        this.eventsList = events;
        this.eventsTree = getIntervalTree(eventsList);
        final List<CalledCopyRatioSegment> overlappingSegments = getOverlappingSegments(eventsList, copyRatioSegmentOverlapDetector, dictionary);
        this.copyNumberInfo = getCopyNumberInfo(overlappingSegments, eventsTree, dictionary);
        this.nearestCallDistances = getNearestCallDistances(eventsList, copyRatioSegmentOverlapDetector, dictionary);
        this.type = events.isEmpty() ? null : events.iterator().next().getEvent().getType();
    }

    private static List<ReadDepthModel.OverlapInfo> getOverlapInfo(final CalledCopyRatioSegment segment, final List<ReadDepthEvent> overlappers) {

        final Map<Integer, LargeSimpleSV> overlapperMap = overlappers.stream().collect(Collectors.toMap(event -> event.getId(), event -> event.getEvent()));
        final SimpleInterval interval = segment.getInterval();
        final List<ReadDepthEvent> startSortedOverlappers = new ArrayList<>(overlappers);
        Collections.sort(startSortedOverlappers, Comparator.comparingInt(event -> event.getEvent().getStart()));

        final List<ReadDepthEvent> endSortedOverlappers = new ArrayList<>(overlappers);
        Collections.sort(endSortedOverlappers, Comparator.comparingInt(event -> event.getEvent().getEnd()));

        int lastStart = interval.getStart();
        int startIndex = 0;
        int endIndex = 0;
        final ArrayList<ReadDepthModel.OverlapInfo> overlapInfoList = new ArrayList<>();
        final List<Integer> openEvents = new ArrayList<>();
        while (startIndex < overlappers.size() || endIndex < overlappers.size()) {
            int nextStart;
            if (startIndex < overlappers.size()) {
                nextStart = startSortedOverlappers.get(startIndex).getEvent().getStart();
                if (nextStart <= interval.getStart()) {
                    openEvents.add(startSortedOverlappers.get(startIndex).getId());
                    startIndex++;
                    continue;
                }
            } else {
                nextStart = interval.getEnd();
            }
            int nextEnd;
            if (endIndex < overlappers.size()) {
                nextEnd = endSortedOverlappers.get(endIndex).getEvent().getEnd();
            } else {
                nextEnd = interval.getEnd();
            }
            final List<Tuple2<Integer, Integer>> idsAndCoefficients = openEvents.stream()
                    .map(id -> new Tuple2<>(id, overlapperMap.get(id).getType() == SimpleSVType.TYPES.DEL ? -1 : 1))
                    .collect(Collectors.toList());
            final int subIntervalLength;
            if (startIndex < overlappers.size() && nextStart < nextEnd) {
                subIntervalLength = nextStart - lastStart;
                lastStart = nextStart;
                openEvents.add(startSortedOverlappers.get(startIndex).getId());
                startIndex++;
            } else if (nextEnd >= interval.getEnd()) {
                subIntervalLength = interval.getEnd() - lastStart;
                lastStart = interval.getEnd();
            } else {
                subIntervalLength = nextEnd - lastStart;
                lastStart = nextEnd;
                openEvents.remove(Integer.valueOf(endSortedOverlappers.get(endIndex).getId()));
                endIndex++;
            }
            if (!idsAndCoefficients.isEmpty()) {
                final double scalingFactor = subIntervalLength;
                overlapInfoList.add(new ReadDepthModel.OverlapInfo(idsAndCoefficients, scalingFactor));
            }
            if (lastStart == interval.getEnd()) {
                break;
            }
        }
        overlapInfoList.trimToSize();
        return overlapInfoList;
    }

    private List<CalledCopyRatioSegment> getOverlappingSegments(final List<ReadDepthEvent> eventsList, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        return eventsList.stream()
                .flatMap(event -> copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(event.getEvent().getInterval(), dictionary)).stream())
                .distinct()
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .collect(Collectors.toList());
    }

    private SVIntervalTree<ReadDepthEvent> getIntervalTree(final List<ReadDepthEvent> eventsList) {
        final SVIntervalTree<ReadDepthEvent> eventsTree = new SVIntervalTree<>();
        for (final ReadDepthEvent call : eventsList) {
            eventsTree.put(call.getEvent().getInterval(), call);
        }
        return eventsTree;
    }

    private List<Tuple2<List<ReadDepthModel.OverlapInfo>, Double>> getCopyNumberInfo(final List<CalledCopyRatioSegment> overlappingSegments, final SVIntervalTree<ReadDepthEvent> eventsTree, final SAMSequenceDictionary dictionary) {
        final ArrayList<Tuple2<List<ReadDepthModel.OverlapInfo>, Double>> result = new ArrayList<>(overlappingSegments.size());
        for (final CalledCopyRatioSegment copyRatioSegment : overlappingSegments) {
            final SVInterval interval = SVIntervalUtils.convertToSVInterval(copyRatioSegment.getInterval(), dictionary);
            final Iterator<SVIntervalTree.Entry<ReadDepthEvent>> iter = eventsTree.overlappers(interval);
            if (iter.hasNext()) {
                final double copyNumberCall = Math.min(Math.pow(2.0, copyRatioSegment.getMeanLog2CopyRatio()) * 2, 4.0); //TODO capping segment copy number at 4
                final List<ReadDepthModel.OverlapInfo> overlapInfo = getOverlapInfo(copyRatioSegment, Utils.stream(iter).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList()));
                result.add(new Tuple2<>(overlapInfo, copyNumberCall));
            }
        }
        result.trimToSize();
        return result;
    }

    private List<Tuple2<Integer, Integer>> getNearestCallDistances(final List<ReadDepthEvent> eventsList, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        final List<Tuple2<Integer, Integer>> result = new ArrayList<>(eventsList.size());
        for (final ReadDepthEvent event : eventsList) {
            final SimpleInterval interval = SVIntervalUtils.convertToSimpleInterval(event.getEvent().getInterval(), dictionary);
            final CalledCopyRatioSegment.Call expectedCall = event.getEvent().getType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION;
            final List<CalledCopyRatioSegment> overlappingCalls = copyRatioSegmentOverlapDetector.getOverlaps(interval).stream().filter(segment -> segment.getCall() == expectedCall).collect(Collectors.toList());
            int leftDistance;
            int rightDistance;
            if (overlappingCalls.isEmpty()) {
                result.add(new Tuple2<>((int) -1e6, (int) 1e6));
            } else {
                leftDistance = overlappingCalls.get(0).getStart() - interval.getStart();
                rightDistance = overlappingCalls.get(0).getEnd() - interval.getEnd();
                for (int i = 1; i < overlappingCalls.size(); i++) {
                    final int leftDistanceTest = overlappingCalls.get(i).getStart() - interval.getStart();
                    final int rightDistanceTest = overlappingCalls.get(i).getEnd() - interval.getEnd();
                    if (Math.abs(leftDistanceTest) < Math.abs(leftDistance)) {
                        leftDistance = leftDistanceTest;
                    }
                    if (Math.abs(rightDistanceTest) < Math.abs(rightDistance)) {
                        rightDistance = rightDistanceTest;
                    }
                }
                event.leftDistance = leftDistance;
                event.rightDistance = rightDistance;
                result.add(new Tuple2<>(leftDistance, rightDistance));
            }
        }
        return result;
    }

    public double[] getStates() {
        return eventsList.stream().mapToDouble(event -> event.getState()).toArray();
    }

    public List<ReadDepthEvent> getEventsList() {
        return eventsList;
    }

    public SVIntervalTree<ReadDepthEvent> getEventsTree() {
        return eventsTree;
    }

    public List<Tuple2<List<ReadDepthModel.OverlapInfo>, Double>> getCopyNumberInfo() {
        return copyNumberInfo;
    }

    public List<Tuple2<Integer, Integer>> getNearestCallDistances() {
        return nearestCallDistances;
    }

    public SimpleSVType.TYPES getType() {
        return type;
    }
}
