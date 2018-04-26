package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

final class ReadDepthCluster implements Serializable {

    public static final long serialVersionUID = 1L;
    private final List<ReadDepthEvent> eventsList;
    private final SVIntervalTree<ReadDepthEvent> eventsTree;
    private final SimpleSVType.TYPES type;

    public ReadDepthCluster(final List<ReadDepthEvent> events, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        this.eventsList = events;
        this.eventsTree = getIntervalTree(eventsList);
        setCopyNumberInfo(copyRatioSegmentOverlapDetector, dictionary);
        final List<CalledCopyRatioSegment> segments = events.stream()
                .flatMap(event -> copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(event.getEvent().getInterval(), dictionary)).stream())
                .collect(Collectors.toList());
        setOverlapInfo(eventsTree, segments, dictionary);
        setNearestCallDistances(copyRatioSegmentOverlapDetector, dictionary);
        this.type = events.isEmpty() ? null : events.iterator().next().getEvent().getType();
    }

    private static void optimizeOverlapInfo(final List<ReadDepthEvent> events) {
        for (final ReadDepthEvent event : events) {
            final Map<Integer, Double> optimizedOverlapInfo = new HashMap<>(SVUtils.hashMapCapacity(event.overlapInfoList.size()));
            double newCopyNumber = 0;
            double totalWeight = 0;
            for (final ReadDepthModel.OverlapInfo info : event.overlapInfoList) {
                newCopyNumber += info.segmentCopyNumber * info.weight;
                totalWeight += info.weight;
                for (final Tuple2<Integer,Integer> idAndSign : info.overlappingIdsAndSigns) {
                    final int id = idAndSign._1;
                    final int sign = idAndSign._2;
                    optimizedOverlapInfo.putIfAbsent(id, 0.);
                    optimizedOverlapInfo.put(id, optimizedOverlapInfo.get(id) + info.weight * sign);
                }
            }
            newCopyNumber /= totalWeight;
            event.optimizedOverlapInfoList = optimizedOverlapInfo.entrySet().stream().map(entry -> new Tuple2<>(entry.getKey(),entry.getValue())).collect(Collectors.toList());
            event.observedCopyNumber = newCopyNumber;
        }
    }

    private static void setOverlapInfo(final SVIntervalTree<ReadDepthEvent> eventsTree, final Collection<CalledCopyRatioSegment> segments, final SAMSequenceDictionary dictionary) {
        final Iterator<SVIntervalTree.Entry<ReadDepthEvent>> iterator = eventsTree.iterator();
        while (iterator.hasNext()) {
            iterator.next().getValue().overlapInfoList = new ArrayList<>();
        }
        for (final CalledCopyRatioSegment segment : segments) {
            final SVInterval segmentInterval = SVIntervalUtils.convertToSVInterval(segment.getInterval(), dictionary);
            final List<ReadDepthEvent> overlappingEvents = Utils.stream(eventsTree.overlappers(segmentInterval)).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList());
            setOverlapInfoOnSegment(segment, overlappingEvents);
        }
        optimizeOverlapInfo(Utils.stream(eventsTree.iterator()).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList()));
    }

    private static void setOverlapInfoOnSegment(final CalledCopyRatioSegment segment, final List<ReadDepthEvent> overlappers) {

        final double copyNumber = Math.pow(2.0, segment.getMeanLog2CopyRatio()) * 2;
        final Map<Integer, ReadDepthEvent> overlapperMap = overlappers.stream().collect(Collectors.toMap(event -> event.getId(), event -> event));
        final SimpleInterval interval = segment.getInterval();
        final List<ReadDepthEvent> startSortedOverlappers = new ArrayList<>(overlappers);
        Collections.sort(startSortedOverlappers, Comparator.comparingInt(event -> event.getEvent().getStart()));

        final List<ReadDepthEvent> endSortedOverlappers = new ArrayList<>(overlappers);
        Collections.sort(endSortedOverlappers, Comparator.comparingInt(event -> event.getEvent().getEnd()));

        int lastStart = interval.getStart();
        int startIndex = 0;
        int endIndex = 0;
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
            for (final Integer id : openEvents) {
                final double weight = subIntervalLength / (double) overlapperMap.get(id).getEvent().getSize();
                final List<Tuple2<Integer,Integer>> overlappingIdsAndSigns = openEvents.stream()
                        .filter(eventId -> !eventId.equals(id))
                        .map(eventId -> new Tuple2<>(id, overlapperMap.get(id).getEvent().getType() == SimpleSVType.TYPES.DEL ? -1 : 1))
                        .collect(Collectors.toList());
                overlapperMap.get(id).overlapInfoList.add(new ReadDepthModel.OverlapInfo(overlappingIdsAndSigns, weight, copyNumber));
            }
            if (lastStart == interval.getEnd()) {
                break;
            }
        }
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

    private void setCopyNumberInfo(final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        for (final ReadDepthEvent event : eventsList) {
            double observedCopyNumber = 0;
            final SVInterval eventInterval = event.getEvent().getInterval();
            final Set<CalledCopyRatioSegment> overlappingSegments = copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(eventInterval, dictionary));
            int copyNumberOverlap = 0;
            final CalledCopyRatioSegment.Call expectedCall = event.getEvent().getType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION;
            for (final CalledCopyRatioSegment copyRatioSegment : overlappingSegments) {
                final SVInterval segmentInterval = SVIntervalUtils.convertToSVInterval(copyRatioSegment.getInterval(), dictionary);
                final double copyNumberCall = Math.pow(2.0, copyRatioSegment.getMeanLog2CopyRatio()) * 2;
                observedCopyNumber += copyNumberCall * segmentInterval.overlapLen(eventInterval);
                if (copyRatioSegment.getCall() == expectedCall) {
                    copyNumberOverlap += segmentInterval.overlapLen(eventInterval);
                }
            }
            event.observedCopyNumber = observedCopyNumber / (double) eventInterval.getLength();
            event.copyNumberCallOverlap = copyNumberOverlap / (double) eventInterval.getLength();
        }
    }

    private void setNearestCallDistances(final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        for (final ReadDepthEvent event : eventsList) {
            final SimpleInterval interval = SVIntervalUtils.convertToSimpleInterval(event.getEvent().getInterval(), dictionary);
            final CalledCopyRatioSegment.Call expectedCall = event.getEvent().getType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION;
            final List<CalledCopyRatioSegment> overlappingCalls = copyRatioSegmentOverlapDetector.getOverlaps(interval).stream().filter(segment -> segment.getCall() == expectedCall).collect(Collectors.toList());
            int leftDistance;
            int rightDistance;
            if (overlappingCalls.isEmpty()) {
                event.leftDistance = 100;
                event.rightDistance = 100;
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
                event.leftDistance = Math.log(1 + Math.abs(leftDistance));
                event.rightDistance = Math.log(1 + Math.abs(rightDistance));
            }
        }
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

    public SimpleSVType.TYPES getType() {
        return type;
    }
}
