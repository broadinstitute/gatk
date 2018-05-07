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
import java.util.stream.Stream;

final class ReadDepthCluster implements Serializable {

    public static final long serialVersionUID = 1L;
    private final List<ReadDepthEvent> eventsList;
    private final SVIntervalTree<ReadDepthEvent> eventsTree;
    private final SimpleSVType.TYPES type;

    public ReadDepthCluster(final List<ReadDepthEvent> events, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        this.eventsList = events;
        this.eventsTree = getIntervalTree(eventsList);
        setCopyNumberInfo(copyRatioSegmentOverlapDetector, dictionary);
        final Set<CalledCopyRatioSegment> segments = events.stream()
                .flatMap(event -> copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(event.getEvent().getInterval(), dictionary)).stream())
                .collect(Collectors.toSet());
        setOverlapInfo(eventsTree, segments, dictionary);
        setNearestCallDistances(copyRatioSegmentOverlapDetector, dictionary);
        this.type = events.isEmpty() ? null : events.iterator().next().getEvent().getEventType();
    }

    private static void optimizeOverlapInfo(final List<ReadDepthEvent> events) {
        for (final ReadDepthEvent event : events) {
            if (!event.overlapInfoList.isEmpty()) {
                final Map<Integer, Double> optimizedOverlapInfo = new HashMap<>(SVUtils.hashMapCapacity(event.overlapInfoList.size()));
                double newCopyNumber = 0;
                double totalWeight = 0;
                for (final ReadDepthModel.OverlapInfo info : event.overlapInfoList) {
                    newCopyNumber += info.segmentCopyNumber * info.weight;
                    totalWeight += info.weight;
                    for (final Tuple2<Integer, Integer> idAndSign : info.overlappingIdsAndSigns) {
                        final int id = idAndSign._1;
                        final int sign = idAndSign._2;
                        optimizedOverlapInfo.putIfAbsent(id, 0.);
                        optimizedOverlapInfo.put(id, optimizedOverlapInfo.get(id) + info.weight * sign);
                    }
                }
                if (totalWeight == 0) {
                    throw new IllegalStateException("Total weight in overlap info was 0");
                }
                newCopyNumber /= totalWeight;
                event.optimizedOverlapInfoList = optimizedOverlapInfo.entrySet().stream().map(entry -> new Tuple2<>(entry.getKey(), entry.getValue())).collect(Collectors.toList());
                event.observedCopyNumber = newCopyNumber;
            } else {
                event.optimizedOverlapInfoList = Collections.emptyList();
                event.observedCopyNumber = 2;
            }
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

    private static final class OverlappingEventBoundary {

        public int pos;
        public int index;

        public OverlappingEventBoundary(final int pos, final int index) {
            this.pos = pos;
            this.index = index;
        }
    }

    private static void setOverlapInfoOnSegment(final CalledCopyRatioSegment segment, final List<ReadDepthEvent> overlappers) {

        if (overlappers.isEmpty()) {
            return;
        }

        final double copyNumber = Math.pow(2.0, segment.getMeanLog2CopyRatio()) * 2;
        final Map<Integer, ReadDepthEvent> overlapperMap = overlappers.stream().collect(Collectors.toMap(event -> event.getId(), event -> event));
        final SimpleInterval interval = segment.getInterval();
        final List<OverlappingEventBoundary> sortedOverlappers = overlappers.stream()
                .flatMap(event -> Stream.of(new OverlappingEventBoundary(event.getEvent().getStart(), event.getId()), new OverlappingEventBoundary(event.getEvent().getEnd(), event.getId())))
                .sorted(Comparator.comparingInt(boundary -> boundary.pos))
                .collect(Collectors.toList());

        final Set<Integer> openEvents = new HashSet<>(SVUtils.hashMapCapacity(overlappers.size()));
        int overlapperIndex = 0;
        while (overlapperIndex < sortedOverlappers.size() && sortedOverlappers.get(overlapperIndex).pos < interval.getStart()) {
            openEvents.add(sortedOverlappers.get(overlapperIndex).index);
            overlapperIndex++;
        }
        int lastPos = interval.getStart();
        boolean reachedEndOfInterval = false;
        while (!reachedEndOfInterval && overlapperIndex < sortedOverlappers.size()) {
            int currentPos = sortedOverlappers.get(overlapperIndex).pos;
            if (currentPos >= interval.getEnd()) {
                currentPos = interval.getEnd();
                reachedEndOfInterval = true;
            }
            final int currentEventId = sortedOverlappers.get(overlapperIndex).index;
            int length = currentPos - lastPos;
            for (final Integer openEventId : openEvents) {
                final double weight = length / (double) overlapperMap.get(openEventId).getEvent().getSize();
                final List<Tuple2<Integer,Integer>> overlappingIdsAndSigns = openEvents.stream()
                        .filter(eventId -> !eventId.equals(openEventId))
                        .map(eventId -> new Tuple2<>(openEventId, overlapperMap.get(openEventId).getEvent().getEventType() == SimpleSVType.TYPES.DEL ? 1 : -1))
                        .collect(Collectors.toList());
                overlapperMap.get(openEventId).overlapInfoList.add(new ReadDepthModel.OverlapInfo(overlappingIdsAndSigns, weight, copyNumber));
            }
            if (openEvents.contains(currentEventId)) {
                openEvents.remove(currentEventId);
            } else {
                openEvents.add(currentEventId);
            }
            lastPos = currentPos;
            overlapperIndex++;
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
            final CalledCopyRatioSegment.Call expectedCall = event.getEvent().getEventType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION;
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
            final CalledCopyRatioSegment.Call expectedCall = event.getEvent().getEventType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION;
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
