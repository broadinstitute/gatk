package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import scala.Tuple2;

import java.util.*;

class BarcodeSetByIntervalIterator implements Iterator<Tuple2<SVInterval, Set<String>>> {

    private final List<Tuple2<Integer, List<String>>> starts;
    private final List<Tuple2<Integer, List<String>>> stops;
    private final int contig;
    private int startIdx;
    private int stopIdx;
    private int prevBoundary;
    private final Set<String> currentBarcodes;

    public BarcodeSetByIntervalIterator(final int contig, final SVIntervalTree<List<String>> intervals) {
        this.contig = contig;
        starts = new ArrayList<>(intervals.size());
        intervals.iterator().forEachRemaining(e -> starts.add(new Tuple2<>(e.getInterval().getStart(), e.getValue())));
        stops = new ArrayList<>(intervals.size());
        intervals.iterator().forEachRemaining(e -> stops.add(new Tuple2<>(e.getInterval().getEnd(), e.getValue())));
        stops.sort(Comparator.comparingInt(Tuple2::_1));

        startIdx = 0;
        stopIdx = 0;
        prevBoundary = 0;

        currentBarcodes = new HashSet<>();
    }

    @Override
    public boolean hasNext() {
        return stopIdx < stops.size();
    }

    @Override
    public Tuple2<SVInterval, Set<String>> next() {
        int nextStart = startIdx < starts.size() ? starts.get(startIdx)._1 : -1;
        int nextStop = stops.get(stopIdx)._1;

        if (currentBarcodes.isEmpty()) {
            currentBarcodes.addAll(starts.get(startIdx)._2);
            prevBoundary = starts.get(startIdx)._1;
            startIdx++;


            while (startIdx < starts.size() && starts.get(startIdx)._1() == prevBoundary) {
                currentBarcodes.addAll(starts.get(startIdx)._2);
                startIdx++;
            }

            nextStart = startIdx < starts.size() ? starts.get(startIdx)._1 : -1;

        }

        final Tuple2<SVInterval, Set<String>> result;
        if (nextStart != -1 && nextStart <= nextStop) {
            // process a start
            //System.out.println("Process a start at " + nextStart + ": " + starts.get(startIdx)._2);

            final SVInterval newInterval = new SVInterval(contig, prevBoundary, nextStart);
            final Set<String> barcodes = new HashSet<>(currentBarcodes);
            result = new Tuple2<>(newInterval, barcodes);

            currentBarcodes.addAll(starts.get(startIdx)._2);
            prevBoundary = starts.get(startIdx)._1;
            startIdx++;

            while (startIdx < starts.size() && starts.get(startIdx)._1() == prevBoundary) {
                currentBarcodes.addAll(starts.get(startIdx)._2);
                startIdx++;
            }

        } else {
            // process a stop
            //System.out.println("Process a stop at " + nextStop + ": " + stops.get(stopIdx)._2);
            final SVInterval newInterval = new SVInterval(contig, prevBoundary, nextStop);
            final Set<String> barcodes = new HashSet<>(currentBarcodes);
            result = new Tuple2<>(newInterval, barcodes);

            currentBarcodes.removeAll(stops.get(stopIdx)._2);
            prevBoundary = stops.get(stopIdx)._1;
            stopIdx++;

            while(stopIdx < stops.size() && stops.get(stopIdx)._1() == prevBoundary) {
                currentBarcodes.removeAll(stops.get(stopIdx)._2);
                stopIdx++;
            }
        }
        return result;
    }
}
