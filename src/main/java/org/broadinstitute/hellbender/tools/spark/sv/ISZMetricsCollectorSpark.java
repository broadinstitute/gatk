package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Worker class to collect insert size metrics.
 */
public class ISZMetricsCollectorSpark implements Serializable{
    private static final long serialVersionUID = 1L;

    // uses TreeMap because it is serializable and sorted
    private TreeMap<Integer, Long> hist;

    public ISZMetricsCollectorSpark(final JavaRDD<GATKRead> filteredReads,
                                    final Set<MetricAccumulationLevel> accumLevels,
                                    final List<SAMReadGroupRecord> samRgRecords,
                                    final double devTolerance) {

        final Map<Integer, Long> toBeCasted = filteredReads.map(read -> Math.abs(read.getFragmentLength()) ).countByValue();
        hist = new TreeMap<>(toBeCasted);
    }

    public void produceMetricsFile(MetricsFile<InsertSizeMetrics, ?> metricsFile) {

        InsertSizeMetrics metric = new InsertSizeMetrics();
        metric.PAIR_ORIENTATION = SamPairUtil.PairOrientation.FR; // what does this mean?

        metric.MIN_INSERT_SIZE = hist.firstKey();
        metric.MAX_INSERT_SIZE = hist.lastKey();

        collectMedianAndMAD(metric);

        collectPercentages(metric);

        metricsFile.addMetric(metric);
    }

    private void collectMedianAndMAD(InsertSizeMetrics metric){

        // first round get total # of valid read pairs and mean size
        Long pairsCount = 0L;
        Double totalLength = 0.0;
        Double totalSquaredLenth = 0.0;
        for(Map.Entry<Integer, Long> pair : hist.entrySet()) {
            final Integer size = pair.getKey();
            final Long cnt = pair.getValue();
            pairsCount += cnt;
            totalLength += ((double) size) * ((double) cnt);
            totalSquaredLenth += ((double) size) * ((double) size) * ((double) cnt);
        }
        metric.READ_PAIRS = pairsCount;
        final double n = (double) pairsCount;
        metric.MEAN_INSERT_SIZE = totalLength/n;
        metric.STANDARD_DEVIATION = Math.sqrt(totalSquaredLenth/n -
                metric.MEAN_INSERT_SIZE*metric.MEAN_INSERT_SIZE);
        metric.STANDARD_DEVIATION *= Math.sqrt(n/(n-1));

        // second round get median done
        final Long midPoint = pairsCount/2;
        pairsCount = 0L; // reset
        for(Map.Entry<Integer, Long> pair : hist.entrySet()){
            final Integer size = pair.getKey();
            final Long cnt = pair.getValue();
            pairsCount += cnt;
            if(pairsCount > midPoint) {
                metric.MEDIAN_INSERT_SIZE = size;
                break;
            }else if(pairsCount.equals(midPoint)){ // borderline case
                final Integer nextSize = hist.higherKey(size);
                metric.MEDIAN_INSERT_SIZE = ((double)(size + nextSize))/2.0;
                break;
            }
        }

        // third round to get MAD done (is there any functional way/map to do this?)
        TreeMap<Double, Long> medianDiff = new TreeMap<Double, Long>();
        for(Map.Entry<Integer, Long> pair : hist.entrySet()){
            Integer size = pair.getKey();
            Long cnt = pair.getValue();
            double diff =  Math.abs((double)size - metric.MEDIAN_INSERT_SIZE);
            if(medianDiff.get(diff)!=null){
                medianDiff.put(diff, medianDiff.get(diff)+cnt);
            }else{
                medianDiff.put(diff, cnt);
            }
        }
        pairsCount = 0L;
        for(Map.Entry<Double, Long> pair : medianDiff.entrySet()){
            final Double diff = pair.getKey();
            final Long cnt = pair.getValue();
            pairsCount += cnt;
            if(pairsCount > midPoint) {
                metric.MEDIAN_ABSOLUTE_DEVIATION = diff;
                break;
            }else if(pairsCount.equals(midPoint)){ // borderline case
                final double nextDiff = medianDiff.higherKey(diff);
                metric.MEDIAN_ABSOLUTE_DEVIATION = ((diff + nextDiff))/2.0;
                break;
            }
        }
    }

    private void collectPercentages(InsertSizeMetrics metric){

        double bins[][] = new double[10][2];
        final double offsets[] = {0.005, 0.05, 0.10, 0.15, 0.2, 0.25, 0.30, 0.35, 0.40, 0.45};
        // final round, percentiles
        for(int i=0; i<10; ++i){
            final double left = (double)metric.MIN_INSERT_SIZE + offsets[i]*(double)metric.READ_PAIRS;
            final Integer leftOne = hist.higherKey((int) Math.round(left));
            final Integer leftTwo = hist.lowerKey((int) Math.round(left));

            if( (left - leftOne) < (leftTwo - left) )
                bins[i][0] = leftOne;
            else
                bins[i][0] = leftTwo;

            final double right = (double)metric.MIN_INSERT_SIZE + (1.0-offsets[i])*(double)metric.READ_PAIRS;
            final Integer rightOne = hist.higherKey((int) Math.round(right));
            final Integer rightTwo = hist.lowerKey((int) Math.round(right));
            if( (right - rightOne) < (rightTwo - right) )
                bins[i][1] = rightOne;
            else
                bins[i][1] = rightTwo;
        }

        metric.WIDTH_OF_10_PERCENT = (int) (bins[0][1]-bins[0][0]);
        metric.WIDTH_OF_20_PERCENT = (int) (bins[1][1]-bins[1][0]);
        metric.WIDTH_OF_30_PERCENT = (int) (bins[2][1]-bins[2][0]);
        metric.WIDTH_OF_40_PERCENT = (int) (bins[3][1]-bins[3][0]);
        metric.WIDTH_OF_50_PERCENT = (int) (bins[4][1]-bins[4][0]);
        metric.WIDTH_OF_60_PERCENT = (int) (bins[5][1]-bins[5][0]);
        metric.WIDTH_OF_70_PERCENT = (int) (bins[6][1]-bins[6][0]);
        metric.WIDTH_OF_80_PERCENT = (int) (bins[7][1]-bins[7][0]);
        metric.WIDTH_OF_90_PERCENT = (int) (bins[8][1]-bins[8][0]);
        metric.WIDTH_OF_99_PERCENT = (int) (bins[9][1]-bins[9][0]);
    }
}