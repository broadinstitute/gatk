package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;

import org.apache.spark.api.java.JavaRDD;

import scala.Tuple2;

import java.util.*;
import java.util.function.Function;
import java.util.stream.StreamSupport;
import java.util.stream.Collectors;
import java.io.Serializable;

// TODO: is it desired that mapping quality is collected as well?
// TODO: consider refactoring code for other metric collection usages
/**
 * Worker class to collect insert size metrics, add metrics to file, and provides accessors to stats of groups of different level.
 */
public final class InsertSizeMetricsCollectorSpark implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Constructor; also acts as driver for organizing workflow.
     *
     * @param filteredReads         reads that pass filters
     * @param header                header in the input
     * @param accumLevels           accumulation level {ALL_READS, SAMPLE, LIBRARY, READ_GROUP}
     * @param histogramMADTolerance MAD tolerance when producing histogram plot
     * @param metricsFile           metrics file to write InsertSizeMetrics to
     */
    public InsertSizeMetricsCollectorSpark(final JavaRDD<GATKRead> filteredReads,
                                           final SAMFileHeader header,
                                           final Set<MetricAccumulationLevel> accumLevels,
                                           final double histogramMADTolerance,
                                           final MetricsFile<InsertSizeMetrics, Integer> metricsFile) {

        /* General strategy:
           construct untrimmed hand rolled "histogram" (SortedMap) of all read groups in three steps,
             because htsjdk Histogram does not play well with Serialization
           so first hand roll a histogram at the read-group level during the 1-pass all reads traversal on Spark.
           when computing InsertSizeMetrics, first convert the SortedMap version histogram to htsjdk Histogram, then compute metrics.*/

        final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> histogramsOfReadGroups = filteredReads.mapToPair(read -> traverseReadsToExtractInfo(read, header))
                                                                                                                                   .groupByKey()
                                                                                                                                   .mapToPair(InsertSizeMetricsCollectorSpark::gatherSizesByOrientation)
                                                                                                                                   .mapToPair(InsertSizeMetricsCollectorSpark::constructHistogramFromList)
                                                                                                                                   .collectAsMap();

        // accumulate for higher levels, if so desired
        // returns a list, of the same size as the number of distinct levels requested,
        // each entry in the list represents a level, and the entry is a map from the group's meta info to its metrics info
        final ArrayList<Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> listOfHistograms = aggregateHistograms(histogramsOfReadGroups, accumLevels);
        final ArrayList<Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>>>> listOfStats = new ArrayList<>();
        for(int i=0; i<accumLevels.size(); ++i){ listOfStats.add(new HashMap<>()); }

        // convert to htsjdk Histogram and compute metrics
        for(int i=0; i<accumLevels.size(); ++i){ convertSortedMapToHTSHistogram(listOfStats.get(i), listOfHistograms.get(i), histogramMADTolerance); }

        for(final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>>> statsInfo : listOfStats) { dumpToFile(metricsFile, statsInfo); }
    }

    /**
     * Using histograms at the read group level to construct histograms at higher levels.
     * @param histogramsAtRGLevel   histograms of read groups
     * @param accumLevels           accumulation level {ALL_READS, SAMPLE, LIBRARY, READ_GROUP}
     * @return                      a list of histograms, where each entry in the list--a map from group into to its histograms--represents one level
     */
    @VisibleForTesting
    static ArrayList<Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> aggregateHistograms(final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> histogramsAtRGLevel,
                                                                                                                         final Set<MetricAccumulationLevel> accumLevels){

        final List<Tuple2<ReadGroupParentExtractor,
                          Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>>> extractors = new ArrayList<>();

        if(accumLevels.contains(MetricAccumulationLevel.LIBRARY))   { extractors.add( new Tuple2<>(new ReadGroupLibraryExtractor(),  new HashMap<>())); }
        if(accumLevels.contains(MetricAccumulationLevel.SAMPLE))    { extractors.add( new Tuple2<>(new ReadGroupSampleExtractor(),   new HashMap<>())); }
        if(accumLevels.contains(MetricAccumulationLevel.ALL_READS)) { extractors.add( new Tuple2<>(new ReadGroupAllReadsExtractor(), new HashMap<>())); }

        for(final GroupMetaInfo groupMetaInfo : histogramsAtRGLevel.keySet()){
            final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> readGroupHistograms = histogramsAtRGLevel.get(groupMetaInfo);
            for(final Tuple2<ReadGroupParentExtractor,
                             Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> extractor : extractors){
                distributeRGHistogramsToAppropriateLevel(groupMetaInfo, readGroupHistograms, extractor);
            }
        }

        final ArrayList<Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> listOfHistograms = new ArrayList<>();

        if(accumLevels.contains(MetricAccumulationLevel.READ_GROUP)){ listOfHistograms.add(histogramsAtRGLevel); }

        for(final Tuple2<ReadGroupParentExtractor,
                         Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> extractor : extractors){
            listOfHistograms.add(extractor._2());
        }
        return listOfHistograms;
    }

    /** Distributes a particular read group's histogram to its "parent" library/sample or all reads if requested.
     * @param readGroupMetaInfo    meta-information of the read group to be distributed
     * @param readGroupHistograms  histogram of the read group to be distributed
     * @param extractor            a tuple, where first is a functor that extracts meta-information of the read group and
     *                             second is the destination of where the read should be distributed to.
     */
    @VisibleForTesting
    static void distributeRGHistogramsToAppropriateLevel(final GroupMetaInfo readGroupMetaInfo,
                                                         final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> readGroupHistograms,
                                                         final Tuple2<ReadGroupParentExtractor, Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> extractor){

        final GroupMetaInfo correspondingHigherLevelGroup = extractor._1().extractParentGroupMetaInfo(readGroupMetaInfo);
        final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> destination = extractor._2();

        // three checks: first check if this higher level group has been seen yet.
        destination.putIfAbsent(correspondingHigherLevelGroup, new HashMap<>());
        final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> higherLevelHistograms = destination.get(correspondingHigherLevelGroup);

        for(final SamPairUtil.PairOrientation orientation : readGroupHistograms.keySet()){
            higherLevelHistograms.putIfAbsent(orientation, new TreeMap<>()); // second check if this orientation has been seen yet.
            final SortedMap<Integer, Long> higherLevelHistogramOfThisOrientation = higherLevelHistograms.get(orientation);
            final SortedMap<Integer, Long> readGroupHistogramOfThisOrientation = readGroupHistograms.get(orientation);
            for(final Integer bin : readGroupHistogramOfThisOrientation.keySet()){
                higherLevelHistogramOfThisOrientation.putIfAbsent(bin, 0L); // third check if this bin has been seen yet.
                Long count = higherLevelHistogramOfThisOrientation.get(bin);
                count += readGroupHistogramOfThisOrientation.get(bin);
                higherLevelHistogramOfThisOrientation.put(bin, count);
            }
        }
    }

    // functor to extract meta-information on a read group's "parent" (eg. library, sample) group
    private interface ReadGroupParentExtractor{
        GroupMetaInfo extractParentGroupMetaInfo(final GroupMetaInfo groupMetaInfo);
    }

    private static final class ReadGroupAllReadsExtractor implements ReadGroupParentExtractor{

        private static final GroupMetaInfo ALLREADS_META_INFO = new GroupMetaInfo(null, null, null, MetricAccumulationLevel.ALL_READS);

        public GroupMetaInfo extractParentGroupMetaInfo(final GroupMetaInfo readGroupMetaInfo){
            return ALLREADS_META_INFO;
        }
    }

    private static final class ReadGroupSampleExtractor implements ReadGroupParentExtractor{
        public GroupMetaInfo extractParentGroupMetaInfo(final GroupMetaInfo readGroupMetaInfo){
            return new GroupMetaInfo(readGroupMetaInfo.sample, null, null, MetricAccumulationLevel.SAMPLE);
        }
    }

    private static final class ReadGroupLibraryExtractor implements ReadGroupParentExtractor{
        public GroupMetaInfo extractParentGroupMetaInfo(final GroupMetaInfo readGroupMetaInfo){
            return new GroupMetaInfo(readGroupMetaInfo.sample, readGroupMetaInfo.library, null, MetricAccumulationLevel.LIBRARY);
        }
    }

    /**
     * Worker function where hand-rolled histograms (SortedMap) is converted to htsjdk Histograms, and metrics information is collected
     * @param htsjdkHistogramsAndMetrics  htsjdk Histogram and InsertSizeMetrics bundled together under particular pair orientations,
     *                                    for the same grouping that's fed in; converted results to be put here
     * @param rawHistograms               hand-rolled histogram to be converted
     * @param histogramMADTolerance       tolerance to trim histogram so "outliers" don't ruin mean and SD values
     */
    @VisibleForTesting
    static void convertSortedMapToHTSHistogram(Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>>> htsjdkHistogramsAndMetrics,
                                               final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> rawHistograms,
                                               final double histogramMADTolerance){

        for(final GroupMetaInfo groupMetaInfo : rawHistograms.keySet()){
            final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> rawHistogramsOfAGroup = rawHistograms.get(groupMetaInfo);
            for(final SamPairUtil.PairOrientation orientation : rawHistogramsOfAGroup.keySet()){
                // convert to htsjdk Histogram
                final Histogram<Integer> htsHist = new Histogram<>("insert_size", groupMetaInfo.getGroupName() + "." + orientationToString(orientation) + "_count");
                final SortedMap<Integer, Long> hist = rawHistogramsOfAGroup.get(orientation);
                for(final int size : hist.keySet()){
                    htsHist.increment(size, hist.get(size));
                }

                final InsertSizeMetrics metrics = new InsertSizeMetrics();
                metrics.PAIR_ORIENTATION = orientation;

                collectMetricsBaseInfo(metrics, groupMetaInfo);
                collectSimpleStats(metrics, htsHist);
                collectSymmetricBinWidth(metrics, htsHist);
                trimHTSHistogramAndSetMean(metrics, htsHist, histogramMADTolerance);

                // save result
                final Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>> mapToStats = new HashMap<>();
                mapToStats.put(orientation, new Tuple2<>(htsHist, metrics));
                htsjdkHistogramsAndMetrics.put(groupMetaInfo, mapToStats);
            }
        }
    }

    // small utility function to collect MetricsBase information
    private static void collectMetricsBaseInfo(InsertSizeMetrics metrics, final GroupMetaInfo groupMetaInfo){
        metrics.SAMPLE     = groupMetaInfo.sample;
        metrics.LIBRARY    = groupMetaInfo.library;
        metrics.READ_GROUP = groupMetaInfo.readGroup;
    }

    // small utility function to collect simple stats information
    private static void collectSimpleStats(InsertSizeMetrics metrics, final Histogram<Integer> htsHist){
        metrics.READ_PAIRS                = (long) htsHist.getSumOfValues();
        metrics.MIN_INSERT_SIZE           = (int) htsHist.getMin();
        metrics.MAX_INSERT_SIZE           = (int) htsHist.getMax();
        metrics.MEDIAN_INSERT_SIZE        = htsHist.getMedian();
        metrics.MEDIAN_ABSOLUTE_DEVIATION = htsHist.getMedianAbsoluteDeviation();
    }

    // small utility function to collect bin width on the untrimmed Histogram, but actual work delegated to computeRanges
    private static void collectSymmetricBinWidth(InsertSizeMetrics metrics, final Histogram<Integer> hist){
        final long bin_widths[] = computeRanges(hist, (int) hist.getMedian(), metrics.READ_PAIRS); // metrics.REAR_PAIRS is assumed to be set properly already
        metrics.WIDTH_OF_10_PERCENT = (int) bin_widths[0];
        metrics.WIDTH_OF_20_PERCENT = (int) bin_widths[1];
        metrics.WIDTH_OF_30_PERCENT = (int) bin_widths[2];
        metrics.WIDTH_OF_40_PERCENT = (int) bin_widths[3];
        metrics.WIDTH_OF_50_PERCENT = (int) bin_widths[4];
        metrics.WIDTH_OF_60_PERCENT = (int) bin_widths[5];
        metrics.WIDTH_OF_70_PERCENT = (int) bin_widths[6];
        metrics.WIDTH_OF_80_PERCENT = (int) bin_widths[7];
        metrics.WIDTH_OF_90_PERCENT = (int) bin_widths[8];
        metrics.WIDTH_OF_99_PERCENT = (int) bin_widths[9];
    }

    // Computes width of symmetrical bins around the histogram's median
    @VisibleForTesting
    @SuppressWarnings("unchecked") // suppress warning on type inference when calling hist.get(int)
    static long[] computeRanges(final Histogram<Integer> hist, final int start, final double totalCount){

        double sum = 0.0;  // for calculating coverage, stored as sum to avoid frequent casting

        int left = start;  // left and right boundaries of histogram bins
        int right = left;  //      start from median, and gradually open up

        final long bin_widths[] = new long[10];   // for storing distance between left and right boundaries of histogram bins
        // dimension is 10 because metrics requires 10 histogram bin width values.
        int i = 0;
        int j = 0;                          // represent lowest and highest indices of bin_widths that needs to be updated

        while (i < 10) {                        // until all width values are computed
            final Histogram<Integer>.Bin leftBin = hist.get(left);
            final Histogram<Integer>.Bin rightBin = (left != right) ? hist.get(right) : null;
            if (null != leftBin) {// since left and right are incremented/decremented by 1, they may end up not in Histogram's bins.
                sum += leftBin.getValue();
            }
            if (null != rightBin) {
                sum += rightBin.getValue();
            }

            j = (int) (10. * sum / totalCount); // if coverage increased by enough value, update necessary ones
            for (int k = i; k < j; ++k) {
                bin_widths[k] = right - left + 1;
            }
            i = j;                          // and update pointers

            --left;
            ++right;
        }

        return bin_widths;
    }

    // small utility function to trim htsjdk Histogram and set corresponding metric's mean and SD
    private static void trimHTSHistogramAndSetMean(InsertSizeMetrics metrics, Histogram<Integer> htsHist, final double histogramMADTolerance){
        htsHist.trimByWidth( (int)(metrics.MEDIAN_INSERT_SIZE + histogramMADTolerance*metrics.MEDIAN_ABSOLUTE_DEVIATION) );
        metrics.MEAN_INSERT_SIZE   = htsHist.getMean();
        if(1==htsHist.getCount()){ // extremely unlikely in reality, but may be true when running tests
            metrics.STANDARD_DEVIATION = 0.0;
        }else{
            metrics.STANDARD_DEVIATION = htsHist.getStandardDeviation();
        }
    }

    private static String orientationToString(final SamPairUtil.PairOrientation orientation){
        return orientation.equals(SamPairUtil.PairOrientation.FR) ? "fr" : (orientation.equals(SamPairUtil.PairOrientation.RF) ? "rf" : "tandem");
    }

    // utility functions to do mapping on RDDs; broken into three steps for easier comprehension
    // this is the first step to traverse all valid reads and extract relevant information
    private static Tuple2<GroupMetaInfo, Tuple2<SamPairUtil.PairOrientation, Integer>> traverseReadsToExtractInfo(final GATKRead read, final SAMFileHeader header){
        final GroupMetaInfo readsGroupMetaInfo = new GroupMetaInfo(read, header, MetricAccumulationLevel.READ_GROUP);
        final Tuple2<SamPairUtil.PairOrientation, Integer> readsPairInfo = new Tuple2<>(SamPairUtil.getPairOrientation(read.convertToSAMRecord(header)), Math.abs(read.getFragmentLength()));
        return new Tuple2<>(readsGroupMetaInfo, readsPairInfo);
    }

    // Maps an iterable of pairs of (orientation, length) to a map where length values are grouped by orientations
    private static Tuple2<GroupMetaInfo, Map<SamPairUtil.PairOrientation, List<Integer>>> gatherSizesByOrientation(final Tuple2<GroupMetaInfo, Iterable<Tuple2<SamPairUtil.PairOrientation, Integer>>> entry){
        return new Tuple2<>(entry._1(), StreamSupport.stream(entry._2().spliterator(), false)
                                                     .collect(Collectors.groupingBy(Tuple2::_1,
                                                                                    Collectors.mapping(Tuple2::_2, Collectors.toList()))));
    }

    // Maps a list of fragment size length values to a histogram, implemented as SortedMap
    private static Tuple2<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> constructHistogramFromList(final Tuple2<GroupMetaInfo, Map<SamPairUtil.PairOrientation, List<Integer>>> entry){

        final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> orientationToSortedMap = new HashMap<>();

        for(final Map.Entry<SamPairUtil.PairOrientation, List<Integer>> e : entry._2().entrySet()){
            orientationToSortedMap.put(e.getKey(),
                                       new TreeMap<>( e.getValue().stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting())) ));
        }
        return new Tuple2<>(entry._1(), orientationToSortedMap);
    }

    /**
     * Write metrics and histograms to flat file, with an order such that coarser level information appear in the file before finer levels.
     * @param metricsFile File to write information to
     * @param stats       stat information for all the groups whose insert size stats has been computed
     */
    private static void dumpToFile(final MetricsFile<InsertSizeMetrics, Integer> metricsFile,
                                   final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>>> stats){
        for(final GroupMetaInfo groupMetaInfo : stats.keySet()){
            for(final SamPairUtil.PairOrientation orientation : stats.get(groupMetaInfo).keySet()){
                metricsFile.addMetric(stats.get(groupMetaInfo).get(orientation)._2());
                metricsFile.addHistogram(stats.get(groupMetaInfo).get(orientation)._1());
            }
        }
    }

    /**
     * A struct containing relevant information for a particular group of reads, where a "group" could be a read group,
     *   a library, a sample, or all reads (where there should be only one such group in the input).
     */
    @VisibleForTesting
    static final class GroupMetaInfo implements Serializable{
        private static final long serialVersionUID = 1L;

        public final String sample;
        public final String library;
        public final String readGroup;
        public final MetricAccumulationLevel level;

        // Intellij-generated code
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            GroupMetaInfo groupMetaInfo = (GroupMetaInfo) o;

            if (sample != null ? !sample.equals(groupMetaInfo.sample) : groupMetaInfo.sample != null) return false;
            if (library != null ? !library.equals(groupMetaInfo.library) : groupMetaInfo.library != null) return false;
            if (readGroup != null ? !readGroup.equals(groupMetaInfo.readGroup) : groupMetaInfo.readGroup != null) return false;
            return level == groupMetaInfo.level;

        }

        // Intellij-generated code
        @Override
        public int hashCode() {
            int result = sample != null ? sample.hashCode() : 0;
            result = 31 * result + (library != null ? library.hashCode() : 0);
            result = 31 * result + (readGroup != null ? readGroup.hashCode() : 0);
            result = 31 * result + (level != null ? level.hashCode() : 0);
            return result;
        }

        public GroupMetaInfo(final GATKRead read, final SAMFileHeader header, final MetricAccumulationLevel level){
            this.sample    = header.getReadGroup(read.getReadGroup()).getSample();
            this.library   = header.getReadGroup(read.getReadGroup()).getLibrary();
            this.readGroup = read.getReadGroup();
            this.level     = level;
        }

        public GroupMetaInfo(final String sample, final String library, final String readGroup, final MetricAccumulationLevel level){
            this.sample    = sample;
            this.library   = library;
            this.readGroup = readGroup;
            this.level     = level;
        }

        /**
         * small utility function to decide what group name to use in the corresponding htsjdk Histogram title/ctor.
         */
        public String getGroupName(){
            String groupName = null;
            switch (this.level){
                case ALL_READS:
                    groupName = "All_reads";
                    break;
                case SAMPLE:
                    groupName = this.sample;
                    break;
                case LIBRARY:
                    groupName = this.library;
                    break;
                case READ_GROUP:
                    groupName = this.readGroup;
                    break;
            }
            return groupName;
        }
    }
}