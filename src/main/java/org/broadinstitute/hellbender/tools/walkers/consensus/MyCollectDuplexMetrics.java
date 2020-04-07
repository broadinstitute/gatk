package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts;
import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "asdf",
        oneLineSummary = "asdf",
        programGroup = ReadDataManipulationProgramGroup.class
)
public class MyCollectDuplexMetrics extends DuplicateSetWalker {
    List<SimpleInterval> intervals;
    private String currentContig = null;

    /**
     * Overlap detector used to determine when read starts overlap with input intervals.
     */
    private CollectReadCounts.CachedOverlapDetector<SimpleInterval> intervalCachedOverlapDetector;

    private Map<SimpleInterval, MyMetrics> intervalMyMetricsMap;

    private Pair<SimpleInterval, MyMetrics> cachedIntervalMetricsPair;

    @Argument(fullName = "targets")
    public File targets;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "")
    public File output;

    public static final String READ_NUMBER_LIMIT_NAME = "read-number-limit";
    @Argument(fullName = READ_NUMBER_LIMIT_NAME, doc="read number limit", optional = true)
    public int readNumberLimit = Integer.MAX_VALUE;

    int count = 0;

    @Override
    public int getReadNumberLimit(){
        return readNumberLimit;
    }

    @Override
    public boolean requiresReference() { return true; }

    // turn this to true for a one-time annotation of the target file with reference context
    private boolean annotateReferenceContext = true;

    @Override
    public void onTraversalStart(){
        intervalMyMetricsMap = new TreeMap<>(new SimpleIntervalComparator(getHeaderForReads()));
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<GenomeLoc> yup = IntervalUtils.parseIntervalArguments(new GenomeLocParser(sequenceDictionary), Arrays.asList(targets.getAbsolutePath()));
        intervals = IntervalUtils.convertGenomeLocsToSimpleIntervals(yup);
        final ReferenceDataSource refDataSource = iDontCareGimmeReferenceDataSource();

        for (SimpleInterval interval : intervals) {
            final MyMetrics metrics = new MyMetrics();
            // see HaplotypeCallerGenotypingEngine::makeAnnotatedCall()

            if (annotateReferenceContext){
                final ReferenceContext referenceContext = new ReferenceContext(refDataSource, interval);
                final String targetReferenceBases = new String(referenceContext.getBases());

                metrics.setReferenceBases(targetReferenceBases);
            }

            intervalMyMetricsMap.put(interval, metrics);
        }


        cachedIntervalMetricsPair = new ImmutablePair<>(intervals.get(0), intervalMyMetricsMap.get(intervals.get(0)));
    }

    @Override
    public void apply(DuplicateSet duplicateSet, ReferenceContext referenceContext, FeatureContext featureContext) {
        SimpleInterval duplicateSetInterval = duplicateSet.getDuplicateSetInterval();
        if (currentContig == null || ! duplicateSetInterval.getContig().equals(currentContig)) {
            //if we are on a new contig, create an OverlapDetector covering the contig
            currentContig = duplicateSetInterval.getContig();
            final List<SimpleInterval> intervalsOnCurrentContig = intervals.stream()
                    .filter(i -> i.getContig().equals(currentContig))
                    .collect(Collectors.toList());
            intervalCachedOverlapDetector = new CollectReadCounts.CachedOverlapDetector<>(intervalsOnCurrentContig);
        }

        // ts: overlap is by the start position...is that a problem?
        final SimpleInterval overlappingInterval = intervalCachedOverlapDetector.getOverlap(
                new SimpleInterval(duplicateSetInterval.getContig(), duplicateSetInterval.getStart(), duplicateSetInterval.getStart()));

        //if read doesn't overlap any of the provided intervals, do nothing
        if (overlappingInterval == null) {
            return;
        }

        if (overlappingInterval.equals(cachedIntervalMetricsPair.getLeft())){
            final MyMetrics metrics = cachedIntervalMetricsPair.getRight();
            update(duplicateSet, metrics);
        } else {
            final MyMetrics metrics = intervalMyMetricsMap.get(overlappingInterval);
            update(duplicateSet, metrics);
            cachedIntervalMetricsPair = new ImmutablePair<>(overlappingInterval, metrics);
        }
    }

    @Override
    public Object onTraversalSuccess(){
        logger.info("Writing metrics to file");
        try {
            final PrintWriter printWriter = new PrintWriter(output);
            StringBuilder header = new StringBuilder()
                    .append("contig,")
                    .append("start,")
                    .append("end,")
                    .append("simplex,")
                    .append("duplex,")
                    .append("reads");
            if (annotateReferenceContext){
                header.append(",reference");
            }
            printWriter.println(header.toString());

            for (Map.Entry<SimpleInterval, MyMetrics> entry : intervalMyMetricsMap.entrySet()){
                final SimpleInterval interval = entry.getKey();
                final MyMetrics metrics = entry.getValue();
                final StringBuilder content = new StringBuilder()
                        .append(interval.getContig() + ",")
                        .append(interval.getStart() + ",")
                        .append(interval.getEnd() + ",")
                        .append(metrics.simplexCount + ",")
                        .append(metrics.duplexCount + ",")
                        .append(metrics.numReads);
                if (annotateReferenceContext){
                    content.append("," + metrics.referenceBases);
                }

                printWriter.println(content.toString());
            }



        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        return "SUCCESS";
    }

    private void update(final DuplicateSet duplicateSet, final MyMetrics metrics){
        if (duplicateSet.isDuplex()){
            metrics.incrementDuplexCount();
        } else {
            metrics.incrementSimplexCount();
        }
        metrics.addReads(duplicateSet.getReads().size());
    }

    // More to come...
    public class MyMetrics {
        int duplexCount;
        int simplexCount;
        int numReads;
        String referenceBases;

        public MyMetrics(){
            this.duplexCount = 0;
            this.simplexCount = 0;
            this.numReads = 0;
        }

        public void incrementDuplexCount(){
            this.duplexCount += 1;
        }

        public void incrementSimplexCount(){
            this.simplexCount += 1;
        }

        public void addReads(int numReads){
            this.numReads += numReads;
        }

        public void setReferenceBases(String referenceBases){
            this.referenceBases = referenceBases;
        }

    }

    private class SimpleIntervalComparator implements Comparator<SimpleInterval> {
        SAMFileHeader header;
        public SimpleIntervalComparator(final SAMFileHeader header){
            this.header = header;
        }

        @Override
        public int compare(SimpleInterval interval1, SimpleInterval interval2) {
            return compareCoordinates(interval1, interval2, header);
        }

        // Copied from ReadCoordinateComparator --- could we combine the two? GATKRead is also Locatable,
        // but there seems to be a reason why (maybe unmapped reads?) they use "getAssignedStart" instead of
        // "getStart()"
        public int compareCoordinates(final Locatable first, final Locatable second, final SAMFileHeader header ) {
            final int firstRefIndex = header.getSequenceIndex(first.getContig());
            final int secondRefIndex = header.getSequenceIndex(second.getContig());

            // First compare contigs
            if ( firstRefIndex == -1 ) {
                return (secondRefIndex == -1 ? 0 : 1);
            }
            else if ( secondRefIndex == -1 ) {
                return -1;
            }

            final int refIndexDifference = firstRefIndex - secondRefIndex;
            if ( refIndexDifference != 0 ) {
                return refIndexDifference;
            }

            // Then the start position
            if (first.getStart() != second.getStart()){
                return Integer.compare(first.getStart(), second.getStart());
            }

            // Lastly, the end position
            return Integer.compare(first.getEnd(), second.getEnd());
        }
    }

}
