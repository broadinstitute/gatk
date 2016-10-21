package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.*;

/**
 * MultiLevelCollector handles accumulating Metrics at different MetricAccumulationLevels(ALL_READS, SAMPLE, LIBRARY, READ_GROUP).
 * Based on the accumulationLevels and readGroup records passed to its constructor, MultiLevelCollector
 * will instantiate the number of PerUnitMetricCollector's needed to generate metrics for each of the levels provided.
 *
 * To Use:
 *
 * Instantiate a MultiLevelCollector and call setup(see thoughts about extending MultiLevelCollector below)
 * setup will create the underlying classes that will handle the accumulation level logic.
 * Pass all reads you wish to collect data against to MultiLevelCollector via the acceptRecord method
 * Call finish and use addAllLevelsToFile to add all of the metrics at each accumulation level to the given file.
 *
 * Extend MultiLevelCollector and implement makeArg and makeChildCollector
 * You will most likely want to make a class that extends PerUnitMetricCollector.  This class should do the work of keeping
 * track of values for one specific "accumulation unit" (e.g. for one library, or for one read group depending on what levels
 * you are accumulating at).
 *
 * If a record has any expensive calculations to be done (that don't need to be done differently depending
 * on what sample/library/read group the read is for) then create a container class for the results of these calculations and pass
 * this class as the ARGTYPE of both the PerUnitMetricCollector and MultiLevelCollector.  You can then do these calculations in the makeArg
 * method and they will only be done once per record.
 *
 * @param <METRIC_TYPE> The type of metrics being collected
 * @param <HISTOGRAM_KEY> If there is are Histograms related to metrics of type <BEAN> then <HKEY> is the key value to these Histograms
 * @param <ARGTYPE> The type of argument passed to individual PerUnitMetricCollector (see SAMRecordMultilevelCollector and PerUnitMetricCollector)
 */
public abstract class MultiLevelCollector<METRIC_TYPE extends MetricBase, HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>, ARGTYPE>  {

    public static final String UNKNOWN = "unknown";
    //The collector that will accept all records (allReads is NULL if !calculateAll)
    private PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> allReadCollector;

    //A list of Distributor that is at most length 4, 1 for each (ALL_READS, SAMPLE, LIBRARY, READ_GROUP) accumulation levels
    //these will be listed in the order in which their children would be added to a metric file
    private List<Distributor> outputOrderedDistributors;

    //Convert the current SAMRecord and the ReferenceSequence for that record into an ARGTYPE object
    //see accept record for use
    protected abstract ARGTYPE makeArg(final SAMRecord samRec, final ReferenceSequence refSeq);

    /**
     * Construct a PerUnitMetricCollector with the given arguments.
     * @param sample If aggregating by ALL_READS this will be null, otherwise the sample that will be used to identify
     *               this collector
     * @param library If aggregating by SAMPLE this will be null, otherwise the library that will be used to identify
     *               this collector
     * @param readGroup If aggregating by LIBRARY this will be null, otherwise the readGroup that will be used to identify
     *                  this collector
     * @return A PerUnitMetricCollector parameterized by the given arguments
     */
    protected abstract PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeChildCollector(final String sample, final String library, final String readGroup);

    //These are exposed here (rather than being encapsulated in the Distributor subclasses below in order
    //to provide subclasses with an explicit point to add initialization (specific to accumulation level) for
    //a PerUnitMetricCollector it is creating
    protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeAllReadCollector() {
        return makeChildCollector(null, null, null);
    }
    protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeSampleCollector(final SAMReadGroupRecord rg) {
        return makeChildCollector(rg.getSample(), null, null);
    }
    protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeLibraryCollector(final SAMReadGroupRecord rg) {
        return makeChildCollector(rg.getSample(), rg.getLibrary(), null);
    }
    protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeReadGroupCollector(final SAMReadGroupRecord rg) {
        return makeChildCollector(rg.getSample(), rg.getLibrary(), rg.getPlatformUnit());
    }

    /**
     * Distributors group PerUnitMetricCollectors based on a AccumulationLevel.  Their structure mimics
     * PerUnitMetricCollectors but instead of adding records to metrics they identify which
     * PerUnitMetricCollector should receive a specific record and distribute records on to the that collector
     *
     * There were will be 0 or 1 Distributors for each of the following MetriAcummulationLevels:
     * ALL_READS, SAMPLE, LIBRARY, READ_GROUP
     */
    private abstract class Distributor {
        //A Map mapping the key for a specific record (as determined by getKey) to the appropriate collector
        private final Map<String, PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>> collectors;

        //Given a SAMReadGroupRecord, return the key that identifies the collector for the corresponding SAMRecord
        protected abstract String getKey(final SAMReadGroupRecord rg);

        //Make a PerUnitMetricCollector for this given Distributor
        protected abstract PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeCollector(final SAMReadGroupRecord rg);

        protected abstract PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeUnknownCollector();

        public Distributor(final List<SAMReadGroupRecord> rgRecs) {
            collectors = new LinkedHashMap<>();
            for(final SAMReadGroupRecord rg : rgRecs) {
                final String key = getKey(rg);
                if(!collectors.containsKey(key)) {
                    collectors.put(key, makeCollector(rg));
                }
            }
        }

        /** Call finish on each PerUnitMetricCollector in this Aggregate Collector */
        public void finish() {
            for(final PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> collector : collectors.values()) {
                collector.finish();
            }
        }

        /** Call acceptRecord(args) on the record collector identified by getKey */
        public void acceptRecord(final ARGTYPE args, final SAMReadGroupRecord rg) {

            String key = UNKNOWN;
            if(rg != null) {
                final String computedKey = getKey(rg);
                if(computedKey != null) {
                    key = computedKey;
                }
            }
            PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> collector = collectors.get(key);
            if (collector == null) {
                if (!UNKNOWN.equals(key)) {
                    throw new GATKException("Could not find collector for " + key);
                }
                collector = makeUnknownCollector();
                collectors.put(key, collector);
            }
            collector.acceptRecord(args);
        }

        /** Add all records to the MetricsFile passed in, this will happen in the order they were
         * found in the input ReadGroup records */
        public void addToFile(final MetricsFile<METRIC_TYPE, HISTOGRAM_KEY> file) {
            for(final PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> collector : collectors.values()) {
                collector.addMetricsToFile(file);
            }
        }
    }

    /** A dummy Distributor to handle the ALL_READS accumulation level.  No distribution is required
     * since there should only ever be one PerUnitMetricCollector for ALL_READS.
     */
    private class AllReadsDistributor extends Distributor {

        public AllReadsDistributor(final List<SAMReadGroupRecord> rgRecs) {
            super(new ArrayList<>());
            makeCollector(null);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return null;
        }

        @Override
        public void acceptRecord(final ARGTYPE args, final SAMReadGroupRecord rg) {
            allReadCollector.acceptRecord(args);
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeCollector(final SAMReadGroupRecord rg) {
            allReadCollector = makeAllReadCollector();
            return allReadCollector;
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeUnknownCollector() {
            throw new UnsupportedOperationException("Should not happen");
        }

        @Override
        public void finish() {
            allReadCollector.finish();
        }

        @Override
        public void addToFile(final MetricsFile<METRIC_TYPE, HISTOGRAM_KEY> file) {
            allReadCollector.addMetricsToFile(file);
        }
    }

    //Discriminates between records based on sample name, and calls acceptRecord on the appropriate PerUnitMetricCollectors
    private class SampleDistributor extends Distributor {
        public SampleDistributor(final List<SAMReadGroupRecord> rgRecs) {
            super(rgRecs);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return rg.getSample();
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeCollector(SAMReadGroupRecord rg) {
            return makeSampleCollector(rg);
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeUnknownCollector() {
            return makeChildCollector(UNKNOWN, null, null);
        }
    }

    //Discriminates between records based on library name, and calls acceptRecord on the appropriate PerUnitMetricCollectors
    private class LibraryDistributor extends Distributor {
        public LibraryDistributor(final List<SAMReadGroupRecord> rgRecs) {
            super(rgRecs);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return rg.getLibrary();
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeCollector(SAMReadGroupRecord rg) {
            return makeLibraryCollector(rg);
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeUnknownCollector() {
            return makeChildCollector(UNKNOWN, UNKNOWN, null);
        }
    }

    //Discriminates between records based on read group name, and calls acceptRecord on the appropriate PerUnitMetricCollectors
    private class ReadGroupCollector extends Distributor {
        public ReadGroupCollector(final List<SAMReadGroupRecord> rgRecs) {
            super(rgRecs);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return rg.getPlatformUnit();
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeCollector(SAMReadGroupRecord rg) {
            return makeReadGroupCollector(rg);
        }

        @Override
        protected PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> makeUnknownCollector() {
            return makeChildCollector(UNKNOWN, UNKNOWN, UNKNOWN);
        }
    }

    /**
     * Use an init method so that overloaded methods in subclasses can pass use parameters that are initialized in their constructor
     * @param accumulationLevels PerUnitMetricCollectors will only be created for the levels identified by accumulationLevels
     * @param samRgRecords PerUnitMetricCollectors will be created for each of the different samples, libraries, and
     *                     readGroups found in the records depending on the accumulationLevels provided
     */
    protected void setup(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords) {
        outputOrderedDistributors = new ArrayList<>(4);
        if(accumulationLevels.contains(MetricAccumulationLevel.ALL_READS)) {
            outputOrderedDistributors.add(new AllReadsDistributor(samRgRecords));
        }
        if (accumulationLevels.contains(MetricAccumulationLevel.SAMPLE)) {
            outputOrderedDistributors.add(new SampleDistributor(samRgRecords));
        }

        if(accumulationLevels.contains(MetricAccumulationLevel.LIBRARY)) {
            outputOrderedDistributors.add(new LibraryDistributor(samRgRecords));
        }

        if(accumulationLevels.contains(MetricAccumulationLevel.READ_GROUP)) {
            outputOrderedDistributors.add(new ReadGroupCollector(samRgRecords));
        }
    }

    /**
     * Construct a argument of ARGTYPE using the given SAMRecord and ReferenceSequence then pass
     * this value to all collectors that should include this record
     */
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        final ARGTYPE arg = makeArg(record, refSeq);

        for(final Distributor collector : outputOrderedDistributors) {
            collector.acceptRecord(arg, record.getReadGroup());
        }
    }

    /**
     * Call finish on all PerUnitMetricCollectors
     */
    public void finish() {
        for(final Distributor collector : outputOrderedDistributors) {
            collector.finish();
        }
    }

    /** Get the PerUnitMetricCollector that collects reads for all levels */
    public PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE> getAllReadsCollector() {
        return allReadCollector;
    }

    /** Add all metrics to the given file in the following MetricAccumulationLevel order
     *  ALL_READS, SAMPLE, LIBRARY, READ_GROUP.
     */
    public void addAllLevelsToFile(final MetricsFile<METRIC_TYPE, HISTOGRAM_KEY> file) {
        for(final Distributor collector : outputOrderedDistributors) {
            collector.addToFile(file);
        }
    }
}
