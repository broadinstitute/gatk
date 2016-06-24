package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;
import java.util.function.BiFunction;

/**
 * Abstract base class for reducible multi-level metrics collectors. Handles accumulating Metrics at
 * different MetricAccumulationLevels (ALL_READS, SAMPLE, LIBRARY, READ_GROUP), and can be reduced/
 * combined for when used to collect metrics parallel on multiple (Spark) partitions.
 *
 * Based on the accumulationLevels and readGroup records passed to its constructor, MultiLevelReducibleCollector
 * will instantiate the correct number of PerUnitMetricCollector's needed to generate metrics for each of
 * the levels provided. After collection is complete, the combine method can be used to reduce a series
 * of collectors into a single aggregate MultiLevelReducibleCollector-derived object that represents the result
 * metrics for an entire set of input reads.
 *
 * MultiLevelReducibleCollector requires a type parameter that represents the collector's UNIT_COLLECTOR
 * class (which must be a subclass of PerUnitMetricsCollector), and must also have a combineUnit method
 * that combines two UNIT_COLLECTOR objects. The default implementation of the combine method in
 * MultiLevelReducibleCollector combines another MultiLevelReducibleCollector by matching
 * up the the various per-unit collectors for each level of collection, and then combining them using the
 * combineUnit method. Generally the combineUnit method just delegates to the UNIT_COLLECTOR's combine method.
 *
 * See {@link org.broadinstitute.hellbender.metrics.InsertSizeMetricsCollector} and
 * {@link org.broadinstitute.hellbender.metrics.PerUnitInsertSizeMetricsCollector} as an example of
 * a pair of classes that follow this pattern. For more information and a higher level description
 * of the component classes involved in a Spark-aware collector, see
 * {@link org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSpark}.
 *
 * To use this class:
 *
 * Subclass UNIT_COLLECTOR and implement a combine method that combines two UNIT_COLLECTORs. Subclass
 * MultiLevelReducibleCollector, providing the UNIT_COLLECTOR-derived class as a type parameter, and implement
 * a combineUnit method that delegates to the UNIT_COLLECTOR's combine method.
 *
 * Instantiate a MultiLevelReducibleCollector and call setup, which will will create the underlying classes
 * that will handle the accumulation level logic. Pass all reads you wish to collect data against to
 * MultiLevelReducibleCollector via the {@link #acceptRecord method}. Call finish when all reads have been processed.
 * When multiple MultiLevelReducibleCollector have been created in parallel, the set of collectors can be
 * reduced using the combine method (which will aggregate like per-unit collectors and call the combineUnit
 * method).
 *
 * Extend MultiLevelReducibleCollector and implement makeArg and makeChildCollector
 * You will need to make a class that extends PerUnitMetricCollector.  This class should do the work of keeping
 * track of values for one specific "accumulation unit" (e.g. for one library, or for one read group depending on
 * what levels you are accumulating at).
 *
 * If a record has any expensive calculations to be done (that don't need to be done differently depending
 * on what sample/library/read group the read is for) then create a container class for the results of these
 * calculations and pass this class as the ARGTYPE of both the PerUnitMetricCollector and MultiLevelReducibleCollector.
 * You can then do these calculations in the makeArg method and they will only be done once per record.
 *
 * @param <METRIC_TYPE> The metrics being collected
 * @param <HISTOGRAM_KEY> If there is are Histograms related to metrics of type <BEAN> then <HKEY> is the key value
 *                       to these Histograms
 * @param <ARGTYPE> The type argument passed to individual PerUnitMetricCollector (see SAMRecordMultilevelCollector
 *                 and PerUnitMetricCollector)
 * @param <UNIT_COLLECTOR> Class to represent metrics collected at a single level
 */
public abstract class MultiLevelReducibleCollector<
        METRIC_TYPE extends MetricBase,
        HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>,
        ARGTYPE,
        UNIT_COLLECTOR extends PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>
        >
    implements Serializable
{
    private static final long serialVersionUID = 1L;

    private static final String UNKNOWN = "unknown";

    /**
     * A list of Distributor that is at most length 4, 1 for each (ALL_READS, SAMPLE, LIBRARY, READ_GROUP)
     * accumulation levels these will be listed in the order in which their children would be added to a
     * metric file
     */
    private List<Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR>> outputOrderedDistributors;

    // Keep track of the all reads distributor and collector key since their
    private AllReadsDistributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> allReadsDistributor;
    private static final String ALL_READS_COLLECTOR_KEY = "ALL_READS_COLLECTOR";

    /**
     * Convert the current SAMRecord and the ReferenceSequence for that record into an ARGTYPE object. See
     * {@link #acceptRecord method} for use.
     */
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
    protected abstract UNIT_COLLECTOR makeChildCollector(final String sample, final String library, final String readGroup);

    /**
     * These are exposed here (rather than being encapsulated in the Distributor subclasses below in order
     * to provide subclasses with an explicit point to add initialization (specific to accumulation level) for
     * a PerUnitMetricCollector it is creating.
     */
    protected UNIT_COLLECTOR makeAllReadCollector() {
        return makeChildCollector(null, null, null);
    }
    protected UNIT_COLLECTOR makeSampleCollector(final SAMReadGroupRecord rg) {
        return makeChildCollector(rg.getSample(), null, null);
    }
    protected UNIT_COLLECTOR makeLibraryCollector(final SAMReadGroupRecord rg) {
        return makeChildCollector(rg.getSample(), rg.getLibrary(), null);
    }
    protected UNIT_COLLECTOR makeReadGroupCollector(final SAMReadGroupRecord rg) {
        return makeChildCollector(rg.getSample(), rg.getLibrary(), rg.getPlatformUnit());
    }

    /**
     * Combine a source MultiLevelReducibleCollector object into this MultiLevelReducibleCollector.
     *
     * Most derived classes will want to provide a type-specific version of this method that delegates
     * to this method, but is typed suitably for use as a JavaRDD.reduce combiner method. Given a
     * derived class of type class "D", the combine method looks like this:
     *
     * public D combine(D target, D source) {
     *     target.combine(source);  // calls this implementation
     *     return target;           // return an object of type D
     * }
     *
     * Such a wrapper can be passed directly to JavaRDD.reduce.
     *
     * @param source source MultiLevelReducibleCollector to be combined with this object. May not be null.
     */
    public void combine(MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> source)
    {
        Utils.nonNull(source);
        if ((this.outputOrderedDistributors.size() != source.outputOrderedDistributors.size())) {
            throw new IllegalArgumentException("MultiLevelCollectors must the same size and level structure to be combined");
        }

        for (int i = 0; i < this.outputOrderedDistributors.size(); i++) {
            this.outputOrderedDistributors.get(i).combine(
                source.outputOrderedDistributors.get(i), this::combineUnit);
        }
    }

    /**
     * Combined two UNIT_COLLECTOR objects into a single object representing the combined collectors.
     *
     * Subclasses must implement this method. Note that it is acceptable for implementations to be
     * destructive and combine objects by mutating and returning one of the input objects.
     *
     * @param c1 first unit collector
     * @param c2 second unit collector
     * @return UNIT_COLLECTOR resulting from the combination of the two input collectors
     */
    public abstract UNIT_COLLECTOR combineUnit(UNIT_COLLECTOR c1, UNIT_COLLECTOR c2);

    /**
     * Distributors group PerUnitMetricCollectors based on a AccumulationLevel.  Their structure mimics
     * PerUnitMetricCollectors but instead of adding records to metrics they identify which
     * PerUnitMetricCollector should receive a specific record and distribute records on to the that collector
     *
     * There were will be 0 or 1 Distributors for each of the following MetricAcummulationLevels:
     * ALL_READS, SAMPLE, LIBRARY, READ_GROUP
     */
    private static abstract class Distributor <
                METRIC_TYPE extends MetricBase,
                HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>,
                ARGTYPE,
                UNIT_COLLECTOR extends PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>
            > implements Serializable {
        private static final long serialVersionUID = 1L;

        final protected MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> multiCollector;

        /**
         * A Map mapping the key for a specific record (as determined by getKey) to the appropriate collector
         */
        protected Map<String, UNIT_COLLECTOR> collectors;

        public Distributor(
                final MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> multiCollector) {
            this.multiCollector = multiCollector;
        }

        /**
         * Initialize the collectors from the given list of read groups.
         * @param rgRecs read groups to use to initialize collectors
         */
        public void initializeFromReadGroups(final List<SAMReadGroupRecord> rgRecs) {
            collectors = new LinkedHashMap<>();
            for(final SAMReadGroupRecord rg : rgRecs) {
                final String key = getKey(rg);
                if (!collectors.containsKey(key)) {
                    collectors.put(key, makeCollector(rg));
                }
            }
        }

        /**
         * Given a SAMReadGroupRecord, return a key that identifies the collector for the corresponding
         * records in this distributor.
         * @param rg the SAMReadGroupRecord
         * @return a key to be used to find the appropriate receiving collector
         */
        protected abstract String getKey(final SAMReadGroupRecord rg);

        /**
         * Make a UNIT_COLLECTOR for this Distributor given a read group record.
         *
         * @param rg The SAMReadGroupRecord to use to create the UNIT_COLLECTOR.
         * @return a UNIT_COLLECTOR for this level
         */
        protected abstract UNIT_COLLECTOR makeCollector(final SAMReadGroupRecord rg);

        /**
         * Create a collector to use as a catch-all for handling records that don't otherwise fit
         * the known level criteria.
         * @return UNIT_COLLECTOR for metrics for unknown records.
         */
        protected abstract UNIT_COLLECTOR makeUnknownCollector();

        /** Call acceptRecord(args) on the record collector identified by getKey */
        public void acceptRecord(final ARGTYPE args, final SAMReadGroupRecord rg) {

            String key = UNKNOWN;
            if(rg != null) {
                final String computedKey = getKey(rg);
                if(computedKey != null) {
                    key = computedKey;
                }
            }
            UNIT_COLLECTOR collector = collectors.get(key);
            if (collector == null) {
                if (!UNKNOWN.equals(key)) {
                    throw new GATKException("Could not find collector for " + key);
                }
                collector = makeUnknownCollector();
                collectors.put(key, collector);
            }
            collector.acceptRecord(args);
        }

        /** Call finish on each PerUnitMetricCollector in this Aggregate Collector */
        public void finish() {
            for(final UNIT_COLLECTOR collector : collectors.values()) {
                collector.finish();
            }
        }

        /**
         * Combine the source distributor into this distributor by combining the collector maps.
         * @param source source Distributor. May not be null.
         * @param reMap map entry merge (reMap) function for merging distributor map entries. This
         *              is usually the combineUnit method of MultiLevelReducibleCollector.
         */
        public void combine(
                Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> source,
                BiFunction<UNIT_COLLECTOR, UNIT_COLLECTOR, UNIT_COLLECTOR> reMap) {
            Utils.nonNull(source);
            Utils.nonNull(reMap);
            collectors.keySet().forEach(k -> collectors.merge(k, source.collectors.get(k), reMap));
        }

        /** Add all records to the MetricsFile passed in, this will happen in the order they were
         * found in the input ReadGroup records */
        public void addToFile(final MetricsFile<METRIC_TYPE, HISTOGRAM_KEY> file) {
            for(final UNIT_COLLECTOR collector : collectors.values()) {
                collector.addMetricsToFile(file);
            }
        }
    }

    /** A dummy Distributor to handle the ALL_READS accumulation level.  No distribution is required
     * since there should only ever be one PerUnitMetricCollector for ALL_READS.
     */
    private static class AllReadsDistributor<
                METRIC_TYPE extends MetricBase,
                HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>,
                ARGTYPE,
                UNIT_COLLECTOR extends PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>
            > extends Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> implements Serializable {

        private static final long serialVersionUID = 1L;

        public AllReadsDistributor(
                final MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> multiCollector) {
            super(multiCollector);
        }

        /**
         * Initialize the collectors from the given list of read groups. Delegate to the base class and
         * and pass in the empty list, and manually create the single all reads collector.
         * @param rgRecs read groups to use to initialize collectors
         */
        @Override
        public void initializeFromReadGroups(final List<SAMReadGroupRecord> rgRecs) {
            super.initializeFromReadGroups(new ArrayList<>());
            collectors.put(ALL_READS_COLLECTOR_KEY, multiCollector.makeAllReadCollector());
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return ALL_READS_COLLECTOR_KEY;
        }

        @Override
        protected UNIT_COLLECTOR makeCollector(final SAMReadGroupRecord rg) {
            return multiCollector.makeAllReadCollector();
        }

        @Override
        protected UNIT_COLLECTOR makeUnknownCollector() {
            throw new UnsupportedOperationException("Should not happen");
        }

    }

    //Discriminates between records based on sample name, and calls acceptRecord on the appropriate PerUnitMetricCollectors
    private static class SampleDistributor<
                METRIC_TYPE extends MetricBase,
                HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>,
                ARGTYPE,
                UNIT_COLLECTOR extends PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>
            > extends Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> implements Serializable {
        private static final long serialVersionUID = 1L;

        public SampleDistributor(
                final MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> multiCollector) {
            super(multiCollector);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return rg.getSample();
        }

        @Override
        protected UNIT_COLLECTOR makeCollector(SAMReadGroupRecord rg) {
            return multiCollector.makeSampleCollector(rg);
        }

        @Override
        protected UNIT_COLLECTOR makeUnknownCollector() {
            return multiCollector.makeChildCollector(UNKNOWN, null, null);
        }
    }

    //Discriminates between records based on library name, and calls acceptRecord on the appropriate PerUnitMetricCollectors
    private static class LibraryDistributor<
                METRIC_TYPE extends MetricBase,
                HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>,
                ARGTYPE,
                UNIT_COLLECTOR extends PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>
            > extends Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> implements Serializable  {
        private static final long serialVersionUID = 1L;

        public LibraryDistributor(
                final MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> multiCollector) {
            super(multiCollector);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return rg.getLibrary();
        }

        @Override
        protected UNIT_COLLECTOR makeCollector(SAMReadGroupRecord rg) {
            return multiCollector.makeLibraryCollector(rg);
        }

        @Override
        protected UNIT_COLLECTOR makeUnknownCollector() {
            return multiCollector.makeChildCollector(UNKNOWN, UNKNOWN, null);
        }
    }

    /**
     * Discriminates between records based on read group name, and calls acceptRecord on the appropriate
     * PerUnitMetricCollectors.
     */
    private class ReadGroupDistributor<
                METRIC_TYPE extends MetricBase,
                HISTOGRAM_KEY extends Comparable<HISTOGRAM_KEY>,
                ARGTYPE,
                UNIT_COLLECTOR extends PerUnitMetricCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE>
            > extends Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> implements Serializable {
        private static final long serialVersionUID = 1L;

        public ReadGroupDistributor(
                final MultiLevelReducibleCollector<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> multiCollector) {
            super(multiCollector);
        }

        @Override
        protected String getKey(SAMReadGroupRecord rg) {
            return rg.getPlatformUnit();
        }

        @Override
        protected UNIT_COLLECTOR makeCollector(SAMReadGroupRecord rg) {
            return multiCollector.makeReadGroupCollector(rg);
        }

        @Override
        protected UNIT_COLLECTOR makeUnknownCollector() {
            return multiCollector.makeChildCollector(UNKNOWN, UNKNOWN, UNKNOWN);
        }
    }

    /**
     * Use an init method so that overloaded methods in subclasses can pass use parameters that are initialized in their constructor
     * @param accumulationLevels PerUnitMetricCollectors will only be created for the levels identified by accumulationLevels
     * @param samRgRecords PerUnitMetricCollectors will be created for each of the different samples, libraries, and
     *                     readGroups found in the records depending on the accumulationLevels provided
     */
    protected void setup(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords) {
        outputOrderedDistributors = new ArrayList<>();
        if(accumulationLevels.contains(MetricAccumulationLevel.ALL_READS)) {
            // cache the allReadsDistributor so we can find it to hand out later
            allReadsDistributor = new AllReadsDistributor<>(this);
            allReadsDistributor.initializeFromReadGroups(samRgRecords);
            outputOrderedDistributors.add(allReadsDistributor);
        }
        if (accumulationLevels.contains(MetricAccumulationLevel.SAMPLE)) {
            SampleDistributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> sampleDistributor =
                    new SampleDistributor<>(this);
            sampleDistributor.initializeFromReadGroups(samRgRecords);
            outputOrderedDistributors.add(sampleDistributor);
        }

        if(accumulationLevels.contains(MetricAccumulationLevel.LIBRARY)) {
            LibraryDistributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> libraryDistributor =
                    new LibraryDistributor<>(this);
            libraryDistributor.initializeFromReadGroups(samRgRecords);
            outputOrderedDistributors.add(libraryDistributor);
        }

        if(accumulationLevels.contains(MetricAccumulationLevel.READ_GROUP)) {
            ReadGroupDistributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> readGroupDistributor =
                    new ReadGroupDistributor<>(this);
            readGroupDistributor.initializeFromReadGroups(samRgRecords);
            outputOrderedDistributors.add(readGroupDistributor);
        }
    }

    /**
     * Construct an argument of ARGTYPE using the given SAMRecord and ReferenceSequence, then pass
     * this value to all collectors that should include this record
     */
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        final ARGTYPE arg = makeArg(record, refSeq);

        for(final Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> distributor : outputOrderedDistributors) {
            distributor.acceptRecord(arg, record.getReadGroup());
        }
    }

    /**
     * Call finish on all PerUnitMetricCollectors
     */
    public void finish() {
        for(final Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> distributor : outputOrderedDistributors) {
            distributor.finish();
        }
    }

    /** Get the PerUnitMetricCollector that collects reads for all levels */
    public UNIT_COLLECTOR getAllReadsCollector() {
        return allReadsDistributor != null ?
                allReadsDistributor.collectors.get(ALL_READS_COLLECTOR_KEY) :
                null;
    }

    /**
     * Add all metrics to the given file in the following MetricAccumulationLevel order
     * ALL_READS, SAMPLE, LIBRARY, READ_GROUP.
     */
    public void addAllLevelsToFile(final MetricsFile<METRIC_TYPE, HISTOGRAM_KEY> file) {
        for(final Distributor<METRIC_TYPE, HISTOGRAM_KEY, ARGTYPE, UNIT_COLLECTOR> distributor : outputOrderedDistributors) {
            distributor.addToFile(file);
        }
    }
}
