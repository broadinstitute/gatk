package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ClusterData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProvider;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProviderFactory;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicBoolean;

import static htsjdk.samtools.util.Log.getInstance;
import static htsjdk.samtools.util.SortingCollection.newInstance;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Runtime.getRuntime;
import static java.lang.String.format;
import static java.lang.System.gc;
import static java.lang.System.runFinalization;
import static java.lang.Thread.currentThread;
import static java.util.Arrays.asList;
import static java.util.Arrays.copyOf;
import static java.util.Collections.sort;
import static java.util.concurrent.TimeUnit.MILLISECONDS;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.*;

/**
 * Manages the conversion of Illumina basecalls into some output format.  Creates multiple threads to manage reading,
 * sorting and writing efficiently.  Output is written in queryname output.  Optionally demultiplexes indexed reads
 * into separate outputs by barcode.
 *
 * @param <CLUSTER_OUTPUT_RECORD> The class to which a ClusterData is converted in preparation for writing.
 */
public class IlluminaBasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    /**
     * Describes the state of a barcode's data's processing in the context of a tile.  It is either not available in
     * that tile, has been read, has been queued to be written to file, or has been written to file.  A barcode only
     * takes on a state once the tile (which is serving as the context of this state) has been read.
     */
    private enum TileBarcodeProcessingState {
        NA, READ, QUEUED_FOR_WRITE, WRITTEN
    }

    /**
     * Describes the state of a tile being processed.  It is either not yet completely read, or read.
     */
    private enum TileProcessingState {
        NOT_DONE_READING, DONE_READING
    }

    private static final Log log = getInstance(IlluminaBasecallsConverter.class);

    public static final IlluminaDataType[] DATA_TYPES_NO_BARCODE =
            {BaseCalls, QualityScores, Position, PF};
    private static final IlluminaDataType[] DATA_TYPES_WITH_BARCODE = copyOf(DATA_TYPES_NO_BARCODE, DATA_TYPES_NO_BARCODE.length + 1);

    static {
        DATA_TYPES_WITH_BARCODE[DATA_TYPES_WITH_BARCODE.length - 1] = Barcodes;
    }

    /**
     * A comparator for tile numbers, which are not necessarily ordered by the number's value.
     */
    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = new Comparator<Integer>() {
        @Override
        public int compare(final Integer integer1, final Integer integer2) {
            final String s1 = integer1.toString();
            final String s2 = integer2.toString();
            // Because a the tile number is followed by a colon, a tile number that
            // is a prefix of another tile number should sort after. (e.g. 10 sorts after 100).
            if (s1.length() < s2.length()) {
                if (s2.startsWith(s1)) {
                    return 1;
                }
            } else if (s2.length() < s1.length()) {
                if (s1.startsWith(s2)) {
                    return -1;
                }
            }
            return s1.compareTo(s2);
        }
    };

    private final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;

    private final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    private final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    private final int maxReadsInRamPerTile;
    private final boolean demultiplex;
    private final List<File> tmpDirs;
    private final IlluminaDataProviderFactory factory;
    private ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter = null;
    private final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    private final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    private int numThreads;
    // If FORCE_GC, this is non-null.  For production this is not necessary because it will run until the JVM
    // ends, but for unit testing it is desirable to stop the task when done with this instance.
    private final TimerTask gcTimerTask;
    private List<Integer> tiles;
    private final boolean includeNonPfReads;
    private final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    // Annoying that we need this.
    private final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;

    /**
     * @param basecallsDir           Where to read basecalls from.
     * @param lane                   What lane to process.
     * @param readStructure          How to interpret each cluster.
     * @param barcodeRecordWriterMap Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                               one writer stored with key=null.
     * @param demultiplex            If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param maxReadsInRamPerTile   Configures number of reads each tile will store in RAM before spilling to disk.
     * @param tmpDirs                For SortingCollection spilling.
     * @param numProcessors          Controls number of threads.  If <= 0, the number of threads allocated is
     *                               available cores - numProcessors.
     * @param forceGc                Force explicit GC periodically.  This is good for causing memory maps to be released.
     * @param firstTile              (For debugging) If non-null, start processing at this tile.
     * @param tileLimit              (For debugging) If non-null, process no more than this many tiles.
     * @param outputRecordComparator For sorting output records within a single tile.
     * @param codecPrototype         For spilling output records to disk.
     * @param outputRecordClass      Inconveniently needed to create SortingCollections.
     * @param includeNonPfReads      If true, will include ALL reads (including those which do not have PF set)
     */
    public IlluminaBasecallsConverter(final File basecallsDir, final int lane, final ReadStructure readStructure,
                                      final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
                                      final boolean demultiplex,
                                      final int maxReadsInRamPerTile,
                                      final List<File> tmpDirs,
                                      final int numProcessors, final boolean forceGc,
                                      final Integer firstTile, final Integer tileLimit,
                                      final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
                                      final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
                                      final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
                                      final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                                      final boolean applyEamssFiltering,
                                      final boolean includeNonPfReads
    ) {
        this(basecallsDir, null, lane, readStructure,
                barcodeRecordWriterMap, demultiplex, maxReadsInRamPerTile,
                tmpDirs, numProcessors, forceGc, firstTile, tileLimit,
                outputRecordComparator, codecPrototype, outputRecordClass,
                bclQualityEvaluationStrategy, applyEamssFiltering,
                includeNonPfReads);
    }

    /**
     * @param basecallsDir           Where to read basecalls from.
     * @param barcodesDir            Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lane                   What lane to process.
     * @param readStructure          How to interpret each cluster.
     * @param barcodeRecordWriterMap Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                               one writer stored with key=null.
     * @param demultiplex            If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param maxReadsInRamPerTile   Configures number of reads each tile will store in RAM before spilling to disk.
     * @param tmpDirs                For SortingCollection spilling.
     * @param numProcessors          Controls number of threads.  If <= 0, the number of threads allocated is
     *                               available cores - numProcessors.
     * @param forceGc                Force explicit GC periodically.  This is good for causing memory maps to be released.
     * @param firstTile              (For debugging) If non-null, start processing at this tile.
     * @param tileLimit              (For debugging) If non-null, process no more than this many tiles.
     * @param outputRecordComparator For sorting output records within a single tile.
     * @param codecPrototype         For spilling output records to disk.
     * @param outputRecordClass      Inconveniently needed to create SortingCollections.
     * @param includeNonPfReads      If true, will include ALL reads (including those which do not have PF set)
     */
    public IlluminaBasecallsConverter(final File basecallsDir, File barcodesDir, final int lane,
                                      final ReadStructure readStructure,
                                      final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
                                      final boolean demultiplex,
                                      final int maxReadsInRamPerTile,
                                      final List<File> tmpDirs, final int numProcessors,
                                      final boolean forceGc, final Integer firstTile,
                                      final Integer tileLimit,
                                      final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
                                      final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
                                      final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
                                      final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                                      final boolean applyEamssFiltering, final boolean includeNonPfReads
    ) {
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        this.demultiplex = demultiplex;
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
        this.tmpDirs = tmpDirs;
        this.outputRecordComparator = outputRecordComparator;
        this.codecPrototype = codecPrototype;
        this.outputRecordClass = outputRecordClass;
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
        this.includeNonPfReads = includeNonPfReads;

        // If we're forcing garbage collection, collect every 5 minutes in a daemon thread.
        if (forceGc) {
            final Timer gcTimer = new Timer(true);
            final long delay = 5 * 1000 * 60;
            gcTimerTask = new TimerTask() {
                @Override
                public void run() {
                    log.info("Before explicit GC, Runtime.totalMemory()=" + getRuntime().totalMemory());
                    gc();
                    runFinalization();
                    log.info("After explicit GC, Runtime.totalMemory()=" + getRuntime().totalMemory());
                }
            };
            gcTimer.scheduleAtFixedRate(gcTimerTask, delay, delay);
        } else {
            gcTimerTask = null;
        }

        this.factory = new IlluminaDataProviderFactory(basecallsDir, barcodesDir, lane, readStructure,
                bclQualityEvaluationStrategy, getDataTypesFromReadStructure(readStructure, demultiplex));
        this.factory.setApplyEamssFiltering(applyEamssFiltering);

        if (numProcessors == 0) {
            this.numThreads = getRuntime().availableProcessors();
        } else if (numProcessors < 0) {
            this.numThreads = getRuntime().availableProcessors() + numProcessors;
        } else {
            this.numThreads = numProcessors;
        }
        this.tiles = new ArrayList<>(factory.getAvailableTiles());
        // Since the first non-fixed part of the read name is the tile number, without preceding zeroes,
        // and the output is sorted by read name, process the tiles in this order.
        sort(tiles, TILE_NUMBER_COMPARATOR);
        if (firstTile != null) {
            int i;
            for (i = 0; i < tiles.size(); ++i) {
                if (tiles.get(i).intValue() == firstTile.intValue()) {
                    tiles = tiles.subList(i, tiles.size());
                    break;
                }
            }
            if (tiles.get(0).intValue() != firstTile.intValue()) {
                throw new UserException("firstTile=" + firstTile + ", but that tile was not found.");
            }
        }
        if (tileLimit != null && tiles.size() > tileLimit) {
            tiles = tiles.subList(0, tileLimit);
        }

        this.numThreads = max(1, min(this.numThreads, tiles.size()));
    }

    /**
     * Must be called before doTileProcessing.  This is not passed in the ctor because often the
     * IlluminaDataProviderFactory is needed in order to construct the converter.
     *
     * @param converter Converts ClusterData to CLUSTER_OUTPUT_RECORD
     */
    public void setConverter(final ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter) {
        this.converter = converter;
    }

    /**
     * In case caller needs to get some info from factory.
     */
    public IlluminaDataProviderFactory getFactory() {
        return factory;
    }

    /**
     * Do the work, i.e. create a bunch of threads to read, sort and write.
     * setConverter() must be called before calling this method.
     */
    public void doTileProcessing() {
        try {
            // Generate the list of tiles that will be processed
            final List<Tile> tiles = new ArrayList<>();
            for (final Integer tileNumber : this.tiles) {
                tiles.add(new Tile(tileNumber));
            }

            final TileReadAggregator tileReadAggregator = new TileReadAggregator(tiles);
            tileReadAggregator.submit();
            try {
                tileReadAggregator.awaitWorkComplete();
            } catch (final InterruptedException e) {
                log.error(e, "Failure encountered in worker thread; attempting to shut down remaining worker threads and terminate ...");
                throw new GATKException("Failure encountered in worker thread; see log for details.");
            } finally {
                tileReadAggregator.shutdown();
            }

            for (final Map.Entry<Byte, Integer> entry : bclQualityEvaluationStrategy.getPoorQualityFrequencies().entrySet()) {
                log.warn(format("Observed low quality of %s %s times.", entry.getKey(), entry.getValue()));
            }
            bclQualityEvaluationStrategy.assertMinimumQualities();

        } finally {
            try {
                if (gcTimerTask != null) gcTimerTask.cancel();
            } catch (final Throwable ex) {
                log.warn(ex, "Ignoring exception stopping background GC thread.");
            }
            // Close the writers
            for (final Map.Entry<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> entry : barcodeRecordWriterMap.entrySet()) {
                final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = entry.getValue();
                log.debug(format("Closing file for barcode %s.", entry.getKey()));
                writer.close();
            }
        }
    }

    /**
     * Simple representation of a tile
     */
    private static class Tile implements Comparable<Tile> {
        private final int tileNumber;

        public Tile(final int i) {
            tileNumber = i;
        }

        public int getNumber() {
            return tileNumber;
        }

        @Override
        public boolean equals(final Object o) {
            return o instanceof Tile && this.getNumber() == ((Tile) o).getNumber();
        }

        @Override
        public int compareTo(final Tile o) {
            return TILE_NUMBER_COMPARATOR.compare(this.getNumber(), o.getNumber());
        }
    }


    /**
     * A Runnable that carries a priority which is used to compare and order other PriorityRunnables in a task queue.
     */
    private abstract class PriorityRunnable implements Runnable {
        private final int priority;

        /**
         * Create a new priority runnable with a default priority of 1.
         */
        public PriorityRunnable() {
            this(1);
        }

        public PriorityRunnable(final int priority) {
            this.priority = priority;
        }

        /**
         * Returns the priority level.  Higher priorities are run earlier.
         *
         * @return
         */
        int getPriority() {
            return this.priority;
        }
    }


    /**
     * Represents the state of a tile's processing and encapsulates the data collected from that tile.
     * <p>
     * TileProcessingRecords are accessed from each worker thread to assess the progress of the run, so its methods
     * are synchronized.
     */
    private class TileProcessingRecord {
        final private Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> barcodeToRecordCollection =
                new HashMap<>();
        final private Map<String, TileBarcodeProcessingState> barcodeToProcessingState = new HashMap<>();
        private TileProcessingState state = TileProcessingState.NOT_DONE_READING;
        private long recordCount = 0;

        /**
         * Returns the state of this tile's processing.
         */
        public synchronized TileProcessingState getState() {
            return this.state;
        }

        /**
         * Sets the state of this tile's processing.
         */
        public synchronized void setState(final TileProcessingState state) {
            this.state = state;
        }

        /**
         * Adds the provided record to this tile.
         */
        public synchronized void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
            this.recordCount += 1;

            // Grab the existing collection, or initialize it if it doesn't yet exist
            SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = this.barcodeToRecordCollection.get(barcode);
            if (recordCollection == null) {
                if (!barcodeRecordWriterMap.containsKey(barcode))
                    throw new UserException(format("Read records with barcode %s, but this barcode was not expected. " +
                            "(Is it referenced in the parameters file?)", barcode));
                recordCollection = this.newSortingCollection();
                this.barcodeToRecordCollection.put(barcode, recordCollection);
                this.barcodeToProcessingState.put(barcode, null);
            }
            recordCollection.add(record);
        }

        private synchronized SortingCollection<CLUSTER_OUTPUT_RECORD> newSortingCollection() {
            final int maxRecordsInRam =
                    maxReadsInRamPerTile /
                            barcodeRecordWriterMap.size();
            return newInstance(
                    outputRecordClass,
                    codecPrototype.clone(),
                    outputRecordComparator,
                    maxRecordsInRam,
                    tmpDirs);
        }

        /**
         * Returns the number of unique barcodes read.
         */
        public synchronized long getBarcodeCount() {
            return this.barcodeToRecordCollection.size();
        }

        /**
         * Returns the number of records read.
         */
        public synchronized long getRecordCount() {
            return recordCount;
        }

        /**
         * Returns the mapping of barcodes to records associated with them.
         */
        public synchronized Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> getBarcodeRecords() {
            return barcodeToRecordCollection;
        }

        /**
         * Gets the state of the provided barcode's data's processing progress.  Only invoke this query if this tile
         * is in a DONE_READING state.
         *
         * @throws IllegalStateException When a barcode is queried before the tile is in the DONE_READING state
         */
        public synchronized TileBarcodeProcessingState getBarcodeState(final String barcode) {
            if (this.getState() == TileProcessingState.NOT_DONE_READING) {
                throw new IllegalStateException(
                        "A tile's barcode data's state cannot be queried until the tile has been completely read.");
            }

            if (this.barcodeToProcessingState.containsKey(barcode)) {
                return this.barcodeToProcessingState.get(barcode);
            } else {
                return TileBarcodeProcessingState.NA;
            }
        }

        public synchronized Map<String, TileBarcodeProcessingState> getBarcodeProcessingStates() {
            return this.barcodeToProcessingState;
        }

        /**
         * Sets the processing state of the provided barcode in this record.
         *
         * @throws java.util.NoSuchElementException When the provided barcode is not one associated with this record.
         */
        public synchronized void setBarcodeState(final String barcode, final TileBarcodeProcessingState state) {
            if (this.barcodeToProcessingState.containsKey(barcode)) {
                this.barcodeToProcessingState.put(barcode, state);
            } else {
                throw new NoSuchElementException(format("No record of the provided barcode, %s.", barcode));
            }
        }

        /**
         * Returns the distinct set of barcodes for which data has been collected in this record.
         *
         * @return
         */
        public synchronized Set<String> getBarcodes() {
            return this.getBarcodeRecords().keySet();
        }
    }

    /**
     * Reads the information from a tile via an IlluminaDataProvider and feeds red information into a processingRecord
     * managed by the TileReadAggregator.
     */
    private class TileReader {
        private final Tile tile;
        private final TileReadAggregator handler;
        private final TileProcessingRecord processingRecord;

        public TileReader(final Tile tile, final TileReadAggregator handler, final TileProcessingRecord processingRecord) {
            this.tile = tile;
            this.handler = handler;
            this.processingRecord = processingRecord;
        }

        /**
         * Reads the data from the appropriate IlluminaDataProvider and feeds it into the TileProcessingRecord for
         * this tile.
         */
        public void process() {
            final IlluminaDataProvider dataProvider = factory.makeDataProvider(asList(this.tile.getNumber()));
            log.debug(format("Reading data from tile %s ...", tile.getNumber()));

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                // If this cluster is passing, or we do NOT want to ONLY emit passing reads, then add it to the next
                if (cluster.isPf() || includeNonPfReads) {
                    final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                    this.processingRecord.addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
                }
            }

            this.handler.completeTile(this.tile);
            dataProvider.close();
        }
    }


    /**
     * Aggregates data collected from tiles and writes them to file. Accepts records from TileReaders and maps
     * them to the appropriate BAM writers.
     */
    private class TileReadAggregator {
        /**
         * The collection of records associated with a particular tile.
         * <p>
         * Implemented as a TreeMap to guarantee tiles are iterated over in natural order.
         */
        private final Map<Tile, TileProcessingRecord> tileRecords = new TreeMap<>();

        /**
         * The executor responsible for doing work.
         * <p>
         * Implemented as a ThreadPoolExecutor with a PriorityBlockingQueue which orders submitted Runnables by their
         * priority.
         */
        private final ExecutorService prioritizingThreadPool = new ThreadPoolExecutor(
                numThreads,
                numThreads,
                0L,
                MILLISECONDS,
                new PriorityBlockingQueue<>(5, new Comparator<Runnable>() {
                    @Override
                    /**
                     * Compare the two Runnables, and assume they are PriorityRunnable; if not something strange is
                     * going on, so allow a ClassCastException be thrown.
                     */
                    public int compare(final Runnable o1, final Runnable o2) {
                        // Higher priority items go earlier in the queue, so reverse the "natural" comparison.
                        return ((PriorityRunnable) o2).getPriority() - ((PriorityRunnable) o1).getPriority();
                    }
                })
        );

        /**
         * The object acting as a latch to notify when the aggregator completes its work.
         */
        private final Object completionLatch = new Object();

        /**
         * Stores the thread that is executing this work so that it can be interrupted upon failure.
         */
        private Thread parentThread;
        private final Object workEnqueueMonitor = new Object();
        private final AtomicBoolean submitted = new AtomicBoolean(false);


        /**
         * Creates a TileReadAggregator that reads from the provided tiles.
         *
         * @param tiles
         */
        public TileReadAggregator(final Collection<Tile> tiles) {
            for (final Tile t : tiles) {
                tileRecords.put(t, new TileProcessingRecord());
            }
        }

        /**
         * Execute the tile aggregator's work.  Creates a thread pool to read data from tiles and write them to file.
         * Invoke this method only once.
         *
         * @throws IllegalStateException If submit was called more than once.
         */
        public void submit() {
            // Ensure the aggregator as not yet been submitted
            if (!this.submitted.compareAndSet(false, true)) {
                throw new IllegalStateException("The submit() method may not be called more than once.");
            }

            // Set the thread that is executing this work
            this.parentThread = currentThread();

            /**
             * For each tile, create and submit a tile processor.  Give it a negative execution priority (so that
             * prioritized tasks with a positive execution priority execute first), and give later tiles a lesser
             * (more negative) priority.
             */
            int priority = 0;
            for (final Tile tile : this.tileRecords.keySet()) {
                final TileReader reader = new TileReader(tile, this, this.tileRecords.get(tile));
                this.prioritizingThreadPool.execute(new PriorityRunnable(--priority) {
                    @Override
                    public void run() {
                        try {
                            reader.process();
                        } catch (final RuntimeException e) {
                            /**
                             * In the event of an internal failure, signal to the parent thread that something has gone
                             * wrong.  This is necessary because if an item of work fails to complete, the aggregator will
                             * will never reach its completed state, and it will never terminate.
                             */
                            parentThread.interrupt();
                            throw e;
                        } catch (final Error e) {
                            parentThread.interrupt();
                            throw e;
                        }
                    }
                });
            }
        }

        /**
         * Signals that a tile's processing is complete.  This must be invoked exactly once per tile, and only after
         * all of that tile has been processed.
         *
         * @throws IllegalStateException When the tile is already in the completed state.
         */
        private void completeTile(final Tile tile) {
            final TileProcessingRecord tileRecord = this.tileRecords.get(tile);

            if (tileRecord.getState() == TileProcessingState.DONE_READING) {
                throw new IllegalStateException("This tile is already in the completed state.");
            }

            // Update all of the barcodes and the tile to be marked as read
            for (final String barcode : tileRecord.getBarcodes()) {
                tileRecord.setBarcodeState(barcode, TileBarcodeProcessingState.READ);
                tileRecord.barcodeToRecordCollection.get(barcode).doneAdding();
            }
            tileRecord.setState(TileProcessingState.DONE_READING);

            log.debug(format("Completed reading tile %s; collected %s reads spanning %s barcodes.",
                    tile.getNumber(), tileRecord.getRecordCount(), tileRecord.getBarcodeCount()));

            //noinspection SynchronizationOnLocalVariableOrMethodParameter
            this.findAndEnqueueWorkOrSignalCompletion();
        }

        /**
         * Blocks until this aggregator completes its work.
         *
         * @throws InterruptedException
         */
        public void awaitWorkComplete() throws InterruptedException {
            synchronized (this.completionLatch) {
                this.completionLatch.wait();
            }
        }

        /**
         * Signals to any thread awaiting via awaitWorkComplete() that no work remains. Called
         * when this aggregator has reached its completed state.
         */
        private void signalWorkComplete() {
            synchronized (this.completionLatch) {
                this.completionLatch.notifyAll();
            }
        }

        /**
         * Poll the aggregator to find more tasks for it to enqueue.  Specifically, searches for un-written data
         * read from tiles for each barcode and enqueues it for writing.
         */
        private void findAndEnqueueWorkOrSignalCompletion() {
            synchronized (this.workEnqueueMonitor) {
                /**
                 * If there is work remaining to be done in this aggregator, walk through all of the barcodes and find
                 * tiles which have not yet written their barcode data but are in a state where they are able to.
                 */
                if (this.isWorkCompleted()) {
                    this.signalWorkComplete();
                } else {
                    final Queue<Runnable> tasks = new LinkedList<>();
                    for (final String barcode : barcodeRecordWriterMap.keySet()) {
                        NEXT_BARCODE:
                        for (final Map.Entry<Tile, TileProcessingRecord> entry : this.tileRecords.entrySet()) {
                            final Tile tile = entry.getKey();
                            final TileProcessingRecord tileRecord = entry.getValue();

                            /**
                             * If this tile has not been read, we cannot write this or later tiles' barcode data;
                             * move to the next barcode.
                             */
                            if (tileRecord.getState() != TileProcessingState.DONE_READING) {
                                break;
                            }
                            switch (tileRecord.getBarcodeState(barcode)) {
                                case NA:
                                case WRITTEN:
                                    /**
                                     * There is no data for this barcode for this tile, or it is already written; in
                                     * either scenario, this barcode will not be processed further for this tile, so
                                     * move onto the next tile as a possible candidate.
                                     */
                                    continue;
                                case QUEUED_FOR_WRITE:
                                    /**
                                     * The write for this barcode is in progress for this tile, so skip to the next
                                     * barcode.
                                     */
                                    break NEXT_BARCODE;
                                case READ:
                                    /**
                                     * This barcode has been read, and all of the earlier tiles have been written
                                     * for this barcode, so queue its writing.
                                     */
                                    tileRecord.setBarcodeState(barcode, TileBarcodeProcessingState.QUEUED_FOR_WRITE);
                                    log.debug(format("Enqueuing work for tile %s and barcode %s.", tile.getNumber(), barcode));
                                    tasks.add(this.newBarcodeWorkInstance(tile, tileRecord, barcode));
                                    break NEXT_BARCODE;
                            }
                        }
                    }

                    for (final Runnable task : tasks) {
                        this.prioritizingThreadPool.execute(task);
                    }
                }
            }
        }

        /**
         * Returns a PriorityRunnable that encapsulates the work involved with writing the provided tileRecord's data
         * for the given barcode to disk.
         *
         * @param tile       The tile from which the record was read
         * @param tileRecord The processing record associated with the tile
         * @param barcode    The barcode whose data within the tileRecord is to be written
         * @return The runnable that upon invocation writes the barcode's data from the tileRecord to disk
         */
        private PriorityRunnable newBarcodeWorkInstance(final Tile tile, final TileProcessingRecord tileRecord, final String barcode) {
            return new PriorityRunnable() {
                @Override
                public void run() {
                    try {
                        final SortingCollection<CLUSTER_OUTPUT_RECORD> records = tileRecord.getBarcodeRecords().get(barcode);
                        final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);

                        log.debug(format("Writing records from tile %s with barcode %s ...", tile.getNumber(), barcode));

                        final PeekIterator<CLUSTER_OUTPUT_RECORD> it = new PeekIterator<>(records.iterator());
                        while (it.hasNext()) {
                            final CLUSTER_OUTPUT_RECORD rec = it.next();

                            /**
                             * PIC-330 Sometimes there are two reads with the same cluster coordinates, and thus
                             * the same read name.  Discard both of them.  This code assumes that the two first of pairs
                             * will come before the two second of pairs, so it isn't necessary to look ahead a different
                             * distance for paired end.  It also assumes that for paired ends there will be duplicates
                             * for both ends, so there is no need to be PE-aware.
                             */
                            if (it.hasNext()) {
                                final CLUSTER_OUTPUT_RECORD lookAhead = it.peek();

                                /* TODO: Put this in SAMFileWriter wrapper
                                if (!rec.getReadUnmappedFlag() || !lookAhead.getReadUnmappedFlag()) {
                                    throw new IllegalStateException("Should not have mapped reads.");
                                }
                                */

                                if (outputRecordComparator.compare(rec, lookAhead) == 0) {
                                    it.next();
                                    log.info("Skipping reads with identical read names: " + rec.toString());
                                    continue;
                                }
                            }

                            writer.write(rec);
                            writeProgressLogger.record(null, 0);
                        }

                        tileRecord.setBarcodeState(barcode, TileBarcodeProcessingState.WRITTEN);
                        findAndEnqueueWorkOrSignalCompletion();

                    } catch (final RuntimeException e) {
                        /**
                         * In the event of an internal failure, signal to the parent thread that something has gone
                         * wrong.  This is necessary because if an item of work fails to complete, the aggregator will
                         * will never reach its completed state, and it will never terminate.
                         */
                        parentThread.interrupt();
                        throw e;
                    } catch (final Error e) {
                        parentThread.interrupt();
                        throw e;
                    }
                }

            };
        }

        /**
         * Returns true if this aggregator has not completed its work.  Specifically, returns false iff
         * any tile's barcode data yas not yet been written.
         *
         * @return True if more work remains to be done, false otherwise
         */
        public boolean isWorkCompleted() {
            for (final Map.Entry<Tile, TileProcessingRecord> entry : this.tileRecords.entrySet()) {
                final TileProcessingRecord tileProcessingRecord = entry.getValue();

                if (tileProcessingRecord.getState() != TileProcessingState.DONE_READING) {
                    log.debug(format("Work is not completed because a tile isn't done being read: %s.", entry.getKey().getNumber()));
                    return false;
                } else {
                    for (final Map.Entry<String, TileBarcodeProcessingState> barcodeStateEntry : tileProcessingRecord.getBarcodeProcessingStates().entrySet()) {
                        final TileBarcodeProcessingState barcodeProcessingState = barcodeStateEntry.getValue();
                        if (barcodeProcessingState != TileBarcodeProcessingState.WRITTEN) {
                            log.debug(format("Work is not completed because a tile isn't done being read: Tile %s, Barcode %s, Processing State %s.",
                                    entry.getKey().getNumber(), barcodeStateEntry.getKey(), barcodeProcessingState));
                            return false;
                        }
                    }
                }
            }
            log.info("All work is complete.");
            return true;
        }

        /**
         * Terminates the threads currently exiting in the thread pool abruptly via ThreadPoolExecutor.shutdownNow().
         */
        public void shutdown() {
            this.prioritizingThreadPool.shutdownNow();
        }
    }

    /**
     * Given a read structure return the data types that need to be parsed for this run
     */
    private static IlluminaDataType[] getDataTypesFromReadStructure(final ReadStructure readStructure,
                                                                    final boolean demultiplex) {
        if (readStructure.barcodes.isEmpty() || !demultiplex) {
            return DATA_TYPES_NO_BARCODE;
        } else {
            return DATA_TYPES_WITH_BARCODE;
        }
    }

    public static interface ClusterDataConverter<OUTPUT_RECORD> {

        /**
         * Creates the OUTPUT_RECORDs from the cluster
         */
        public OUTPUT_RECORD convertClusterToOutputRecord(final ClusterData cluster);
    }

    public static interface ConvertedClusterDataWriter<OUTPUT_RECORD> {
        void write(final OUTPUT_RECORD rec);

        void close();
    }
}
