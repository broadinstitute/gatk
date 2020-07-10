package org.broadinstitute.hellbender.engine;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.*;

import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.DelegatingIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.htsgetreader.*;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.SAMRecordToReadIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import javax.annotation.Nonnull;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Manages traversals and queries over sources of reads which are accessible via {@link GATKPath}s pointing to a file
 * behind an htsget server
 * (for now, BAM/CRAM files only).
 * <p>
 * Two basic operations are available:
 * <p>
 * -Iteration over all reads, optionally restricted to reads that overlap a set of intervals
 * -Targeted queries by one interval at a time
 */
public final class ReadsHtsgetDataSource implements ReadsDataSource {
    private static final Logger logger = LogManager.getLogger(ReadsHtsgetDataSource.class);

    private static final int numThreads = 8;

    /**
     * The sources provided to this data source
     */
    private final List<GATKPath> sources;

    /**
     * Only reads that overlap these intervals (and unmapped reads, if {@link #traverseUnmapped} is set) will be returned
     * during a full iteration. Null if iteration is unbounded.
     * <p>
     * Individual queries are unaffected by these intervals -- only traversals initiated via {@link #iterator} are affected.
     */
    private List<SimpleInterval> intervals;

    /**
     * If true, restrict traversals to unmapped reads (and reads overlapping any {@link #intervals}, if set).
     * False if iteration is unbounded or bounded only by our {@link #intervals}.
     * <p>
     * Note that this setting covers only unmapped reads that have no position -- unmapped reads that are assigned the
     * position of their mates will be returned by queries overlapping that position.
     * <p>
     * Individual queries are unaffected by this setting  -- only traversals initiated via {@link #iterator} are affected.
     */
    private boolean traverseUnmapped;

    /**
     * Mapping from the paths provided to this data source and their SAM headers
     */
    private final Map<GATKPath, SAMFileHeader> headers;

    /**
     * Hang onto all the iterators we've provided so we can close them when this source is closed
     */
    private final List<CloseableIterator<SAMRecord>> iterators;

    /**
     * The merged header for this data source
     */
    private final SAMFileHeader header;

    /**
     * SamReaderFactory used to create the SamReaders used by this data source
     */
    private final SamReaderFactory readerFactory;

    /**
     * Executor service used to spawn threads to initialize iterators in parallel
     */
    private final ExecutorService executorService;

    /**
     * Initialize this data source with a single SAM/BAM file and validation stringency SILENT.
     *
     * @param source SAM/BAM file, not null.
     */
    public ReadsHtsgetDataSource(final GATKPath source) {
        this(source != null ? Collections.singletonList(source) : null, null);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files and validation stringency SILENT.
     *
     * @param sources SAM/BAM files, not null.
     */
    public ReadsHtsgetDataSource(final List<GATKPath> sources) {
        this(sources, null);
    }

    /**
     * Initialize this data source with a single SAM/BAM file and a custom SamReaderFactory
     *
     * @param source                 path to SAM/BAM file, not null.
     * @param customSamReaderFactory SamReaderFactory to use, if null a default factory with no reference and validation
     *                               stringency SILENT is used.
     */
    public ReadsHtsgetDataSource(final GATKPath source, final SamReaderFactory customSamReaderFactory) {
        this(Collections.singletonList(source), customSamReaderFactory);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files and a custom SamReaderFactory
     *
     * @param sources                path to SAM/BAM file, not null.
     * @param customSamReaderFactory SamReaderFactory to use, if null a default factory with no reference and validation
     *                               stringency SILENT is used.
     */
    public ReadsHtsgetDataSource(final List<GATKPath> sources, final SamReaderFactory customSamReaderFactory) {
        Utils.nonNull(sources);
        Utils.nonEmpty(sources, "ReadsHtsgetDataSource cannot be created from empty source list");

        final String nonHtsgetSources = sources.stream()
            .filter(source -> !source.getScheme().equals(GATKPath.HTSGET_SCHEME))
            .map(GATKPath::toString)
            .collect(Collectors.joining(", "));

        if (!nonHtsgetSources.isEmpty()) {
            throw new UserException("This source can only be instantiated from htsget paths: " + nonHtsgetSources);
        }

        this.readerFactory =
            customSamReaderFactory == null
                ? SamReaderFactory.makeDefault().validationStringency(ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY)
                : customSamReaderFactory;

        this.sources = sources;
        this.headers = new LinkedHashMap<>(sources.size());
        this.iterators = new ArrayList<>();

        this.executorService = Executors.newFixedThreadPool(ReadsHtsgetDataSource.numThreads, new ThreadFactoryBuilder()
            .setNameFormat("htsget-reader-thread-%d")
            .setDaemon(true)
            .build());

        final List<Future<Map.Entry<GATKPath, SAMFileHeader>>> futures = new ArrayList<>(sources.size());

        logger.info("Downloading headers from htsget.");

        for (final GATKPath source : sources) {
            futures.add(this.executorService.submit(() -> {
                try {
                    // Ensure each source can be obtained from htsget server
                    final HtsgetRequest req = new HtsgetRequest(source).withDataClass(HtsgetClass.header);
                    // Request only the headers and use them to construct SAMFileHeader for each source
                    try (final InputStream headerStream = req.getResponse().getDataStream()) {
                        final SAMFileHeader header = readerFactory.open(SamInputResource.of(headerStream)).getFileHeader();
                        return new AbstractMap.SimpleImmutableEntry<>(source, header);
                    }
                } catch (final UserException e) {
                    throw new UserException("Failed to load header from htsget source " + source.toString(), e);
                }
            }));
        }

        try {
            for (final Future<Map.Entry<GATKPath, SAMFileHeader>> future : futures) {
                final Map.Entry<GATKPath, SAMFileHeader> entry = future.get();
                this.headers.put(entry.getKey(), entry.getValue());
            }
        } catch (final ExecutionException | InterruptedException e) {
            throw new UserException("Interrupted while initializing iterator", e);
        }

        logger.info("Finished loading headers from htsget.");

        if (sources.size() > 1) {
            this.header = createHeaderMerger(this.headers.values()).getMergedHeader();
        } else {
            this.header = this.headers.values().iterator().next();
        }

        if (this.header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            logger.warn("Files not in coordinate sorted order");
        }
    }

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads that overlap the given intervals,
     * and to unmapped reads if specified.
     * <p>
     * Calls to {@link #query} are not affected by this method.
     *
     * @param intervals        Our next full traversal will return reads overlapping these intervals
     * @param traverseUnmapped Our next full traversal will return unmapped reads (this affects only unmapped reads that
     *                         have no position -- unmapped reads that have the position of their mapped mates will be
     *                         included if the interval overlapping that position is included).
     */
    @Override
    public void setTraversalBounds(final List<SimpleInterval> intervals, final boolean traverseUnmapped) {
        this.intervals = intervals != null && !intervals.isEmpty() ? intervals : null;
        this.traverseUnmapped = traverseUnmapped;
    }

    /**
     * @return True if traversals initiated via {@link #iterator} will be restricted to reads that overlap intervals
     * as configured via {@link #setTraversalBounds}, otherwise false
     */
    @Override
    public boolean traversalIsBounded() {
        return this.intervals != null || this.traverseUnmapped;
    }

    /**
     * This data source can be queried even without index files
     *
     * @return always true
     */
    @Override
    public boolean isQueryableByInterval() {
        return true;
    }

    /**
     * Iterate over all reads in this data source. If intervals were provided via {@link #setTraversalBounds},
     * iteration is limited to reads that overlap that set of intervals.
     *
     * @return An iterator over the reads in this data source, limited to reads that overlap the intervals supplied
     * via {@link #setTraversalBounds} (if intervals were provided)
     */
    @Override
    @Nonnull
    public Iterator<GATKRead> iterator() {
        return this.prepareIteratorsForTraversal(this.intervals, this.traverseUnmapped);
    }

    /**
     * Query reads over a specific interval. This operation is not affected by prior calls to
     * {@link #setTraversalBounds}
     *
     * @param interval The interval over which to query
     * @return Iterator over reads overlapping the query interval
     */
    @Override
    public Iterator<GATKRead> query(final SimpleInterval interval) {
        return this.prepareIteratorsForTraversal(Collections.singletonList(interval), false);
    }

    /**
     * @return An iterator over just the unmapped reads with no assigned position. This operation is not affected
     * by prior calls to {@link #setTraversalBounds}.
     */
    @Override
    public Iterator<GATKRead> queryUnmapped() {
        return this.prepareIteratorsForTraversal(null, true);
    }

    private Iterator<GATKRead> prepareIteratorsForTraversal(final List<SimpleInterval> intervals, final boolean traverseUnmapped) {
        if (this.executorService.isShutdown()) {
            throw new UserException("This data source has already been shut down.");
        }

        final Map<SamReader, CloseableIterator<SAMRecord>> mapping = intervals == null
            ? this.getReaders(traverseUnmapped)
            : this.getReadersWithIntervals(intervals, traverseUnmapped);

        if (mapping.isEmpty()) {
            return Collections.emptyIterator();
        }
        final CloseableIterator<SAMRecord> samRecordIterator;
        if (mapping.size() == 1) {
            samRecordIterator = mapping.values().iterator().next();
        } else {
            final Set<SAMFileHeader> headers = mapping.keySet().stream()
                .map(SamReader::getFileHeader)
                .collect(Collectors.toCollection(LinkedHashSet::new));
            samRecordIterator = new MergingSamRecordIterator(createHeaderMerger(headers), mapping, true);
        }
        this.iterators.add(samRecordIterator);
        return new SAMRecordToReadIterator(samRecordIterator);
    }

    private Map<SamReader, CloseableIterator<SAMRecord>> getReaders(final boolean onlyUnplacedUnmapped) {
        final Map<SamReader, CloseableIterator<SAMRecord>> mapping = new LinkedHashMap<>(this.sources.size());
        final List<Future<Map.Entry<SamReader, CloseableIterator<SAMRecord>>>> futures = new ArrayList<>(this.sources.size());
        this.sources.forEach(source -> futures.add(this.executorService.submit(() -> {
            final HtsgetRequest req = new HtsgetRequest(source);
            if (onlyUnplacedUnmapped) {
                req.setInterval(HtsgetRequest.UNMAPPED_UNPLACED_INTERVAL);
            }
            final SamReader reader = this.readerFactory.open(SamInputResource.of(req.getResponse().getDataStream()));
            return new AbstractMap.SimpleImmutableEntry<>(reader, wrapIteratorWithClose(reader.iterator(), reader));
        })));
        try {
            for (final Future<Map.Entry<SamReader, CloseableIterator<SAMRecord>>> future : futures) {
                final Map.Entry<SamReader, CloseableIterator<SAMRecord>> entry = future.get();
                mapping.put(entry.getKey(), entry.getValue());
            }
        } catch (final ExecutionException | InterruptedException e) {
            throw new UserException("Interrupted while initializing iterator", e);
        }
        return mapping;
    }

    private Map<SamReader, CloseableIterator<SAMRecord>> getReadersWithIntervals(final List<SimpleInterval> intervals, final boolean includeUnplacedUnmapped) {
        final Map<SamReader, CloseableIterator<SAMRecord>> mapping = new LinkedHashMap<>(this.sources.size() * intervals.size());
        SimpleInterval prevInterval = null;
        final List<Future<Map.Entry<SamReader, CloseableIterator<SAMRecord>>>> futures = new ArrayList<>(this.sources.size() * intervals.size());
        for (final SimpleInterval interval : intervals) {
            final SimpleInterval finalPrevInterval = prevInterval;
            this.sources.parallelStream()
                .filter(source -> IntervalUtils.intervalIsOnDictionaryContig(
                    interval,
                    this.headers.get(source).getSequenceDictionary()))
                .forEach(source -> futures.add(this.executorService.submit(() -> {
                    final HtsgetRequest req = new HtsgetRequest(source).withInterval(interval);
                    final SamReader reader = this.readerFactory.open(SamInputResource.of(req.getResponse().getDataStream()));
                    return new AbstractMap.SimpleEntry<>(reader, getIterWithInterval(reader, interval, finalPrevInterval));
                })));
            prevInterval = interval;
        }
        try {
            for (final Future<Map.Entry<SamReader, CloseableIterator<SAMRecord>>> future : futures) {
                final Map.Entry<SamReader, CloseableIterator<SAMRecord>> entry = future.get();
                mapping.put(entry.getKey(), entry.getValue());
            }
        } catch (final ExecutionException | InterruptedException e) {
            throw new UserException("Interrupted while initializing iterator", e);
        }
        if (includeUnplacedUnmapped) {
            mapping.putAll(this.getReaders(true));
        }
        return mapping;
    }

    @Override
    public boolean supportsSerialIteration() {
        return true;
    }

    /**
     * Shut down this data source permanently, closing all iterations and readers.
     */
    @Override
    public void close() {
        this.iterators.forEach(CloseableIterator::close);
        this.executorService.shutdownNow();
    }

    /**
     * Returns the SAM header for this data source. Will be a merged header if there are multiple readers.
     * If there is only a single reader, returns its header directly.
     *
     * @return SAM header for this data source
     */
    @Override
    public SAMFileHeader getHeader() {
        return this.header;
    }

    private static CloseableIterator<SAMRecord> getIterWithInterval(final SamReader samReader, final SimpleInterval currInterval, final SimpleInterval prevInterval) {
        /*
        To remove reads duplicated across two subsequent intervals, we take any read which is in
        the current interval but NOT in the previous interval, unless the current interval is the first,
        in which case we take any read overlapping it
        */
        final CloseableIterator<SAMRecord> filteredSamRecords = new FilteringSamIterator(samReader.iterator(), new SamRecordFilter() {
            @Override
            public boolean filterOut(final SAMRecord record) {
                return !currInterval.overlaps(record) || (prevInterval != null && prevInterval.overlaps(record));
            }

            @Override
            public boolean filterOut(final SAMRecord first, final SAMRecord second) {
                throw new UnsupportedOperationException();
            }
        });
        return wrapIteratorWithClose(filteredSamRecords, samReader);
    }

    /**
     * Wrap an iterator to allow us to free the backing SamReader without holding onto an explicit reference to it
     *
     * @param iterator  the iterator to wrap
     * @param samReader the SamReader to close once the iterator has been used up
     * @return a wrapped CloseableIterator
     */
    private static CloseableIterator<SAMRecord> wrapIteratorWithClose(final CloseableIterator<SAMRecord> iterator, final SamReader samReader) {
        return new DelegatingIterator<SAMRecord>(iterator) {
            @Override
            public void close() {
                try {
                    iterator.close();
                    samReader.close();
                } catch (final IOException e) {
                    throw new GATKException("Error closing SAMReader", e);
                }
            }
        };
    }

    // TODO: Push this and following method down into htsjdk
    private static SamFileHeaderMerger createHeaderMerger(final Collection<SAMFileHeader> headers) {
        return new SamFileHeaderMerger(identifySortOrder(headers), headers, true);
    }

    private static SAMFileHeader.SortOrder identifySortOrder(final Collection<SAMFileHeader> headers) {
        final Set<SAMFileHeader.SortOrder> sortOrders = headers.stream()
            .map(SAMFileHeader::getSortOrder)
            .collect(Collectors.toSet());
        final SAMFileHeader.SortOrder order;
        if (sortOrders.size() == 1) {
            order = sortOrders.iterator().next();
        } else {
            order = SAMFileHeader.SortOrder.unsorted;
            logger.warn("Inputs have different sort orders. Assuming {} sorted reads for all of them.", order);
        }
        return order;
    }
}
