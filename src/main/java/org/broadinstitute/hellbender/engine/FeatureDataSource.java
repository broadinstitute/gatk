package org.broadinstitute.hellbender.engine;

import com.intel.genomicsdb.model.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.reader.GenomicsDBFeatureReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.*;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBConstants;
import org.broadinstitute.hellbender.utils.IndexUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.SeekableByteChannelPrefetcher;

import java.io.File;
import java.io.IOException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

/**
 * Enables traversals and queries over sources of Features, which are metadata associated with a location
 * on the genome in a format supported by our file parsing framework, Tribble. Examples of Features are
 * VCF records and hapmap records.
 *
 * Two basic operations are available on this data source:
 *
 * -Iteration over all Features in this data source, optionally restricted to Features overlapping
 *  a set of intervals if intervals are provided via {@link #setIntervalsForTraversal(List)}. Traversal
 *  by a set of intervals requires the file to have been indexed using the bundled tool IndexFeatureFile.
 *  The set of intervals provided MUST be non-overlapping and sorted in increasing order of start position.
 *
 * -Targeted queries by one interval at a time. This also requires the file to have been indexed using
 *  the bundled tool IndexFeatureFile. Targeted queries by one interval at a time are unaffected by
 *  any intervals for full traversal set via {@link #setIntervalsForTraversal(List)}.
 *
 * To improve performance in the case of targeted queries by one interval at a time, this class uses a caching
 * scheme that is optimized for the common access pattern of multiple separate queries over intervals with
 * gradually increasing start positions. It optimizes for this use case by pre-fetching records immediately
 * following each interval during a query and caching them. Performance will suffer if the access pattern is
 * random, involves queries over intervals with DECREASING start positions instead of INCREASING start positions,
 * or involves lots of very large jumps forward on the genome or lots of contig switches. Query caching
 * can be disabled, if desired.
 *
 * @param <T> The type of Feature returned by this data source
 */
public final class FeatureDataSource<T extends Feature> implements GATKDataSource<T>, AutoCloseable {
    private static final Logger logger = LogManager.getLogger(FeatureDataSource.class);

    /**
     * identifies a path as a GenomicsDB URI
     */
    public static final String GENOMIC_DB_URI_SCHEME = "gendb://";

    /**
     * Feature reader used to retrieve records from our file
     */
    private final FeatureReader<T> featureReader;

    /**
     * Iterator representing an open traversal over this data source initiated via a call to {@link #iterator}
     * (null if there is no open traversal). We need this to ensure that each iterator is properly closed,
     * and to enforce the constraint (required by Tribble) that we never have more than one iterator open
     * over our feature reader.
     */
    private CloseableTribbleIterator<T> currentIterator;

    /**
     * Our intervals for traversal. If set, restricts full traversals initiated via {@link #iterator} to
     * return only Features overlapping this set of intervals. Does not affect individual queries
     * initiated via {@link #query(SimpleInterval)} and/or {@link #queryAndPrefetch(SimpleInterval)}.
     */
    private List<SimpleInterval> intervalsForTraversal;

    /**
     * Cache containing Features from recent queries initiated via {@link #query(SimpleInterval)} and/or
     * {@link #queryAndPrefetch(SimpleInterval)}. This is guaranteed to start at the start position of the
     * most recent query, but will typically end well after the end of the most recent query. Designed to
     * improve performance of the common access pattern involving multiple queries across nearby intervals
     * with gradually increasing start positions.
     */
    private final FeatureCache<T> queryCache;

    /**
     * When we experience a cache miss (ie., a query interval not fully contained within our cache) and need
     * to re-populate the Feature cache from disk to satisfy a query, this controls the number of extra bases
     * AFTER the end of our interval to fetch. Should be sufficiently large so that typically a significant number
     * of subsequent queries will be cache hits (ie., query intervals fully contained within our cache) before
     * we have another cache miss and need to go to disk again.
     */
    private final int queryLookaheadBases;

    /**
     * Holds information about the path this datasource reads from.
     */
    private final FeatureInput<T> featureInput;

    /**
     * True if this datasource is backed by a file that has an associated index file, false if it doesn't
     */
    private final boolean hasIndex;

    /**
     * True if this datasource supports efficient random access queries.
     *
     * For a file, this is the same as {@link #hasIndex}, but there are non-file data sources (eg., GenomicsDB)
     * that don't have a separate index file but do support random access.
     */
    private final boolean supportsRandomAccess;

    /**
     * Default value for queryLookaheadBases, if none is specified. This is designed to be large enough
     * so that in typical usage (ie., query intervals with gradually increasing start locations) there will
     * be a substantial number of cache hits between cache misses, reducing the number of times we need to
     * repopulate the cache from disk.
     */
    public static final int DEFAULT_QUERY_LOOKAHEAD_BASES = 1000;

    /**
     * Creates a FeatureDataSource backed by the provided File. The data source will have an automatically
     * generated name, and will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES})
     * during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     */
    public FeatureDataSource(final File featureFile) {
        this(featureFile, null);
    }

    /**
     * Creates a FeatureDataSource backed by the provided path. The data source will have an automatically
     * generated name, and will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES})
     * during queries that produce cache misses.
     *
     * @param featurePath path or URI to source of Features
     */
    public FeatureDataSource(final String featurePath) {
        this(featurePath, null, DEFAULT_QUERY_LOOKAHEAD_BASES, null);
    }

    /**
     * Creates a FeatureDataSource backed by the provided File and assigns this data source the specified logical
     * name. We will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES}) during queries
     * that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param name logical name for this data source (may be null)
     */
    public FeatureDataSource(final File featureFile, final String name) {
        this(featureFile, name, DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    /**
     * Creates a FeatureDataSource backed by the provided File and assigns this data source the specified logical
     * name. We will look ahead the specified number of bases during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param name logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     */
    public FeatureDataSource(final File featureFile, final String name, final int queryLookaheadBases){
        this(Utils.nonNull(featureFile).getAbsolutePath(), name, queryLookaheadBases, null);
    }

    /**
     * Creates a FeatureDataSource backed by the resource at the provided path.
     *
     * @param featurePath path to file or GenomicsDB url containing features
     * @param name logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     */
    public FeatureDataSource(final String featurePath, final String name, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType) {
        this(new FeatureInput<>(featurePath, name != null ? name : featurePath), queryLookaheadBases, targetFeatureType);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput a FeatureInput specifying a source of Features
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType) {
        this(featureInput, queryLookaheadBases, targetFeatureType, 0, 0);
    }

    /**
     * Creates a FeatureDataSource backed by the resource at the provided path.
     *
     * @param featurePath path to file or GenomicsDB url containing features
     * @param name logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer  MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     */
    public FeatureDataSource(final String featurePath, final String name, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer ) {
        this(new FeatureInput<>(featurePath, name != null ? name : featurePath), queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput a FeatureInput specifying a source of Features
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer  MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer) {
        this(featureInput, queryLookaheadBases, targetFeatureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
             null);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput a FeatureInput specifying a source of Features
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     * @param cloudPrefetchBuffer  MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param reference Path to a reference. May be null. Needed only for reading from GenomicsDB.
     */
    public FeatureDataSource(final FeatureInput<T> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType,
                             final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final Path reference) {
        Utils.validateArg( queryLookaheadBases >= 0, "Query lookahead bases must be >= 0");
        this.featureInput = Utils.nonNull(featureInput, "featureInput must not be null");

        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = (cloudPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudPrefetchBuffer, is) : Function.identity());
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = (cloudIndexPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudIndexPrefetchBuffer, is) : Function.identity());

        // Create a feature reader without requiring an index.  We will require one ourselves as soon as
        // a query by interval is attempted.
        this.featureReader = getFeatureReader(featureInput, targetFeatureType, cloudWrapper, cloudIndexWrapper, reference);

        if (isGenomicsDBPath(featureInput.getFeaturePath())) {
            //genomics db uri's have no associated index file to read from, but they do support random access
            this.hasIndex = false;
            this.supportsRandomAccess = true;
        } else if (featureReader instanceof AbstractFeatureReader) {
            this.hasIndex = ((AbstractFeatureReader<T, ?>) featureReader).hasIndex();
            this.supportsRandomAccess = hasIndex;
        } else {
            throw new GATKException("Found a feature input that was neither GenomicsDB or a Tribble AbstractFeatureReader.  Input was " + featureInput.toString() + ".");
        }
        // Due to a bug in HTSJDK, unindexed block compressed input files may fail to parse completely. For safety,
        // these files have been disabled. See https://github.com/broadinstitute/gatk/issues/4224 for discussion
        if (!hasIndex && IOUtil.hasBlockCompressedExtension(featureInput.getFeaturePath())) {
            throw new UserException.MissingIndex(featureInput.toString(),"Support for unindexed block-compressed files has been temporarily disabled. Try running IndexFeatureFile on the input.");
        }

        this.currentIterator = null;
        this.intervalsForTraversal = null;
        this.queryCache = new FeatureCache<>();
        this.queryLookaheadBases = queryLookaheadBases;
    }

    /**
     * @param path String containing the path to test
     * @return true if path represent a GenomicsDB URI, otherwise false
     */
    public static boolean isGenomicsDBPath(final String path) {
        return path != null && path.startsWith(GENOMIC_DB_URI_SCHEME);
    }

    @SuppressWarnings("unchecked")
    private static <T extends Feature> FeatureReader<T> getFeatureReader(final FeatureInput<T> featureInput, final Class<? extends Feature> targetFeatureType,
                                                                         final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper,
                                                                         final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper,
                                                                         final Path reference) {
        if (isGenomicsDBPath(featureInput.getFeaturePath())) {
            try {
                if (reference == null) {
                    throw new UserException.MissingReference("You must provide a reference if you want to load from GenomicsDB");
                }
                try {
                    final File referenceAsFile = reference.toFile();
                    return (FeatureReader<T>)getGenomicsDBFeatureReader(featureInput.getFeaturePath(), referenceAsFile);
                } catch (final UnsupportedOperationException e){
                    throw new UserException.BadInput("GenomicsDB requires that the reference be a local file.", e);
                }
            } catch (final ClassCastException e) {
                throw new UserException("GenomicsDB inputs can only be used to provide VariantContexts.", e);
            }
        } else {
            FeatureCodec<T, ?> codec = getCodecForFeatureInput(featureInput, targetFeatureType);
            return getTribbleFeatureReader(featureInput, codec, cloudWrapper, cloudIndexWrapper);
        }
    }

    /**
     * Get a new FeatureCodec instance to use for a FeatureInput. Avoid re-discovering which codec class to
     * use by checking to see if the FeatureInput already has a cached codec class. It not, discover the codec class
     * and cache it for next time.
     *
     * @return A new FeatureCodec instance to use for the FeatureInput.
     */
    @SuppressWarnings("unchecked")
    private static <T extends Feature> FeatureCodec<T, ?> getCodecForFeatureInput(final FeatureInput<T> featureInput,
                                                                                  final Class<? extends Feature> targetFeatureType) {
        FeatureCodec<T, ?> codec;
        final Class<FeatureCodec<T, ?>> codecClass = featureInput.getFeatureCodecClass();
        if (codecClass == null) {
            final Path featurePath = IOUtils.getPath(featureInput.getFeaturePath());
            IOUtils.assertFileIsReadable(featurePath);
            codec = (FeatureCodec<T, ?>) FeatureManager.getCodecForFile(featurePath, targetFeatureType);
            featureInput.setFeatureCodecClass((Class<FeatureCodec<T, ?>>)codec.getClass());
        } else {
            try {
                codec = codecClass.newInstance();
            }
            catch ( InstantiationException | IllegalAccessException e ) {
                throw new GATKException("Unable to automatically instantiate codec " + codecClass.getName());
            }
        }
        return codec;
    }

    private static <T extends Feature> AbstractFeatureReader<T, ?> getTribbleFeatureReader(final FeatureInput<T> featureInput, final FeatureCodec<T, ?> codec, Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper, Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper) {
        Utils.nonNull(codec);
        try {
            final String absolutePath = IOUtils.getPath(featureInput.getFeaturePath()).toAbsolutePath().toUri().toString();

            // Instruct the reader factory to not require an index. We will require one ourselves as soon as
            // a query by interval is attempted.
            final boolean requireIndex = false;

            // Only apply the wrappers if the feature input is on Google Cloud Storage
            if ( BucketUtils.isCloudStorageUrl(absolutePath) ) {
                return AbstractFeatureReader.getFeatureReader(absolutePath, null, codec, requireIndex, cloudWrapper, cloudIndexWrapper);
            } else {
                return AbstractFeatureReader.getFeatureReader(absolutePath, null, codec, requireIndex, Function.identity(), Function.identity());
            }
        }
        catch ( final TribbleException e ) {
            throw new GATKException("Error initializing feature reader for path " +  featureInput.getFeaturePath(), e);
        }
    }

    private static FeatureReader<VariantContext> getGenomicsDBFeatureReader(final String path, final File reference) {
        if( !isGenomicsDBPath(path) ) {
            throw new IllegalArgumentException("Trying to create a GenomicsDBReader from a non-GenomicsDB input");
        }

        final String noheader = path.replace(GENOMIC_DB_URI_SCHEME, "");
        final File workspace = new File(noheader);
        final File callsetJson = new File(noheader, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME);
        final File vidmapJson = new File(noheader, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME);
        final File vcfHeader = new File(noheader, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME);

        if ( ! workspace.exists() || ! workspace.canRead() || ! workspace.isDirectory() ) {
            throw new UserException("GenomicsDB workspace " + workspace.getAbsolutePath() + " does not exist, " +
                    " is not readable, or is not a directory");
        }

        try {
            IOUtils.canReadFile(callsetJson);
            IOUtils.canReadFile(vidmapJson);
            IOUtils.canReadFile(vcfHeader);
        } catch ( UserException.CouldNotReadInputFile e ) {
            throw new UserException("Couldn't connect to GenomicsDB because the vidmap, callset JSON files, or gVCF Header (" +
                    GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME + "," + GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME + "," +
                    GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME + ") could not be read from GenomicsDB workspace " + workspace.getAbsolutePath(), e);
        }

        final GenomicsDBExportConfiguration.ExportConfiguration exportConfigurationBuilder =
                createExportConfiguration(reference, workspace, callsetJson, vidmapJson, vcfHeader);

        try {
            return new GenomicsDBFeatureReader<>(exportConfigurationBuilder, new BCF2Codec(), Optional.empty());
        } catch (final IOException e) {
            throw new UserException("Couldn't create GenomicsDBFeatureReader", e);
        }
    }

    private static GenomicsDBExportConfiguration.ExportConfiguration createExportConfiguration(final File reference, final File workspace,
                                                                                               final File callsetJson, final File vidmapJson,
                                                                                               final File vcfHeader) {
        GenomicsDBExportConfiguration.ExportConfiguration.Builder exportConfigurationBuilder =
                GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                        .setWorkspace(workspace.getAbsolutePath())
                        .setReferenceGenome(reference.getAbsolutePath())
                        .setVidMappingFile(vidmapJson.getAbsolutePath())
                        .setCallsetMappingFile(callsetJson.getAbsolutePath())
                        .setVcfHeaderFilename(vcfHeader.getAbsolutePath())
                        .setProduceGTField(false)
                        .setProduceGTWithMinPLValueForSpanningDeletions(false)
                        .setSitesOnlyQuery(false)
                        .setMaxDiploidAltAllelesThatCanBeGenotyped(GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
        Path arrayFolder = Paths.get(workspace.getAbsolutePath(), GenomicsDBConstants.DEFAULT_ARRAY_NAME).toAbsolutePath();

        // For the multi-interval support, we create multiple arrays (directories) in a single workspace -
        // one per interval. So, if you wish to import intervals ("chr1", [ 1, 100M ]) and ("chr2", [ 1, 100M ]),
        // you end up with 2 directories named chr1$1$100M and chr2$1$100M. So, the array names depend on the
        // partition bounds.

        // During the read phase, the user only supplies the workspace. The array names are obtained by scanning
        // the entries in the workspace and reading the right arrays. For example, if you wish to read ("chr2",
        // 50, 50M), then only the second array is queried.

        // In the previous version of the tool, the array name was a constant - genomicsdb_array. The new version
        // will be backward compatible with respect to reads. Hence, if a directory named genomicsdb_array is found,
        // the array name is passed to the GenomicsDBFeatureReader otherwise the array names are generated from the
        // directory entries.
        if (Files.exists(arrayFolder)) {
            exportConfigurationBuilder.setArrayName(GenomicsDBConstants.DEFAULT_ARRAY_NAME);
        } else {
            exportConfigurationBuilder.setGenerateArrayNameFromPartitionBounds(true);
        }

        return exportConfigurationBuilder.build();
    }

    /**
     * Returns the sequence dictionary for this source of Features.
     * Uses the dictionary from the VCF header (if present) for variant inputs,
     * otherwise attempts to create a sequence dictionary from the index file (if present).
     * Returns null if no dictionary could be created from either the header or the index.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        SAMSequenceDictionary dict = null;
        final Object header = getHeader();
        if (header instanceof VCFHeader) {
            dict = ((VCFHeader) header).getSequenceDictionary();
        }
        if (dict != null && !dict.isEmpty()) {
            return dict;
        }
        if (hasIndex) {
            return IndexUtils.createSequenceDictionaryFromFeatureIndex(new File(featureInput.getFeaturePath()));
        }
        return null;
    }

    /**
     * Restricts traversals of this data source via {@link #iterator} to only return Features that overlap the provided
     * intervals. Calls to {@link #query(SimpleInterval)} and/or {@link #queryAndPrefetch(SimpleInterval)} are not
     * affected by these intervals.
     *
     * Intervals MUST be non-overlapping and sorted in order of increasing start position, otherwise traversal
     * results will be incorrect.
     *
     * Passing in a null or empty interval List clears the intervals for traversal, making future iterations
     * over this data source unrestricted by intervals.
     *
     * @param intervals Our next full traversal will return only Features overlapping these intervals
     */
    public void setIntervalsForTraversal( final List<SimpleInterval> intervals ) {
        // Treat null and empty interval lists the same
        intervalsForTraversal = (intervals != null && !intervals.isEmpty()) ? intervals : null;

        if ( intervalsForTraversal != null && ! supportsRandomAccess ) {
            throw new UserException("Input " + featureInput.getFeaturePath() + " must support random access to enable traversal by intervals. " +
                                    "If it's a file, please index it using the bundled tool " + IndexFeatureFile.class.getSimpleName());
        }
    }


    /**
     * Gets an iterator over all Features in this data source, restricting traversal to Features
     * overlapping our intervals if intervals were provided via {@link #setIntervalsForTraversal(List)}
     *
     * Calling this method invalidates (closes) any previous iterator obtained from this method.
     *
     * @return an iterator over all Features in this data source, limited to Features that overlap the intervals supplied via {@link #setIntervalsForTraversal(List)} (if intervals were provided)
     */
    @Override
    public Iterator<T> iterator() {
        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        try {
            // Save the iterator returned so that we can close it properly later
            currentIterator = intervalsForTraversal != null ? new FeatureIntervalIterator<T>(intervalsForTraversal, featureReader, featureInput.getFeaturePath())
                                                            : featureReader.iterator();
            return currentIterator;
        }
        catch ( IOException e ) {
            throw new GATKException("Error creating iterator over file " + featureInput.getFeaturePath(), e);
        }
    }

    /**
     * Gets an iterator over all Features in this data source that overlap the provided interval.
     *
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     *
     * Requires the backing file to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     *
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     *
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Features overlapping this interval
     * @return an iterator over all Features in this data source that overlap the provided interval
     */
    @Override
    public Iterator<T> query( final SimpleInterval interval ) {
        return queryAndPrefetch(interval).iterator();
    }

    /**
     * Returns a List of all Features in this data source that overlap the provided interval.
     *
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     *
     * Requires the backing file to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     *
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     *
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Features overlapping this interval
     * @return a List of all Features in this data source that overlap the provided interval
     */
    public List<T> queryAndPrefetch( final SimpleInterval interval ) {
        if ( ! supportsRandomAccess ) {
            throw new UserException("Input " + featureInput.getFeaturePath() + " must support random access to enable queries by interval. " +
                                    "If it's a file, please index it using the bundled tool " + IndexFeatureFile.class.getSimpleName());
        }

        // If the query can be satisfied using existing cache contents, prepare for retrieval
        // by discarding all Features at the beginning of the cache that end before the start
        // of our query interval.
        if ( queryCache.cacheHit(interval) ) {
            queryCache.trimToNewStartPosition(interval.getStart());
        }
        // Otherwise, we have a cache miss, so go to disk to refill our cache.
        else {
            refillQueryCache(interval);
        }

        // Return the subset of our cache that overlaps our query interval
        return queryCache.getCachedFeaturesUpToStopPosition(interval.getEnd());
    }

    /**
     * Refill our cache from disk after a cache miss. Will prefetch Features overlapping an additional
     * queryLookaheadBases bases after the end of the provided interval, in addition to those overlapping
     * the interval itself.
     *
     * Calling this has the side effect of invalidating (closing) any currently-open iteration over
     * this data source.
     *
     * @param interval the query interval that produced a cache miss
     */
    private void refillQueryCache( final SimpleInterval interval ) {
        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        // Expand the end of our query by the configured number of bases, in anticipation of probable future
        // queries with slightly larger start/stop positions.
        //
        // Note that it doesn't matter if we go off the end of the contig in the process, since
        // our reader's query operation is not aware of (and does not care about) contig boundaries.
        // Note: we use addExact to blow up on overflow rather than propagate negative results downstream
        final SimpleInterval queryInterval = new SimpleInterval(interval.getContig(), interval.getStart(), Math.addExact(interval.getEnd(), queryLookaheadBases));

        // Query iterator over our reader will be immediately closed after re-populating our cache
        try ( CloseableTribbleIterator<T> queryIter = featureReader.query(queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd()) ) {
            queryCache.fill(queryIter, queryInterval);
        }
        catch ( IOException e ) {
            throw new GATKException("Error querying file " + featureInput + " over interval " + interval, e);
        }
    }

    /**
     * Get the logical name of this data source.
     *
     * @return the logical name of this data source
     */
    public String getName() {
        return featureInput.getName();
    }

    /**
     * Gets the header associated with this data source
     *
     * @return header associated with this data source as an Object
     */
    public Object getHeader() {
        return featureReader.getHeader();
    }

    /**
     * Permanently close this data source, invalidating any open iteration over it, and making it invalid for future
     * iterations and queries.
     */
    @Override
    public void close() {
        closeOpenIterationIfNecessary();

        logger.debug(String.format("Cache statistics for FeatureInput %s:", featureInput));
        queryCache.printCacheStatistics();

        try {
            if ( featureReader != null ) {
                featureReader.close();
            }
        }
        catch ( IOException e ) {
            throw new GATKException("Error closing Feature reader for input " + featureInput);
        }
    }

    /**
     * Close the iterator currently open over this data source, if there is one.
     */
    private void closeOpenIterationIfNecessary() {
        if ( currentIterator != null ) {
            currentIterator.close();
            currentIterator = null;
        }
    }
}
